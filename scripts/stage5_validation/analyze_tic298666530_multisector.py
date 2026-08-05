#!/usr/bin/env python3
"""Exploratory multi-sector fold for TIC 298666530.

This is a candidate-specific quick-look, not the frozen Stage 4 multi-sector
search contract. It consumes sector-level A2v1 FITS products, applies the same
per-sector internal-quality and upper-tail cleaning used by the transparent
BLS baseline, refines the supplied period in the combined data, and writes
provenance, diagnostic tables, and publication-style PNG/PDF figures.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from astropy.io import fits  # noqa: E402
from astropy.table import Table  # noqa: E402
from astropy.timeseries import BoxLeastSquares  # noqa: E402

from twirl.io.hlsp import BJDREFI, HLSPLightCurve, read_hlsp  # noqa: E402
from twirl.plotting.style import apply_twirl_style, get_ordered_palette  # noqa: E402
from twirl.search.bls import BLSConfig, prepare_bls_inputs_from_arrays  # noqa: E402


TIC = 298666530
GAIA_DR3 = 2154363086096593536
SMALL_APERTURE = "DET_FLUX_ADP_SML"
PRIMARY_APERTURE = "DET_FLUX_ADP"
APERTURES = (SMALL_APERTURE, PRIMARY_APERTURE)
EDGE_WARN_SECTORS = {79}


@dataclass(frozen=True)
class PreparedSector:
    sector: int
    camera: int
    ccd: int
    time: np.ndarray
    flux: np.ndarray
    error: np.ndarray
    orbitid: np.ndarray
    n_total: int
    n_quality: int
    n_kept: int
    robust_sigma: float


@dataclass(frozen=True)
class RefinedBLS:
    aperture: str
    period_d: float
    transit_time_btjd: float
    duration_d: float
    depth: float
    depth_err: float
    depth_snr: float
    power: float
    period_grid: np.ndarray
    power_grid: np.ndarray


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fits-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--combined-fits", type=Path, required=True)
    parser.add_argument("--period-d", type=float, default=0.5733)
    parser.add_argument("--period-window-d", type=float, default=0.003)
    parser.add_argument("--n-periods", type=int, default=60_001)
    parser.add_argument(
        "--sector-n-periods",
        type=int,
        default=12_001,
        help="Period samples for each sector's independent local-prior BLS check.",
    )
    parser.add_argument(
        "--durations-min",
        type=float,
        nargs="+",
        default=[8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0],
    )
    return parser.parse_args(argv)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def robust_sigma(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return float("nan")
    median = np.median(finite)
    return float(1.4826 * np.median(np.abs(finite - median)))


def load_sector_light_curves(fits_root: Path) -> list[HLSPLightCurve]:
    paths = sorted(Path(fits_root).rglob("*.fits"))
    if not paths:
        raise ValueError(f"no FITS files found under {fits_root}")

    light_curves: list[HLSPLightCurve] = []
    seen_sectors: set[int] = set()
    for path in paths:
        lc = read_hlsp(path, columns=APERTURES)
        if lc is None:
            raise ValueError(f"failed to read {path}")
        if lc.tic != TIC:
            continue
        missing = sorted(set(APERTURES) - set(lc.flux))
        if missing:
            raise ValueError(f"{path} is missing {missing}")
        if lc.sector in seen_sectors:
            raise ValueError(f"duplicate sector {lc.sector} under {fits_root}")
        seen_sectors.add(lc.sector)
        light_curves.append(lc)
    if not light_curves:
        raise ValueError(f"no TIC {TIC} FITS files found under {fits_root}")
    return sorted(light_curves, key=lambda lc: lc.sector)


def prepare_sectors(
    light_curves: list[HLSPLightCurve],
    aperture: str,
) -> list[PreparedSector]:
    config = BLSConfig(apertures=(aperture,))
    prepared: list[PreparedSector] = []
    for lc in light_curves:
        result = prepare_bls_inputs_from_arrays(
            time=lc.time,
            flux=lc.flux[aperture],
            quality=lc.quality,
            orbitid=lc.orbitid,
            cfg=config,
        )
        if not result.ready:
            raise ValueError(
                f"S{lc.sector} {aperture} BLS preparation failed: {result.status}"
            )
        sigma = robust_sigma(result.normalized_flux)
        if not np.isfinite(sigma) or sigma <= 0:
            raise ValueError(f"S{lc.sector} {aperture} has invalid robust scatter")
        prepared.append(
            PreparedSector(
                sector=lc.sector,
                camera=lc.cam,
                ccd=lc.ccd,
                time=result.time,
                flux=result.normalized_flux,
                error=np.full(result.time.shape, sigma, dtype=float),
                orbitid=result.orbitid,
                n_total=result.n_total,
                n_quality=result.n_cad_quality,
                n_kept=result.n_cad_kept,
                robust_sigma=sigma,
            )
        )
    return prepared


def concatenate_prepared(
    sectors: list[PreparedSector],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    time = np.concatenate([item.time for item in sectors])
    flux = np.concatenate([item.flux for item in sectors])
    error = np.concatenate([item.error for item in sectors])
    sector = np.concatenate(
        [np.full(item.time.shape, item.sector, dtype=np.int16) for item in sectors]
    )
    orbitid = np.concatenate([item.orbitid for item in sectors])
    order = np.argsort(time, kind="stable")
    return time[order], flux[order], error[order], sector[order], orbitid[order]


def refine_bls(
    sectors: list[PreparedSector],
    *,
    aperture: str,
    period_center_d: float,
    period_window_d: float,
    n_periods: int,
    durations_min: list[float],
) -> RefinedBLS:
    time, flux, error, _, _ = concatenate_prepared(sectors)
    period_grid = np.linspace(
        period_center_d - period_window_d,
        period_center_d + period_window_d,
        n_periods,
        dtype=float,
    )
    duration_grid = np.asarray(durations_min, dtype=float) / 1440.0
    result = BoxLeastSquares(time, flux, dy=error).power(
        period_grid,
        duration_grid,
        oversample=10,
    )
    index = int(np.nanargmax(result.power))
    return RefinedBLS(
        aperture=aperture,
        period_d=float(result.period[index]),
        transit_time_btjd=float(result.transit_time[index]),
        duration_d=float(result.duration[index]),
        depth=float(result.depth[index]),
        depth_err=float(result.depth_err[index]),
        depth_snr=float(result.depth_snr[index]),
        power=float(result.power[index]),
        period_grid=np.asarray(result.period, dtype=float),
        power_grid=np.asarray(result.power, dtype=float),
    )


def central_epoch(transit_time_btjd: float, period_d: float, time: np.ndarray) -> float:
    cycles = np.rint((np.median(time) - transit_time_btjd) / period_d)
    return float(transit_time_btjd + cycles * period_d)


def phase_days(time: np.ndarray, epoch_btjd: float, period_d: float) -> np.ndarray:
    return (time - epoch_btjd + 0.5 * period_d) % period_d - 0.5 * period_d


def binned_curve(
    x: np.ndarray,
    y: np.ndarray,
    *,
    edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    centers = 0.5 * (edges[:-1] + edges[1:])
    median = np.full(centers.shape, np.nan, dtype=float)
    error = np.full(centers.shape, np.nan, dtype=float)
    count = np.zeros(centers.shape, dtype=int)
    index = np.digitize(x, edges) - 1
    for bin_index in range(centers.size):
        values = y[index == bin_index]
        values = values[np.isfinite(values)]
        count[bin_index] = values.size
        if values.size == 0:
            continue
        median[bin_index] = np.median(values)
        scatter = robust_sigma(values)
        error[bin_index] = scatter / np.sqrt(values.size)
    return centers, median, error, count


def sector_metrics(
    by_aperture: dict[str, list[PreparedSector]],
    *,
    period_d: float,
    epoch_btjd: float,
    duration_d: float,
    local_period_center_d: float,
    local_period_window_d: float,
    local_n_periods: int,
    local_durations_min: list[float],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    sectors = [item.sector for item in by_aperture[SMALL_APERTURE]]
    local_period_grid = np.linspace(
        local_period_center_d - local_period_window_d,
        local_period_center_d + local_period_window_d,
        local_n_periods,
        dtype=float,
    )
    local_duration_grid = np.asarray(local_durations_min, dtype=float) / 1440.0
    for sector in sectors:
        row: dict[str, object] = {
            "sector": sector,
            "edge_warn": sector in EDGE_WARN_SECTORS,
        }
        for aperture, tag in (
            (SMALL_APERTURE, "small"),
            (PRIMARY_APERTURE, "primary"),
        ):
            item = next(value for value in by_aperture[aperture] if value.sector == sector)
            result = BoxLeastSquares(
                item.time,
                item.flux,
                dy=item.error,
            ).power(
                np.asarray([period_d]),
                np.asarray([duration_d]),
                oversample=20,
            )
            local_result = BoxLeastSquares(
                item.time,
                item.flux,
                dy=item.error,
            ).power(
                local_period_grid,
                local_duration_grid,
                oversample=10,
            )
            local_index = int(np.nanargmax(local_result.power))
            transit_time = float(result.transit_time[0])
            oc_d = float(phase_days(
                np.asarray([transit_time]),
                epoch_btjd,
                period_d,
            )[0])
            cadence_cycle = np.rint((item.time - epoch_btjd) / period_d).astype(int)
            local_phase = item.time - (epoch_btjd + cadence_cycle * period_d)
            observed_cycles = np.unique(
                cadence_cycle[np.abs(local_phase) <= 0.5 * duration_d]
            )
            row.update(
                {
                    f"{tag}_n_total": item.n_total,
                    f"{tag}_n_quality": item.n_quality,
                    f"{tag}_n_kept": item.n_kept,
                    f"{tag}_scatter": item.robust_sigma,
                    f"{tag}_depth": float(result.depth[0]),
                    f"{tag}_depth_err": float(result.depth_err[0]),
                    f"{tag}_depth_snr": float(result.depth_snr[0]),
                    f"{tag}_transit_time_btjd": transit_time,
                    f"{tag}_oc_min": oc_d * 1440.0,
                    f"{tag}_observed_transits": int(observed_cycles.size),
                    f"{tag}_local_period_d": float(
                        local_result.period[local_index]
                    ),
                    f"{tag}_local_duration_min": float(
                        local_result.duration[local_index] * 1440.0
                    ),
                    f"{tag}_local_depth": float(local_result.depth[local_index]),
                    f"{tag}_local_depth_snr": float(
                        local_result.depth_snr[local_index]
                    ),
                }
            )
            row["camera"] = item.camera
            row["ccd"] = item.ccd
        small_depth = float(row["small_depth"])
        row["primary_to_small_depth_ratio"] = (
            float(row["primary_depth"]) / small_depth
            if np.isfinite(small_depth) and small_depth != 0
            else np.nan
        )
        rows.append(row)
    return pd.DataFrame(rows).sort_values("sector").reset_index(drop=True)


def leave_one_sector_out_periods(
    sectors: list[PreparedSector],
    *,
    best_period_d: float,
    duration_d: float,
) -> np.ndarray:
    trial_periods = np.linspace(best_period_d - 0.0003, best_period_d + 0.0003, 3001)
    estimates: list[float] = []
    for omitted in [item.sector for item in sectors]:
        retained = [item for item in sectors if item.sector != omitted]
        time, flux, error, _, _ = concatenate_prepared(retained)
        result = BoxLeastSquares(time, flux, dy=error).power(
            trial_periods,
            np.asarray([duration_d]),
            oversample=10,
        )
        estimates.append(float(result.period[int(np.nanargmax(result.power))]))
    return np.asarray(estimates)


def write_combined_fits(
    light_curves: list[HLSPLightCurve],
    path: Path,
) -> None:
    rows: list[Table] = []
    for lc in light_curves:
        table = Table()
        table["TIME"] = lc.time
        table["SECTOR"] = np.full(len(lc.time), lc.sector, dtype=np.int16)
        table["ORBITID"] = lc.orbitid.astype(np.int16)
        table["QUALITY"] = lc.quality.astype(np.int64)
        table[SMALL_APERTURE] = lc.flux[SMALL_APERTURE].astype(np.float32)
        table[PRIMARY_APERTURE] = lc.flux[PRIMARY_APERTURE].astype(np.float32)
        rows.append(table)
    combined = Table(np.hstack([])) if not rows else Table(
        np.concatenate([row.as_array() for row in rows])
    )
    combined.sort("TIME")
    combined.meta.update(
        {
            "TICID": TIC,
            "GAIADR3": GAIA_DR3,
            "BJDREFI": BJDREFI,
            "PRODUCT": "A2v1 target quick-look",
            "QUALITY": "Internal TGLC quality only; no inferred external overlay",
        }
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    combined.write(path, overwrite=True)


def plot_summary(
    *,
    by_aperture: dict[str, list[PreparedSector]],
    refined: dict[str, RefinedBLS],
    metrics: pd.DataFrame,
    period_d: float,
    epoch_btjd: float,
    duration_d: float,
    out_dir: Path,
) -> None:
    template = apply_twirl_style("full_page")
    colors = {
        SMALL_APERTURE: get_ordered_palette(2, "viridis")[0],
        PRIMARY_APERTURE: get_ordered_palette(2, "viridis")[1],
    }
    fig, axes = plt.subplots(3, 2, figsize=(template["figsize"][0], 7.2))
    ax_period, ax_small, ax_primary, ax_depth, ax_oc, ax_ratio = axes.flat

    for aperture, label in (
        (SMALL_APERTURE, "1×1 ADP"),
        (PRIMARY_APERTURE, "3×3 ADP"),
    ):
        result = refined[aperture]
        power = result.power_grid
        relative = (power - np.nanmin(power)) / (
            np.nanmax(power) - np.nanmin(power)
        )
        ax_period.plot(
            result.period_grid,
            relative,
            color=colors[aperture],
            lw=1.0,
            label=label,
        )
        ax_period.axvline(result.period_d, color=colors[aperture], lw=0.8, ls="--")
    ax_period.set_xlabel("Period (d)")
    ax_period.set_ylabel("Relative BLS power")
    ax_period.legend(loc="best")
    ax_period.text(0.02, 0.95, "(a)", transform=ax_period.transAxes, va="top")

    phase_edges = np.linspace(-60.0, 60.0, 61)
    sector_palette = get_ordered_palette(len(metrics), "viridis")
    for aperture, axis, panel, label in (
        (SMALL_APERTURE, ax_small, "(b)", "1×1 ADP relative flux"),
        (PRIMARY_APERTURE, ax_primary, "(c)", "3×3 ADP relative flux"),
    ):
        all_phase: list[np.ndarray] = []
        all_flux: list[np.ndarray] = []
        for color, item in zip(sector_palette, by_aperture[aperture], strict=True):
            phase_min = phase_days(item.time, epoch_btjd, period_d) * 1440.0
            visible = np.abs(phase_min) <= 60.0
            axis.scatter(
                phase_min[visible],
                item.flux[visible],
                s=1.0,
                color=color,
                alpha=0.12,
                linewidths=0,
                rasterized=True,
            )
            centers, medians, _, counts = binned_curve(
                phase_min[visible],
                item.flux[visible],
                edges=phase_edges,
            )
            valid = counts >= 2
            axis.plot(centers[valid], medians[valid], color=color, lw=0.5, alpha=0.7)
            all_phase.append(phase_min[visible])
            all_flux.append(item.flux[visible])
        centers, medians, errors, counts = binned_curve(
            np.concatenate(all_phase),
            np.concatenate(all_flux),
            edges=phase_edges,
        )
        valid = counts >= 5
        axis.errorbar(
            centers[valid],
            medians[valid],
            yerr=errors[valid],
            color="black",
            lw=1.1,
            marker="o",
            ms=2.0,
            capsize=1.0,
            label="all-sector median",
        )
        axis.axvspan(
            -0.5 * duration_d * 1440.0,
            0.5 * duration_d * 1440.0,
            color="0.75",
            alpha=0.25,
            zorder=0,
        )
        axis.axhline(1.0, color="0.35", lw=0.6, ls=":")
        axis.set_xlim(-60.0, 60.0)
        axis.set_xlabel("Minutes from refined mid-transit")
        axis.set_ylabel(label)
        axis.text(0.02, 0.95, panel, transform=axis.transAxes, va="top")

    sectors = metrics["sector"].to_numpy()
    for tag, label, color, marker in (
        ("small", "1×1 ADP", colors[SMALL_APERTURE], "o"),
        ("primary", "3×3 ADP", colors[PRIMARY_APERTURE], "s"),
    ):
        ax_depth.errorbar(
            sectors,
            100.0 * metrics[f"{tag}_depth"],
            yerr=100.0 * metrics[f"{tag}_depth_err"],
            color=color,
            marker=marker,
            ms=3.0,
            lw=0.8,
            capsize=1.5,
            label=label,
        )
        ax_oc.plot(
            sectors,
            metrics[f"{tag}_oc_min"],
            color=color,
            marker=marker,
            ms=3.0,
            lw=0.8,
            label=label,
        )
    ax_depth.set_xlabel("TESS sector")
    ax_depth.set_ylabel("Fixed-ephemeris depth (%)")
    ax_depth.legend(loc="best")
    ax_depth.text(0.02, 0.95, "(d)", transform=ax_depth.transAxes, va="top")

    ax_oc.axhline(0.0, color="0.35", lw=0.7, ls=":")
    ax_oc.set_xlabel("TESS sector")
    ax_oc.set_ylabel("Sector phase offset (min)")
    ax_oc.text(0.02, 0.95, "(e)", transform=ax_oc.transAxes, va="top")

    ax_ratio.plot(
        sectors,
        metrics["primary_to_small_depth_ratio"],
        color="0.15",
        marker="o",
        ms=3.0,
        lw=0.8,
    )
    ax_ratio.axhline(1.0, color="0.35", lw=0.7, ls=":")
    ax_ratio.set_xlabel("TESS sector")
    ax_ratio.set_ylabel("3×3 / 1×1 depth ratio")
    ax_ratio.text(0.02, 0.95, "(f)", transform=ax_ratio.transAxes, va="top")

    fig.tight_layout(h_pad=1.1, w_pad=1.0)
    for suffix in ("png", "pdf"):
        fig.savefig(out_dir / f"tic298666530_multisector_summary.{suffix}", dpi=300)
    plt.close(fig)


def plot_sector_grid(
    *,
    sectors: list[PreparedSector],
    metrics: pd.DataFrame,
    period_d: float,
    epoch_btjd: float,
    duration_d: float,
    out_dir: Path,
) -> None:
    apply_twirl_style("full_page")
    n_columns = 4
    n_rows = int(np.ceil(len(sectors) / n_columns))
    fig, axes = plt.subplots(n_rows, n_columns, figsize=(7.1, 1.65 * n_rows), sharex=True)
    axes_array = np.atleast_1d(axes).ravel()
    edges = np.linspace(-45.0, 45.0, 46)
    palette = get_ordered_palette(len(sectors), "viridis")
    metrics_by_sector = metrics.set_index("sector")
    for axis, item, color in zip(axes_array, sectors, palette, strict=False):
        local_phase = phase_days(item.time, epoch_btjd, period_d) * 1440.0
        visible = np.abs(local_phase) <= 45.0
        axis.scatter(
            local_phase[visible],
            item.flux[visible],
            s=1.0,
            color="0.55",
            alpha=0.22,
            linewidths=0,
            rasterized=True,
        )
        centers, medians, errors, counts = binned_curve(
            local_phase[visible],
            item.flux[visible],
            edges=edges,
        )
        valid = counts >= 2
        axis.errorbar(
            centers[valid],
            medians[valid],
            yerr=errors[valid],
            color=color,
            lw=0.9,
            marker="o",
            ms=1.8,
            capsize=1.0,
        )
        axis.axvspan(
            -0.5 * duration_d * 1440.0,
            0.5 * duration_d * 1440.0,
            color=color,
            alpha=0.12,
            zorder=0,
        )
        axis.axhline(1.0, color="0.35", lw=0.5, ls=":")
        depth = 100.0 * metrics_by_sector.loc[item.sector, "small_depth"]
        axis.set_title(f"S{item.sector}  depth={depth:.1f}%")
        axis.set_xlim(-45.0, 45.0)
    for axis in axes_array[len(sectors):]:
        axis.set_visible(False)
    for axis in axes_array[-n_columns:]:
        if axis.get_visible():
            axis.set_xlabel("Minutes from mid-transit")
    for row in range(n_rows):
        axis = axes_array[row * n_columns]
        axis.set_ylabel("1×1 ADP flux")
    fig.tight_layout(h_pad=0.7, w_pad=0.6)
    for suffix in ("png", "pdf"):
        fig.savefig(out_dir / f"tic298666530_sector_folds.{suffix}", dpi=300)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if (
        args.period_window_d <= 0
        or args.n_periods < 3
        or args.sector_n_periods < 3
    ):
        raise ValueError("period refinement requires a positive window and >=3 periods")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    light_curves = load_sector_light_curves(args.fits_root)
    by_aperture = {
        aperture: prepare_sectors(light_curves, aperture)
        for aperture in APERTURES
    }
    refined = {
        aperture: refine_bls(
            by_aperture[aperture],
            aperture=aperture,
            period_center_d=args.period_d,
            period_window_d=args.period_window_d,
            n_periods=args.n_periods,
            durations_min=args.durations_min,
        )
        for aperture in APERTURES
    }
    anchor = refined[SMALL_APERTURE]
    small_time, _, _, _, _ = concatenate_prepared(by_aperture[SMALL_APERTURE])
    epoch_btjd = central_epoch(anchor.transit_time_btjd, anchor.period_d, small_time)

    metrics = sector_metrics(
        by_aperture,
        period_d=anchor.period_d,
        epoch_btjd=epoch_btjd,
        duration_d=anchor.duration_d,
        local_period_center_d=args.period_d,
        local_period_window_d=args.period_window_d,
        local_n_periods=args.sector_n_periods,
        local_durations_min=args.durations_min,
    )
    metrics.to_csv(args.out_dir / "sector_metrics.csv", index=False)

    jackknife = leave_one_sector_out_periods(
        by_aperture[SMALL_APERTURE],
        best_period_d=anchor.period_d,
        duration_d=anchor.duration_d,
    )
    pd.DataFrame(
        {
            "omitted_sector": metrics["sector"],
            "period_d": jackknife,
        }
    ).to_csv(args.out_dir / "leave_one_sector_out_periods.csv", index=False)

    provenance_rows = []
    for lc in light_curves:
        provenance_rows.append(
            {
                "sector": lc.sector,
                "camera": lc.cam,
                "ccd": lc.ccd,
                "n_rows": len(lc.time),
                "path": str(lc.path.resolve()),
                "bytes": lc.path.stat().st_size,
                "sha256": sha256_file(lc.path),
            }
        )
    pd.DataFrame(provenance_rows).to_csv(
        args.out_dir / "input_provenance.csv",
        index=False,
    )

    periodogram = pd.DataFrame(
        {
            "period_d": refined[SMALL_APERTURE].period_grid,
            "small_power": refined[SMALL_APERTURE].power_grid,
            "primary_power": refined[PRIMARY_APERTURE].power_grid,
        }
    )
    periodogram.to_csv(args.out_dir / "refined_periodogram.csv", index=False)

    write_combined_fits(light_curves, args.combined_fits)
    plot_summary(
        by_aperture=by_aperture,
        refined=refined,
        metrics=metrics,
        period_d=anchor.period_d,
        epoch_btjd=epoch_btjd,
        duration_d=anchor.duration_d,
        out_dir=args.out_dir,
    )
    plot_sector_grid(
        sectors=by_aperture[SMALL_APERTURE],
        metrics=metrics,
        period_d=anchor.period_d,
        epoch_btjd=epoch_btjd,
        duration_d=anchor.duration_d,
        out_dir=args.out_dir,
    )

    period_step = float(np.diff(anchor.period_grid[:2])[0])
    summary = {
        "tic": TIC,
        "gaia_dr3_source_id": GAIA_DR3,
        "status": "exploratory_multisector_quicklook",
        "sectors": metrics["sector"].astype(int).tolist(),
        "n_sectors": int(len(metrics)),
        "n_sector_fits": int(len(light_curves)),
        "baseline_d": float(np.ptp(small_time)),
        "period_prior_d": float(args.period_d),
        "period_refinement_window_d": float(args.period_window_d),
        "period_grid_step_d": period_step,
        "period_d": anchor.period_d,
        "period_jackknife_median_d": float(np.median(jackknife)),
        "period_jackknife_mad_d": float(1.4826 * np.median(
            np.abs(jackknife - np.median(jackknife))
        )),
        "period_jackknife_min_d": float(np.min(jackknife)),
        "period_jackknife_max_d": float(np.max(jackknife)),
        "epoch_btjd": epoch_btjd,
        "epoch_bjd": epoch_btjd + BJDREFI,
        "duration_min": anchor.duration_d * 1440.0,
        "small_depth": anchor.depth,
        "small_depth_err": anchor.depth_err,
        "small_depth_snr": anchor.depth_snr,
        "primary_period_d": refined[PRIMARY_APERTURE].period_d,
        "primary_period_delta_d": refined[PRIMARY_APERTURE].period_d - anchor.period_d,
        "primary_depth": refined[PRIMARY_APERTURE].depth,
        "primary_depth_err": refined[PRIMARY_APERTURE].depth_err,
        "primary_depth_snr": refined[PRIMARY_APERTURE].depth_snr,
        "small_sector_depth_min": float(metrics["small_depth"].min()),
        "small_sector_depth_max": float(metrics["small_depth"].max()),
        "primary_sector_depth_min": float(metrics["primary_depth"].min()),
        "primary_sector_depth_max": float(metrics["primary_depth"].max()),
        "edge_warn_sectors": sorted(EDGE_WARN_SECTORS.intersection(
            set(metrics["sector"].astype(int))
        )),
        "combined_fits": str(args.combined_fits.resolve()),
        "quality_scope": (
            "FITS-internal TGLC quality only; this candidate quick-look does not "
            "apply a separately hash-bound external cadence-quality overlay"
        ),
        "claim_boundary": (
            "Algorithmic multi-sector recurrence diagnostic, not astrophysical "
            "confirmation and not the frozen Stage 4 multi-sector search contract"
        ),
    }
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (args.out_dir / "README.md").write_text(
        "\n".join(
            [
                "# TIC 298666530 multi-sector quick-look",
                "",
                (
                    f"- Inputs: `{len(light_curves)}` A2v1 sector FITS spanning "
                    f"S{min(summary['sectors'])}--S{max(summary['sectors'])} "
                    f"with gaps; exact paths and hashes are in [input_provenance.csv]"
                    "(input_provenance.csv)."
                ),
                (
                    f"- Refined 1×1-ADP ephemeris: `P = {anchor.period_d:.9f} d`, "
                    f"`T0 = {summary['epoch_bjd']:.8f} BJD`, "
                    f"`duration = {summary['duration_min']:.2f} min`."
                ),
                (
                    "- Independent wide-window weighted BLS optima: "
                    f"`P = {anchor.period_d:.9f} d`, "
                    f"`depth = {100.0 * anchor.depth:.2f}%` (1×1 ADP); "
                    f"`P = {refined[PRIMARY_APERTURE].period_d:.9f} d`, "
                    f"`depth = {100.0 * refined[PRIMARY_APERTURE].depth:.2f}%` "
                    "(3×3 ADP). A period disagreement is diagnostic rather than "
                    "a consensus ephemeris."
                ),
                (
                    "- Sector-level depths vary substantially and are tabulated "
                    "rather than interpreted as one physical depth; S79 is also "
                    "retained with its detector-edge warning explicit."
                ),
                (
                    "- The sector table also records an independent local-prior "
                    "BLS period, duration, depth, and depth S/N for each aperture; "
                    "these columns expose sectors that do not reproduce the "
                    "combined ephemeris."
                ),
                (
                    "- Cleaning matches the transparent per-sector BLS baseline "
                    "(finite, internal `QUALITY==0`, per-sector median normalization, "
                    "upper-tail clipping); the combined refinement additionally uses "
                    "per-sector robust-MAD weights."
                ),
                (
                    "- This is an exploratory candidate diagnostic. It does not "
                    "replace the pending frozen multi-sector search/false-alarm "
                    "contract and is not an astrophysical confirmation."
                ),
                "",
                "Outputs: [summary figure PNG](tic298666530_multisector_summary.png), "
                "[summary figure PDF](tic298666530_multisector_summary.pdf), "
                "[sector folds PNG](tic298666530_sector_folds.png), "
                "[sector metrics](sector_metrics.csv), and "
                "[summary JSON](summary.json).",
                "",
            ]
        ),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
