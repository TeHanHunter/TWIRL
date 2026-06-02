#!/usr/bin/env python
"""Plot a TGLC saved-aperture pixel light-curve map for one target.

The input must be a single-target FITS file produced by
``tglc.quick_lc.tglc_lc(..., save_aper=True)``. The primary HDU is expected to
hold a ``(cadence, y, x)`` aperture cube.
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

try:
    from scipy.ndimage import median_filter
except ModuleNotFoundError:
    median_filter = None

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from twirl.plotting.style import apply_twirl_style  # noqa: E402

BJD_MINUS_BTJD = 2_457_000.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fits",
        type=Path,
        required=True,
        help="TGLC FITS file produced with save_aper=True.",
    )
    parser.add_argument(
        "--out-prefix",
        type=Path,
        default=Path("reports/stage1_lightcurves/pixel_maps/wd1856_s56_tglc_pixel_map"),
        help="Output path without extension. PNG, PDF, and JSON sidecar are written.",
    )
    parser.add_argument("--period-d", type=float, default=1.407983, help="Transit period in days.")
    parser.add_argument(
        "--t0-bjd",
        type=float,
        default=2459825.473398,
        help="Transit epoch as BJD. Converted to BTJD automatically if the time column is BTJD.",
    )
    parser.add_argument(
        "--duration-min",
        type=float,
        default=6.0,
        help="Nominal event duration in minutes for the shaded window and depth map.",
    )
    parser.add_argument("--n-bins", type=int, default=80, help="Number of phase bins.")
    parser.add_argument(
        "--phase-window-min",
        type=float,
        default=0.0,
        help="If positive, zoom to +/- this many minutes around the event instead of plotting phase 0-1.",
    )
    parser.add_argument(
        "--pixel-detrend-window-d",
        type=float,
        default=0.3,
        help="Running-median detrend window for pixel and aperture panels. Set <=0 to disable.",
    )
    parser.add_argument(
        "--clip-sigma",
        type=float,
        default=7.0,
        help="MAD sigma clip applied to each pixel before binning and measuring event depth.",
    )
    parser.add_argument(
        "--pixel-normalization",
        choices=("aperture", "pixel"),
        default="aperture",
        help=(
            "Normalize per-pixel panels by the total aperture median after subtracting each "
            "pixel median, or by each pixel's own median."
        ),
    )
    parser.add_argument(
        "--source-pickle",
        type=Path,
        default=None,
        help=(
            "Optional TGLC source pickle from the same quick_lc run. If given, "
            "nearby Gaia DR3 source positions are overlaid on the aperture map."
        ),
    )
    parser.add_argument(
        "--neighbor-gmag-max",
        type=float,
        default=21.0,
        help="Faint Gaia G-band limit for source markers on the aperture map.",
    )
    parser.add_argument(
        "--label-neighbor-gmag-max",
        type=float,
        default=18.5,
        help="Faint Gaia G-band limit for text labels on non-target source markers.",
    )
    parser.add_argument(
        "--show-model-contours",
        action="store_true",
        help="Also draw TGLC target and field-star model contours from the FITS image extension.",
    )
    parser.add_argument(
        "--right-panel",
        choices=("pixel-grid", "aperture-fold"),
        default="pixel-grid",
        help=(
            "Right-hand diagnostic: a 5x5 grid of folded pixel light curves, "
            "or one larger phase-folded aperture/PSF transit panel."
        ),
    )
    parser.add_argument(
        "--show-aperture-overlay",
        action="store_true",
        help="In pixel-grid mode, overlay the total calibrated aperture flux in each pixel panel.",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional figure title. By default a compact target label is inferred from the FITS header.",
    )
    return parser.parse_args()


def output_base(path: Path) -> Path:
    return path.with_suffix("") if path.suffix else path


def find_column(table, *names: str) -> np.ndarray:
    lower_to_name = {name.lower(): name for name in table.names}
    for name in names:
        actual = lower_to_name.get(name.lower())
        if actual is not None:
            return np.asarray(table[actual])
    raise KeyError(f"None of the requested columns are present: {names}")


def optional_column(table, *names: str) -> np.ndarray | None:
    try:
        return find_column(table, *names)
    except KeyError:
        return None


def robust_mad(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    if finite.sum() == 0:
        return float("nan")
    med = np.nanmedian(values[finite])
    return float(1.4826 * np.nanmedian(np.abs(values[finite] - med)))


def sigma_clip_good(values: np.ndarray, good: np.ndarray, sigma: float) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    use = good & np.isfinite(values)
    if sigma <= 0 or use.sum() < 10:
        return use
    med = float(np.nanmedian(values[use]))
    scatter = robust_mad(values[use])
    if not np.isfinite(scatter) or scatter <= 0:
        return use
    return use & (np.abs(values - med) <= sigma * scatter)


def normalize_relative(flux: np.ndarray, good: np.ndarray) -> tuple[np.ndarray, str]:
    flux = np.asarray(flux, dtype=float)
    ref = flux[good & np.isfinite(flux)]
    if ref.size == 0:
        return np.full_like(flux, np.nan, dtype=float), "none"
    med = float(np.nanmedian(ref))
    if np.isfinite(med) and abs(med) > 1e-12:
        return flux / med - 1.0, "median"
    center = float(np.nanmedian(ref))
    scale = robust_mad(ref)
    if not np.isfinite(scale) or scale <= 0:
        scale = float(np.nanpercentile(np.abs(ref - center), 95))
    if not np.isfinite(scale) or scale <= 0:
        return flux - center, "centered"
    return (flux - center) / scale, "mad"


def normalize_pixel_flux(
    flux: np.ndarray,
    good: np.ndarray,
    mode: str,
    aperture_scale: float,
) -> tuple[np.ndarray, str]:
    if mode == "pixel":
        return normalize_relative(flux, good)
    flux = np.asarray(flux, dtype=float)
    ref = flux[good & np.isfinite(flux)]
    if ref.size == 0:
        return np.full_like(flux, np.nan, dtype=float), "none"
    center = float(np.nanmedian(ref))
    if np.isfinite(aperture_scale) and abs(aperture_scale) > 1e-12:
        return (flux - center) / aperture_scale, "aperture"
    return normalize_relative(flux, good)


def infer_epoch_in_time_system(time: np.ndarray, t0_bjd: float) -> tuple[float, str]:
    med_time = float(np.nanmedian(time[np.isfinite(time)]))
    if med_time > 100_000:
        return t0_bjd, "BJD"
    return t0_bjd - BJD_MINUS_BTJD, "BTJD"


def event_phase(time: np.ndarray, period_d: float, t0_time: float) -> np.ndarray:
    # Shift by half a cycle so the event sits at phase 0.5 in every panel.
    return ((time - t0_time) / period_d + 0.5) % 1.0


def phase_distance_from_event(phase: np.ndarray) -> np.ndarray:
    return np.abs(phase - 0.5)


def binned_median(phase: np.ndarray, flux: np.ndarray, good: np.ndarray, n_bins: int):
    return binned_median_edges(phase, flux, good, np.linspace(0.0, 1.0, n_bins + 1))


def binned_median_edges(x_value: np.ndarray, flux: np.ndarray, good: np.ndarray, edges: np.ndarray):
    centers = 0.5 * (edges[:-1] + edges[1:])
    n_bins = len(centers)
    med = np.full(n_bins, np.nan, dtype=float)
    err = np.full(n_bins, np.nan, dtype=float)
    count = np.zeros(n_bins, dtype=int)
    use = good & np.isfinite(x_value) & np.isfinite(flux)
    which = np.searchsorted(edges, x_value[use], side="right") - 1
    which = np.clip(which, 0, n_bins - 1)
    values = flux[use]
    for idx in range(n_bins):
        sub = values[which == idx]
        count[idx] = sub.size
        if sub.size:
            med[idx] = np.nanmedian(sub)
            scatter = robust_mad(sub)
            if np.isfinite(scatter):
                err[idx] = scatter / np.sqrt(sub.size)
    return centers, med, err, count


def running_median_detrend(
    flux: np.ndarray,
    time: np.ndarray,
    good: np.ndarray,
    window_d: float,
) -> np.ndarray:
    flux = np.asarray(flux, dtype=float)
    if window_d <= 0 or median_filter is None:
        return np.array(flux, copy=True)
    finite = np.isfinite(flux)
    if finite.sum() < 5:
        return np.full_like(flux, np.nan, dtype=float)
    time_good = time[good & np.isfinite(time)]
    if time_good.size < 2:
        return np.array(flux, copy=True)
    cadence_d = float(np.nanmedian(np.diff(np.sort(time_good))))
    if not np.isfinite(cadence_d) or cadence_d <= 0:
        return np.array(flux, copy=True)
    window_cadences = max(5, int(round(window_d / cadence_d)))
    if window_cadences % 2 == 0:
        window_cadences += 1
    filled = np.array(flux, copy=True)
    idx = np.arange(len(filled))
    filled[~finite] = np.interp(idx[~finite], idx[finite], filled[finite])
    trend = median_filter(filled, size=window_cadences, mode="nearest")
    return flux - trend


def read_tglc_saveaper(path: Path):
    with fits.open(path, memmap=False) as hdul:
        if hdul[0].data is None:
            raise RuntimeError("Primary HDU has no aperture cube. Re-run TGLC with save_aper=True.")
        cube = np.asarray(hdul[0].data, dtype=float)
        table = hdul[1].data
        header = hdul[0].header.copy()
        model_image = None
        if len(hdul) > 2 and hdul[2].data is not None:
            model_image = np.asarray(hdul[2].data, dtype=float)
    if cube.ndim != 3:
        raise RuntimeError(f"Expected aperture cube with shape (time, y, x), got {cube.shape}")
    return cube, table, header, model_image


def make_quality_mask(table, time: np.ndarray) -> np.ndarray:
    tess_flags = optional_column(table, "TESS_flags", "TESS_FLAG", "QUALITY")
    tglc_flags = optional_column(table, "TGLC_flags", "TGLC_FLAG")
    quality = np.zeros(len(time), dtype=int)
    if tess_flags is not None:
        quality |= np.asarray(tess_flags, dtype=int)
    if tglc_flags is not None:
        quality |= np.asarray(tglc_flags, dtype=int)
    return (quality == 0) & np.isfinite(time)


def measure_depth(
    rel_flux: np.ndarray,
    good: np.ndarray,
    distance: np.ndarray,
    half_width_phase: float,
    n_bins: int,
) -> tuple[float, float]:
    in_event = good & np.isfinite(rel_flux) & (distance <= half_width_phase)
    oot = good & np.isfinite(rel_flux) & (distance >= max(4.0 * half_width_phase, 2.0 / n_bins))
    if in_event.sum() < 3 or oot.sum() < 10:
        return float("nan"), float("nan")
    depth = float(np.nanmedian(rel_flux[oot]) - np.nanmedian(rel_flux[in_event]))
    scatter = robust_mad(rel_flux[oot])
    sigma = depth / scatter if np.isfinite(scatter) and scatter > 0 else float("nan")
    return depth, sigma


def add_model_contours(ax, model_image: np.ndarray | None) -> None:
    if model_image is None or model_image.ndim < 3 or model_image.shape[1:] != (5, 5):
        return
    target = model_image[0]
    field = model_image[1]
    yy, xx = np.mgrid[: target.shape[0], : target.shape[1]]
    if np.isfinite(target).any() and np.nanmax(target) > 0:
        levels = np.nanmax(target) * np.array([0.15, 0.4, 0.7])
        ax.contour(xx, yy, target, levels=levels, colors="black", linewidths=0.9)
    if np.isfinite(field).any() and np.nanmax(field) > 0:
        levels = np.nanmax(field) * np.array([0.2, 0.5, 0.8])
        ax.contour(xx, yy, field, levels=levels, colors="#b2182b", linestyles="--", linewidths=0.8)


def _designation_to_source_id(value) -> int | None:
    if value is None:
        return None
    text = value.decode() if isinstance(value, bytes) else str(value)
    try:
        return int(text.split()[-1])
    except (TypeError, ValueError):
        return None


def load_gaia_neighbors(
    source_pickle: Path | None,
    header,
    star_x: float,
    star_y: float,
    ny: int,
    nx: int,
    gmag_max: float,
) -> tuple[list[dict], dict]:
    if source_pickle is None:
        return [], {}
    with source_pickle.open("rb") as fh:
        source = pickle.load(fh)

    gaia = source.gaia
    sector = int(header.get("SECTOR", getattr(source, "sector", 0)))
    x_col = f"sector_{sector}_x"
    y_col = f"sector_{sector}_y"
    if x_col not in gaia.colnames or y_col not in gaia.colnames:
        raise RuntimeError(f"Gaia source table lacks {x_col}/{y_col} columns.")

    target_id = _designation_to_source_id(f"Gaia DR3 {header.get('GAIADR3')}")
    ids = np.array([_designation_to_source_id(value) for value in gaia["DESIGNATION"]], dtype=object)
    target_mask = ids == target_id
    if not np.any(target_mask):
        raise RuntimeError(f"Target Gaia DR3 source {target_id} is not in {source_pickle}.")

    target_idx = int(np.flatnonzero(target_mask)[0])
    target_sector_x = float(gaia[x_col][target_idx])
    target_sector_y = float(gaia[y_col][target_idx])
    origin_x = target_sector_x - float(star_x)
    origin_y = target_sector_y - float(star_y)

    x_cut = np.asarray(gaia[x_col], dtype=float) - origin_x
    y_cut = np.asarray(gaia[y_col], dtype=float) - origin_y
    gmag = np.asarray(gaia["phot_g_mean_mag"], dtype=float)
    tess_mag = np.asarray(gaia["tess_mag"], dtype=float) if "tess_mag" in gaia.colnames else np.full(len(gaia), np.nan)
    in_cut = (
        np.isfinite(x_cut)
        & np.isfinite(y_cut)
        & np.isfinite(gmag)
        & (x_cut >= -0.5)
        & (x_cut <= nx - 0.5)
        & (y_cut >= -0.5)
        & (y_cut <= ny - 0.5)
        & (gmag <= gmag_max)
    )

    neighbors = []
    for idx in np.flatnonzero(in_cut):
        source_id = ids[idx]
        neighbors.append(
            {
                "designation": str(gaia["DESIGNATION"][idx]),
                "source_id": None if source_id is None else int(source_id),
                "x": float(x_cut[idx]),
                "y": float(y_cut[idx]),
                "phot_g_mean_mag": float(gmag[idx]),
                "tess_mag": None if not np.isfinite(tess_mag[idx]) else float(tess_mag[idx]),
                "is_target": bool(source_id == target_id),
            }
        )
    neighbors.sort(key=lambda item: (not item["is_target"], item["phot_g_mean_mag"]))

    geometry = {
        "source_pickle": str(source_pickle),
        "target_gaia_dr3": None if target_id is None else int(target_id),
        "target_sector_x": target_sector_x,
        "target_sector_y": target_sector_y,
        "crop_origin_x": origin_x,
        "crop_origin_y": origin_y,
        "star_x_from_header": float(star_x),
        "star_y_from_header": float(star_y),
    }
    return neighbors, geometry


def add_gaia_neighbor_markers(
    ax,
    neighbors: list[dict],
    label_neighbor_gmag_max: float,
) -> None:
    if not neighbors:
        return

    target = [item for item in neighbors if item["is_target"]]
    background = [item for item in neighbors if not item["is_target"]]

    if background:
        x = np.array([item["x"] for item in background], dtype=float)
        y = np.array([item["y"] for item in background], dtype=float)
        gmag = np.array([item["phot_g_mean_mag"] for item in background], dtype=float)
        sizes = np.clip(220.0 * 10 ** (-0.4 * (gmag - 12.0)), 24.0, 145.0)
        ax.scatter(
            x,
            y,
            marker="o",
            s=sizes,
            facecolors="none",
            edgecolors="#2166ac",
            linewidths=1.25,
            label="Gaia DR3 background",
            zorder=5,
        )
        for item in background:
            if item["phot_g_mean_mag"] <= label_neighbor_gmag_max:
                x0, x1 = ax.get_xlim()
                y0, y1 = ax.get_ylim()
                near_right = item["x"] > x0 + 0.72 * (x1 - x0)
                near_top = item["y"] > y0 + 0.72 * (y1 - y0)
                ax.text(
                    item["x"] + (-0.06 if near_right else 0.06),
                    item["y"] + (-0.06 if near_top else 0.06),
                    f"G={item['phot_g_mean_mag']:.1f}",
                    color="#2166ac",
                    fontsize=7,
                    ha="right" if near_right else "left",
                    va="top" if near_top else "bottom",
                    zorder=6,
                )

    if target:
        item = target[0]
        ax.scatter(
            [item["x"]],
            [item["y"]],
            marker="x",
            s=105,
            color="#b2182b",
            linewidths=2.0,
            label="WD 1856 Gaia DR3",
            zorder=7,
        )
        ax.text(
            item["x"] + 0.12,
            item["y"] + 0.12,
            f"WD 1856\nG={item['phot_g_mean_mag']:.1f}",
            color="#b2182b",
            fontsize=8,
            ha="left",
            va="bottom",
            zorder=8,
        )


def normalize_to_oot(
    flux: np.ndarray,
    good: np.ndarray,
    phase_min: np.ndarray,
    duration_min: float,
    phase_window_min: float,
) -> tuple[np.ndarray, float]:
    flux = np.asarray(flux, dtype=float)
    event_half_width = 0.5 * duration_min
    oot_inner = max(2.0 * event_half_width, 10.0)
    in_window = good & np.isfinite(flux) & np.isfinite(phase_min) & (np.abs(phase_min) <= phase_window_min)
    oot = in_window & (np.abs(phase_min) >= oot_inner)
    if oot.sum() < 10:
        oot = good & np.isfinite(flux)
    baseline = float(np.nanmedian(flux[oot])) if oot.sum() else float("nan")
    if np.isfinite(baseline) and abs(baseline) > 1e-12:
        return flux / baseline, baseline
    rel, _ = normalize_relative(flux, good)
    return rel + 1.0, baseline


def plot_folded_transit_panel(
    ax,
    table,
    time: np.ndarray,
    good: np.ndarray,
    phase_min: np.ndarray,
    x_edges: np.ndarray,
    phase_window_min: float,
    duration_min: float,
    clip_sigma: float,
) -> dict:
    series_specs = [
        ("cal_aper_flux", "calibrated aperture", "black", "o", 0.13),
        ("cal_psf_flux", "calibrated PSF", "#b2182b", "s", 0.08),
    ]
    plotted = []
    binned_values = []
    for column, label, color, marker, raw_alpha in series_specs:
        flux = optional_column(table, column)
        if flux is None:
            continue
        flux = np.asarray(flux, dtype=float)
        rel_flux, baseline = normalize_to_oot(flux, good, phase_min, duration_min, phase_window_min)
        panel_good = (
            good
            & np.isfinite(time)
            & np.isfinite(phase_min)
            & np.isfinite(rel_flux)
            & (np.abs(phase_min) <= phase_window_min)
        )
        panel_good = sigma_clip_good(rel_flux, panel_good, clip_sigma)
        ax.scatter(
            phase_min[panel_good],
            rel_flux[panel_good],
            s=4.5,
            color=color,
            alpha=raw_alpha,
            linewidths=0,
            rasterized=True,
        )
        centers, med, err, count = binned_median_edges(phase_min, rel_flux, panel_good, x_edges)
        use = np.isfinite(med) & (count >= 3)
        binned_values.append(med[use])
        ax.errorbar(
            centers[use],
            med[use],
            yerr=err[use],
            fmt=marker + "-",
            color=color,
            ecolor=color,
            elinewidth=0.7,
            capsize=0,
            markersize=4.2,
            linewidth=1.25,
            label=label,
            zorder=4,
        )
        in_event = panel_good & (np.abs(phase_min) <= 0.5 * duration_min)
        oot = panel_good & (np.abs(phase_min) >= max(duration_min, 10.0))
        depth = float(np.nanmedian(rel_flux[oot]) - np.nanmedian(rel_flux[in_event])) if in_event.sum() and oot.sum() else float("nan")
        plotted.append(
            {
                "column": column,
                "baseline": baseline,
                "depth": depth,
                "raw_points": int(panel_good.sum()),
                "binned_points": int(use.sum()),
            }
        )

    ax.axvspan(-0.5 * duration_min, 0.5 * duration_min, color="0.88", lw=0, zorder=0)
    ax.axvline(0.0, color="0.25", lw=0.9, ls=":")
    ax.axhline(1.0, color="0.35", lw=0.8, ls=":")
    ax.set_xlim(-phase_window_min, phase_window_min)
    finite_bins = np.concatenate([values for values in binned_values if values.size]) if binned_values else np.array([])
    if finite_bins.size:
        ylo, yhi = np.nanpercentile(finite_bins[np.isfinite(finite_bins)], [1, 99])
        ypad = max(0.04, 0.15 * (yhi - ylo))
        ax.set_ylim(max(0.0, ylo - ypad), min(1.25, yhi + ypad))
    else:
        ax.set_ylim(0.55, 1.08)
    ax.set_xlabel("phase from T0 (min)")
    ax.set_ylabel("normalized flux")
    ax.text(
        0.02,
        0.98,
        "phase folded around WD 1856 b",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        color="0.2",
    )
    ax.legend(loc="lower right", fontsize=8, frameon=True)
    ax.grid(color="0.86", linewidth=0.7)
    return {"series": plotted}


def main() -> None:
    args = parse_args()
    if args.right_panel == "aperture-fold" and args.phase_window_min <= 0:
        args.phase_window_min = 90.0
    base = output_base(args.out_prefix)
    base.parent.mkdir(parents=True, exist_ok=True)

    cube, table, header, model_image = read_tglc_saveaper(args.fits)
    time = np.asarray(find_column(table, "time"), dtype=float)
    good = make_quality_mask(table, time)
    if cube.shape[0] != len(time):
        raise RuntimeError(f"Cube has {cube.shape[0]} cadences but table has {len(time)} rows.")

    t0_time, time_system = infer_epoch_in_time_system(time, args.t0_bjd)
    phase = event_phase(time, args.period_d, t0_time)
    phase_min = (((time - t0_time + 0.5 * args.period_d) % args.period_d) - 0.5 * args.period_d) * 1440.0
    distance = phase_distance_from_event(phase)
    half_width_phase = (args.duration_min / 1440.0) / (2.0 * args.period_d)
    if args.phase_window_min > 0:
        x_value = phase_min
        x_good = good & (np.abs(phase_min) <= args.phase_window_min)
        x_edges = np.linspace(-args.phase_window_min, args.phase_window_min, args.n_bins + 1)
        xlim = (-args.phase_window_min, args.phase_window_min)
        event_span = (-0.5 * args.duration_min, 0.5 * args.duration_min)
        xlabel = "phase from T0 (min)"
    else:
        x_value = phase
        x_good = good
        x_edges = np.linspace(0.0, 1.0, args.n_bins + 1)
        xlim = (0.0, 1.0)
        event_span = (0.5 - half_width_phase, 0.5 + half_width_phase)
        xlabel = "phase"

    cal_aper = optional_column(table, "cal_aper_flux")
    if cal_aper is None:
        cal_aper = np.nansum(cube, axis=(1, 2))
    cal_aper = np.asarray(cal_aper, dtype=float)
    aper_scale = float(np.nanmedian(cal_aper[good & np.isfinite(cal_aper)]))
    ny, nx = cube.shape[1:]
    raw_aperture = optional_column(table, "aperture_flux_raw")
    if raw_aperture is not None:
        pixel_aperture_flux = np.asarray(raw_aperture, dtype=float)
        pixel_aperture_scale_source = "aperture_flux_raw"
    else:
        star_x_int = int(round(float(header.get("STAR_X", nx // 2))))
        star_y_int = int(round(float(header.get("STAR_Y", ny // 2))))
        pixel_aperture_flux = np.nansum(
            cube[
                :,
                max(0, star_y_int - 1) : min(ny, star_y_int + 2),
                max(0, star_x_int - 1) : min(nx, star_x_int + 2),
            ],
            axis=(1, 2),
        )
        pixel_aperture_scale_source = "target_centered_3x3_sum"
    pixel_aperture_scale = float(
        np.nanmedian(pixel_aperture_flux[good & np.isfinite(pixel_aperture_flux)])
    )
    if not np.isfinite(pixel_aperture_scale) or abs(pixel_aperture_scale) <= 1e-12:
        pixel_aperture_scale = aper_scale
    if args.pixel_detrend_window_d > 0:
        aper_rel = running_median_detrend(
            cal_aper, time, good, args.pixel_detrend_window_d
        ) / aper_scale
    else:
        aper_rel, _ = normalize_relative(cal_aper, good)
    aper_phase, aper_bin, _, _ = binned_median_edges(x_value, aper_rel, x_good, x_edges)

    depths = np.full((ny, nx), np.nan, dtype=float)
    sigmas = np.full((ny, nx), np.nan, dtype=float)
    rel_pixels: dict[tuple[int, int], np.ndarray] = {}
    norm_modes: dict[tuple[int, int], str] = {}
    panel_good_masks: dict[tuple[int, int], np.ndarray] = {}

    for y in range(ny):
        for x in range(nx):
            pixel_flux = running_median_detrend(
                cube[:, y, x], time, good, args.pixel_detrend_window_d
            )
            rel, mode = normalize_pixel_flux(
                pixel_flux, good, args.pixel_normalization, pixel_aperture_scale
            )
            rel_pixels[(y, x)] = rel
            norm_modes[(y, x)] = mode
            pixel_good = sigma_clip_good(rel, good, args.clip_sigma)
            panel_good_masks[(y, x)] = pixel_good
            depths[y, x], sigmas[y, x] = measure_depth(
                rel, pixel_good, distance, half_width_phase, args.n_bins
            )

    all_rel = np.concatenate(
        [
            rel_pixels[(y, x)][x_good & panel_good_masks[(y, x)]]
            for y in range(ny)
            for x in range(nx)
        ]
    )
    y_abs = np.nanpercentile(np.abs(all_rel[np.isfinite(all_rel)]), 98.5)
    if np.isfinite(y_abs):
        if args.pixel_normalization == "aperture":
            y_abs = float(np.clip(1.3 * y_abs, 0.0015, 0.025))
        else:
            y_abs = float(np.clip(y_abs, 0.025, 0.18))
    else:
        y_abs = 0.004 if args.pixel_normalization == "aperture" else 0.10

    positive_depths = np.where(np.isfinite(depths), np.maximum(depths, 0.0), np.nan)
    depth_max = float(np.nanpercentile(positive_depths, 95)) if np.isfinite(positive_depths).any() else 1.0
    depth_max = max(depth_max, float(np.nanmax(positive_depths)) if np.isfinite(positive_depths).any() else 0.0)
    depth_max = depth_max if np.isfinite(depth_max) and depth_max > 0 else 1.0

    apply_twirl_style("full_page")
    fig = plt.figure(figsize=(15.5, 7.6), constrained_layout=False)
    outer = gridspec.GridSpec(1, 2, width_ratios=[0.86, 1.55], wspace=0.13, figure=fig)

    ax_map = fig.add_subplot(outer[0, 0])
    image = np.nanmedian(cube[good], axis=0)
    finite_image = image[np.isfinite(image)]
    vmin, vmax = np.nanpercentile(finite_image, [5, 98]) if finite_image.size else (None, None)
    im = ax_map.imshow(image, origin="lower", cmap="Greys", vmin=vmin, vmax=vmax)
    ax_map.set_xlabel("cube column (x)")
    ax_map.set_ylabel("cube row (y)")
    ax_map.set_xticks(np.r_[-0.5, np.arange(nx), nx - 0.5])
    ax_map.set_yticks(np.r_[-0.5, np.arange(ny), ny - 0.5])
    ax_map.set_xticks(np.arange(-0.5, nx, 1.0), minor=True)
    ax_map.set_yticks(np.arange(-0.5, ny, 1.0), minor=True)
    ax_map.set_xlim(-0.5, nx - 0.5)
    ax_map.set_ylim(-0.5, ny - 0.5)
    ax_map.grid(False)
    ax_map.grid(which="minor", color="0.75", linewidth=0.6)
    ax_map.tick_params(which="minor", length=0)
    star_x = header.get("STAR_X", np.nan)
    star_y = header.get("STAR_Y", np.nan)
    neighbors, gaia_geometry = load_gaia_neighbors(
        args.source_pickle,
        header,
        star_x,
        star_y,
        ny,
        nx,
        args.neighbor_gmag_max,
    )
    if args.show_model_contours:
        add_model_contours(ax_map, model_image)
    if neighbors:
        add_gaia_neighbor_markers(ax_map, neighbors, args.label_neighbor_gmag_max)
    elif np.isfinite(star_x) and np.isfinite(star_y):
        ax_map.scatter([star_x], [star_y], marker="x", s=95, color="#b2182b", linewidths=2.0)
        ax_map.text(star_x + 0.12, star_y + 0.12, "target", color="#b2182b", fontsize=8)
    if neighbors:
        ax_map.legend(loc="lower left", fontsize=7, frameon=True, borderpad=0.35, handletextpad=0.4)
    if args.right_panel == "pixel-grid":
        cbar_image = fig.colorbar(im, ax=ax_map, fraction=0.047, pad=0.03)
        cbar_image.set_label("")

    right_panel_summary = {}
    if args.right_panel == "pixel-grid":
        grid = gridspec.GridSpecFromSubplotSpec(ny, nx, subplot_spec=outer[0, 1], wspace=0.06, hspace=0.06)
        cmap = plt.get_cmap("viridis")
        norm = mpl.colors.Normalize(vmin=0.0, vmax=depth_max)
        mini_axes = []
        for y in range(ny):
            for x in range(nx):
                ax = fig.add_subplot(grid[y, x])
                mini_axes.append(ax)
                rel = rel_pixels[(y, x)]
                panel_good = x_good & panel_good_masks[(y, x)]
                ax.axvspan(event_span[0], event_span[1], color="0.90", lw=0, zorder=0)
                ax.scatter(
                    x_value[panel_good],
                    rel[panel_good],
                    s=2.2,
                    color="0.35",
                    alpha=0.11,
                    linewidths=0,
                    rasterized=True,
                    zorder=1,
                )
                centers, med, _, count = binned_median_edges(x_value, rel, panel_good, x_edges)
                use_bins = np.isfinite(med) & (count > 0)
                ax.plot(centers[use_bins], med[use_bins], color="black", lw=1.2, zorder=3)
                if args.show_aperture_overlay:
                    ap_use = np.isfinite(aper_bin)
                    ax.plot(aper_phase[ap_use], aper_bin[ap_use], color="#b2182b", lw=0.9, alpha=0.85, zorder=2)
                ax.axhline(0.0, color="0.5", ls=":", lw=0.7)
                ax.set_xlim(*xlim)
                ax.set_ylim(-y_abs, y_abs)
                depth = depths[y, x]
                edge_color = cmap(norm(max(depth, 0.0))) if np.isfinite(depth) else "0.55"
                for spine in ax.spines.values():
                    spine.set_color(edge_color)
                    spine.set_linewidth(2.0)
                label = f"({y},{x})"
                if np.isfinite(depth):
                    label += f"  d={depth:.3f}"
                ax.text(0.03, 0.93, label, transform=ax.transAxes, ha="left", va="top", fontsize=6)
                if y < ny - 1:
                    ax.set_xticklabels([])
                else:
                    ax.set_xlabel(xlabel, labelpad=0)
                if x > 0:
                    ax.set_yticklabels([])
                else:
                    if y == ny // 2:
                        ylabel = "frac. aperture flux" if args.pixel_normalization == "aperture" else "rel. flux"
                        ax.set_ylabel(ylabel, labelpad=-1)
                    else:
                        ax.set_ylabel("")
                ax.tick_params(length=2, pad=1)

        sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cax = fig.add_axes([0.952, 0.19, 0.012, 0.56])
        cbar = fig.colorbar(sm, cax=cax)
        cbar.set_label("pixel event depth")
    else:
        ax_fold = fig.add_subplot(outer[0, 1])
        right_panel_summary = plot_folded_transit_panel(
            ax_fold,
            table,
            time,
            good,
            phase_min,
            x_edges,
            args.phase_window_min,
            args.duration_min,
            args.clip_sigma,
        )

    tic = header.get("TICID", "unknown")
    sector = header.get("SECTOR", "unknown")
    tmag = header.get("TESSMAG", np.nan)
    title = args.title
    if title is None:
        tmag_label = f"T={tmag:.2f}" if isinstance(tmag, (int, float)) and np.isfinite(tmag) else "T=?"
        title = f"TIC {tic}  S{sector}  TGLC save_aper=True  {tmag_label}"
    subtitle = (
        f"P={args.period_d:.6f} d, T0={args.t0_bjd:.6f} BJD; "
        f"{good.sum():,}/{len(good):,} quality-zero cadences"
    )
    fig.text(0.5, 0.985, title, ha="center", va="top", fontsize=13, color="#b2182b", weight="bold")
    fig.text(0.5, 0.958, subtitle, ha="center", va="top", fontsize=9, color="0.2")
    fig.subplots_adjust(top=0.92, bottom=0.08, left=0.055, right=0.935)

    png_path = base.with_suffix(".png")
    pdf_path = base.with_suffix(".pdf")
    json_path = base.with_suffix(".json")
    fig.savefig(png_path, dpi=220)
    fig.savefig(pdf_path)
    plt.close(fig)

    best_flat = int(np.nanargmax(positive_depths)) if np.isfinite(positive_depths).any() else -1
    best_yx = divmod(best_flat, nx) if best_flat >= 0 else (None, None)
    summary = {
        "input_fits": str(args.fits),
        "output_png": str(png_path),
        "output_pdf": str(pdf_path),
        "period_d": args.period_d,
        "t0_bjd": args.t0_bjd,
        "t0_in_time_column": t0_time,
        "time_system_inferred": time_system,
        "duration_min": args.duration_min,
        "pixel_normalization": args.pixel_normalization,
        "phase_window_min": args.phase_window_min,
        "pixel_detrend_window_d": args.pixel_detrend_window_d,
        "clip_sigma": args.clip_sigma,
        "right_panel": args.right_panel,
        "right_panel_summary": right_panel_summary,
        "show_aperture_overlay": bool(args.show_aperture_overlay),
        "pixel_aperture_scale": pixel_aperture_scale,
        "pixel_aperture_scale_source": pixel_aperture_scale_source,
        "source_pickle": None if args.source_pickle is None else str(args.source_pickle),
        "show_model_contours": bool(args.show_model_contours),
        "gaia_cutout_geometry": gaia_geometry,
        "gaia_neighbors": neighbors,
        "cube_shape": list(cube.shape),
        "quality_zero_cadences": int(good.sum()),
        "total_cadences": int(len(good)),
        "best_depth_pixel_yx": list(best_yx),
        "best_depth": None if best_yx[0] is None else float(depths[best_yx]),
        "best_depth_sigma": None if best_yx[0] is None else float(sigmas[best_yx]),
        "normalization_modes": {f"{y},{x}": norm_modes[(y, x)] for y in range(ny) for x in range(nx)},
        "header": {
            "TICID": str(header.get("TICID", "")),
            "GAIADR3": str(header.get("GAIADR3", "")),
            "SECTOR": str(header.get("SECTOR", "")),
            "CAMERA": str(header.get("CAMERA", "")),
            "CCD": str(header.get("CCD", "")),
            "FFIVER": str(header.get("FFIVER", "")),
            "TGLCVER": str(header.get("TGLCVER", "")),
            "CONTAMRT": str(header.get("CONTAMRT", "")),
        },
    }
    json_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote {png_path}")
    print(f"Wrote {pdf_path}")
    print(f"Wrote {json_path}")


if __name__ == "__main__":
    main()
