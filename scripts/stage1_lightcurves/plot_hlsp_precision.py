#!/usr/bin/env python3
"""Plot photometric precision versus TESS magnitude.

This is the precision-panel-only counterpart to ``qc_sector_pdf.py``.  It is
useful for sharded HLSP trees where the full QC PDF's unrestricted recursive
glob is unnecessarily slow.  It can also read TGLC HDF5 files and reproduce
the TGLC raw-flux MAD normalization.
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from astropy.io import fits

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.lightcurves.precision import (  # noqa: E402
    point_to_point_mad_precision,
    quality_good_mask,
)
from twirl.plotting.style import apply_twirl_style  # noqa: E402


def discover_hlsp_paths(root: Path, sector: int, pattern: str | None) -> list[Path]:
    if pattern:
        return sorted(root.glob(pattern))
    sector_token = f"s{sector:04d}"
    sharded = sorted(root.glob(f"*/*/*/*/hlsp_*_tess_ffi_{sector_token}-*.fits"))
    if sharded:
        return sharded
    return sorted(root.rglob(f"hlsp_*_tess_ffi_{sector_token}-*.fits"))


def discover_tglc_h5_paths(root: Path, pattern: str | None) -> list[Path]:
    if pattern:
        return sorted(root.glob(pattern))
    return sorted(root.rglob("LC/*.h5"))


def discover_tglc_h5_orbit_paths(root: Path, orbits: list[int]) -> list[Path]:
    paths: list[Path] = []
    for orbit in orbits:
        orbit_root = root / f"orbit-{orbit}" / "ffi"
        paths.extend(sorted(orbit_root.glob("cam*/ccd*/LC/*.h5")))
    return paths


def tic_from_h5_path(path: Path) -> int | None:
    try:
        return int(path.stem)
    except ValueError:
        return None


def group_h5_paths_by_tic(paths: list[Path]) -> dict[int, list[Path]]:
    groups: dict[int, list[Path]] = {}
    for path in paths:
        tic = tic_from_h5_path(path)
        if tic is None:
            continue
        groups.setdefault(tic, []).append(path)
    return {tic: sorted(group) for tic, group in sorted(groups.items())}


def good_mask(quality: np.ndarray, drop_bits: int) -> np.ndarray:
    return quality_good_mask(quality, drop_bits)


def per_lc_precision_dmad(
    flux: np.ndarray,
    quality: np.ndarray,
    *,
    bin_minutes: float,
    cadence_s: float,
    drop_bits: int,
) -> float:
    return point_to_point_mad_precision(
        flux,
        quality,
        bin_minutes=bin_minutes,
        cadence_s=cadence_s,
        drop_bits=drop_bits,
    )


def tess_flux_scale_from_tmag(tmag: float) -> float:
    """TGLC/TESS flux scale in e-/s for a given TESS magnitude."""
    return 1.5e4 * 10.0 ** ((10.0 - np.asarray(tmag, dtype=float)) / 2.5)


def infer_h5_aperture_flux_scale(
    flux: np.ndarray,
    raw_magnitude: np.ndarray,
    quality: np.ndarray,
    *,
    tmag: float,
    cadence_s: float,
    drop_bits: int,
    min_points: int = 30,
) -> float:
    """Infer the expected aperture flux per cadence for a TGLC HDF5 aperture.

    TGLC HDF5 ``RawFlux`` is aperture flux in electrons per cadence, while
    ``RawMagnitude`` is computed after dividing by the aperture flux fraction
    and exposure time.  The original TGLC precision plots work in total-flux
    units, so for HDF5 raw flux we infer the aperture fraction and normalize by
    the corresponding expected aperture flux.
    """
    flux = np.asarray(flux, dtype=float)
    raw_magnitude = np.asarray(raw_magnitude, dtype=float)
    quality = np.asarray(quality, dtype=np.int64)
    mask = (
        good_mask(quality, drop_bits)
        & np.isfinite(flux)
        & (flux > 0)
        & np.isfinite(raw_magnitude)
    )
    if np.count_nonzero(mask) < min_points:
        return np.nan

    inferred_total_per_cad = tess_flux_scale_from_tmag(raw_magnitude[mask]) * float(cadence_s)
    aperture_fraction = flux[mask] / inferred_total_per_cad
    aperture_fraction = aperture_fraction[np.isfinite(aperture_fraction) & (aperture_fraction > 0)]
    if aperture_fraction.size < min_points:
        return np.nan

    expected_total_per_cad = tess_flux_scale_from_tmag(tmag) * float(cadence_s)
    return float(expected_total_per_cad * np.nanmedian(aperture_fraction))


def tglc_binned_rawflux_precision(
    flux: np.ndarray,
    raw_magnitude: np.ndarray,
    quality: np.ndarray,
    *,
    tmag: float,
    bin_minutes: float,
    cadence_s: float,
    drop_bits: int,
    min_points: int = 30,
    mad_scale: float = 1.48,
) -> float:
    """TGLC-style MAD precision for raw aperture flux.

    This mirrors the TGLC plotting helper for the HDF5 schema: optionally
    average adjacent cadences into the requested bin size, take the
    point-to-point MAD, and normalize by the expected aperture flux scale
    implied by ``TESSMAG`` and the aperture fraction inferred from
    ``RawMagnitude``.
    """
    flux = np.asarray(flux, dtype=float)
    raw_magnitude = np.asarray(raw_magnitude, dtype=float)
    quality = np.asarray(quality, dtype=np.int64)
    clean = flux[good_mask(quality, drop_bits) & np.isfinite(flux)]
    if clean.size < min_points:
        return np.nan

    n_per_bin = max(1, int(round(float(bin_minutes) * 60.0 / float(cadence_s))))
    if n_per_bin > 1:
        n_keep = clean.size // n_per_bin * n_per_bin
        if n_keep < max(2 * n_per_bin, min_points):
            return np.nan
        clean = np.mean(clean[:n_keep].reshape(-1, n_per_bin), axis=1)

    if clean.size < 2:
        return np.nan
    scale = infer_h5_aperture_flux_scale(
        flux,
        raw_magnitude,
        quality,
        tmag=tmag,
        cadence_s=cadence_s,
        drop_bits=drop_bits,
        min_points=min_points,
    )
    if not (np.isfinite(scale) and scale > 0):
        return np.nan
    mad = np.nanmedian(np.abs(np.diff(clean)))
    return float(mad_scale * mad / (np.sqrt(2.0) * scale))


def tglc_binned_rawflux_absdiff(
    flux: np.ndarray,
    raw_magnitude: np.ndarray,
    quality: np.ndarray,
    *,
    tmag: float,
    bin_minutes: float,
    cadence_s: float,
    drop_bits: int,
    min_points: int = 30,
) -> np.ndarray:
    flux = np.asarray(flux, dtype=float)
    quality = np.asarray(quality, dtype=np.int64)
    clean = flux[good_mask(quality, drop_bits) & np.isfinite(flux)]
    if clean.size < min_points:
        return np.array([], dtype=float)

    n_per_bin = max(1, int(round(float(bin_minutes) * 60.0 / float(cadence_s))))
    if n_per_bin > 1:
        n_keep = clean.size // n_per_bin * n_per_bin
        if n_keep < max(2 * n_per_bin, min_points):
            return np.array([], dtype=float)
        clean = np.mean(clean[:n_keep].reshape(-1, n_per_bin), axis=1)

    if clean.size < 2:
        return np.array([], dtype=float)
    scale = infer_h5_aperture_flux_scale(
        flux,
        raw_magnitude,
        quality,
        tmag=tmag,
        cadence_s=cadence_s,
        drop_bits=drop_bits,
        min_points=min_points,
    )
    if not (np.isfinite(scale) and scale > 0):
        return np.array([], dtype=float)
    return np.abs(np.diff(clean)) / scale


def read_precision_row(
    path: Path,
    *,
    flux_column: str,
    bin_minutes: float,
    cadence_s: float,
    drop_bits: int,
) -> tuple[float, float] | None:
    try:
        with fits.open(path, memmap=False) as hdul:
            hdr = hdul[0].header
            tab = hdul[1].data
            names = set(tab.columns.names)
            if flux_column not in names:
                return None
            quality = np.asarray(tab["QUALITY"] if "QUALITY" in names else np.zeros(len(tab), dtype=np.int32))
            flux = np.asarray(tab[flux_column], dtype=float)
            precision = per_lc_precision_dmad(
                flux,
                quality,
                bin_minutes=bin_minutes,
                cadence_s=cadence_s,
                drop_bits=drop_bits,
            )
            tmag = float(hdr.get("TESSMAG", np.nan))
    except Exception:
        return None
    if not (np.isfinite(tmag) and np.isfinite(precision) and precision > 0):
        return None
    return tmag, precision


def read_tglc_h5_precision_row(
    path: Path,
    *,
    aperture: str,
    bin_minutes: float,
    cadence_s: float,
    drop_bits: int,
) -> tuple[float, float] | None:
    try:
        import h5py

        with h5py.File(path, "r") as h5:
            tmag = float(h5.attrs.get("TessMag", np.nan))
            quality = np.asarray(h5["LightCurve/QualityFlag"][:], dtype=np.int64)
            flux = np.asarray(
                h5[f"LightCurve/AperturePhotometry/{aperture}Aperture/RawFlux"][:],
                dtype=float,
            )
            raw_magnitude = np.asarray(
                h5[f"LightCurve/AperturePhotometry/{aperture}Aperture/RawMagnitude"][:],
                dtype=float,
            )
            precision = tglc_binned_rawflux_precision(
                flux,
                raw_magnitude,
                quality,
                tmag=tmag,
                bin_minutes=bin_minutes,
                cadence_s=cadence_s,
                drop_bits=drop_bits,
            )
    except Exception:
        return None
    if not (np.isfinite(tmag) and np.isfinite(precision) and precision > 0):
        return None
    return tmag, precision


def read_tglc_h5_group_precision_row(
    paths: list[Path],
    *,
    aperture: str,
    bin_minutes: float,
    cadence_s: float,
    drop_bits: int,
) -> tuple[float, float] | None:
    try:
        import h5py

        tmags: list[float] = []
        diffs: list[np.ndarray] = []
        for path in sorted(paths):
            with h5py.File(path, "r") as h5:
                tmag = float(h5.attrs.get("TessMag", np.nan))
                quality = np.asarray(h5["LightCurve/QualityFlag"][:], dtype=np.int64)
                flux = np.asarray(
                    h5[f"LightCurve/AperturePhotometry/{aperture}Aperture/RawFlux"][:],
                    dtype=float,
                )
                raw_magnitude = np.asarray(
                    h5[f"LightCurve/AperturePhotometry/{aperture}Aperture/RawMagnitude"][:],
                    dtype=float,
                )
            if np.isfinite(tmag):
                tmags.append(tmag)
            diff = tglc_binned_rawflux_absdiff(
                flux,
                raw_magnitude,
                quality,
                tmag=tmag,
                bin_minutes=bin_minutes,
                cadence_s=cadence_s,
                drop_bits=drop_bits,
            )
            if diff.size:
                diffs.append(diff)
    except Exception:
        return None

    if not tmags or not diffs:
        return None
    tmag = float(np.nanmedian(tmags))
    precision = float(1.48 * np.nanmedian(np.concatenate(diffs)) / np.sqrt(2.0))
    if not (np.isfinite(tmag) and np.isfinite(precision) and precision > 0):
        return None
    return tmag, precision


def running_median(y: np.ndarray, window: int) -> np.ndarray:
    if y.size == 0:
        return y
    half = max(1, window // 2)
    out = np.empty_like(y, dtype=float)
    for i in range(y.size):
        lo = max(0, i - half)
        hi = min(y.size, i + half + 1)
        segment = y[lo:hi]
        finite = np.isfinite(segment)
        out[i] = np.nanmedian(segment[finite]) if finite.any() else np.nan
    return out


def finite_interp(x: np.ndarray, xp: np.ndarray, fp: np.ndarray) -> np.ndarray:
    out = np.interp(x, xp, fp)
    out[(x < np.nanmin(xp)) | (x > np.nanmax(xp))] = np.nan
    return out


def load_noisemodel(path: Path | None) -> tuple[np.ndarray, np.ndarray] | None:
    if path is None or not path.is_file():
        return None
    data = np.loadtxt(path)
    return np.asarray(data[:, 0], dtype=float), np.asarray(data[:, 1], dtype=float)


def plot_precision(
    *,
    tmag: np.ndarray,
    precision: np.ndarray,
    noisemodel: tuple[np.ndarray, np.ndarray] | None,
    sector: int,
    bin_minutes: float,
    label: str,
    output: Path,
    style: str,
    tmag_min: float,
    tmag_max: float,
    running_median_window: int,
) -> None:
    template_name = "tglc_precision" if style == "tglc" else "full_page"
    template = apply_twirl_style(template_name)
    order = np.argsort(tmag)
    tmag_s = tmag[order]
    precision_s = precision[order]
    win = max(1, int(running_median_window))
    rmed = running_median(precision_s, win)

    sigma_base = None
    if noisemodel is not None:
        sigma_base = finite_interp(tmag_s, noisemodel[0], noisemodel[1])

    fig, axes = plt.subplots(
        2,
        1,
        figsize=template["figsize"],
        sharex=True,
        gridspec_kw={"height_ratios": [3, 2], "hspace": 0.08},
    )
    ax_t, ax_b = axes
    scatter_color = "#d95f02"
    line_color = "#d95f02" if style == "tglc" else "#9b3f13"
    alpha = 0.22 if style == "tglc" else 0.16
    marker_size = max(template["dense_marker_size"], 4.0 if style == "tglc" else 3.0)
    median_lw = 2.8 if style == "tglc" else 1.8
    baseline_lw = 2.0 if style == "tglc" else 1.3
    median_path_effects = [
        pe.Stroke(linewidth=median_lw + 4.0, foreground="white"),
        pe.Stroke(linewidth=median_lw + 2.0, foreground="black"),
        pe.Normal(),
    ]

    ax_t.scatter(
        tmag_s,
        precision_s,
        s=marker_size,
        alpha=alpha,
        color=scatter_color,
        rasterized=True,
        label=label,
    )
    ax_t.plot(
        tmag_s,
        rmed,
        "-",
        color=line_color,
        lw=median_lw,
        path_effects=median_path_effects,
        zorder=8,
    )
    if noisemodel is not None:
        ax_t.plot(noisemodel[0], noisemodel[1], "k-", lw=baseline_lw, label=r"$\sigma_{\rm base}(T)$")
    ax_t.set_yscale("log")
    ylabel = "30-min precision" if style == "tglc" else "Estimated Photometric Precision"
    if style != "tglc":
        ylabel += f"\n(frac, {int(bin_minutes)}-min bin)"
    ax_t.set_ylabel(ylabel)
    ax_t.set_ylim(1e-4, 1.0)
    ax_t.set_xlim(tmag_min, tmag_max)
    ax_t.set_title(f"S{sector} {int(bin_minutes)}-min bin")
    ax_t.legend(loc="lower right" if style == "tglc" else "upper left", markerscale=4)

    if sigma_base is not None:
        ratio = precision_s / sigma_base
        finite_ratio = np.isfinite(ratio)
        if style != "tglc":
            ax_b.scatter(
                tmag_s,
                ratio,
                s=marker_size,
                alpha=alpha,
                color=scatter_color,
                rasterized=True,
            )
        rmed_r = running_median(ratio, win)
        ax_b.plot(
            tmag_s,
            rmed_r,
            "-",
            color=line_color,
            lw=median_lw,
            path_effects=median_path_effects,
            label=label,
            zorder=8,
        )
        ax_b.axhline(1.0, color="k", lw=baseline_lw, label=r"$\sigma_{\rm base}(T)$")
        ax_b.set_yscale("log")
        ratio_floor = float(np.nanpercentile(ratio[finite_ratio], 0.5)) if finite_ratio.any() else 0.5
        ymin = max(0.05, min(0.5, 0.8 * ratio_floor))
        ax_b.set_ylim(ymin, 5.0)
        if ymin < 0.15:
            tick_values = [0.1, 0.2, 0.5, 1, 2, 5]
        elif ymin < 0.35:
            tick_values = [0.2, 0.5, 1, 2, 5]
        else:
            tick_values = [0.5, 1, 2, 5]
        ax_b.set_yticks(tick_values)
        ax_b.set_yticklabels([f"{v:g}" for v in tick_values])
        ax_b.legend(loc="lower right" if style == "tglc" else "upper left", ncol=2 if style == "tglc" else 1)
    ax_b.set_xlabel("TESS magnitude")
    ax_b.set_ylabel("Precision ratio" if style == "tglc" else r"Precision / $\sigma_{\rm base}(T)$")
    for ax in axes:
        for spine in ax.spines.values():
            spine.set_linewidth(1.4 if style == "tglc" else 0.8)
            spine.set_color("0.2")

    if style == "tglc":
        fig.subplots_adjust(left=0.24, right=0.98, top=0.93, bottom=0.09, hspace=0.08)
    else:
        fig.subplots_adjust(left=0.12, right=0.98, top=0.93, bottom=0.12, hspace=0.08)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hlsp-root", type=Path, required=True)
    ap.add_argument("--sector", type=int, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--input-kind", choices=["hlsp", "tglc-h5"], default="hlsp")
    ap.add_argument("--n-precision", type=int, default=5000)
    ap.add_argument("--all", action="store_true", help="Use every discovered HLSP instead of a random sample.")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--flux-column", default="DET_FLUX")
    ap.add_argument("--h5-aperture", choices=["Small", "Primary", "Large"], default="Primary")
    ap.add_argument("--h5-orbits", default=None, help="Comma-separated orbit numbers to read under --hlsp-root.")
    ap.add_argument(
        "--stitch-h5-by-tic",
        action="store_true",
        help="For TGLC HDF5 inputs, combine orbit-level files into one sector-level row per TIC.",
    )
    ap.add_argument("--label", default=None)
    ap.add_argument("--style", choices=["twirl", "tglc"], default="twirl")
    ap.add_argument("--glob-pattern", default=None)
    ap.add_argument("--quality-drop-bits", type=lambda s: int(s, 0), default=0xFFFFFFFF)
    ap.add_argument("--bin-minutes", type=float, default=30.0)
    ap.add_argument("--cadence-s", type=float, default=200.0)
    ap.add_argument(
        "--running-median-window",
        type=int,
        default=5,
        help="Number of magnitude-sorted light curves in the running-median center line.",
    )
    ap.add_argument("--noisemodel", type=Path, default=Path("data_local/refs/noisemodel.dat"))
    ap.add_argument("--tmag-min", type=float, default=7.0)
    ap.add_argument("--tmag-max", type=float, default=20.5)
    args = ap.parse_args()

    h5_orbits = None
    if args.h5_orbits:
        h5_orbits = [int(item) for item in args.h5_orbits.split(",") if item.strip()]

    if args.input_kind == "tglc-h5":
        if h5_orbits is not None:
            paths = discover_tglc_h5_orbit_paths(args.hlsp_root, h5_orbits)
        else:
            paths = discover_tglc_h5_paths(args.hlsp_root, args.glob_pattern)
        print(f"[precision] discovered {len(paths)} TGLC HDF5 files under {args.hlsp_root}", flush=True)
    else:
        paths = discover_hlsp_paths(args.hlsp_root, args.sector, args.glob_pattern)
        print(f"[precision] discovered {len(paths)} HLSP FITS under {args.hlsp_root}", flush=True)
    if not paths:
        return 1

    h5_groups: dict[int, list[Path]] | None = None
    if args.input_kind == "tglc-h5" and args.stitch_h5_by_tic:
        h5_groups = group_h5_paths_by_tic(paths)
        print(f"[precision] grouped into {len(h5_groups)} TIC-level HDF5 rows", flush=True)
        items = list(h5_groups.values())
    else:
        items = paths

    if args.all or args.n_precision <= 0:
        sample = items
    else:
        rng = random.Random(args.seed)
        sample = rng.sample(items, min(args.n_precision, len(items)))

    rows: list[tuple[float, float]] = []
    for i, item in enumerate(sample, start=1):
        if args.input_kind == "tglc-h5" and args.stitch_h5_by_tic:
            row = read_tglc_h5_group_precision_row(
                item,
                aperture=args.h5_aperture,
                bin_minutes=args.bin_minutes,
                cadence_s=args.cadence_s,
                drop_bits=args.quality_drop_bits,
            )
        elif args.input_kind == "tglc-h5":
            row = read_tglc_h5_precision_row(
                item,
                aperture=args.h5_aperture,
                bin_minutes=args.bin_minutes,
                cadence_s=args.cadence_s,
                drop_bits=args.quality_drop_bits,
            )
        else:
            row = read_precision_row(
                item,
                flux_column=args.flux_column,
                bin_minutes=args.bin_minutes,
                cadence_s=args.cadence_s,
                drop_bits=args.quality_drop_bits,
            )
        if row is not None and args.tmag_min <= row[0] <= args.tmag_max:
            rows.append(row)
        if i % 500 == 0:
            print(f"[precision] read {i}/{len(sample)} sampled LCs", flush=True)
    if not rows:
        print("[precision] no valid precision rows", file=sys.stderr, flush=True)
        return 1

    arr = np.asarray(rows, dtype=float)
    tmag = arr[:, 0]
    precision = arr[:, 1]
    noisemodel = load_noisemodel(args.noisemodel if args.noisemodel.is_file() else None)
    if args.label is not None:
        label = args.label
    elif args.input_kind == "tglc-h5":
        label = f"TGLC RawFlux {args.h5_aperture} aperture"
    else:
        label = "TGLC Aperture" if args.style == "tglc" else f"{args.flux_column}"

    plot_precision(
        tmag=tmag,
        precision=precision,
        noisemodel=noisemodel,
        sector=args.sector,
        bin_minutes=args.bin_minutes,
        label=label,
        output=args.output,
        style=args.style,
        tmag_min=args.tmag_min,
        tmag_max=args.tmag_max,
        running_median_window=args.running_median_window,
    )
    np.savez_compressed(
        args.output.with_suffix(".npz"),
        sector=args.sector,
        input_kind=args.input_kind,
        flux_column=args.flux_column,
        h5_aperture=args.h5_aperture,
        h5_orbits=np.asarray(h5_orbits if h5_orbits is not None else [], dtype=int),
        stitch_h5_by_tic=bool(args.stitch_h5_by_tic),
        h5_normalization="inferred_aperture_fraction_from_raw_magnitude",
        bin_minutes=args.bin_minutes,
        running_median_window=args.running_median_window,
        n_total_hlsp=len(paths),
        n_total_h5_groups=(len(h5_groups) if h5_groups is not None else 0),
        n_sampled=len(sample),
        n_valid=len(rows),
        used_all=bool(args.all or args.n_precision <= 0),
        tmag_min=args.tmag_min,
        tmag_max=args.tmag_max,
        style=args.style,
        tmag_prec_sample=tmag,
        prec_30min=precision,
        sigma_base_tmag=(noisemodel[0] if noisemodel is not None else np.array([])),
        sigma_base_value=(noisemodel[1] if noisemodel is not None else np.array([])),
    )
    print(f"[precision] wrote {args.output}", flush=True)
    print(f"[precision] wrote {args.output.with_suffix('.pdf')}", flush=True)
    print(f"[precision] wrote {args.output.with_suffix('.npz')}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
