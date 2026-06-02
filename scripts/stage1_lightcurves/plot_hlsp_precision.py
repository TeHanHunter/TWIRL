#!/usr/bin/env python3
"""Plot HLSP photometric precision versus TESS magnitude.

This is the precision-panel-only counterpart to ``qc_sector_pdf.py``.  It is
useful for sharded HLSP trees where the full QC PDF's unrestricted recursive
glob is unnecessarily slow.
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
            good = good_mask(quality, drop_bits) & np.isfinite(flux)
            scale = np.nanmedian(flux[good])
            if not (np.isfinite(scale) and scale > 0):
                return None
            precision = per_lc_precision_dmad(
                flux / scale,
                quality,
                bin_minutes=bin_minutes,
                cadence_s=cadence_s,
                drop_bits=drop_bits,
            )
    except Exception:
        return None
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
        out[i] = np.nanmedian(y[lo:hi])
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
) -> None:
    template_name = "tglc_precision" if style == "tglc" else "full_page"
    template = apply_twirl_style(template_name)
    order = np.argsort(tmag)
    tmag_s = tmag[order]
    precision_s = precision[order]
    win = max(50, int(len(precision_s) // 60))
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
    alpha = 0.06 if style == "tglc" else 0.055
    median_lw = 2.2 if style == "tglc" else 1.8
    baseline_lw = 2.0 if style == "tglc" else 1.3
    stroke_width = 0.0 if style == "tglc" else 3.0

    ax_t.scatter(
        tmag_s,
        precision_s,
        s=template["dense_marker_size"],
        alpha=alpha,
        color=scatter_color,
        rasterized=True,
        label=label,
    )
    path_effects = []
    if stroke_width > 0:
        path_effects = [pe.Stroke(linewidth=stroke_width, foreground="k"), pe.Normal()]
    ax_t.plot(
        tmag_s,
        rmed,
        "-",
        color=line_color,
        lw=median_lw,
        path_effects=path_effects,
    )
    if noisemodel is not None:
        ax_t.plot(noisemodel[0], noisemodel[1], "k-", lw=baseline_lw, label=r"$\sigma_{\rm base}(T)$")
    ax_t.set_yscale("log")
    ylabel = "Estimated Photometric Precision"
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
                s=template["dense_marker_size"],
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
            path_effects=path_effects,
            label=label,
        )
        ax_b.axhline(1.0, color="k", lw=baseline_lw, label=r"$\sigma_{\rm base}(T)$")
        if style == "tglc":
            ax_b.set_ylim(0.5, 2.5)
            ax_b.set_yticks([0.5, 1, 1.5, 2])
            ax_b.set_yticklabels(["0.5", "1", "1.5", "2"])
        else:
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
    ax_b.set_ylabel("Precision Ratio" if style == "tglc" else r"Precision / $\sigma_{\rm base}(T)$")
    for ax in axes:
        for spine in ax.spines.values():
            spine.set_linewidth(1.4 if style == "tglc" else 0.8)
            spine.set_color("0.2")

    if style == "tglc":
        fig.subplots_adjust(left=0.16, right=0.98, top=0.93, bottom=0.09, hspace=0.08)
    else:
        fig.subplots_adjust(left=0.12, right=0.98, top=0.93, bottom=0.12, hspace=0.08)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    fig.savefig(output.with_suffix(".pdf"))
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
    ap.add_argument("--label", default=None)
    ap.add_argument("--style", choices=["twirl", "tglc"], default="twirl")
    ap.add_argument("--glob-pattern", default=None)
    ap.add_argument("--quality-drop-bits", type=lambda s: int(s, 0), default=0xFFFFFFFF)
    ap.add_argument("--bin-minutes", type=float, default=30.0)
    ap.add_argument("--cadence-s", type=float, default=200.0)
    ap.add_argument("--noisemodel", type=Path, default=Path("data_local/refs/noisemodel.dat"))
    ap.add_argument("--tmag-min", type=float, default=7.0)
    ap.add_argument("--tmag-max", type=float, default=20.5)
    args = ap.parse_args()

    if args.input_kind == "tglc-h5":
        paths = discover_tglc_h5_paths(args.hlsp_root, args.glob_pattern)
        print(f"[precision] discovered {len(paths)} TGLC HDF5 files under {args.hlsp_root}", flush=True)
    else:
        paths = discover_hlsp_paths(args.hlsp_root, args.sector, args.glob_pattern)
        print(f"[precision] discovered {len(paths)} HLSP FITS under {args.hlsp_root}", flush=True)
    if not paths:
        return 1
    if args.all or args.n_precision <= 0:
        sample = paths
    else:
        rng = random.Random(args.seed)
        sample = rng.sample(paths, min(args.n_precision, len(paths)))

    rows: list[tuple[float, float]] = []
    for i, path in enumerate(sample, start=1):
        if args.input_kind == "tglc-h5":
            row = read_tglc_h5_precision_row(
                path,
                aperture=args.h5_aperture,
                bin_minutes=args.bin_minutes,
                cadence_s=args.cadence_s,
                drop_bits=args.quality_drop_bits,
            )
        else:
            row = read_precision_row(
                path,
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
    label = args.label or ("TGLC Aperture" if args.style == "tglc" else f"{args.flux_column}")

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
    )
    np.savez_compressed(
        args.output.with_suffix(".npz"),
        sector=args.sector,
        input_kind=args.input_kind,
        flux_column=args.flux_column,
        h5_aperture=args.h5_aperture,
        bin_minutes=args.bin_minutes,
        n_total_hlsp=len(paths),
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
