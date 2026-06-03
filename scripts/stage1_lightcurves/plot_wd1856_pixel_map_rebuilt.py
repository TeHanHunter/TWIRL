#!/usr/bin/env python
"""Rebuild the WD 1856 S56 pixel-map diagnostic from the save-aperture cube.

This plot is intentionally separate from ``plot_tglc_pixel_lightcurve_map.py``:
it uses local per-transit baselines and lays out the 5x5 folded pixel panels in
the same orientation as the left ``imshow(origin="lower")`` aperture map.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import pickle
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from twirl.plotting.style import apply_twirl_style  # noqa: E402

BJD_MINUS_BTJD = 2_457_000.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fits",
        type=Path,
        default=Path(
            "data_local/stage1_lightcurves/tglc_saveaper_wd1856_s56_spoc/lc/"
            "hlsp_tglc_tess_ffi_gaiaid-2146576589564898688-s0056-cam4-ccd1_tess_v2.1_llc.fits"
        ),
        help="TGLC save_aper=True FITS product.",
    )
    parser.add_argument(
        "--source-pickle",
        type=Path,
        default=Path(
            "data_local/stage1_lightcurves/tglc_saveaper_wd1856_s56_spoc/source/"
            "source_SPOC_TIC 267574918_sector_56.pkl"
        ),
        help="Matching TGLC source pickle for Gaia source positions.",
    )
    parser.add_argument(
        "--out-prefix",
        type=Path,
        default=Path("reports/stage1_lightcurves/pixel_maps/wd1856_s56_pixel_map_rebuilt"),
        help="Output prefix; .png, .pdf, and .json are written.",
    )
    parser.add_argument("--period-d", type=float, default=1.407983)
    parser.add_argument("--t0-bjd", type=float, default=2459825.473398)
    parser.add_argument("--duration-min", type=float, default=6.0)
    parser.add_argument("--phase-window-min", type=float, default=60.0)
    parser.add_argument("--baseline-inner-min", type=float, default=14.0)
    parser.add_argument("--baseline-outer-min", type=float, default=56.0)
    parser.add_argument("--bin-width-min", type=float, default=4.0)
    parser.add_argument("--neighbor-gmag-max", type=float, default=21.0)
    parser.add_argument("--label-neighbor-gmag-max", type=float, default=18.5)
    return parser.parse_args()


def robust_mad(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    if finite.sum() == 0:
        return float("nan")
    med = np.nanmedian(values[finite])
    return float(1.4826 * np.nanmedian(np.abs(values[finite] - med)))


def optional_column(table, name: str) -> np.ndarray | None:
    lower_to_name = {col.lower(): col for col in table.names}
    actual = lower_to_name.get(name.lower())
    return None if actual is None else np.asarray(table[actual])


def quality_mask(table, time: np.ndarray) -> np.ndarray:
    good = np.isfinite(time)
    for name in ("TESS_flags", "TGLC_flags", "QUALITY"):
        col = optional_column(table, name)
        if col is not None:
            good &= np.asarray(col, dtype=int) == 0
    return good


def infer_t0_for_time_column(time: np.ndarray, t0_bjd: float) -> tuple[float, str]:
    med = float(np.nanmedian(time[np.isfinite(time)]))
    if med > 100_000:
        return t0_bjd, "BJD"
    return t0_bjd - BJD_MINUS_BTJD, "BTJD"


def phase_minutes(time: np.ndarray, period_d: float, t0: float) -> np.ndarray:
    return (((time - t0 + 0.5 * period_d) % period_d) - 0.5 * period_d) * 1440.0


def event_numbers(time: np.ndarray, period_d: float, t0: float) -> np.ndarray:
    return np.rint((time - t0) / period_d).astype(int)


def binned_median(x: np.ndarray, y: np.ndarray, edges: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    centers = 0.5 * (edges[:-1] + edges[1:])
    med = np.full(len(centers), np.nan)
    err = np.full(len(centers), np.nan)
    use = np.isfinite(x) & np.isfinite(y)
    which = np.searchsorted(edges, x[use], side="right") - 1
    vals = y[use]
    for idx in range(len(centers)):
        sub = vals[which == idx]
        if sub.size:
            med[idx] = np.nanmedian(sub)
            scatter = robust_mad(sub)
            if np.isfinite(scatter):
                err[idx] = scatter / np.sqrt(sub.size)
    return centers, med, err


def local_event_residuals(
    flux: np.ndarray,
    time: np.ndarray,
    good: np.ndarray,
    event_no: np.ndarray,
    phase_min: np.ndarray,
    phase_window_min: float,
    baseline_inner_min: float,
    baseline_outer_min: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return event-folded phase, residual flux, and per-cadence local baseline."""
    flux = np.asarray(flux, dtype=float)
    residual = np.full_like(flux, np.nan, dtype=float)
    baseline = np.full_like(flux, np.nan, dtype=float)
    in_window_all = good & np.isfinite(flux) & np.isfinite(phase_min) & (np.abs(phase_min) <= phase_window_min)
    for evt in np.unique(event_no[in_window_all]):
        evt_mask = in_window_all & (event_no == evt)
        base_mask = evt_mask & (np.abs(phase_min) >= baseline_inner_min) & (np.abs(phase_min) <= baseline_outer_min)
        if base_mask.sum() < 4:
            continue
        base = float(np.nanmedian(flux[base_mask]))
        residual[evt_mask] = flux[evt_mask] - base
        baseline[evt_mask] = base
    keep = np.isfinite(residual)
    return phase_min[keep], residual[keep], baseline[keep]


def measure_local_depth(
    phase: np.ndarray,
    residual: np.ndarray,
    duration_min: float,
    baseline_inner_min: float,
    baseline_outer_min: float,
) -> tuple[float, float, int, int]:
    in_event = np.isfinite(residual) & (np.abs(phase) <= 0.5 * duration_min)
    oot = (
        np.isfinite(residual)
        & (np.abs(phase) >= baseline_inner_min)
        & (np.abs(phase) <= baseline_outer_min)
    )
    if in_event.sum() < 3 or oot.sum() < 12:
        return float("nan"), float("nan"), int(in_event.sum()), int(oot.sum())
    depth = float(np.nanmedian(residual[oot]) - np.nanmedian(residual[in_event]))
    scatter = robust_mad(residual[oot])
    snr = depth / scatter if np.isfinite(scatter) and scatter > 0 else float("nan")
    return depth, snr, int(in_event.sum()), int(oot.sum())


def designation_to_source_id(value) -> int | None:
    if value is None:
        return None
    text = value.decode() if isinstance(value, bytes) else str(value)
    try:
        return int(text.split()[-1])
    except ValueError:
        return None


def load_gaia_neighbors(source_pickle: Path | None, header, nx: int, ny: int, gmag_max: float) -> tuple[list[dict], dict]:
    if source_pickle is None or not source_pickle.exists():
        return [], {}
    with source_pickle.open("rb") as fh:
        with contextlib.redirect_stdout(io.StringIO()):
            source = pickle.load(fh)
    gaia = source.gaia
    sector = int(header.get("SECTOR", getattr(source, "sector", 0)))
    x_col = f"sector_{sector}_x"
    y_col = f"sector_{sector}_y"
    if x_col not in gaia.colnames or y_col not in gaia.colnames:
        return [], {"source_pickle": str(source_pickle), "error": f"missing {x_col}/{y_col}"}

    target_id = designation_to_source_id(f"Gaia DR3 {header.get('GAIADR3')}")
    ids = np.array([designation_to_source_id(value) for value in gaia["DESIGNATION"]], dtype=object)
    target_idx = np.flatnonzero(ids == target_id)
    if target_idx.size == 0:
        return [], {"source_pickle": str(source_pickle), "target_gaia_dr3": target_id, "error": "target missing"}
    target_idx = int(target_idx[0])

    star_x = float(header.get("STAR_X", nx // 2))
    star_y = float(header.get("STAR_Y", ny // 2))
    target_sector_x = float(gaia[x_col][target_idx])
    target_sector_y = float(gaia[y_col][target_idx])
    origin_x = target_sector_x - star_x
    origin_y = target_sector_y - star_y

    x_cut = np.asarray(gaia[x_col], dtype=float) - origin_x
    y_cut = np.asarray(gaia[y_col], dtype=float) - origin_y
    gmag = np.asarray(gaia["phot_g_mean_mag"], dtype=float)
    tmag = np.asarray(gaia["tess_mag"], dtype=float) if "tess_mag" in gaia.colnames else np.full(len(gaia), np.nan)
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
        sid = ids[idx]
        neighbors.append(
            {
                "source_id": None if sid is None else int(sid),
                "x": float(x_cut[idx]),
                "y": float(y_cut[idx]),
                "phot_g_mean_mag": float(gmag[idx]),
                "tess_mag": None if not np.isfinite(tmag[idx]) else float(tmag[idx]),
                "is_target": bool(sid == target_id),
            }
        )
    neighbors.sort(key=lambda item: (not item["is_target"], item["phot_g_mean_mag"]))
    geometry = {
        "source_pickle": str(source_pickle),
        "target_gaia_dr3": None if target_id is None else int(target_id),
        "star_x_from_header": star_x,
        "star_y_from_header": star_y,
        "target_sector_x": target_sector_x,
        "target_sector_y": target_sector_y,
        "crop_origin_x": origin_x,
        "crop_origin_y": origin_y,
    }
    return neighbors, geometry


def draw_gaia_neighbors(ax, neighbors: list[dict], label_gmag_max: float) -> None:
    background = [item for item in neighbors if not item["is_target"]]
    target = [item for item in neighbors if item["is_target"]]
    if background:
        x = np.array([item["x"] for item in background])
        y = np.array([item["y"] for item in background])
        g = np.array([item["phot_g_mean_mag"] for item in background])
        sizes = np.clip(220.0 * 10 ** (-0.4 * (g - 12.0)), 24.0, 145.0)
        ax.scatter(x, y, marker="o", s=sizes, facecolors="none", edgecolors="#2166ac", lw=1.2, label="Gaia DR3 background", zorder=6)
        for item in background:
            if item["phot_g_mean_mag"] <= label_gmag_max:
                ax.text(item["x"] + 0.06, item["y"] + 0.06, f"G={item['phot_g_mean_mag']:.1f}", color="#2166ac", fontsize=7, zorder=7)
    if target:
        item = target[0]
        ax.scatter([item["x"]], [item["y"]], marker="x", s=115, color="#b2182b", lw=2.1, label="WD 1856", zorder=8)
        ax.text(item["x"] + 0.11, item["y"] + 0.10, f"WD 1856\nG={item['phot_g_mean_mag']:.1f}", color="#b2182b", fontsize=8, ha="left", va="bottom", zorder=9)


def main() -> None:
    args = parse_args()
    out_prefix = args.out_prefix.with_suffix("") if args.out_prefix.suffix else args.out_prefix
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    with fits.open(args.fits, memmap=False) as hdul:
        cube = np.asarray(hdul[0].data, dtype=float)
        table = hdul[1].data
        header = hdul[0].header.copy()
    if cube.ndim != 3:
        raise RuntimeError(f"Expected (cadence, y, x) cube, got {cube.shape}")

    time = np.asarray(table["time"], dtype=float)
    good = quality_mask(table, time)
    t0_time, time_system = infer_t0_for_time_column(time, args.t0_bjd)
    pmin = phase_minutes(time, args.period_d, t0_time)
    evt_no = event_numbers(time, args.period_d, t0_time)
    ny, nx = cube.shape[1:]

    bins = np.arange(-args.phase_window_min, args.phase_window_min + args.bin_width_min, args.bin_width_min)
    if bins[-1] < args.phase_window_min:
        bins = np.append(bins, args.phase_window_min)

    pixel_profiles: dict[tuple[int, int], dict] = {}
    depth = np.full((ny, nx), np.nan)
    snr = np.full((ny, nx), np.nan)
    event_counts = {}
    all_binned = []

    for y in range(ny):
        for x in range(nx):
            phase, residual, baseline = local_event_residuals(
                cube[:, y, x],
                time,
                good,
                evt_no,
                pmin,
                args.phase_window_min,
                args.baseline_inner_min,
                args.baseline_outer_min,
            )
            centers, med, err = binned_median(phase, residual, bins)
            dep, dep_snr, n_in, n_oot = measure_local_depth(
                phase,
                residual,
                args.duration_min,
                args.baseline_inner_min,
                args.baseline_outer_min,
            )
            depth[y, x] = dep
            snr[y, x] = dep_snr
            event_counts[(y, x)] = {"n_in_event": n_in, "n_oot": n_oot, "n_folded": int(residual.size)}
            pixel_profiles[(y, x)] = {
                "phase": phase,
                "residual": residual,
                "baseline": baseline,
                "centers": centers,
                "med": med,
                "err": err,
            }
            finite_med = med[np.isfinite(med)]
            if finite_med.size:
                all_binned.append(finite_med)

    all_binned_arr = np.concatenate(all_binned) if all_binned else np.array([])
    if all_binned_arr.size:
        ylo, yhi = np.nanpercentile(all_binned_arr, [1, 99])
        ypad = max(0.7, 0.25 * (yhi - ylo))
        ylim = (min(-0.8, ylo - ypad), max(0.8, yhi + ypad))
    else:
        ylim = (-5.0, 4.0)
    # Keep the diagnostic visually focused on minute-scale transit depths.
    ylim = (max(ylim[0], -8.0), min(ylim[1], 6.0))

    positive_depth = np.where(np.isfinite(depth), np.maximum(depth, 0.0), np.nan)
    vmax_depth = float(np.nanmax(positive_depth)) if np.isfinite(positive_depth).any() else 1.0
    vmax_depth = max(vmax_depth, 1.0)
    depth_norm = mpl.colors.Normalize(vmin=0.0, vmax=vmax_depth)
    cmap = mpl.cm.viridis

    apply_twirl_style("full_page")
    fig = plt.figure(figsize=(15.0, 7.0), constrained_layout=False)
    outer = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[0.78, 1.55], wspace=0.10)

    ax_map = fig.add_subplot(outer[0, 0])
    image = np.nanmedian(cube[good], axis=0)
    finite = image[np.isfinite(image)]
    vmin, vmax = np.nanpercentile(finite, [5, 98]) if finite.size else (None, None)
    im = ax_map.imshow(image, origin="lower", cmap="Greys", vmin=vmin, vmax=vmax)
    ax_map.set_xlim(-0.5, nx - 0.5)
    ax_map.set_ylim(-0.5, ny - 0.5)
    ax_map.set_xlabel("cube column (x)")
    ax_map.set_ylabel("cube row (y)")
    ax_map.set_xticks(np.arange(nx))
    ax_map.set_yticks(np.arange(ny))
    ax_map.set_xticks(np.arange(-0.5, nx, 1.0), minor=True)
    ax_map.set_yticks(np.arange(-0.5, ny, 1.0), minor=True)
    ax_map.grid(False)
    ax_map.grid(which="minor", color="0.73", lw=0.7)
    ax_map.tick_params(which="minor", length=0)

    neighbors, gaia_geometry = load_gaia_neighbors(args.source_pickle, header, nx, ny, args.neighbor_gmag_max)
    if neighbors:
        draw_gaia_neighbors(ax_map, neighbors, args.label_neighbor_gmag_max)
        ax_map.legend(loc="lower left", fontsize=7, frameon=True, borderpad=0.35, handletextpad=0.4)
    else:
        star_x = float(header.get("STAR_X", nx // 2))
        star_y = float(header.get("STAR_Y", ny // 2))
        ax_map.scatter([star_x], [star_y], marker="x", s=115, color="#b2182b", lw=2.1, label="target")

    cbar_image = fig.colorbar(im, ax=ax_map, fraction=0.050, pad=0.025)
    cbar_image.ax.tick_params(labelsize=7)
    ax_map.text(
        0.02,
        0.98,
        "gray: median e-/s",
        transform=ax_map.transAxes,
        ha="left",
        va="top",
        fontsize=7,
        color="0.25",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 1.5},
    )
    ax_map.set_title("Gaia positions on the saved-aperture pixels", fontsize=10)

    grid = gridspec.GridSpecFromSubplotSpec(ny, nx, subplot_spec=outer[0, 1], wspace=0.075, hspace=0.075)
    mini_axes = {}
    for grid_row, y in enumerate(range(ny - 1, -1, -1)):
        for x in range(nx):
            ax = fig.add_subplot(grid[grid_row, x])
            mini_axes[(y, x)] = ax
            profile = pixel_profiles[(y, x)]
            phase = profile["phase"]
            residual = profile["residual"]
            med = profile["med"]
            centers = profile["centers"]
            err = profile["err"]

            ax.axvspan(-0.5 * args.duration_min, 0.5 * args.duration_min, color="0.90", lw=0, zorder=0)
            ax.scatter(phase, residual, s=2.5, color="0.35", alpha=0.075, linewidths=0, rasterized=True, zorder=1)
            use = np.isfinite(med)
            ax.plot(centers[use], med[use], color="black", lw=1.15, zorder=3)
            err_use = use & np.isfinite(err)
            ax.fill_between(
                centers[err_use],
                med[err_use] - err[err_use],
                med[err_use] + err[err_use],
                color="#b2182b",
                alpha=0.18,
                lw=0,
                zorder=2,
            )
            ax.axhline(0.0, color="0.45", ls=":", lw=0.75)
            ax.axvline(0.0, color="0.30", ls=":", lw=0.70)
            ax.set_xlim(-args.phase_window_min, args.phase_window_min)
            ax.set_ylim(*ylim)
            dep = depth[y, x]
            edge_color = cmap(depth_norm(max(dep, 0.0))) if np.isfinite(dep) else "0.55"
            for spine in ax.spines.values():
                spine.set_color(edge_color)
                spine.set_linewidth(2.0)
            ax.text(
                0.04,
                0.92,
                f"x={x}, y={y}\nD={dep:.2f}" if np.isfinite(dep) else f"x={x}, y={y}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=6.2,
                color="0.15",
            )
            if grid_row < ny - 1:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel("min from T0", labelpad=0)
            if x > 0:
                ax.set_yticklabels([])
            elif y == ny // 2:
                ax.set_ylabel("pixel residual (e-/s)", labelpad=1)
            else:
                ax.set_ylabel("")
            ax.tick_params(length=2, pad=1)

    sm = mpl.cm.ScalarMappable(norm=depth_norm, cmap=cmap)
    sm.set_array([])
    cax = fig.add_axes([0.948, 0.18, 0.012, 0.58])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label("local dip depth (e-/s)")

    tic = header.get("TICID", "267574918")
    sector = header.get("SECTOR", "56")
    fig.text(
        0.50,
        0.985,
        f"TIC {tic} S{sector}: rebuilt pixel-map transit diagnostic",
        ha="center",
        va="top",
        fontsize=14,
        color="#b2182b",
        weight="bold",
    )
    fig.text(
        0.50,
        0.955,
        (
            f"P={args.period_d:.6f} d, T0={args.t0_bjd:.6f} BJD; "
            f"per-event baseline from {args.baseline_inner_min:.0f}-{args.baseline_outer_min:.0f} min; "
            f"right grid is y=4 top to y=0 bottom"
        ),
        ha="center",
        va="top",
        fontsize=9,
        color="0.20",
    )
    fig.subplots_adjust(top=0.90, bottom=0.08, left=0.055, right=0.928)

    png = out_prefix.with_suffix(".png")
    pdf = out_prefix.with_suffix(".pdf")
    js = out_prefix.with_suffix(".json")
    fig.savefig(png, dpi=220)
    fig.savefig(pdf)
    plt.close(fig)

    best_flat = int(np.nanargmax(positive_depth)) if np.isfinite(positive_depth).any() else -1
    best_y, best_x = divmod(best_flat, nx) if best_flat >= 0 else (None, None)
    summary = {
        "input_fits": str(args.fits),
        "source_pickle": str(args.source_pickle) if args.source_pickle else None,
        "output_png": str(png),
        "output_pdf": str(pdf),
        "period_d": args.period_d,
        "t0_bjd": args.t0_bjd,
        "t0_in_time_column": t0_time,
        "time_system_inferred": time_system,
        "duration_min": args.duration_min,
        "phase_window_min": args.phase_window_min,
        "baseline_inner_min": args.baseline_inner_min,
        "baseline_outer_min": args.baseline_outer_min,
        "bin_width_min": args.bin_width_min,
        "right_panel_orientation": "Grid rows are plotted y=4 top to y=0 bottom, matching imshow(origin='lower') on the left.",
        "cube_shape": list(cube.shape),
        "quality_zero_cadences": int(good.sum()),
        "total_cadences": int(len(good)),
        "depth_e_per_s_yx": depth.tolist(),
        "snr_yx": snr.tolist(),
        "event_counts": {f"{y},{x}": event_counts[(y, x)] for y in range(ny) for x in range(nx)},
        "best_depth_pixel_yx": [best_y, best_x] if best_y is not None else None,
        "best_depth": None if best_y is None else float(depth[best_y, best_x]),
        "best_depth_snr": None if best_y is None else float(snr[best_y, best_x]),
        "gaia_cutout_geometry": gaia_geometry,
        "gaia_neighbors": neighbors,
        "header": {
            "TICID": str(header.get("TICID", "")),
            "GAIADR3": str(header.get("GAIADR3", "")),
            "SECTOR": str(header.get("SECTOR", "")),
            "CAMERA": str(header.get("CAMERA", "")),
            "CCD": str(header.get("CCD", "")),
            "FFIVER": str(header.get("FFIVER", "")),
            "TGLCVER": str(header.get("TGLCVER", "")),
        },
    }
    js.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote {png}")
    print(f"Wrote {pdf}")
    print(f"Wrote {js}")


if __name__ == "__main__":
    main()
