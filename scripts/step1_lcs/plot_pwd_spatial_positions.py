#!/usr/bin/env python

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from matplotlib.lines import Line2D

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.plotting.style import apply_twirl_style, get_ordered_palette


DEFAULT_CATALOG = Path("data_local/catalogs/GaiaEDR3_WD_main.fits")
DEFAULT_OUTDIR = Path("reports/step1_lcs")

PWD_BIN_LABELS = [
    "0.00-0.25",
    "0.25-0.50",
    "0.50-0.75",
    ">0.75",
]

LEGEND_LABELS = [
    r"$0.00 \leq P_{\mathrm{wd}} < 0.25$",
    r"$0.25 \leq P_{\mathrm{wd}} < 0.50$",
    r"$0.50 \leq P_{\mathrm{wd}} \leq 0.75$",
    r"$P_{\mathrm{wd}} > 0.75$",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot the spatial distribution of the Gentile Fusillo Gaia EDR3 WD "
            "catalog in Pwd bins, following the proposal Figure 1 polar style."
        )
    )
    parser.add_argument(
        "--catalog",
        type=Path,
        default=DEFAULT_CATALOG,
        help="Path to the local Gaia EDR3 WD FITS catalog.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help="Directory for generated figure and summary files.",
    )
    parser.add_argument(
        "--max-points-per-bin",
        type=int,
        default=80_000,
        help="Maximum number of plotted points per Pwd bin.",
    )
    parser.add_argument(
        "--distance-clip-pct",
        type=float,
        default=99.5,
        help="Upper percentile used to clip radial distances for plotting.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed used for reproducible sampling.",
    )
    parser.add_argument(
        "--template",
        choices=("column", "full_page"),
        default="full_page",
        help="Named TWIRL plotting template from doc/plotting_style.md.",
    )
    parser.add_argument(
        "--plot-all-points",
        action="store_true",
        help="Plot all valid rows instead of a stratified per-bin subsample.",
    )
    parser.add_argument(
        "--dense-marker-scale",
        type=float,
        default=1.0,
        help="Scale factor applied to the default dense scatter marker size.",
    )
    parser.add_argument(
        "--output-stem",
        type=str,
        default="wd_pwd_spatial_positions",
        help="Base filename for output figure products.",
    )
    return parser.parse_args()


def load_catalog(path: Path) -> dict[str, np.ndarray]:
    table = Table.read(path, hdu="Joined", unit_parse_strict="silent")
    table.convert_bytestring_to_unicode()

    data = {
        "source_id": np.asarray(table["source_id"]),
        "ra": np.asarray(table["ra"], dtype=float),
        "dec": np.asarray(table["dec"], dtype=float),
        "parallax": np.asarray(table["parallax"], dtype=float),
        "pwd": np.asarray(table["Pwd"], dtype=float),
    }
    return data


def build_masks(data: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    finite_mask = (
        np.isfinite(data["ra"])
        & np.isfinite(data["dec"])
        & np.isfinite(data["parallax"])
        & np.isfinite(data["pwd"])
    )
    positive_parallax = data["parallax"] > 0
    base_mask = finite_mask & positive_parallax
    return finite_mask, base_mask


def assign_pwd_bins(pwd: np.ndarray) -> np.ndarray:
    bin_ids = np.full(pwd.shape, -1, dtype=int)
    finite = np.isfinite(pwd)

    bin_ids[finite & (pwd >= 0.0) & (pwd < 0.25)] = 0
    bin_ids[finite & (pwd >= 0.25) & (pwd < 0.5)] = 1
    bin_ids[finite & (pwd >= 0.5) & (pwd <= 0.75)] = 2
    bin_ids[finite & (pwd > 0.75) & (pwd <= 1.0)] = 3

    return bin_ids


def stratified_sample_indices(
    bin_ids: np.ndarray,
    base_mask: np.ndarray,
    max_points_per_bin: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    sampled_indices = []

    for bin_index in range(len(PWD_BIN_LABELS)):
        candidate_idx = np.flatnonzero(base_mask & (bin_ids == bin_index))
        if candidate_idx.size == 0:
            continue
        if candidate_idx.size <= max_points_per_bin:
            sampled_indices.append(candidate_idx)
            continue
        sampled_indices.append(
            np.sort(rng.choice(candidate_idx, size=max_points_per_bin, replace=False))
        )

    if not sampled_indices:
        return np.array([], dtype=int)

    return np.concatenate(sampled_indices)


def all_valid_indices(base_mask: np.ndarray) -> np.ndarray:
    return np.flatnonzero(base_mask)


def write_summary_csv(
    outpath: Path,
    bin_ids: np.ndarray,
    base_mask: np.ndarray,
    sampled_indices: np.ndarray,
    distance_pc: np.ndarray,
) -> None:
    sampled_mask = np.zeros(base_mask.size, dtype=bool)
    sampled_mask[sampled_indices] = True

    with outpath.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "pwd_bin",
                "bin_rule",
                "catalog_count",
                "positive_parallax_count",
                "plotted_count",
                "median_distance_pc",
                "p95_distance_pc",
            ]
        )
        bin_rules = [
            "0.00 <= Pwd < 0.25",
            "0.25 <= Pwd < 0.50",
            "0.50 <= Pwd <= 0.75",
            "0.75 < Pwd <= 1.00",
        ]
        for bin_index, label in enumerate(PWD_BIN_LABELS):
            bin_mask = bin_ids == bin_index
            valid_mask = base_mask & bin_mask
            valid_dist = distance_pc[valid_mask]
            if valid_dist.size:
                median_distance = float(np.nanmedian(valid_dist))
                p95_distance = float(np.nanpercentile(valid_dist, 95))
            else:
                median_distance = np.nan
                p95_distance = np.nan

            writer.writerow(
                [
                    label,
                    bin_rules[bin_index],
                    int(np.count_nonzero(bin_mask)),
                    int(np.count_nonzero(valid_mask)),
                    int(np.count_nonzero(sampled_mask & bin_mask)),
                    median_distance,
                    p95_distance,
                ]
            )


def galactic_plane_segments() -> list[tuple[np.ndarray, np.ndarray]]:
    galactic_longitude = np.linspace(0.0, 360.0, 1441)
    coords = SkyCoord(
        l=galactic_longitude * u.deg,
        b=np.zeros_like(galactic_longitude) * u.deg,
        frame="galactic",
    ).icrs
    ra_deg = coords.ra.wrap_at(360 * u.deg).deg
    dec_deg = coords.dec.deg

    split_idx = np.where(np.abs(np.diff(ra_deg)) > 180.0)[0] + 1
    ra_segments = np.split(ra_deg, split_idx)
    dec_segments = np.split(dec_deg, split_idx)
    return [
        (ra_segment, dec_segment)
        for ra_segment, dec_segment in zip(ra_segments, dec_segments)
        if ra_segment.size > 1
    ]


def plot_pwd_bins(
    ra_deg: np.ndarray,
    dec_deg: np.ndarray,
    distance_pc: np.ndarray,
    pwd_bin_ids: np.ndarray,
    outdir: Path,
    distance_clip_pct: float,
    template_name: str,
    dense_marker_scale: float,
    output_stem: str,
) -> tuple[Path, Path]:
    template = apply_twirl_style(template_name)
    palette = get_ordered_palette(len(PWD_BIN_LABELS))
    marker_size = template["dense_marker_size"] * dense_marker_scale

    clipped_rmax = float(np.nanpercentile(distance_pc, distance_clip_pct))
    clipped_distance = np.clip(distance_pc, 0, clipped_rmax)
    theta = np.deg2rad(dec_deg)

    fig = plt.figure(figsize=template["figsize"])
    grid = fig.add_gridspec(
        1,
        2,
        width_ratios=[1.0, 1.45],
        wspace=max(template["panel_wspace"], 0.12),
    )
    ax = fig.add_subplot(grid[0, 0], projection="polar")
    ax_sky = fig.add_subplot(grid[0, 1])

    ax.set_thetalim(-np.pi / 2, np.pi / 2)
    ax.set_theta_zero_location("W")
    ax.set_theta_direction(-1)

    for bin_index, label in enumerate(PWD_BIN_LABELS):
        mask = pwd_bin_ids == bin_index
        if not np.any(mask):
            continue
        ax.scatter(
            theta[mask],
            clipped_distance[mask],
            s=marker_size,
            alpha=0.09,
            linewidths=0,
            rasterized=True,
            color=palette[bin_index],
            label=f"Pwd {label}",
        )
        ax_sky.scatter(
            ra_deg[mask],
            dec_deg[mask],
            s=max(0.45, marker_size - 0.2),
            alpha=0.18,
            linewidths=0,
            rasterized=True,
            color=palette[bin_index],
        )

    wd_dec = 53 + 30 / 60 + 33.3 / 3600
    wd_ra = 18 * 15 + 57 / 60 * 15 + 39.0 / 3600 * 15
    wd_dist = 1000.0 / 40.3983
    wd_theta = np.deg2rad(wd_dec)
    ax.scatter(
        [wd_theta],
        [wd_dist],
        s=34 if template_name == "column" else 42,
        marker="o",
        facecolor="white",
        edgecolor="k",
        linewidth=0.9,
        alpha=0.82,
        zorder=5,
    )
    ax.text(
        wd_theta ,
        wd_dist + 200,
        "WD 1856",
        ha="right",
        va="center",
        fontsize=template["annotation_size"],
    )
    ax_sky.scatter(
        [wd_ra],
        [wd_dec],
        s=34 if template_name == "column" else 42,
        marker="o",
        facecolor="white",
        edgecolor="k",
        linewidth=0.9,
        alpha=0.82,
        zorder=5,
    )
    ax_sky.text(
        wd_ra - 10,
        wd_dec,
        "WD 1856",
        ha="left",
        va="center",
        fontsize=template["annotation_size"],
    )

    dec_ticks = [90, 60, 30, 0, -30, -60, -90]
    ax.set_thetagrids(dec_ticks, labels=[f"{deg}°" for deg in dec_ticks])
    ax.set_rlim(0, clipped_rmax)
    max_kpc_tick = max(1, int(clipped_rmax // 1000))
    radial_ticks = 1000.0 * np.arange(1, max_kpc_tick + 1)
    ax.set_rticks(radial_ticks)
    ax.set_yticklabels([str(idx) for idx in range(1, max_kpc_tick + 1)])
    ax.set_rlabel_position(135)
    ax.tick_params(axis="x", labelsize=template["tick_size"], pad=-3)
    ax.tick_params(axis="y", labelsize=template["tick_size"], pad=2)
    ax.grid(color="0.72", linewidth=template["grid_linewidth"])
    ax.spines["polar"].set_color("0.4")
    ax.text(
        0.8,
        0.82,
        "Distance (kpc)",
        ha="center",
        va="center",
        rotation=90,
        fontsize=template["label_size"],
        transform=ax.transAxes,
    )
    ax.text(
        0.1,
        0.52,
        "Dec (degree)",
        ha="center",
        va="center",
        rotation=90,
        fontsize=template["label_size"],
        transform=ax.transAxes,
    )

    ax_sky.set_xlim(360, 0)
    ax_sky.set_ylim(-90, 90)
    ax_sky.set_xticks([360, 300, 240, 180, 120, 60, 0])
    ax_sky.set_yticks([-80, -40, 0, 40, 80])
    ax_sky.set_xlabel("RA (degree)", fontsize=template["label_size"])
    ax_sky.set_ylabel("Dec (degree)", fontsize=template["label_size"], labelpad=6)
    ax_sky.yaxis.set_label_position("left")
    ax_sky.yaxis.set_label_coords(-0.095, 0.5)
    ax_sky.grid(color="0.85", linewidth=template["grid_linewidth"])
    ax_sky.set_facecolor("1")
    ax_sky.tick_params(axis="x", labelsize=template["tick_size"], pad=2)
    ax_sky.tick_params(axis="y", labelsize=template["tick_size"], pad=2)

    for ra_segment, dec_segment in galactic_plane_segments():
        ax_sky.plot(
            ra_segment,
            dec_segment,
            color="black",
            linewidth=0.7,
            linestyle="--",
            alpha=0.95,
            zorder=4,
        )

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="None",
            markerfacecolor=palette[idx],
            markeredgecolor="none",
            markersize=4.6,
            alpha=0.78,
            label=LEGEND_LABELS[idx],
        )
        for idx, _ in enumerate(PWD_BIN_LABELS)
    ]
    legend = fig.legend(
        legend_handles,
        LEGEND_LABELS,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.975),
        ncol=4,
        frameon=True,
        title="White dwarf probability",
        fontsize=template["legend_size"],
        title_fontsize=template["legend_title_size"],
        handletextpad=0.35,
        columnspacing=0.9,
        borderpad=0.35,
    )
    legend.get_frame().set_edgecolor("black")
    legend.get_frame().set_facecolor("1")

    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.14, top=0.80)

    pdf_path = outdir / f"{output_stem}.pdf"
    png_path = outdir / f"{output_stem}.png"
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return pdf_path, png_path


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    data = load_catalog(args.catalog)
    _, base_mask = build_masks(data)

    pwd_bin_ids = assign_pwd_bins(data["pwd"])
    distance_pc = np.full_like(data["parallax"], np.nan, dtype=float)
    distance_pc[base_mask] = 1000.0 / data["parallax"][base_mask]

    if args.plot_all_points:
        sampled_indices = all_valid_indices(base_mask)
    else:
        sampled_indices = stratified_sample_indices(
            bin_ids=pwd_bin_ids,
            base_mask=base_mask,
            max_points_per_bin=args.max_points_per_bin,
            seed=args.seed,
        )

    summary_csv = args.outdir / "wd_pwd_bin_summary.csv"
    write_summary_csv(
        outpath=summary_csv,
        bin_ids=pwd_bin_ids,
        base_mask=base_mask,
        sampled_indices=sampled_indices,
        distance_pc=distance_pc,
    )

    pdf_path, png_path = plot_pwd_bins(
        ra_deg=data["ra"][sampled_indices],
        dec_deg=data["dec"][sampled_indices],
        distance_pc=distance_pc[sampled_indices],
        pwd_bin_ids=pwd_bin_ids[sampled_indices],
        outdir=args.outdir,
        distance_clip_pct=args.distance_clip_pct,
        template_name=args.template,
        dense_marker_scale=args.dense_marker_scale,
        output_stem=args.output_stem,
    )

    total_valid = int(np.count_nonzero(base_mask))
    print(f"Catalog path: {args.catalog}")
    print(f"Positive-parallax rows plotted from: {total_valid}")
    for bin_index, label in enumerate(PWD_BIN_LABELS):
        count = int(np.count_nonzero(base_mask & (pwd_bin_ids == bin_index)))
        sampled = int(np.count_nonzero(pwd_bin_ids[sampled_indices] == bin_index))
        print(f"Pwd {label}: {count} valid rows, {sampled} plotted")
    print(f"Summary CSV: {summary_csv}")
    print(f"Figure PDF: {pdf_path}")
    print(f"Figure PNG: {png_path}")


if __name__ == "__main__":
    main()
