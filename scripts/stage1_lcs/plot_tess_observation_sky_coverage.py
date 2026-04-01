#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from matplotlib import patheffects
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm, ListedColormap

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.plotting.style import apply_twirl_style


DEFAULT_CATALOG = Path(
    "data_local/catalogs/twirl_master_catalog/twirl_wd_master_catalog_v0_tesscoverage.fits"
)
DEFAULT_OBSERVATIONS = Path(
    "data_local/catalogs/twirl_master_catalog/twirl_wd_tess_observations_v0.fits"
)
DEFAULT_OUTDIR = Path("reports/stage1_lcs")
DEFAULT_OUTPUT_STEM = "wd_tess_observation_sky_coverage"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot TWIRL TESS sky coverage with targets colored by the number of attached "
            "sectors or detector-hit records, split into observed-now and future-planned sectors."
        )
    )
    parser.add_argument(
        "--catalog",
        type=Path,
        default=DEFAULT_CATALOG,
        help="Coverage-enriched TWIRL master catalog FITS file.",
    )
    parser.add_argument(
        "--observations",
        type=Path,
        default=DEFAULT_OBSERVATIONS,
        help="One-row-per-hit TWIRL TESS observation table FITS file.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help="Directory for generated figure and summary products.",
    )
    parser.add_argument(
        "--output-stem",
        type=str,
        default=DEFAULT_OUTPUT_STEM,
        help="Base filename for output figure products.",
    )
    parser.add_argument(
        "--summary-name",
        type=str,
        default="wd_tess_observation_sky_coverage_summary.csv",
        help="Filename for the CSV summary sidecar.",
    )
    parser.add_argument(
        "--observed-sector-max",
        type=int,
        default=99,
        help="Treat sectors <= this value as currently observed and larger sectors as future-planned.",
    )
    parser.add_argument(
        "--count-mode",
        choices=("observations", "sectors"),
        default="sectors",
        help=(
            "Count attached detector hits or unique sectors per target within each panel's sector range."
        ),
    )
    parser.add_argument(
        "--color-scale",
        choices=("log", "linear"),
        default="log",
        help="Color normalization for nonzero count values.",
    )
    parser.add_argument(
        "--template",
        choices=("column", "full_page"),
        default="full_page",
        help="Named TWIRL plotting template from doc/plotting_style.md.",
    )
    parser.add_argument(
        "--highconf-only",
        action="store_true",
        help="Restrict the plot to the Pwd > 0.75 high-confidence sample.",
    )
    parser.add_argument(
        "--dense-marker-scale",
        type=float,
        default=0.16,
        help="Scale factor applied to the default dense scatter marker size.",
    )
    parser.add_argument(
        "--figure-title",
        type=str,
        default=None,
        help="Optional figure title. Omit to leave the figure untitled.",
    )
    parser.add_argument(
        "--observed-panel-title",
        type=str,
        default=None,
        help="Optional title for the observed-sector panel. Omit to leave it untitled.",
    )
    parser.add_argument(
        "--future-panel-title",
        type=str,
        default=None,
        help="Optional title for the future-sector panel. Omit to leave it untitled.",
    )
    return parser.parse_args()


def progress(message: str) -> None:
    print(message, flush=True)


def load_catalog(path: Path, highconf_only: bool) -> Table:
    table = Table.read(path, hdu=1, unit_parse_strict="silent")
    table.convert_bytestring_to_unicode()

    keep = np.isfinite(np.asarray(table["ra"], dtype=float)) & np.isfinite(np.asarray(table["dec"], dtype=float))
    if highconf_only and "is_highconf_wd" in table.colnames:
        keep &= np.asarray(table["is_highconf_wd"], dtype=bool)
    return table[keep]


def load_observations(path: Path) -> Table:
    table = Table.read(path, hdu=1, unit_parse_strict="silent")
    table.convert_bytestring_to_unicode()
    return table


def galactic_aitoff_coords(ra_deg: np.ndarray, dec_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    coords = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs").galactic
    lon_wrap = coords.l.wrap_at(180 * u.deg).radian
    return -lon_wrap, coords.b.radian


AITOFF_XTICK_DEG = np.array([-180, -120, -60, 0, 60, 120], dtype=float)
AITOFF_YTICK_DEG = np.array([-45, -30, -15, 0, 15, 30, 45, 60, 75], dtype=float)
AITOFF_LONGITUDE_LABELS = ["", "120°", "60°", "0°", "300°", "240°"]


def aggregate_counts(
    source_ids: np.ndarray,
    observations: Table,
    observed_sector_max: int,
    count_mode: str,
) -> tuple[np.ndarray, np.ndarray]:
    source_ids = np.asarray(source_ids, dtype=np.int64)
    index_by_source_id = {int(source_id): idx for idx, source_id in enumerate(source_ids)}

    observed_counts = np.zeros(len(source_ids), dtype=np.int16)
    future_counts = np.zeros(len(source_ids), dtype=np.int16)

    if len(observations) == 0:
        return observed_counts, future_counts

    obs_source_id = np.asarray(observations["source_id"], dtype=np.int64)
    obs_sector = np.asarray(observations["sector"], dtype=np.int16)

    if count_mode == "observations":
        for source_id, sector in zip(obs_source_id, obs_sector):
            idx = index_by_source_id.get(int(source_id))
            if idx is None:
                continue
            if int(sector) <= observed_sector_max:
                observed_counts[idx] += 1
            else:
                future_counts[idx] += 1
        return observed_counts, future_counts

    grouped_observed: dict[int, set[int]] = {}
    grouped_future: dict[int, set[int]] = {}
    for source_id, sector in zip(obs_source_id, obs_sector):
        idx = index_by_source_id.get(int(source_id))
        if idx is None:
            continue
        if int(sector) <= observed_sector_max:
            grouped_observed.setdefault(idx, set()).add(int(sector))
        else:
            grouped_future.setdefault(idx, set()).add(int(sector))

    for idx, sectors in grouped_observed.items():
        observed_counts[idx] = len(sectors)
    for idx, sectors in grouped_future.items():
        future_counts[idx] = len(sectors)
    return observed_counts, future_counts


def build_discrete_color_mapping(
    now_counts: np.ndarray,
    future_counts: np.ndarray,
    color_scale: str,
) -> tuple[ListedColormap, BoundaryNorm, int, np.ndarray]:
    combined_nonzero = np.concatenate(
        [now_counts[now_counts > 0], future_counts[future_counts > 0]]
    )
    if combined_nonzero.size == 0:
        cmap = ListedColormap([plt.get_cmap("viridis")(0.6)])
        norm = BoundaryNorm(np.array([0.5, 1.5]), ncolors=1, clip=True)
        return cmap, norm, 1, np.array([1], dtype=int)

    vmax = int(combined_nonzero.max())
    values = np.arange(1, vmax + 1, dtype=float)

    if vmax == 1:
        color_positions = np.array([0.6], dtype=float)
    elif color_scale == "log":
        color_positions = np.log10(values)
        color_positions = (color_positions - color_positions.min()) / (
            color_positions.max() - color_positions.min()
        )
    else:
        color_positions = np.linspace(0.0, 1.0, num=vmax)

    cmap = ListedColormap(plt.get_cmap("viridis")(color_positions))
    boundaries = np.arange(0.5, vmax + 1.5, 1.0)
    norm = BoundaryNorm(boundaries, ncolors=cmap.N, clip=True)

    if vmax <= 10:
        tick_values = np.arange(1, vmax + 1, dtype=int)
    else:
        if color_scale == "log":
            tick_values = np.rint(np.geomspace(1, vmax, num=8)).astype(int)
        else:
            tick_values = np.rint(np.linspace(1, vmax, num=8)).astype(int)
        tick_values = np.unique(np.clip(tick_values, 1, vmax))
        if tick_values[-1] != vmax:
            tick_values = np.append(tick_values, vmax)

    return cmap, norm, vmax, tick_values


def write_summary_csv(
    path: Path,
    observed_counts: np.ndarray,
    all_counts: np.ndarray,
    observed_sector_max: int,
    min_sector: int | None,
    count_mode: str,
    highconf_only: bool,
) -> None:
    all_panel_label = "all_coverage"
    if min_sector is not None:
        all_panel_label = f"all_ge_s{min_sector}"

    rows = [
        (
            f"observed_lt_s{observed_sector_max + 1}",
            int(np.count_nonzero(observed_counts > 0)),
            float(np.nanmedian(observed_counts[observed_counts > 0])) if np.any(observed_counts > 0) else 0.0,
            int(observed_counts.max()) if observed_counts.size else 0,
        ),
        (
            all_panel_label,
            int(np.count_nonzero(all_counts > 0)),
            float(np.nanmedian(all_counts[all_counts > 0])) if np.any(all_counts > 0) else 0.0,
            int(all_counts.max()) if all_counts.size else 0,
        ),
    ]

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["panel", "nonzero_targets", "median_positive_count", "max_count"])
        writer.writerows(rows)
        writer.writerow([])
        writer.writerow(["setting", "value"])
        writer.writerow(["count_mode", count_mode])
        writer.writerow(["observed_sector_max", observed_sector_max])
        writer.writerow(["highconf_only", highconf_only])


def plot_counts(
    catalog: Table,
    observed_counts: np.ndarray,
    all_counts: np.ndarray,
    outdir: Path,
    output_stem: str,
    summary_name: str,
    observed_sector_max: int,
    min_sector: int | None,
    count_mode: str,
    color_scale: str,
    template_name: str,
    dense_marker_scale: float,
    highconf_only: bool,
    figure_title: str | None,
    observed_panel_title: str | None,
    future_panel_title: str | None,
) -> tuple[Path, Path, Path]:
    template = apply_twirl_style(template_name)
    marker_size = template["dense_marker_size"] * dense_marker_scale
    background_size = max(0.3, marker_size - 0.2)

    ra_deg = np.asarray(catalog["ra"], dtype=float)
    dec_deg = np.asarray(catalog["dec"], dtype=float)
    sky_x, sky_y = galactic_aitoff_coords(ra_deg, dec_deg)

    cmap, norm, vmax, tick_values = build_discrete_color_mapping(
        observed_counts,
        all_counts,
        color_scale,
    )

    fig = plt.figure(figsize=(template["figsize"][0], template["figsize"][1] * 1.55))
    grid = fig.add_gridspec(
        2,
        1,
        hspace=0.28,
    )
    ax_now = fig.add_subplot(grid[0, 0], projection="aitoff")
    ax_future = fig.add_subplot(grid[1, 0], projection="aitoff")

    if observed_panel_title is None:
        observed_panel_title = f"Observed Now ({min_sector} <= Sector < {observed_sector_max + 1})"
    if future_panel_title is None:
        if min_sector is None:
            future_panel_title = "All Coverage"
        else:
            future_panel_title = f"All Coverage (Sector >= {min_sector})"

    panels = [
        (
            ax_now,
            observed_counts,
            observed_panel_title,
        ),
        (
            ax_future,
            all_counts,
            future_panel_title,
        ),
    ]

    for ax, counts, title in panels:
        ax.scatter(
            sky_x,
            sky_y,
            c="0.85",
            s=background_size,
            alpha=0.18,
            linewidths=0,
            rasterized=True,
        )

        positive = counts > 0
        if np.any(positive):
            order = np.argsort(counts[positive], kind="mergesort")
            ax.scatter(
                sky_x[positive][order],
                sky_y[positive][order],
                c=counts[positive][order],
                cmap=cmap,
                norm=norm,
                s=marker_size,
                alpha=0.35,
                linewidths=0,
                rasterized=True,
            )

        if title:
            ax.set_title(title, pad=20)
        ax.set_ylabel("Galactic latitude (degree)", fontsize=template["label_size"], labelpad=4)
        ax.yaxis.set_label_coords(-0.07, 0.5)
        ax.set_xlabel("Galactic longitude (degree)", fontsize=template["label_size"], labelpad=7)
        ax.set_xticks(np.radians(AITOFF_XTICK_DEG))
        ax.set_xticklabels([])
        ax.set_yticks(np.radians(AITOFF_YTICK_DEG))
        ax.tick_params(axis="both", labelsize=template["tick_size"], colors="0.28", pad=4)
        ax.grid(color="0.85", linewidth=template["grid_linewidth"])
        latitudes = np.linspace(-np.pi / 2.0, np.pi / 2.0, 361)
        for x_tick in np.radians(AITOFF_XTICK_DEG):
            ax.plot(
                np.full_like(latitudes, x_tick),
                latitudes,
                color="black",
                linewidth=0.55,
                linestyle=":",
                alpha=0.75,
                zorder=5.2,
            )
        ax.axhline(0.0, color="black", linewidth=0.75, linestyle="--", alpha=0.95, zorder=5.0)
        for x_tick, label in zip(np.radians(AITOFF_XTICK_DEG), AITOFF_LONGITUDE_LABELS):
            if not label:
                continue
            text = ax.text(
                x_tick,
                0.03,
                label,
                ha="center",
                va="center",
                color="black",
                fontsize=max(template["tick_size"] - 1.0, 5.0),
                zorder=5.5,
            )
            text.set_path_effects(
                [patheffects.withStroke(linewidth=1.3, foreground=(1.0, 1.0, 1.0, 0.8))]
            )

    count_label = (
        "Number of sectors per target"
        if count_mode == "sectors"
        else "Number of detector observations per target"
    )
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cax = fig.add_axes([0.28, 0.065, 0.44, 0.020])
    cbar = fig.colorbar(
        sm,
        cax=cax,
        orientation="horizontal",
    )
    cbar.set_label(count_label, labelpad=5)
    cbar.set_ticks(tick_values)

    if figure_title:
        fig.suptitle(figure_title, y=0.98)

    outdir.mkdir(parents=True, exist_ok=True)
    png_path = outdir / f"{output_stem}.png"
    pdf_path = outdir / f"{output_stem}.pdf"
    summary_path = outdir / summary_name

    fig.subplots_adjust(top=0.92 if figure_title else 0.95, bottom=0.16)
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    write_summary_csv(
        path=summary_path,
        observed_counts=observed_counts,
        all_counts=all_counts,
        observed_sector_max=observed_sector_max,
        min_sector=min_sector,
        count_mode=count_mode,
        highconf_only=highconf_only,
    )
    return png_path, pdf_path, summary_path


def main() -> None:
    args = parse_args()

    progress(f"[plot] loading catalog from {args.catalog}")
    catalog = load_catalog(args.catalog, highconf_only=args.highconf_only)
    progress(f"[plot] loaded catalog rows={len(catalog)}")
    progress(f"[plot] loading observation table from {args.observations}")
    observations = load_observations(args.observations)
    progress(f"[plot] loaded observation rows={len(observations)}")

    min_sector: int | None = None
    if len(observations) > 0:
        min_sector = int(np.min(np.asarray(observations["sector"], dtype=np.int16)))

    if args.highconf_only and "is_highconf_wd" in observations.colnames:
        observations = observations[np.asarray(observations["is_highconf_wd"], dtype=bool)]
        progress(f"[plot] filtered observation rows={len(observations)} to the high-confidence sample")
        if len(observations) > 0:
            min_sector = int(np.min(np.asarray(observations["sector"], dtype=np.int16)))
        else:
            min_sector = None

    observed_counts, future_counts = aggregate_counts(
        source_ids=np.asarray(catalog["source_id"], dtype=np.int64),
        observations=observations,
        observed_sector_max=args.observed_sector_max,
        count_mode=args.count_mode,
    )
    all_counts = observed_counts + future_counts
    progress(
        "[plot] counts ready: "
        f"observed_nonzero={int(np.count_nonzero(observed_counts > 0))}, "
        f"all_nonzero={int(np.count_nonzero(all_counts > 0))}, "
        f"count_mode={args.count_mode}"
    )

    png_path, pdf_path, summary_path = plot_counts(
        catalog=catalog,
        observed_counts=observed_counts,
        all_counts=all_counts,
        outdir=args.outdir,
        output_stem=args.output_stem,
        summary_name=args.summary_name,
        observed_sector_max=args.observed_sector_max,
        min_sector=min_sector,
        count_mode=args.count_mode,
        color_scale=args.color_scale,
        template_name=args.template,
        dense_marker_scale=args.dense_marker_scale,
        highconf_only=args.highconf_only,
        figure_title=args.figure_title,
        observed_panel_title=args.observed_panel_title,
        future_panel_title=args.future_panel_title,
    )

    progress(f"[plot] wrote PNG: {png_path}")
    progress(f"[plot] wrote PDF: {pdf_path}")
    progress(f"[plot] wrote summary CSV: {summary_path}")


if __name__ == "__main__":
    main()
