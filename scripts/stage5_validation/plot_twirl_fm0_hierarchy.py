#!/usr/bin/env python3
"""Render the TWIRL-FM0 data hierarchy and downstream model tree."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from twirl.plotting.style import apply_twirl_style, get_ordered_palette


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT_DIR = (
    REPO_ROOT / "reports" / "exploratory" / "twirl_fm0_working_plan_v0_1"
)
DESIGN_VERSION = "TWIRL-FM0.1"
PDF_METADATA_TIMESTAMP = datetime(2026, 8, 14, 12, tzinfo=timezone.utc)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _box(
    ax: plt.Axes,
    *,
    xy: tuple[float, float],
    width: float,
    height: float,
    text: str,
    edgecolor: object,
    face_alpha: float = 0.10,
    fontsize: float = 9.0,
    linewidth: float = 1.4,
    textcolor: object = "0.12",
    align: str = "center",
) -> FancyBboxPatch:
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.012",
        linewidth=linewidth,
        edgecolor=edgecolor,
        facecolor=to_rgba(edgecolor, face_alpha),
        zorder=2,
    )
    ax.add_patch(patch)
    x, y = xy
    if align == "left":
        text_x = x + 0.018
        horizontal_alignment = "left"
    else:
        text_x = x + width / 2.0
        horizontal_alignment = "center"
    ax.text(
        text_x,
        y + height / 2.0,
        text,
        ha=horizontal_alignment,
        va="center",
        fontsize=fontsize,
        color=textcolor,
        linespacing=1.25,
        zorder=3,
    )
    return patch


def _arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: object,
    linewidth: float = 1.5,
    connectionstyle: str = "arc3,rad=0",
) -> None:
    """Draw one directed connector from ``start`` to destination ``end``."""
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=12,
            linewidth=linewidth,
            color=color,
            connectionstyle=connectionstyle,
            shrinkA=1.5,
            shrinkB=1.5,
            zorder=4,
        )
    )


def _connector(
    ax: plt.Axes,
    x_values: list[float],
    y_values: list[float],
    *,
    color: object = "0.55",
    linewidth: float = 1.15,
) -> None:
    ax.plot(x_values, y_values, color=color, linewidth=linewidth, zorder=1)


def _routed_arrow(
    ax: plt.Axes,
    points: list[tuple[float, float]],
    *,
    color: object,
    linewidth: float = 1.35,
) -> None:
    """Draw an orthogonal route with one arrowhead at the final destination."""
    if len(points) < 2:
        raise ValueError("a routed arrow requires at least two points")
    for start, end in zip(points[:-2], points[1:-1], strict=True):
        _connector(
            ax,
            [start[0], end[0]],
            [start[1], end[1]],
            color=color,
            linewidth=linewidth,
        )
    _arrow(
        ax,
        points[-2],
        points[-1],
        color=color,
        linewidth=linewidth,
    )


def build_figure() -> plt.Figure:
    apply_twirl_style("full_page")
    palette = get_ordered_palette(5, "viridis")
    purple, blue, teal, green, gold = palette
    red = "#A23B3B"
    gray = "0.48"

    fig, ax = plt.subplots(figsize=(12.6, 7.4))
    fig.subplots_adjust(left=0.025, right=0.985, top=0.965, bottom=0.075)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis("off")
    ax.grid(False)

    ax.text(
        0.025,
        0.965,
        "A   One physical host, many observations",
        ha="left",
        va="top",
        fontsize=12.2,
        fontweight="bold",
        color="0.12",
    )
    ax.text(
        0.405,
        0.965,
        "B   Shared encoder and three reusable representation levels",
        ha="left",
        va="top",
        fontsize=12.2,
        fontweight="bold",
        color="0.12",
    )

    _box(
        ax,
        xy=(0.025, 0.735),
        width=0.145,
        height=0.145,
        text="Physical source\nGaia DR3 source_id",
        edgecolor=blue,
        fontsize=9.2,
        linewidth=1.7,
    )
    ax.text(
        0.0975,
        0.703,
        "source partition / sampling group\nincludes every linked TIC alias",
        ha="center",
        va="top",
        fontsize=7.7,
        color=blue,
        fontweight="bold",
        linespacing=1.2,
    )

    sector_y = (0.805, 0.655, 0.505)
    sector_text = (
        "Sector visit A\nwindows 1 … k",
        "Sector visit B\nwindows 1 … m",
        "Sector visit N\nwindows 1 … q",
    )
    for y, label in zip(sector_y, sector_text, strict=True):
        _box(
            ax,
            xy=(0.245, y - 0.050),
            width=0.125,
            height=0.100,
            text=label,
            edgecolor=blue,
            fontsize=8.5,
            face_alpha=0.075,
        )

    _connector(ax, [0.170, 0.205], [0.8075, 0.8075], color=blue)
    _connector(ax, [0.205, 0.205], [0.505, 0.805], color=blue)
    for y in sector_y:
        _connector(ax, [0.205, 0.245], [y, y], color=blue)
    ax.text(
        0.205,
        0.900,
        "contains",
        ha="center",
        va="bottom",
        fontsize=7.6,
        color=gray,
    )
    ax.text(
        0.025,
        0.425,
        "The host enters one source partition before windows are sampled.\n"
        "Sector visits and windows remain repeated observations;\n"
        "they are not averaged into one light curve or counted as new stars.",
        ha="left",
        va="top",
        fontsize=8.6,
        color="0.28",
        linespacing=1.35,
    )

    # Every sector visit follows the same shared encoder path.  The right-hand
    # collection bus avoids visually privileging the middle example visit.
    for y in sector_y:
        _connector(ax, [0.370, 0.395], [y, y], color=teal, linewidth=1.35)
    _connector(ax, [0.395, 0.395], [sector_y[-1], sector_y[0]], color=teal, linewidth=1.35)
    _arrow(ax, (0.395, 0.655), (0.425, 0.655), color=teal, linewidth=1.8)

    _box(
        ax,
        xy=(0.425, 0.585),
        width=0.135,
        height=0.140,
        text="Shared cadence +\nwindow encoder\nlocal stem + context",
        edgecolor=teal,
        fontsize=9.0,
        linewidth=1.8,
    )
    _arrow(ax, (0.560, 0.645), (0.585, 0.645), color=teal, linewidth=1.8)
    _box(
        ax,
        xy=(0.585, 0.600),
        width=0.095,
        height=0.090,
        text="Per visit\nz_window[1 … k]",
        edgecolor=teal,
        fontsize=7.8,
        linewidth=1.6,
    )
    _arrow(ax, (0.680, 0.645), (0.705, 0.645), color=green, linewidth=1.8)
    _box(
        ax,
        xy=(0.705, 0.590),
        width=0.120,
        height=0.130,
        text="Sector\naggregator",
        edgecolor=green,
        fontsize=9.2,
        linewidth=1.7,
    )
    _arrow(ax, (0.825, 0.645), (0.842, 0.645), color=green, linewidth=1.8)
    _box(
        ax,
        xy=(0.842, 0.600),
        width=0.105,
        height=0.090,
        text="Across visits\nz_sector[A … N]",
        edgecolor=green,
        fontsize=7.8,
        linewidth=1.6,
    )

    _arrow(
        ax,
        (0.8945, 0.600),
        (0.800, 0.510),
        color=gold,
        linewidth=1.7,
        connectionstyle="arc3,rad=-0.08",
    )
    _box(
        ax,
        xy=(0.735, 0.385),
        width=0.130,
        height=0.120,
        text="Multi-sector\nhost aggregator",
        edgecolor=gold,
        fontsize=9.1,
        linewidth=1.7,
    )
    _arrow(ax, (0.865, 0.445), (0.900, 0.445), color=gold, linewidth=1.7)
    _box(
        ax,
        xy=(0.900, 0.400),
        width=0.072,
        height=0.090,
        text="z_host",
        edgecolor=gold,
        fontsize=8.7,
        linewidth=1.6,
    )

    head_specs = (
        (
            (0.420, 0.120),
            0.165,
            "Window / cadence heads\nevent localization\nmasked reconstruction",
            [(0.6325, 0.600), (0.6325, 0.350), (0.5025, 0.350), (0.5025, 0.265)],
        ),
        (
            (0.610, 0.120),
            0.165,
            "Sector heads\nmorphology classification\ndepth + duration + uncertainty",
            [(0.8945, 0.600), (0.8945, 0.350), (0.6925, 0.350), (0.6925, 0.265)],
        ),
        (
            (0.800, 0.120),
            0.175,
            "Host heads\nmulti-sector classification\nretrieval + anomaly search",
            [(0.936, 0.400), (0.936, 0.310), (0.8875, 0.310), (0.8875, 0.265)],
        ),
    )
    for xy, width, label, route in head_specs:
        _box(
            ax,
            xy=xy,
            width=width,
            height=0.145,
            text=label,
            edgecolor=purple,
            fontsize=8.4,
            face_alpha=0.075,
        )
        _routed_arrow(
            ax,
            route,
            color=purple,
            linewidth=1.35,
        )

    _box(
        ax,
        xy=(0.025, 0.185),
        width=0.335,
        height=0.160,
        text=(
            "Model-visible in TWIRL-FM0.1\n"
            "local + host-relative time in 200 s units • validity masks\n"
            "raw-relative, ADP, ADP015 × 1×1 and 3×3 apertures\n"
            "two aperture-level measurement-error proxies\n"
            "not the A2v1 FITS error columns"
        ),
        edgecolor=teal,
        fontsize=7.8,
        face_alpha=0.065,
        align="left",
    )
    _box(
        ax,
        xy=(0.025, 0.035),
        width=0.335,
        height=0.120,
        text=(
            "Not exposed during pretraining\n"
            "BLS, phase folds, candidates, labels, injections, Teacher scores,\n"
            "Gaia/TIC identity, or explicit sector/detector IDs"
        ),
        edgecolor=red,
        fontsize=8.3,
        face_alpha=0.035,
        align="left",
        linewidth=1.6,
    )

    fig.text(
        0.025,
        0.018,
        (
            "TWIRL-FM0.1 frozen design overview. Arrowheads point toward the receiving "
            "representation or downstream head; branch lines at left denote containment."
        ),
        ha="left",
        va="bottom",
        fontsize=7.2,
        color="0.42",
    )
    return fig


def write_outputs(output_dir: Path) -> tuple[Path, Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = output_dir / "twirl_fm0_hierarchy.png"
    pdf_path = output_dir / "twirl_fm0_hierarchy.pdf"
    provenance_path = output_dir / "twirl_fm0_hierarchy.provenance.json"

    figure = build_figure()
    figure.savefig(png_path, dpi=240)
    figure.savefig(
        pdf_path,
        metadata={
            "Title": "TWIRL-FM0 data hierarchy and downstream model tree",
            "Author": "TWIRL",
            "Subject": DESIGN_VERSION,
            "CreationDate": PDF_METADATA_TIMESTAMP,
            "ModDate": PDF_METADATA_TIMESTAMP,
        },
    )
    plt.close(figure)

    source_documents = [
        "doc/foundation_model_design.md",
        "configs/models/twirl_fm0_1_s56_s67_poc.yaml",
        "reports/stage5_validation/twirl_fm0_1_design_freeze_v1/freeze.json",
        "doc/twirl_plan.md",
        "doc/a2v1_production_protocol.md",
        "doc/plotting_style.md",
    ]
    provenance = {
        "artifact": "twirl_fm0_hierarchy",
        "design_version": DESIGN_VERSION,
        "status": "scientific_design_frozen_implementation_not_started",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "generator": str(Path(__file__).resolve().relative_to(REPO_ROOT)),
        "source_documents": source_documents,
        "source_sha256": {
            source: _sha256(REPO_ROOT / source) for source in source_documents
        },
        "model_visible_summary": [
            "local and host-relative elapsed time in nominal 200 s units",
            "effective-validity and channel masks",
            "raw-relative, ADP, and ADP015 views for 1x1 and 3x3 apertures",
            "two aperture-level measurement-error proxies, not A2v1 FITS error columns",
        ],
        "pretraining_forbidden_summary": [
            "BLS and phase-folded products",
            "candidate selection",
            "labels and injections",
            "Teacher scores",
            "Gaia/TIC identity",
            "explicit sector/detector identifiers",
        ],
        "outputs": {
            png_path.name: {"sha256": _sha256(png_path)},
            pdf_path.name: {"sha256": _sha256(pdf_path)},
        },
    }
    provenance_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return png_path, pdf_path, provenance_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    for path in write_outputs(args.output_dir.resolve()):
        print(path)


if __name__ == "__main__":
    main()
