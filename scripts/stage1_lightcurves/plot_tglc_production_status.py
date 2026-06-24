#!/usr/bin/env python3
"""Plot a sector-by-sector TGLC production status snapshot."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.plotting.style import apply_twirl_style  # noqa: E402


STAGE_ROWS = [
    "Prepared\n(cutouts)",
    "ePSF\nfitted",
    "TGLC HDF5\nLCs",
    "Final FITS\nproduct",
]

STATUS_TO_CODE = {
    "missing": 0,
    "pending": 0,
    "partial": 1,
    "active": 1,
    "complete": 2,
    "current": 4,
    "non_current": 3,
    "legacy": 3,
    "twirl_fs_v1": 3,
    "twirl_fs_v2": 4,
}

CODE_LABELS = {
    0: "missing / pending",
    1: "partial / active",
    2: "complete",
    3: "product exists, not current FSv2",
    4: "current TWIRL-FS v2",
}

CODE_COLORS = {
    0: "#eeeeee",
    1: "#f4a261",
    2: "#2a9d8f",
    3: "#577590",
    4: "#1b5e20",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def status_code(value: str) -> int:
    key = (value or "missing").strip().lower()
    return STATUS_TO_CODE.get(key, 0)


def build_matrix(rows: list[dict[str, str]]) -> np.ndarray:
    matrix = np.zeros((len(STAGE_ROWS), len(rows)), dtype=float)
    for j, row in enumerate(rows):
        matrix[0, j] = status_code(row["prepared_status"])
        matrix[1, j] = status_code(row["epsf_status"])
        matrix[2, j] = status_code(row["hdf5_status"])
        final = row["final_status"]
        product = row.get("hlsp_product", "")
        if final == "complete" and product == "twirl_fs_v2":
            matrix[3, j] = 4
        elif final == "complete":
            matrix[3, j] = 3
        else:
            matrix[3, j] = status_code(final)
    return matrix


def add_cell_labels(ax, rows: list[dict[str, str]], matrix: np.ndarray) -> None:
    for j, row in enumerate(rows):
        sector = int(row["sector"])
        source_count = int(row["source_count"])
        expected = int(row["source_expected"])
        if 0 < source_count < expected:
            ax.text(j, 0, f"{source_count/expected:.0%}", ha="center", va="center", fontsize=6)
        final_code = matrix[3, j]
        product = row.get("hlsp_product", "")
        if final_code == 4:
            label = "v2"
        elif product == "twirl_fs_v1":
            label = "v1"
        elif final_code == 3:
            label = "L"
        else:
            label = ""
        if label:
            color = "white"
            ax.text(j, 3, label, ha="center", va="center", fontsize=6.5, color=color)
        note = row.get("note", "")
        note_lower = note.lower()
        if "active" in note_lower or "stale" in note_lower:
            marker = "*" if "stale" in note_lower else "+"
            ax.text(j, -0.43, marker, ha="center", va="center", fontsize=8, color="0.15")


def plot_status(rows: list[dict[str, str]], out_png: Path, out_pdf: Path) -> None:
    apply_twirl_style("full_page")
    matrix = build_matrix(rows)
    sectors = [int(r["sector"]) for r in rows]
    codes = sorted(CODE_LABELS)
    cmap = ListedColormap([CODE_COLORS[c] for c in codes])

    fig, ax = plt.subplots(figsize=(11.0, 4.7), constrained_layout=False)
    ax.imshow(matrix, cmap=cmap, vmin=-0.5, vmax=4.5, aspect="auto")
    ax.set_xticks(np.arange(len(sectors)))
    ax.set_xticklabels([str(s) for s in sectors], rotation=90)
    ax.set_yticks(np.arange(len(STAGE_ROWS)))
    ax.set_yticklabels(STAGE_ROWS)
    ax.set_xlabel("TESS sector")
    ax.set_title("TGLC production status: S56-S93", pad=8)

    ax.set_xticks(np.arange(-0.5, len(sectors), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(STAGE_ROWS), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)

    if 64 in sectors:
        idx = sectors.index(64)
        ax.axvline(idx - 0.5, color="0.2", linestyle="--", linewidth=0.9)
        ax.text(
            idx - 0.6,
            -0.78,
            "normal finalize paused",
            ha="right",
            va="center",
            fontsize=7,
            color="0.15",
        )

    add_cell_labels(ax, rows, matrix)

    patches = [
        mpatches.Patch(facecolor=CODE_COLORS[c], edgecolor="0.25", label=CODE_LABELS[c])
        for c in codes
    ]
    legend = ax.legend(
        handles=patches,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0.0,
        frameon=True,
    )
    legend.get_frame().set_linewidth(0.8)

    snapshot = rows[0].get("snapshot_remote_time", "")
    note = (
        f"Snapshot: pdogpu1 {snapshot}. "
        "qc_pause.flag present; '+' marks active prep tail, '*' stale lease gap. "
        "Final-product labels: v2=current TWIRL-FS v2, v1=TWIRL-FS v1, L=legacy."
    )
    fig.text(0.08, 0.025, note, ha="left", va="bottom", fontsize=7)
    fig.subplots_adjust(left=0.08, right=0.78, bottom=0.22, top=0.84)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220)
    fig.savefig(out_pdf)
    plt.close(fig)


def write_summary(rows: list[dict[str, str]], path: Path, png_path: Path, pdf_path: Path) -> None:
    n = len(rows)
    complete_prep = sum(r["prepared_status"] == "complete" for r in rows)
    partial_prep = sum(r["prepared_status"] in {"partial", "active"} for r in rows)
    epsf = sum(r["epsf_status"] == "complete" for r in rows)
    hdf5 = sum(r["hdf5_status"] == "complete" for r in rows)
    final = sum(r["final_status"] == "complete" for r in rows)
    current_v2 = sum(r["hlsp_product"] == "twirl_fs_v2" for r in rows)
    active = [r for r in rows if "active" in r.get("note", "").lower()]
    stale = [r for r in rows if "stale" in r.get("note", "").lower()]

    lines = [
        "# TGLC production status S56-S93",
        "",
        f"Snapshot: `pdogpu1 {rows[0].get('snapshot_remote_time', '')}`.",
        "",
        f"- Prepared cutouts complete: `{complete_prep}/{n}` sectors; partial/active prep: `{partial_prep}`.",
        f"- ePSF fitted: `{epsf}/{n}` sectors.",
        f"- TGLC HDF5 light curves complete: `{hdf5}/{n}` sectors.",
        f"- Final FITS products exist: `{final}/{n}` sectors; current TWIRL-FS v2 sectors: `{current_v2}`.",
        "- `qc_pause.flag` is present, so normal GPU finalize is intentionally paused.",
        "",
        "Active / problematic tail:",
    ]
    for r in active + stale:
        lines.append(f"- S{r['sector']}: {r['note']}")
    lines.extend(
        [
            "",
            "Artifacts:",
            f"- PNG: `{png_path.name}`",
            f"- PDF: `{pdf_path.name}`",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", required=True, type=Path, help="Status snapshot CSV.")
    parser.add_argument("--out-dir", required=True, type=Path, help="Output directory.")
    parser.add_argument(
        "--stem",
        default="tglc_s56_s93_production_status",
        help="Output filename stem.",
    )
    args = parser.parse_args()

    rows = read_rows(args.csv)
    out_png = args.out_dir / f"{args.stem}.png"
    out_pdf = args.out_dir / f"{args.stem}.pdf"
    out_json = args.out_dir / f"{args.stem}.json"
    out_md = args.out_dir / "summary.md"

    plot_status(rows, out_png, out_pdf)
    write_summary(rows, out_md, out_png, out_pdf)

    payload = {"rows": rows, "png": str(out_png), "pdf": str(out_pdf)}
    out_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
