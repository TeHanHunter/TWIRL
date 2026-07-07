#!/usr/bin/env python3
"""Plot a compact S56 recovery50 CNN teacher diagnostic summary."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.plotting.style import apply_twirl_style


DEFAULT_ROOT = (
    REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue_2k"
)
DEFAULT_OUT_DIR = DEFAULT_ROOT / "model_summary"

PROFILE_DIRS = {
    "shape_only": "cnn_shape_only",
    "shape_plus_bls": "cnn_shape_plus_bls",
}

PROFILE_LABELS = {
    "shape_only": "Shape only",
    "shape_plus_bls": "Shape + BLS",
}

PROFILE_COLORS = {
    "shape_only": "#4C78A8",
    "shape_plus_bls": "#2F7F6F",
}

DISPLAY_CLASS = {
    "planet_like": "Planet",
    "instrumental_or_systematic": "Negative",
    "stellar_variability": "Variable",
}

LABEL_ROWS = [
    ("uncertain", "Flat/no signal\n(raw uncertain)", "negative target", "#E6E6E6"),
    ("instrumental_or_systematic", "Systematic", "negative target", "#E6E6E6"),
    ("planet_like", "Planet-like", "planet target", "#CFE4D4"),
    ("stellar_variability", "Variability", "variability target", "#D9E3F2"),
    ("eclipsing_binary_or_pceb", "EB/PCEB", "audit only", "#EFE4CC"),
    ("skip", "Skip", "excluded", "#F2F2F2"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
        help="S56 recovery50 teacher queue output root.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help="Output directory for PNG/PDF figure.",
    )
    parser.add_argument(
        "--png-dpi",
        type=int,
        default=320,
        help="PNG resolution.",
    )
    return parser.parse_args()


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def close_panel(ax: plt.Axes) -> None:
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("0.25")
        spine.set_linewidth(0.8)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        0.015,
        0.985,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        fontweight="bold",
        color="0.12",
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.0, "alpha": 0.9},
        clip_on=True,
    )


def draw_label_table(ax: plt.Axes, audit_summary: dict) -> None:
    counts = audit_summary.get("label_counts", {})
    total = int(audit_summary.get("n_labeled", sum(counts.values())))

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    close_panel(ax)
    panel_label(ax, "a")
    ax.set_title("Human labels to teacher roles", pad=8)

    left = 0.04
    col_label = left
    col_n = 0.56
    col_role = 0.70
    top = 0.88
    row_h = 0.115

    ax.add_patch(
        Rectangle(
            (0.02, top - 0.045),
            0.96,
            0.055,
            facecolor="#F0F0F0",
            edgecolor="none",
            zorder=0,
        )
    )
    ax.text(col_label, top - 0.015, "Raw human label", ha="left", va="center", weight="bold")
    ax.text(col_n, top - 0.015, "Rows", ha="right", va="center", weight="bold")
    ax.text(col_role, top - 0.015, "Teacher role", ha="left", va="center", weight="bold")

    y = top - 0.105
    for raw_label, display, role, fill in LABEL_ROWS:
        n = int(counts.get(raw_label, 0))
        ax.hlines(y + row_h * 0.42, 0.02, 0.98, color="0.88", linewidth=0.55)
        ax.text(col_label, y, display, ha="left", va="center", linespacing=1.05)
        ax.text(col_n, y, f"{n:,}", ha="right", va="center")
        ax.add_patch(
            Rectangle(
                (col_role - 0.012, y - 0.025),
                0.245,
                0.05,
                facecolor=fill,
                edgecolor="0.65",
                linewidth=0.45,
                zorder=0,
            )
        )
        ax.text(col_role + 0.004, y, role, ha="left", va="center", fontsize=6.2)
        y -= row_h

    ax.hlines(0.125, 0.02, 0.98, color="0.88", linewidth=0.55)
    ax.text(0.04, 0.075, f"Total labeled rows: {total:,}", ha="left", va="center")
    ax.text(
        0.04,
        0.033,
        "Inputs exclude truth/recovery/source columns.",
        ha="left",
        va="center",
        color="0.35",
        fontsize=6.2,
    )


def confusion_array(summary: dict, split: str = "test") -> tuple[list[str], np.ndarray]:
    classes = list(summary["classes"])
    cm = summary["metrics"][split]["confusion_matrix"]
    arr = np.zeros((len(classes), len(classes)), dtype=int)
    for i, true_class in enumerate(classes):
        row = cm.get(true_class, {})
        for j, pred_class in enumerate(classes):
            arr[i, j] = int(row.get(pred_class, 0))
    return classes, arr


def draw_confusion(ax: plt.Axes, summary: dict) -> None:
    classes, arr = confusion_array(summary, "test")
    labels = [DISPLAY_CLASS.get(name, name) for name in classes]
    max_value = max(int(arr.max()), 1)
    cmap = LinearSegmentedColormap.from_list(
        "twirl_blue_gray",
        ["#FFFFFF", "#DEE6EE", "#8CA1B6", "#23384E"],
    )

    ax.imshow(arr, cmap=cmap, vmin=0, vmax=max_value, aspect="equal")
    ax.set_box_aspect(1)
    ax.set_title("Held-out test confusion matrix", pad=8)
    ax.set_xlabel("Predicted label")
    ax.set_ylabel("Human label")
    ax.set_xticks(np.arange(len(labels)))
    ax.set_yticks(np.arange(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)
    ax.set_xticks(np.arange(-0.5, len(labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(labels), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)
    close_panel(ax)
    panel_label(ax, "b")

    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            value = int(arr[i, j])
            color = "white" if value > max_value * 0.42 else "0.12"
            weight = "bold" if i == j else "regular"
            ax.text(j, i, f"{value}", ha="center", va="center", color=color, weight=weight)

    metrics = summary["metrics"]["test"]
    ax.text(
        0.02,
        0.02,
        f"BA = {metrics['balanced_accuracy']:.3f}\nacc = {metrics['accuracy']:.3f}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        color="0.20",
        fontsize=6.3,
        bbox={"facecolor": "white", "edgecolor": "0.70", "linewidth": 0.4, "pad": 2.0},
    )


def load_histories(root: Path) -> dict[str, pd.DataFrame]:
    histories: dict[str, pd.DataFrame] = {}
    for profile, dirname in PROFILE_DIRS.items():
        path = root / "cnn_teacher" / dirname / "training_history.csv"
        histories[profile] = pd.read_csv(path)
    return histories


def draw_loss(ax: plt.Axes, histories: dict[str, pd.DataFrame]) -> None:
    for profile, history in histories.items():
        ax.plot(
            history["epoch"],
            history["train_loss"],
            color=PROFILE_COLORS[profile],
            linewidth=1.6,
            label=PROFILE_LABELS[profile],
        )
    ax.set_title("Training loss")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Cross-entropy")
    ax.set_xlim(left=1)
    ax.legend(loc="upper right", frameon=True, borderpad=0.35, handlelength=2.2)
    close_panel(ax)
    panel_label(ax, "c")


def draw_accuracy(ax: plt.Axes, histories: dict[str, pd.DataFrame]) -> None:
    for profile, history in histories.items():
        color = PROFILE_COLORS[profile]
        ax.plot(
            history["epoch"],
            history["validation_balanced_accuracy"],
            color=color,
            linewidth=1.6,
            label=f"{PROFILE_LABELS[profile]} validation",
        )
        ax.plot(
            history["epoch"],
            history["test_balanced_accuracy"],
            color=color,
            linewidth=1.35,
            linestyle=(0, (3.0, 1.7)),
            label=f"{PROFILE_LABELS[profile]} test",
        )
    ax.set_title("Validation and test accuracy")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Balanced accuracy")
    ax.set_xlim(left=1)
    ax.set_ylim(0.35, 1.02)
    ax.legend(
        loc="lower right",
        frameon=True,
        borderpad=0.35,
        handlelength=2.3,
        fontsize=6.2,
    )
    close_panel(ax)
    panel_label(ax, "d")


def make_figure(root: Path, out_dir: Path, png_dpi: int) -> tuple[Path, Path]:
    apply_twirl_style("full_page")
    audit_summary = load_json(root / "human_label_audit" / "summary.json")
    combined_summary = load_json(
        root / "cnn_teacher" / "cnn_shape_plus_bls" / "summary.json"
    )
    histories = load_histories(root)

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(7.1, 5.2),
        gridspec_kw={"height_ratios": [1.0, 1.0], "width_ratios": [1.0, 1.0]},
    )

    draw_label_table(axes[0, 0], audit_summary)
    draw_confusion(axes[0, 1], combined_summary)
    draw_loss(axes[1, 0], histories)
    draw_accuracy(axes[1, 1], histories)

    fig.subplots_adjust(
        left=0.075,
        right=0.985,
        top=0.93,
        bottom=0.095,
        wspace=0.32,
        hspace=0.55,
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / "s56_recovery50_2k_model_summary.png"
    pdf_path = out_dir / "s56_recovery50_2k_model_summary.pdf"
    fig.savefig(png_path, dpi=png_dpi)
    fig.savefig(pdf_path)
    plt.close(fig)
    return png_path, pdf_path


def main() -> None:
    args = parse_args()
    png_path, pdf_path = make_figure(args.root, args.out_dir, args.png_dpi)
    print(f"Wrote {png_path}")
    print(f"Wrote {pdf_path}")


if __name__ == "__main__":
    main()
