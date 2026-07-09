#!/usr/bin/env python3
"""Plot a triage-focused S56 recovery50 CNN teacher performance summary."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
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

CLASS_ORDER = [
    "planet_like",
    "instrumental_or_systematic",
    "stellar_variability",
]

DISPLAY_CLASS = {
    "planet_like": "Planet",
    "instrumental_or_systematic": "Negative",
    "stellar_variability": "Variable",
}

EXAMPLE_COLORS = {
    "TP": "#2F7F6F",
    "FP": "#B95C4A",
    "FN": "#7A5AA6",
    "TN": "#666666",
    "VAR": "#4C78A8",
    "LOW": "#9A7B30",
}


@dataclass(frozen=True)
class ExampleSpec:
    label: str
    title: str
    color_key: str
    row: pd.Series | None


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
        "--teacher-dirname",
        default="cnn_teacher",
        help="Teacher output subdirectory under --root.",
    )
    parser.add_argument(
        "--tensor-dirname",
        default="cnn_tensors",
        help="Tensor output subdirectory under --root.",
    )
    parser.add_argument(
        "--output-stem",
        default="s56_recovery50_2k_model_summary",
        help="Output filename stem.",
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
        zorder=10,
    )


def load_predictions(root: Path, teacher_dirname: str) -> dict[str, pd.DataFrame]:
    predictions: dict[str, pd.DataFrame] = {}
    for profile, dirname in PROFILE_DIRS.items():
        path = root / teacher_dirname / dirname / "cnn_predictions.parquet"
        frame = pd.read_parquet(path)
        predictions[profile] = frame.loc[
            frame["cnn_training_split"].fillna("").astype(str).eq("test")
        ].copy()
    return predictions


def average_precision_step(y_true: np.ndarray, score: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    y = np.asarray(y_true, dtype=bool)
    s = np.asarray(score, dtype=float)
    finite = np.isfinite(s)
    y = y[finite]
    s = s[finite]
    n_positive = int(y.sum())
    if y.size == 0 or n_positive == 0:
        return np.array([0.0]), np.array([0.0]), float("nan")
    order = np.argsort(-s, kind="mergesort")
    y_sorted = y[order]
    tp = np.cumsum(y_sorted)
    rank = np.arange(1, len(y_sorted) + 1)
    precision = tp / rank
    recall = tp / n_positive
    precision_plot = np.r_[1.0, precision]
    recall_plot = np.r_[0.0, recall]
    ap = float(np.sum((recall_plot[1:] - recall_plot[:-1]) * precision_plot[1:]))
    return recall_plot, precision_plot, ap


def review_budget_curve(y_true: np.ndarray, score: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    y = np.asarray(y_true, dtype=bool)
    s = np.asarray(score, dtype=float)
    finite = np.isfinite(s)
    y = y[finite]
    s = s[finite]
    n_positive = int(y.sum())
    if y.size == 0 or n_positive == 0:
        return np.array([0.0]), np.array([0.0])
    order = np.argsort(-s, kind="mergesort")
    y_sorted = y[order]
    x = np.arange(1, len(y_sorted) + 1, dtype=float) / float(len(y_sorted))
    recovered = np.cumsum(y_sorted) / float(n_positive)
    return np.r_[0.0, x], np.r_[0.0, recovered]


def y_planet_and_score(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    y = frame["main_teacher_target"].fillna("").astype(str).eq("planet_like").to_numpy()
    score = frame["cnn_p_planet_like"].to_numpy(dtype=float)
    return y, score


def metric_at_budget(x: np.ndarray, y: np.ndarray, budget: float) -> float:
    idx = int(np.searchsorted(x, budget, side="left"))
    idx = min(max(idx, 0), len(y) - 1)
    return float(y[idx])


def metric_at_threshold(frame: pd.DataFrame, threshold: float) -> tuple[float, float, int]:
    y_true, score = y_planet_and_score(frame)
    selected = score >= threshold
    n_selected = int(selected.sum())
    n_positive = int(y_true.sum())
    if n_selected == 0:
        return float("nan"), 0.0, 0
    tp = int(np.sum(selected & y_true))
    precision = tp / n_selected
    recall = tp / n_positive if n_positive else float("nan")
    return float(precision), float(recall), n_selected


def argmax_planet_point(frame: pd.DataFrame) -> tuple[float, float, int]:
    y_true = frame["main_teacher_target"].fillna("").astype(str).eq("planet_like").to_numpy()
    predicted = frame["cnn_label"].fillna("").astype(str).eq("planet_like").to_numpy()
    n_predicted = int(predicted.sum())
    n_positive = int(y_true.sum())
    if n_predicted == 0:
        return float("nan"), 0.0, 0
    tp = int(np.sum(predicted & y_true))
    precision = tp / n_predicted
    recall = tp / n_positive if n_positive else float("nan")
    return float(precision), float(recall), n_predicted


def draw_review_budget(ax: plt.Axes, predictions: dict[str, pd.DataFrame]) -> None:
    budget = 0.20
    for profile, frame in predictions.items():
        y_true, score = y_planet_and_score(frame)
        x, recovered = review_budget_curve(y_true, score)
        color = PROFILE_COLORS[profile]
        linewidth = 1.8 if profile == "shape_plus_bls" else 1.35
        zorder = 4 if profile == "shape_plus_bls" else 3
        label = f"{PROFILE_LABELS[profile]} ({metric_at_budget(x, recovered, budget):.0%} @ 20%)"
        if profile == "shape_plus_bls":
            ax.fill_between(x, 0.0, recovered, color=color, alpha=0.10, linewidth=0)
        ax.plot(x, recovered, color=color, linewidth=linewidth, label=label, zorder=zorder)

    ax.plot([0, 1], [0, 1], color="0.55", linewidth=1.0, linestyle=(0, (2, 2)), label="Random")
    combined_x, combined_y = review_budget_curve(*y_planet_and_score(predictions["shape_plus_bls"]))
    recovered_20 = metric_at_budget(combined_x, combined_y, budget)
    ax.scatter([budget], [recovered_20], s=24, color=PROFILE_COLORS["shape_plus_bls"], zorder=5)
    ax.annotate(
        f"{recovered_20:.0%} recovered\nat 20% review",
        xy=(budget, recovered_20),
        xytext=(0.34, 0.52),
        textcoords="axes fraction",
        arrowprops={"arrowstyle": "-", "color": "0.25", "linewidth": 0.7},
        ha="left",
        va="center",
        fontsize=7,
        bbox={"facecolor": "white", "edgecolor": "0.70", "linewidth": 0.4, "pad": 2.0},
    )
    ax.set_title("Human-review efficiency")
    ax.set_xlabel("Fraction of held-out rows reviewed")
    ax.set_ylabel("Planet-like labels recovered")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.03)
    ax.legend(loc="lower right", frameon=True, borderpad=0.35, handlelength=2.2)
    close_panel(ax)
    panel_label(ax, "a")


def draw_precision_recall(ax: plt.Axes, predictions: dict[str, pd.DataFrame]) -> None:
    combined_frame = predictions["shape_plus_bls"]
    y_combined, _ = y_planet_and_score(combined_frame)
    prevalence = float(np.mean(y_combined)) if len(y_combined) else float("nan")

    for profile, frame in predictions.items():
        y_true, score = y_planet_and_score(frame)
        recall, precision, ap = average_precision_step(y_true, score)
        color = PROFILE_COLORS[profile]
        linewidth = 1.8 if profile == "shape_plus_bls" else 1.35
        ax.plot(
            recall,
            precision,
            color=color,
            linewidth=linewidth,
            label=f"{PROFILE_LABELS[profile]} AP={ap:.3f}",
            zorder=4 if profile == "shape_plus_bls" else 3,
        )

    if np.isfinite(prevalence):
        ax.axhline(
            prevalence,
            color="0.55",
            linestyle=(0, (2, 2)),
            linewidth=1.0,
            label=f"Prevalence={prevalence:.2f}",
        )

    arg_precision, arg_recall, arg_n = argmax_planet_point(combined_frame)
    if np.isfinite(arg_precision):
        ax.scatter(
            [arg_recall],
            [arg_precision],
            s=28,
            color=PROFILE_COLORS["shape_plus_bls"],
            edgecolor="white",
            linewidth=0.6,
            zorder=5,
        )
        ax.annotate(
            f"argmax\nn={arg_n}",
            xy=(arg_recall, arg_precision),
            xytext=(0.58, 0.74),
            textcoords="axes fraction",
            arrowprops={"arrowstyle": "-", "color": "0.25", "linewidth": 0.7},
            fontsize=7,
            ha="left",
            va="center",
            bbox={"facecolor": "white", "edgecolor": "0.70", "linewidth": 0.4, "pad": 2.0},
        )

    pseudo_precision, pseudo_recall, pseudo_n = metric_at_threshold(combined_frame, 0.98)
    if pseudo_n > 0 and np.isfinite(pseudo_precision):
        ax.scatter(
            [pseudo_recall],
            [pseudo_precision],
            s=28,
            marker="D",
            color="#9A7B30",
            edgecolor="white",
            linewidth=0.6,
            zorder=5,
        )
        ax.annotate(
            f"p>=0.98\nn={pseudo_n}",
            xy=(pseudo_recall, pseudo_precision),
            xytext=(0.13, 0.88),
            textcoords="axes fraction",
            arrowprops={"arrowstyle": "-", "color": "0.25", "linewidth": 0.7},
            fontsize=7,
            ha="left",
            va="center",
            bbox={"facecolor": "white", "edgecolor": "0.70", "linewidth": 0.4, "pad": 2.0},
        )

    ax.set_title("Planet-like precision-recall")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.03)
    ax.legend(loc="lower left", frameon=True, borderpad=0.35, handlelength=2.1)
    close_panel(ax)
    panel_label(ax, "b")


def confusion_from_predictions(frame: pd.DataFrame) -> np.ndarray:
    arr = np.zeros((len(CLASS_ORDER), len(CLASS_ORDER)), dtype=int)
    true_labels = frame["main_teacher_target"].fillna("").astype(str)
    pred_labels = frame["cnn_label"].fillna("").astype(str)
    for i, true_class in enumerate(CLASS_ORDER):
        mask = true_labels.eq(true_class)
        for j, pred_class in enumerate(CLASS_ORDER):
            arr[i, j] = int(np.sum(mask & pred_labels.eq(pred_class)))
    return arr


def draw_confusion(ax: plt.Axes, frame: pd.DataFrame) -> None:
    counts = confusion_from_predictions(frame)
    row_sums = counts.sum(axis=1, keepdims=True)
    norm = np.divide(counts, row_sums, out=np.zeros_like(counts, dtype=float), where=row_sums > 0)
    labels = [DISPLAY_CLASS.get(name, name) for name in CLASS_ORDER]
    cmap = LinearSegmentedColormap.from_list(
        "twirl_confusion",
        ["#FFFFFF", "#E2E8EC", "#8DA0AE", "#263B4E"],
    )
    ax.imshow(norm, cmap=cmap, vmin=0.0, vmax=1.0, aspect="equal")
    ax.set_box_aspect(1)
    ax.set_title("Held-out confusion matrix")
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
    for i in range(counts.shape[0]):
        for j in range(counts.shape[1]):
            value = counts[i, j]
            pct = norm[i, j]
            color = "white" if pct > 0.55 else "0.12"
            weight = "bold" if i == j else "regular"
            ax.text(
                j,
                i,
                f"{value}\n{pct:.0%}",
                ha="center",
                va="center",
                color=color,
                weight=weight,
                fontsize=6.4,
                linespacing=0.95,
            )
    close_panel(ax)
    panel_label(ax, "c")


def pick_first(frame: pd.DataFrame, mask: pd.Series, sort_col: str, ascending: bool = False) -> pd.Series | None:
    subset = frame.loc[mask].copy()
    if subset.empty:
        return None
    return subset.sort_values(sort_col, ascending=ascending).iloc[0]


def select_examples(frame: pd.DataFrame) -> list[ExampleSpec]:
    true_label = frame["main_teacher_target"].fillna("").astype(str)
    pred_label = frame["cnn_label"].fillna("").astype(str)
    is_planet = true_label.eq("planet_like")
    pred_planet = pred_label.eq("planet_like")
    is_negative = true_label.eq("instrumental_or_systematic")
    is_variable = true_label.eq("stellar_variability")

    correct = true_label.eq(pred_label)
    margin_sorted = frame.assign(_abs_margin=frame["cnn_margin"].abs())

    return [
        ExampleSpec(
            "TP",
            "Planet recovered",
            "TP",
            pick_first(frame, is_planet & pred_planet, "cnn_p_planet_like", ascending=False),
        ),
        ExampleSpec(
            "FN",
            "Planet missed",
            "FN",
            pick_first(frame, is_planet & ~pred_planet, "cnn_p_planet_like", ascending=True),
        ),
        ExampleSpec(
            "FP",
            "False alarm",
            "FP",
            pick_first(frame, ~is_planet & pred_planet, "cnn_p_planet_like", ascending=False),
        ),
        ExampleSpec(
            "TN",
            "Flat/systematic rejected",
            "TN",
            pick_first(frame, is_negative & pred_label.eq("instrumental_or_systematic"), "cnn_confidence", ascending=False),
        ),
        ExampleSpec(
            "VAR",
            "Variable recovered",
            "VAR",
            pick_first(frame, is_variable & pred_label.eq("stellar_variability"), "cnn_confidence", ascending=False),
        ),
        ExampleSpec(
            "LOW",
            "Low-margin decision",
            "LOW",
            pick_first(margin_sorted, correct, "_abs_margin", ascending=True),
        ),
    ]


def load_folded_tensors(root: Path, tensor_dirname: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    path = root / tensor_dirname / "recovery50_cnn_tensors.npz"
    tensors = np.load(path, allow_pickle=False)
    return tensors["folded"], tensors["folded_mask"], tensors["folded_x_duration"]


def robust_ylim(y: np.ndarray) -> tuple[float, float]:
    finite = y[np.isfinite(y)]
    if finite.size == 0:
        return -1.0, 1.0
    lo, hi = np.nanpercentile(finite, [8, 92])
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        center = float(np.nanmedian(finite))
        return center - 1.0, center + 1.0
    center = float(np.nanmedian(finite))
    half_span = max(abs(hi - center), abs(lo - center)) * 1.35
    half_span = min(max(half_span, 0.35), 65.0)
    return center - half_span, center + half_span


def draw_one_example(
    ax: plt.Axes,
    spec: ExampleSpec,
    folded: np.ndarray,
    folded_mask: np.ndarray,
    x_grid: np.ndarray,
) -> None:
    close_panel(ax)
    ax.axvspan(-0.5, 0.5, color="0.93", zorder=0)
    ax.axvline(0.0, color="0.72", linewidth=0.7, zorder=1)
    ax.axhline(0.0, color="0.80", linewidth=0.7, zorder=1)
    color = EXAMPLE_COLORS[spec.color_key]
    if spec.row is None:
        ax.text(0.5, 0.5, "No example", transform=ax.transAxes, ha="center", va="center")
        ax.set_xticks([])
        ax.set_yticks([])
        return

    tensor_index = int(spec.row["tensor_index"])
    y = np.asarray(folded[tensor_index, 0, :], dtype=float) * 100.0
    mask = np.asarray(folded_mask[tensor_index, 0, :], dtype=bool) & np.isfinite(y)
    ymin, ymax = robust_ylim(y[mask])
    y_plot = np.clip(y[mask], ymin, ymax)
    smooth = (
        pd.Series(y_plot)
        .rolling(window=9, center=True, min_periods=1)
        .median()
        .to_numpy(dtype=float)
    )
    ax.scatter(x_grid[mask], y_plot, s=3.0, color=color, alpha=0.25, linewidths=0)
    ax.plot(x_grid[mask], smooth, color=color, linewidth=1.45)
    ax.set_xlim(float(np.nanmin(x_grid)), float(np.nanmax(x_grid)))
    ax.set_ylim(ymin, ymax)
    human = DISPLAY_CLASS.get(str(spec.row["main_teacher_target"]), str(spec.row["main_teacher_target"]))
    pred = DISPLAY_CLASS.get(str(spec.row["cnn_label"]), str(spec.row["cnn_label"]))
    p_planet = float(spec.row["cnn_p_planet_like"])
    short_title = {
        "TP": "TP: recovered",
        "FN": "FN: missed",
        "FP": "FP: false alarm",
        "TN": "TN: rejected",
        "VAR": "VAR: recovered",
        "LOW": "Low margin",
    }.get(spec.label, spec.title)
    ax.set_title(short_title, fontsize=6.1, pad=2)
    ax.text(
        0.03,
        0.06,
        f"H:{human}  M:{pred}\np_planet={p_planet:.2f}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=5.4,
        color="0.15",
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.2, "alpha": 0.78},
    )
    ax.tick_params(labelsize=5.5, length=2.0)


def draw_examples(
    fig: plt.Figure,
    parent_spec,
    frame: pd.DataFrame,
    folded: np.ndarray,
    folded_mask: np.ndarray,
    x_grid: np.ndarray,
) -> None:
    subgrid = parent_spec.subgridspec(
        3,
        3,
        height_ratios=[0.14, 1.0, 1.0],
        wspace=0.22,
        hspace=0.42,
    )
    header_ax = fig.add_subplot(subgrid[0, :])
    header_ax.set_axis_off()
    header_ax.text(
        0.015,
        0.98,
        "d",
        transform=header_ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        fontweight="bold",
        color="0.12",
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.0, "alpha": 0.9},
    )
    header_ax.text(
        0.5,
        0.65,
        "Held-out folded examples",
        transform=header_ax.transAxes,
        ha="center",
        va="center",
        fontsize=9,
        color="0.15",
    )
    examples = select_examples(frame)
    axes: list[plt.Axes] = []
    for idx, spec in enumerate(examples):
        ax = fig.add_subplot(subgrid[1 + idx // 3, idx % 3])
        draw_one_example(ax, spec, folded, folded_mask, x_grid)
        axes.append(ax)
        if idx // 3 == 1:
            ax.set_xlabel("phase / duration", fontsize=6)
        else:
            ax.set_xticklabels([])
        if idx % 3 == 0:
            ax.set_ylabel("rel. flux (%)", fontsize=6)
        else:
            ax.set_yticklabels([])


def make_figure(
    root: Path,
    out_dir: Path,
    png_dpi: int,
    *,
    teacher_dirname: str,
    tensor_dirname: str,
    output_stem: str,
) -> tuple[Path, Path]:
    apply_twirl_style("full_page")
    predictions = load_predictions(root, teacher_dirname)
    combined = predictions["shape_plus_bls"]
    folded, folded_mask, folded_x = load_folded_tensors(root, tensor_dirname)
    summary = load_json(root / teacher_dirname / "cnn_shape_plus_bls" / "summary.json")

    fig = plt.figure(figsize=(7.1, 6.15))
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.0, 1.05],
        width_ratios=[1.0, 1.0],
        left=0.075,
        right=0.985,
        top=0.875,
        bottom=0.085,
        wspace=0.30,
        hspace=0.42,
    )

    draw_review_budget(fig.add_subplot(grid[0, 0]), predictions)
    draw_precision_recall(fig.add_subplot(grid[0, 1]), predictions)
    draw_confusion(fig.add_subplot(grid[1, 0]), combined)
    draw_examples(fig, grid[1, 1], combined, folded, folded_mask, folded_x)

    fig.text(
        0.075,
        0.965,
        "S56 recovery50 CNN teacher: held-out triage performance",
        ha="left",
        va="top",
        fontsize=9.5,
        color="0.12",
    )
    fig.text(
        0.985,
        0.965,
        (
            f"test n={summary['metrics']['test']['n']}; "
            f"balanced accuracy={summary['metrics']['test']['balanced_accuracy']:.3f}"
        ),
        ha="right",
        va="top",
        fontsize=7,
        color="0.35",
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / f"{output_stem}.png"
    pdf_path = out_dir / f"{output_stem}.pdf"
    fig.savefig(png_path, dpi=png_dpi)
    fig.savefig(pdf_path)
    plt.close(fig)
    return png_path, pdf_path


def main() -> None:
    args = parse_args()
    png_path, pdf_path = make_figure(
        args.root,
        args.out_dir,
        args.png_dpi,
        teacher_dirname=args.teacher_dirname,
        tensor_dirname=args.tensor_dirname,
        output_stem=args.output_stem,
    )
    print(f"Wrote {png_path}")
    print(f"Wrote {pdf_path}")


if __name__ == "__main__":
    main()
