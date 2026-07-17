#!/usr/bin/env python3
"""Plot the selected S56 harmonic-CNN ensemble performance diagnostics."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.plotting.style import apply_twirl_style, get_ordered_palette


DEFAULT_ROOT = (
    REPO_ROOT
    / "reports/stage5_validation/s56_label_adjudication_real343/retrained/s56_harmonic_cnn_v1"
)
CLASS_KEYS = ("planet_like", "eclipse_contact", "smooth_variable", "other")
CLASS_LABELS = ("Planet-like", "Eclipse/\ncontact", "Smooth\nvariable", "Other")
PROFILE_LABELS = {
    "metadata_only": "Metadata only",
    "single_period_native_fold": "Single-period native fold",
    "seven_harmonic_shape": "Seven-harmonic shape",
    "shape_plus_raw_chronology": "Shape + raw chronology",
    "shape_plus_periodogram_bls": "Shape + periodogram/BLS",
    "full_combined": "Full combined",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument(
        "--output-stem",
        default="s56_harmonic_cnn_v1_performance",
    )
    parser.add_argument("--png-dpi", type=int, default=320)
    return parser.parse_args()


def _load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _close_panel(ax: plt.Axes) -> None:
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("0.25")
        spine.set_linewidth(0.8)


def _panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.14,
        1.04,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        fontweight="bold",
        color="0.12",
        clip_on=False,
    )


def _confusion_matrix(frame: pd.DataFrame) -> np.ndarray:
    truth = pd.to_numeric(frame["morphology_target_index"], errors="raise").to_numpy(int)
    predicted = pd.to_numeric(
        frame["morphology_prediction_index"], errors="raise"
    ).to_numpy(int)
    n_classes = len(CLASS_KEYS)
    active = (truth >= 0) & (truth < n_classes)
    truth = truth[active]
    predicted = predicted[active]
    if not len(truth):
        raise ValueError("fixed test contains no active four-class morphology targets")
    if np.any(truth >= n_classes):
        raise ValueError("fixed-test morphology targets are outside the four-class contract")
    if np.any((predicted < 0) | (predicted >= n_classes)):
        raise ValueError("fixed-test predictions are outside the four-class contract")
    matrix = np.zeros((n_classes, n_classes), dtype=int)
    np.add.at(matrix, (truth, predicted), 1)
    return matrix


def _draw_confusion(ax: plt.Axes, frame: pd.DataFrame, color: tuple[float, ...]) -> None:
    matrix = _confusion_matrix(frame)
    support = matrix.sum(axis=1)
    normalized = np.divide(
        matrix,
        support[:, None],
        out=np.zeros_like(matrix, dtype=float),
        where=support[:, None] > 0,
    )
    cmap = LinearSegmentedColormap.from_list(
        "twirl_confusion",
        [(1.0, 1.0, 1.0), (*color[:3], 1.0)],
    )
    ax.imshow(normalized, vmin=0.0, vmax=1.0, cmap=cmap, interpolation="nearest")
    for row in range(len(CLASS_KEYS)):
        for column in range(len(CLASS_KEYS)):
            value = normalized[row, column]
            ax.text(
                column,
                row,
                f"{100.0 * value:.0f}%\n({matrix[row, column]})",
                ha="center",
                va="center",
                fontsize=6.5,
                color="white" if value >= 0.55 else "0.16",
            )
    ax.set_xticks(np.arange(len(CLASS_KEYS)), CLASS_LABELS)
    ax.set_yticks(
        np.arange(len(CLASS_KEYS)),
        [f"{label.replace(chr(10), ' ')}  (n={n})" for label, n in zip(CLASS_LABELS, support)],
    )
    ax.tick_params(axis="x", rotation=0, pad=2)
    ax.set_xlabel("Predicted human morphology")
    ax.set_ylabel("Human morphology")
    ax.set_title("Locked-test confusion")
    ax.grid(False)
    ax.set_box_aspect(1)
    _close_panel(ax)


def _calibration_bins(frame: pd.DataFrame, n_bins: int = 10) -> tuple[pd.DataFrame, float]:
    probability = frame.loc[:, [f"p_{label}" for label in CLASS_KEYS]].to_numpy(float)
    truth = pd.to_numeric(frame["morphology_target_index"], errors="raise").to_numpy(int)
    active = (truth >= 0) & (truth < len(CLASS_KEYS))
    probability = probability[active]
    truth = truth[active]
    if not len(truth):
        raise ValueError("fixed test contains no active four-class morphology targets")
    confidence = probability.max(axis=1)
    correct = probability.argmax(axis=1) == truth
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    rows: list[dict[str, float | int]] = []
    ece = 0.0
    for low, high in zip(edges[:-1], edges[1:]):
        mask = (confidence >= low) & (
            confidence < high if high < 1.0 else confidence <= high
        )
        n = int(mask.sum())
        if not n:
            continue
        accuracy = float(np.mean(correct[mask]))
        mean_confidence = float(np.mean(confidence[mask]))
        ece += n / len(frame) * abs(accuracy - mean_confidence)
        rows.append(
            {
                "confidence": mean_confidence,
                "accuracy": accuracy,
                "n": n,
                "low": float(low),
                "high": float(high),
            }
        )
    return pd.DataFrame(rows), float(ece)


def _draw_calibration(
    ax: plt.Axes,
    frame: pd.DataFrame,
    color: tuple[float, ...],
) -> float:
    bins, ece = _calibration_bins(frame)
    ax.plot([0.0, 1.0], [0.0, 1.0], color="0.35", lw=1.0, ls="--", label="Ideal")
    if not bins.empty:
        sizes = 18.0 + 54.0 * np.sqrt(bins["n"] / bins["n"].max())
        ax.plot(
            bins["confidence"],
            bins["accuracy"],
            color=color,
            lw=1.5,
            zorder=2,
        )
        ax.scatter(
            bins["confidence"],
            bins["accuracy"],
            s=sizes,
            facecolor=color,
            edgecolor="white",
            linewidth=0.7,
            zorder=3,
            label="Test bins",
        )
    ax.text(
        0.04,
        0.95,
        f"ECE = {ece:.3f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7,
        bbox={"facecolor": "white", "edgecolor": "0.75", "pad": 2.0},
    )
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xticks(np.linspace(0.0, 1.0, 6))
    ax.set_yticks(np.linspace(0.0, 1.0, 6))
    ax.set_xlabel("Mean predicted confidence")
    ax.set_ylabel("Observed accuracy")
    ax.set_title("Locked-test reliability")
    ax.set_aspect("equal", adjustable="box")
    ax.set_box_aspect(1)
    legend = ax.legend(loc="lower right", frameon=True)
    legend.get_frame().set_edgecolor("black")
    _close_panel(ax)
    return ece


def _load_histories(root: Path, selected_profile: str) -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    for fold in range(5):
        path = root / selected_profile / f"fold_{fold}" / "history.csv"
        if not path.exists():
            raise FileNotFoundError(path)
        frame = pd.read_csv(path)
        required = {
            "epoch",
            "train_loss",
            "validation_morphology_loss",
            "validation_macro_f1",
            "validation_balanced_accuracy",
        }
        missing = sorted(required - set(frame.columns))
        if missing:
            raise ValueError(f"{path} is missing history columns: {missing}")
        frame["fold"] = fold
        parts.append(frame)
    return pd.concat(parts, ignore_index=True)


def _draw_fold_summary(
    ax: plt.Axes,
    histories: pd.DataFrame,
    column: str,
    *,
    color: tuple[float, ...],
    label: str,
    linestyle: str = "-",
) -> None:
    for _, fold in histories.groupby("fold", sort=True):
        ax.plot(
            fold["epoch"],
            fold[column],
            color=color,
            lw=0.65,
            alpha=0.18,
            linestyle=linestyle,
        )
    summary = (
        histories.groupby("epoch", as_index=False)[column]
        .agg(
            median="median",
            lower=lambda values: values.quantile(0.16),
            upper=lambda values: values.quantile(0.84),
        )
        .sort_values("epoch")
    )
    ax.fill_between(
        summary["epoch"].to_numpy(float),
        summary["lower"].to_numpy(float),
        summary["upper"].to_numpy(float),
        color=color,
        alpha=0.12,
        linewidth=0,
    )
    ax.plot(
        summary["epoch"],
        summary["median"],
        color=color,
        lw=1.7,
        linestyle=linestyle,
        label=label,
    )


def _draw_loss(
    ax: plt.Axes,
    histories: pd.DataFrame,
    train_color: tuple[float, ...],
    validation_color: tuple[float, ...],
) -> None:
    _draw_fold_summary(
        ax,
        histories,
        "train_loss",
        color=train_color,
        label="Train multitask loss",
    )
    _draw_fold_summary(
        ax,
        histories,
        "validation_morphology_loss",
        color=validation_color,
        label="Validation morphology loss",
    )
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Loss")
    ax.set_title("Optimization history (five folds)")
    ax.set_xlim(left=1)
    legend = ax.legend(loc="best", frameon=True)
    legend.get_frame().set_edgecolor("black")
    _close_panel(ax)


def _draw_scores(
    ax: plt.Axes,
    histories: pd.DataFrame,
    summary: dict,
    f1_color: tuple[float, ...],
    accuracy_color: tuple[float, ...],
) -> tuple[float, float]:
    _draw_fold_summary(
        ax,
        histories,
        "validation_macro_f1",
        color=f1_color,
        label="Validation macro F1",
    )
    _draw_fold_summary(
        ax,
        histories,
        "validation_balanced_accuracy",
        color=accuracy_color,
        label="Validation balanced accuracy",
    )
    all_test = summary["test_metrics"]["morphology_by_source"]["all"]
    test_f1 = float(all_test["macro_f1"])
    test_accuracy = float(all_test["balanced_accuracy"])
    ax.axhline(
        test_f1,
        color=f1_color,
        lw=1.1,
        ls=":",
        label=f"Locked test macro F1 = {test_f1:.2f}",
    )
    ax.axhline(
        test_accuracy,
        color=accuracy_color,
        lw=1.1,
        ls=":",
        label=f"Locked test balanced accuracy = {test_accuracy:.2f}",
    )
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Score")
    ax.set_title("Development-fold performance")
    ax.set_xlim(left=1)
    ax.set_ylim(0.0, 1.02)
    legend = ax.legend(loc="lower right", frameon=True, ncol=1)
    legend.get_frame().set_edgecolor("black")
    _close_panel(ax)
    return test_f1, test_accuracy


def _humanize_profile(profile: str) -> str:
    return PROFILE_LABELS.get(profile, profile.replace("_", " ").title())


def _validate_inputs(root: Path, summary: dict, frame: pd.DataFrame) -> str:
    selected = str(summary.get("selected_profile", ""))
    if not selected:
        raise ValueError("summary.json has no selected_profile")
    if frame.empty:
        raise ValueError("fixed_test_predictions.csv is empty")
    expected = {
        "morphology_target_index",
        "morphology_prediction_index",
        "is_injected_row",
        *[f"p_{label}" for label in CLASS_KEYS],
    }
    missing = sorted(expected - set(frame.columns))
    if missing:
        raise ValueError(f"fixed test predictions are missing columns: {missing}")
    for fold in range(5):
        checkpoint = root / selected / f"fold_{fold}" / "teacher.pt"
        if not checkpoint.exists():
            raise FileNotFoundError(checkpoint)
    return selected


def _as_bool(values: Iterable[object]) -> np.ndarray:
    return np.asarray(
        [str(value).strip().lower() in {"1", "true", "t", "yes", "y"} for value in values],
        dtype=bool,
    )


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    out_dir = (args.out_dir or (root / "performance_figure")).resolve()
    summary = _load_json(root / "summary.json")
    predictions = pd.read_csv(root / "fixed_test_predictions.csv", low_memory=False)
    selected_profile = _validate_inputs(root, summary, predictions)
    histories = _load_histories(root, selected_profile)
    morphology_target = pd.to_numeric(
        predictions["morphology_target_index"], errors="raise"
    )
    morphology_predictions = predictions.loc[
        morphology_target.between(0, len(CLASS_KEYS) - 1)
    ].copy()
    if morphology_predictions.empty:
        raise ValueError("fixed test contains no active four-class morphology rows")

    template = apply_twirl_style("full_page")
    palette = get_ordered_palette(5, palette="Dark2")
    confusion_color = palette[0]
    calibration_color = palette[1]
    train_color = palette[2]
    f1_color = palette[3]
    accuracy_color = palette[4]

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(template["figsize"][0], 6.25),
        gridspec_kw={"height_ratios": (1.0, 0.86)},
    )
    _draw_confusion(axes[0, 0], morphology_predictions, confusion_color)
    ece = _draw_calibration(axes[0, 1], morphology_predictions, calibration_color)
    _draw_loss(axes[1, 0], histories, train_color, calibration_color)
    test_f1, test_accuracy = _draw_scores(
        axes[1, 1],
        histories,
        summary,
        f1_color,
        accuracy_color,
    )
    for label, ax in zip(("a", "b", "c", "d"), axes.flat):
        _panel_label(ax, label)

    injected = _as_bool(morphology_predictions["is_injected_row"])
    planet = morphology_predictions["morphology_target_index"].eq(0).to_numpy()
    n_planet_real = int(np.count_nonzero(planet & ~injected))
    n_planet_injected = int(np.count_nonzero(planet & injected))
    footer = (
        f"Selected on development folds: {_humanize_profile(selected_profile)}.  "
        f"Locked test: n={len(predictions)} total, {len(morphology_predictions)} morphology "
        f"({np.count_nonzero(~injected)} real, {np.count_nonzero(injected)} injected); "
        "opened once after selection.\n"
        f"Planet-like support: {n_planet_real} real + {n_planet_injected} injected; "
        "all targets are human labels, with injection truth kept audit-only."
    )
    fig.text(
        0.5,
        0.012,
        footer,
        ha="center",
        va="bottom",
        fontsize=6.3,
        color="0.25",
        linespacing=1.35,
    )
    fig.subplots_adjust(left=0.125, right=0.985, top=0.955, bottom=0.12, wspace=0.34, hspace=0.36)

    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / f"{args.output_stem}.png"
    pdf_path = out_dir / f"{args.output_stem}.pdf"
    fig.savefig(png_path, dpi=args.png_dpi, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    figure_summary = {
        "selected_profile": selected_profile,
        "n_fixed_test_total": int(len(predictions)),
        "n_fixed_test_morphology": int(len(morphology_predictions)),
        "n_fixed_test_morphology_real": int(np.count_nonzero(~injected)),
        "n_fixed_test_morphology_injected": int(np.count_nonzero(injected)),
        "n_fixed_test_planet_real": n_planet_real,
        "n_fixed_test_planet_injected": n_planet_injected,
        "locked_test_macro_f1": test_f1,
        "locked_test_balanced_accuracy": test_accuracy,
        "locked_test_ece_recomputed": ece,
        "inputs": {
            "summary": str(root / "summary.json"),
            "predictions": str(root / "fixed_test_predictions.csv"),
            "selected_profile_histories": str(root / selected_profile / "fold_*/history.csv"),
        },
        "outputs": {"png": str(png_path), "pdf": str(pdf_path)},
    }
    (out_dir / "figure_summary.json").write_text(
        json.dumps(figure_summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(figure_summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
