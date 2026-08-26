#!/usr/bin/env python3
"""Plot the development-only FM0.2.1 step-500 to step-1000 trajectory."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.plotting.style import apply_twirl_style, get_ordered_palette  # noqa: E402


EXPECTED_TRAINING_SHA = "ddf442aafb8f62966e549e2287abad3474dd556a"
EXPECTED_EVALUATOR_SHA = "83816d07975eebe3825d76dfe7096d22b70376f5"
HISTORY_FIELDS = [
    "step",
    "total",
    "reconstruction_first",
    "reconstruction_second",
    "reconstruction_mean",
    "reconstruction",
    "vicreg_invariance",
    "vicreg_variance",
    "vicreg_covariance",
    "vicreg",
    "vicreg_weighted",
    "embedding_projection_gradient_norm",
    "learning_rate",
    "window_draws_seen",
    "masked_views_seen",
]


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step500-health", type=Path, required=True)
    parser.add_argument("--step1000-health", type=Path, required=True)
    parser.add_argument("--step1000-receipt", type=Path, required=True)
    parser.add_argument("--loss-history", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"JSON root is not a mapping: {path}")
    return payload


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _interval(payload: dict[str, Any], key: str) -> tuple[float, float, float]:
    interval = payload[key]
    return tuple(float(interval[name]) for name in ("mean", "lower", "upper"))


def _validate_inputs(
    step500: dict[str, Any],
    step1000: dict[str, Any],
    receipt: dict[str, Any],
    history: dict[str, Any],
    step1000_path: Path,
) -> list[dict[str, float]]:
    for expected_step, health in ((500, step500), (1000, step1000)):
        if health.get("schema_version") != "twirl_fm0_representation_health_v2":
            raise ValueError(f"step {expected_step} representation schema differs")
        if health.get("passed") is not True:
            raise ValueError(f"step {expected_step} evaluation did not pass")
        if int(health["run"]["global_step"]) != expected_step:
            raise ValueError(f"step {expected_step} milestone differs")
        if health["run"]["run_git_sha"] != EXPECTED_TRAINING_SHA:
            raise ValueError(f"step {expected_step} training revision differs")
        if health.get("evaluator_git_sha") != EXPECTED_EVALUATOR_SHA:
            raise ValueError(f"step {expected_step} evaluator revision differs")
        if health["data_access"]["sealed_test_access_count"] != 0:
            raise ValueError(f"step {expected_step} sealed-test access is nonzero")

    population_keys = (
        "selected_leakage_components_sha256",
        "selected_observation_keys_sha256",
    )
    for key in population_keys:
        if step500["evaluation_population"][key] != step1000["evaluation_population"][key]:
            raise ValueError(f"development population differs: {key}")
    authority500 = step500["evaluation_population"]["observation_sector_authority"]
    authority1000 = step1000["evaluation_population"]["observation_sector_authority"]
    if authority500["sha256"] != authority1000["sha256"]:
        raise ValueError("observation-sector authority differs")

    if receipt.get("passed") is not True or receipt.get("stop_after_step") != 1000:
        raise ValueError("step-1000 receipt did not pass")
    if receipt.get("sealed_test_access_count") != 0:
        raise ValueError("step-1000 receipt records sealed-test access")
    if receipt.get("scientific_go_no_go_applied") is not False:
        raise ValueError("formal scientific gate was unexpectedly applied")
    if receipt.get("artifact_git_sha") != EXPECTED_TRAINING_SHA:
        raise ValueError("step-1000 receipt training revision differs")
    if receipt.get("evaluator_git_sha") != EXPECTED_EVALUATOR_SHA:
        raise ValueError("step-1000 receipt evaluator revision differs")
    if receipt.get("evaluation_sha256") != _sha256(step1000_path):
        raise ValueError("step-1000 evaluation hash differs from its receipt")
    if receipt.get("checkpoint_sha256") != history.get("checkpoint_sha256"):
        raise ValueError("loss history is not bound to the evaluated checkpoint")
    if history.get("schema_version") != "twirl_fm0_2_loss_history_extract_v1":
        raise ValueError("loss-history extract schema differs")
    if int(history.get("global_step", -1)) != 1000:
        raise ValueError("loss history does not end at step 1000")

    rows = history.get("loss_history")
    if not isinstance(rows, list) or len(rows) != 1000:
        raise ValueError("loss history must contain exactly 1000 rows")
    expected_steps = list(range(1, 1001))
    actual_steps = [int(row["step"]) for row in rows]
    if actual_steps != expected_steps:
        raise ValueError("loss-history steps are not contiguous 1--1000")
    for row in rows:
        if not set(HISTORY_FIELDS).issubset(row):
            raise ValueError("loss-history fields differ")
        if not all(math.isfinite(float(row[name])) for name in HISTORY_FIELDS):
            raise ValueError("loss history contains nonfinite values")
    return [{name: float(row[name]) for name in HISTORY_FIELDS} for row in rows]


def _write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _panel_label(ax: Any, label: str, *, right: bool = False) -> None:
    ax.text(
        0.98 if right else 0.02,
        0.96,
        label,
        transform=ax.transAxes,
        ha="right" if right else "left",
        va="top",
        fontweight="bold",
    )


def _plot_interval_trajectory(
    ax: Any,
    steps: np.ndarray,
    intervals: list[tuple[float, float, float]],
    *,
    color: Any,
    marker: str,
    label: str,
) -> None:
    means = np.asarray([value[0] for value in intervals])
    lower = np.asarray([value[1] for value in intervals])
    upper = np.asarray([value[2] for value in intervals])
    ax.errorbar(
        steps,
        means,
        yerr=np.vstack((means - lower, upper - means)),
        color=color,
        marker=marker,
        ms=4.0,
        capsize=2.5,
        lw=1.15,
        label=label,
    )


def _representation_rows(health_by_step: dict[int, dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for step, health in health_by_step.items():
        trained = health["trained_encoder"]["representations"]
        deltas = health["source_paired_trained_minus_random"]
        masked = health["masked_reconstruction"]
        projection = health["trained_encoder"]["projection_spectrum"]
        for branch in ("h_window", "z_window"):
            representation = trained[branch]
            paired = _interval(
                representation["safe_mask_pair_separation"],
                "paired_minus_unrelated_source_clustered_95_interval",
            )
            delta = _interval(
                deltas[branch], "trained_minus_random_source_clustered_95_interval"
            )
            retrieval = _interval(
                representation["query_sector_excluded_cross_visit_retrieval"],
                "top1_source_clustered_95_interval",
            )
            rows.append(
                {
                    "step": step,
                    "representation": branch,
                    "effective_rank": representation["embedding_health"]["effective_rank"],
                    "zero_or_constant_dimensions": representation["embedding_health"]["zero_or_constant_dimensions"],
                    "paired_minus_unrelated_mean": paired[0],
                    "paired_minus_unrelated_lower": paired[1],
                    "paired_minus_unrelated_upper": paired[2],
                    "trained_minus_random_mean": delta[0],
                    "trained_minus_random_lower": delta[1],
                    "trained_minus_random_upper": delta[2],
                    "cross_sector_retrieval_mean": retrieval[0],
                    "cross_sector_retrieval_lower": retrieval[1],
                    "cross_sector_retrieval_upper": retrieval[2],
                    "masked_huber_to_zero_ratio": masked["model_to_zero_baseline_ratio"],
                    "projection_effective_rank": projection["effective_rank"],
                    "projection_condition_number": projection["condition_number_on_numerical_support"],
                }
            )
    return rows


def _plot_representation(
    step500: dict[str, Any], step1000: dict[str, Any], output_dir: Path
) -> None:
    healths = [step500, step1000]
    steps = np.asarray([500.0, 1000.0])
    colors = get_ordered_palette(8)
    apply_twirl_style("full_page")
    fig, axes = plt.subplots(3, 2, figsize=(7.1, 7.5))
    ax_rank, ax_paired, ax_delta, ax_recon, ax_retrieval, ax_spectrum = axes.flat

    for branch, color, marker, label in (
        ("h_window", colors[5], "o", r"trained $h$"),
        ("z_window", colors[7], "s", r"trained $z$"),
    ):
        values = [
            float(h["trained_encoder"]["representations"][branch]["embedding_health"]["effective_rank"])
            for h in healths
        ]
        ax_rank.plot(steps, values, color=color, marker=marker, ms=4, lw=1.2, label=label)
    controls = step1000["random_encoder_control"]["representations"]
    pca_rank = float(
        step1000["baselines"]["pca"]["development_evaluation"]["embedding_health"]["effective_rank"]
    )
    ax_rank.axhline(
        float(controls["z_window"]["embedding_health"]["effective_rank"]),
        color="0.55", ls=":", lw=1.0, label="matched random $z$",
    )
    ax_rank.axhline(pca_rank, color="0.3", ls="--", lw=1.0, label="train-fit PCA")
    ax_rank.set(xlabel="Optimizer step", ylabel="Effective rank", xticks=steps, xlim=(440, 1060))
    ax_rank.legend(loc="upper left", ncol=2)
    ax_rank.text(
        0.98, 0.06, r"constant dims: $z=0$, $h=1$ at both milestones",
        transform=ax_rank.transAxes, ha="right", va="bottom",
    )
    _panel_label(ax_rank, "a", right=True)

    for branch, color, marker, label in (
        ("h_window", colors[5], "o", r"trained $h$"),
        ("z_window", colors[7], "s", r"trained $z$"),
    ):
        intervals = [
            _interval(
                h["trained_encoder"]["representations"][branch]["safe_mask_pair_separation"],
                "paired_minus_unrelated_source_clustered_95_interval",
            )
            for h in healths
        ]
        _plot_interval_trajectory(
            ax_paired, steps, intervals, color=color, marker=marker, label=label
        )
    random_z = _interval(
        controls["z_window"]["safe_mask_pair_separation"],
        "paired_minus_unrelated_source_clustered_95_interval",
    )[0]
    scalar = float(
        step1000["baselines"]["robust_scalar"]["development_evaluation"]
        ["safe_mask_pair_separation"]["paired_minus_unrelated_mean"]
    )
    pca = float(
        step1000["baselines"]["pca"]["development_evaluation"]
        ["safe_mask_pair_separation"]["paired_minus_unrelated_mean"]
    )
    ax_paired.axhline(random_z, color="0.55", ls=":", lw=1.0, label="matched random $z$")
    ax_paired.axhline(scalar, color=colors[2], ls="-.", lw=1.0, label="robust scalar")
    ax_paired.axhline(pca, color="0.3", ls="--", lw=1.0, label="train-fit PCA")
    ax_paired.set(
        xlabel="Optimizer step",
        ylabel="Paired − unrelated cosine",
        xticks=steps,
        xlim=(440, 1060),
        yscale="log",
    )
    ax_paired.legend(loc="upper left", ncol=2)
    _panel_label(ax_paired, "b", right=True)

    for branch, color, marker, label in (
        ("h_window", colors[5], "o", r"$h$"),
        ("z_window", colors[7], "s", r"$z$"),
    ):
        intervals = [
            _interval(
                h["source_paired_trained_minus_random"][branch],
                "trained_minus_random_source_clustered_95_interval",
            )
            for h in healths
        ]
        _plot_interval_trajectory(
            ax_delta, steps, intervals, color=color, marker=marker, label=label
        )
    ax_delta.axhline(0.0, color="0.25", ls="--", lw=0.9)
    ax_delta.set(
        xlabel="Optimizer step",
        ylabel="Trained − random separation",
        xticks=steps,
        xlim=(440, 1060),
    )
    ax_delta.legend(loc="upper left")
    _panel_label(ax_delta, "c", right=True)

    excess = np.asarray(
        [100.0 * (float(h["masked_reconstruction"]["model_to_zero_baseline_ratio"]) - 1.0) for h in healths]
    )
    ax_recon.plot(steps, excess, color=colors[6], marker="o", ms=4, lw=1.2)
    ax_recon.axhline(0.0, color="0.25", ls="--", lw=0.9, label="zero predictor")
    ax_recon.set(
        xlabel="Optimizer step",
        ylabel="Masked Huber excess vs zero (%)",
        xticks=steps,
        xlim=(440, 1060),
    )
    ax_recon.annotate(
        f"{excess[-1]:.2f}% worse",
        xy=(steps[-1], excess[-1]),
        xytext=(-70, -20),
        textcoords="offset points",
        arrowprops={"arrowstyle": "-", "color": "0.35", "lw": 0.7},
    )
    ax_recon.legend(loc="lower left")
    _panel_label(ax_recon, "d")

    for branch, color, marker, label in (
        ("h_window", colors[5], "o", r"trained $h$"),
        ("z_window", colors[7], "s", r"trained $z$"),
    ):
        intervals = [
            _interval(
                h["trained_encoder"]["representations"][branch]
                ["query_sector_excluded_cross_visit_retrieval"],
                "top1_source_clustered_95_interval",
            )
            for h in healths
        ]
        _plot_interval_trajectory(
            ax_retrieval, steps, intervals, color=color, marker=marker, label=label
        )
    random_retrieval = _interval(
        controls["z_window"]["query_sector_excluded_cross_visit_retrieval"],
        "top1_source_clustered_95_interval",
    )[0]
    pca_retrieval = _interval(
        step1000["baselines"]["pca"]["development_evaluation"]
        ["query_sector_excluded_cross_visit_retrieval"],
        "top1_source_clustered_95_interval",
    )[0]
    ax_retrieval.axhline(random_retrieval, color="0.55", ls=":", lw=1.0, label="matched random $z$")
    ax_retrieval.axhline(pca_retrieval, color="0.3", ls="--", lw=1.0, label="train-fit PCA")
    ax_retrieval.set(
        xlabel="Optimizer step",
        ylabel="Cross-sector top-1 retrieval",
        xticks=steps,
        xlim=(440, 1060),
        ylim=(0.0, 0.075),
    )
    ax_retrieval.legend(loc="upper left", ncol=2)
    _panel_label(ax_retrieval, "e", right=True)

    indices = np.arange(1, 257)
    for health, color, label in (
        (step500, colors[4], "trained step 500"),
        (step1000, colors[7], "trained step 1000"),
    ):
        singular_values = np.asarray(
            health["trained_encoder"]["projection_spectrum"]["singular_values"],
            dtype=float,
        )
        ax_spectrum.plot(indices, singular_values, color=color, lw=1.05, label=label)
    random_singular = np.asarray(
        step1000["random_encoder_control"]["projection_spectrum"]["singular_values"],
        dtype=float,
    )
    ax_spectrum.plot(indices, random_singular, color="0.55", ls=":", lw=1.0, label="matched random")
    ax_spectrum.set(
        xlabel="Singular-value index",
        ylabel="Projection singular value",
        yscale="log",
        xlim=(1, 256),
    )
    ax_spectrum.legend(loc="upper right")
    _panel_label(ax_spectrum, "f")

    fig.subplots_adjust(left=0.10, right=0.985, bottom=0.075, top=0.985, hspace=0.34, wspace=0.30)
    stem = output_dir / "fm0_2_step1000_representation_trajectory"
    fig.savefig(stem.with_suffix(".png"), dpi=240)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


def _positive_scale(ax: Any, values: list[np.ndarray]) -> None:
    if all(np.all(value > 0.0) for value in values):
        ax.set_yscale("log")
    else:
        positive = np.concatenate([value[value > 0.0] for value in values])
        linthresh = max(float(np.min(positive)) * 0.5, 1.0e-12)
        ax.set_yscale("symlog", linthresh=linthresh)


def _plot_training(rows: list[dict[str, float]], output_dir: Path) -> None:
    colors = get_ordered_palette(8)
    apply_twirl_style("full_page")
    fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.7))
    ax_objective, ax_reconstruction, ax_vicreg, ax_optimization = axes.flat
    values = {name: np.asarray([row[name] for row in rows]) for name in HISTORY_FIELDS}
    steps = values["step"]

    ax_objective.plot(steps, values["total"], color=colors[7], lw=0.95, label="total")
    ax_objective.plot(steps, values["reconstruction_mean"], color=colors[4], lw=0.9, label="reconstruction mean")
    ax_objective.plot(steps, values["vicreg_weighted"], color=colors[1], lw=0.9, label="weighted VICReg")
    ax_objective.axvline(500, color="0.55", ls=":", lw=0.9)
    ax_objective.set(xlabel="Optimizer step", ylabel="Objective contribution", xlim=(1, 1000))
    ax_objective.legend(loc="upper right")
    _panel_label(ax_objective, "a")

    for name, color, linestyle, label in (
        ("reconstruction_first", colors[2], "-", "first view"),
        ("reconstruction_second", colors[6], "--", "second view"),
        ("reconstruction_mean", colors[4], ":", "mean"),
    ):
        ax_reconstruction.plot(steps, values[name], color=color, ls=linestyle, lw=0.9, label=label)
    ax_reconstruction.axvline(500, color="0.55", ls=":", lw=0.9)
    ax_reconstruction.set(xlabel="Optimizer step", ylabel="Reconstruction loss", xlim=(1, 1000))
    ax_reconstruction.legend(loc="upper right")
    _panel_label(ax_reconstruction, "b")

    vicreg_arrays: list[np.ndarray] = []
    for name, color, linestyle, label in (
        ("vicreg_invariance", colors[2], "-", "invariance"),
        ("vicreg_variance", colors[5], "--", "variance"),
        ("vicreg_covariance", colors[7], ":", "covariance"),
    ):
        array = values[name]
        vicreg_arrays.append(array)
        ax_vicreg.plot(steps, array, color=color, ls=linestyle, lw=0.9, label=label)
    _positive_scale(ax_vicreg, vicreg_arrays)
    ax_vicreg.axvline(500, color="0.55", ls=":", lw=0.9)
    ax_vicreg.set(xlabel="Optimizer step", ylabel="Raw VICReg component", xlim=(1, 1000))
    ax_vicreg.legend(loc="upper right")
    _panel_label(ax_vicreg, "c")

    gradient = values["embedding_projection_gradient_norm"]
    line_gradient = ax_optimization.plot(
        steps, gradient, color=colors[6], lw=0.9, label="projection-gradient norm"
    )[0]
    _positive_scale(ax_optimization, [gradient])
    ax_optimization.axvline(500, color="0.55", ls=":", lw=0.9)
    ax_optimization.set(
        xlabel="Optimizer step",
        ylabel="Projection-gradient norm",
        xlim=(1, 1000),
    )
    ax_lr = ax_optimization.twinx()
    line_lr = ax_lr.plot(
        steps, values["learning_rate"], color=colors[1], lw=1.0, label="learning rate"
    )[0]
    ax_lr.set_ylabel("Learning rate")
    ax_lr.grid(False)
    ax_optimization.legend(
        [line_gradient, line_lr],
        [line_gradient.get_label(), line_lr.get_label()],
        loc="center right",
    )
    _panel_label(ax_optimization, "d")

    fig.subplots_adjust(left=0.10, right=0.86, bottom=0.10, top=0.98, hspace=0.32, wspace=0.34)
    stem = output_dir / "fm0_2_step1000_training_components"
    fig.savefig(stem.with_suffix(".png"), dpi=240)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    args = _parser().parse_args()
    step500 = _load_json(args.step500_health)
    step1000 = _load_json(args.step1000_health)
    receipt = _load_json(args.step1000_receipt)
    history = _load_json(args.loss_history)
    rows = _validate_inputs(step500, step1000, receipt, history, args.step1000_health)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(args.output_dir / "training_trace.csv", rows, HISTORY_FIELDS)
    metric_rows = _representation_rows({500: step500, 1000: step1000})
    _write_csv(
        args.output_dir / "step500_step1000_metrics.csv",
        metric_rows,
        list(metric_rows[0]),
    )
    _plot_representation(step500, step1000, args.output_dir)
    _plot_training(rows, args.output_dir)
    print(args.output_dir / "fm0_2_step1000_representation_trajectory.png")
    print(args.output_dir / "fm0_2_step1000_training_components.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
