#!/usr/bin/env python3
"""Plot and tabulate the development-only FM0.2.1 trajectory through step 2000."""
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
EXPECTED_CHECKPOINT_SHA = "976b5053c857c38b9fbf7a35c9d0605f0023318b2c9bd37a88995992e5aa7bd2"
EXPECTED_COMPONENT_HASH = "50e0b59e1bc80ded61e30c1a6c07a56eee6306b03c3bd52e1f30ebdc43f5c72f"
EXPECTED_OBSERVATION_HASH = "9ffdc49e042cf4e8dccbe53e70fd849723abf96625c36a00adea8f4fde20546a"
EXPECTED_AUTHORITY_HASH = "bdd3e8039b312aeb21557662226a8a11b761dd1da031742be42ee9b0d1c6edc5"
Z_RANK_MINIMUM = 26.0
MASKED_HUBER_MAXIMUM = 0.004864579067
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
    for step in (500, 1000, 2000):
        parser.add_argument(f"--step{step}-health", type=Path, required=True)
    parser.add_argument("--step2000-receipt", type=Path, required=True)
    parser.add_argument("--step2000-post-validation", type=Path, required=True)
    parser.add_argument("--step1000-loss-history", type=Path, required=True)
    parser.add_argument("--step2000-loss-history", type=Path, required=True)
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


def _close(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=2.0e-6, abs_tol=2.0e-9)


def _validate_health(step: int, health: dict[str, Any]) -> None:
    if health.get("schema_version") != "twirl_fm0_representation_health_v2":
        raise ValueError(f"step {step} representation schema differs")
    if health.get("passed") is not True:
        raise ValueError(f"step {step} evaluator did not pass")
    if int(health["run"]["global_step"]) != step:
        raise ValueError(f"step {step} milestone differs")
    if health["run"]["run_git_sha"] != EXPECTED_TRAINING_SHA:
        raise ValueError(f"step {step} training revision differs")
    if health.get("evaluator_git_sha") != EXPECTED_EVALUATOR_SHA:
        raise ValueError(f"step {step} evaluator revision differs")
    if int(health["data_access"]["sealed_test_access_count"]) != 0:
        raise ValueError(f"step {step} sealed-test access is nonzero")
    population = health["evaluation_population"]
    if population["selected_leakage_components_sha256"] != EXPECTED_COMPONENT_HASH:
        raise ValueError(f"step {step} development component population differs")
    if population["selected_observation_keys_sha256"] != EXPECTED_OBSERVATION_HASH:
        raise ValueError(f"step {step} development observation population differs")
    if population["observation_sector_authority"]["sha256"] != EXPECTED_AUTHORITY_HASH:
        raise ValueError(f"step {step} observation-sector authority differs")


def _validated_history(payload: dict[str, Any], expected_step: int) -> list[dict[str, float]]:
    if payload.get("schema_version") != "twirl_fm0_2_loss_history_extract_v1":
        raise ValueError("loss-history extract schema differs")
    if int(payload.get("global_step", -1)) != expected_step:
        raise ValueError("loss-history global step differs")
    rows = payload.get("loss_history")
    if not isinstance(rows, list) or len(rows) != expected_step:
        raise ValueError("loss-history row count differs")
    validated: list[dict[str, float]] = []
    for expected_row_step, raw in enumerate(rows, start=1):
        if not isinstance(raw, dict) or not set(HISTORY_FIELDS).issubset(raw):
            raise ValueError(f"loss-history fields differ at step {expected_row_step}")
        row = {field: float(raw[field]) for field in HISTORY_FIELDS}
        if int(row["step"]) != expected_row_step:
            raise ValueError("loss-history steps are not contiguous")
        if not all(math.isfinite(value) for value in row.values()):
            raise ValueError(f"loss history contains nonfinite values at step {expected_row_step}")
        if not _close(row["reconstruction"], row["reconstruction_first"]):
            raise ValueError("optimized reconstruction identity differs")
        if not _close(
            row["reconstruction_mean"],
            0.5 * (row["reconstruction_first"] + row["reconstruction_second"]),
        ):
            raise ValueError("reconstruction mean identity differs")
        if not _close(row["vicreg_weighted"], 0.0002 * row["vicreg"]):
            raise ValueError("weighted VICReg identity differs")
        if not _close(
            row["vicreg"],
            25.0 * row["vicreg_invariance"]
            + 25.0 * row["vicreg_variance"]
            + row["vicreg_covariance"],
        ):
            raise ValueError("VICReg component identity differs")
        if not _close(row["total"], row["reconstruction_first"] + row["vicreg_weighted"]):
            raise ValueError("total-loss identity differs")
        validated.append(row)
    return validated


def _validate_inputs(
    health_by_step: dict[int, dict[str, Any]],
    receipt: dict[str, Any],
    post_validation: dict[str, Any],
    step1000_history: dict[str, Any],
    step2000_history: dict[str, Any],
    step2000_health_path: Path,
) -> list[dict[str, float]]:
    for step, health in health_by_step.items():
        _validate_health(step, health)
    step2000 = health_by_step[2000]
    if step2000["run"]["checkpoint_sha256"] != EXPECTED_CHECKPOINT_SHA:
        raise ValueError("step-2000 checkpoint hash differs")
    if receipt.get("passed") is not True or receipt.get("stop_after_step") != 2000:
        raise ValueError("step-2000 evaluator receipt did not pass")
    if receipt.get("scientific_go_no_go_applied") is not False:
        raise ValueError("formal scientific gate was unexpectedly applied")
    if receipt.get("sealed_test_access_count") != 0:
        raise ValueError("step-2000 evaluator receipt records sealed access")
    if receipt.get("artifact_git_sha") != EXPECTED_TRAINING_SHA:
        raise ValueError("step-2000 receipt training revision differs")
    if receipt.get("evaluator_git_sha") != EXPECTED_EVALUATOR_SHA:
        raise ValueError("step-2000 receipt evaluator revision differs")
    if receipt.get("evaluation_sha256") != _sha256(step2000_health_path):
        raise ValueError("step-2000 evaluation hash differs from its receipt")
    if receipt.get("checkpoint_sha256") != EXPECTED_CHECKPOINT_SHA:
        raise ValueError("step-2000 evaluator checkpoint binding differs")
    if post_validation.get("passed") is not True or post_validation.get("stop_after_step") != 2000:
        raise ValueError("step-2000 post-validation did not pass")
    if post_validation.get("sealed_test_access_count") != 0:
        raise ValueError("step-2000 post-validation records sealed access")
    if post_validation.get("artifacts_sha256", {}).get("checkpoint_step_00002000.pt") != EXPECTED_CHECKPOINT_SHA:
        raise ValueError("post-validation checkpoint binding differs")

    old_rows = _validated_history(step1000_history, 1000)
    rows = _validated_history(step2000_history, 2000)
    if step2000_history.get("checkpoint_sha256") != EXPECTED_CHECKPOINT_SHA:
        raise ValueError("step-2000 loss history checkpoint binding differs")
    if old_rows != rows[:1000]:
        raise ValueError("step-2000 loss history does not preserve the step-1000 prefix")
    return rows


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


def _plot_interval(
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
    for step, health in sorted(health_by_step.items()):
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
            retrieval = representation["query_sector_excluded_cross_visit_retrieval"]
            retrieval_interval = _interval(retrieval, "top1_source_clustered_95_interval")
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
                    "cross_sector_retrieval_raw_top1": retrieval["top1_same_component_retrieval"],
                    "cross_sector_retrieval_clustered_mean": retrieval_interval[0],
                    "cross_sector_retrieval_clustered_lower": retrieval_interval[1],
                    "cross_sector_retrieval_clustered_upper": retrieval_interval[2],
                    "masked_huber": masked["masked_huber_mean"],
                    "zero_prediction_masked_huber": masked["zero_prediction_masked_huber_mean"],
                    "masked_huber_to_zero_ratio": masked["model_to_zero_baseline_ratio"],
                    "projection_numerical_rank": projection["numerical_rank"],
                    "projection_effective_rank": projection["effective_rank"],
                    "projection_condition_number": projection["condition_number_on_numerical_support"],
                    "projection_leading_energy_share": projection["leading_singular_energy_share"],
                    "projection_smallest_singular_value": projection["smallest_singular_value"],
                }
            )
    return rows


def _gate_rows(step2000: dict[str, Any]) -> list[dict[str, Any]]:
    z = step2000["trained_encoder"]["representations"]["z_window"]
    paired_lower = _interval(
        z["safe_mask_pair_separation"],
        "paired_minus_unrelated_source_clustered_95_interval",
    )[1]
    delta_lower = _interval(
        step2000["source_paired_trained_minus_random"]["z_window"],
        "trained_minus_random_source_clustered_95_interval",
    )[1]
    huber = float(step2000["masked_reconstruction"]["masked_huber_mean"])
    observed = [
        ("z_window_effective_rank", float(z["embedding_health"]["effective_rank"]), ">= 26", float(z["embedding_health"]["effective_rank"]) >= Z_RANK_MINIMUM),
        ("z_window_zero_or_constant_dimensions", int(z["embedding_health"]["zero_or_constant_dimensions"]), "<= 0", int(z["embedding_health"]["zero_or_constant_dimensions"]) <= 0),
        ("z_paired_minus_unrelated_95_lower", paired_lower, "> 0", paired_lower > 0.0),
        ("z_trained_minus_random_95_lower", delta_lower, "> 0", delta_lower > 0.0),
        ("development_masked_huber", huber, f"<= {MASKED_HUBER_MAXIMUM:.12f}", huber <= MASKED_HUBER_MAXIMUM),
        ("sealed_test_access_count", int(step2000["data_access"]["sealed_test_access_count"]), "== 0", int(step2000["data_access"]["sealed_test_access_count"]) == 0),
    ]
    return [
        {"criterion": name, "observed": value, "threshold": threshold, "meets_criterion": passed}
        for name, value, threshold, passed in observed
    ]


def _plot_representation(health_by_step: dict[int, dict[str, Any]], output_dir: Path) -> None:
    ordered_steps = sorted(health_by_step)
    healths = [health_by_step[step] for step in ordered_steps]
    steps = np.asarray(ordered_steps, dtype=float)
    final = health_by_step[2000]
    controls = final["random_encoder_control"]["representations"]
    colors = get_ordered_palette(8)
    apply_twirl_style("full_page")
    fig, axes = plt.subplots(3, 2, figsize=(7.1, 7.5))
    ax_rank, ax_paired, ax_delta, ax_recon, ax_retrieval, ax_spectrum = axes.flat

    for branch, color, marker, label in (
        ("h_window", colors[5], "o", r"trained $h$"),
        ("z_window", colors[7], "s", r"trained $z$"),
    ):
        values = [
            float(health["trained_encoder"]["representations"][branch]["embedding_health"]["effective_rank"])
            for health in healths
        ]
        ax_rank.plot(steps, values, color=color, marker=marker, ms=4, lw=1.2, label=label)
    pca_rank = float(final["baselines"]["pca"]["development_evaluation"]["embedding_health"]["effective_rank"])
    ax_rank.axhline(float(controls["z_window"]["embedding_health"]["effective_rank"]), color="0.55", ls=":", lw=1.0, label="matched random $z$")
    ax_rank.axhline(pca_rank, color="0.3", ls="--", lw=1.0, label="train-fit PCA")
    ax_rank.axhline(Z_RANK_MINIMUM, color=colors[2], ls="-.", lw=1.0, label="$z$ gate")
    ax_rank.set(xlabel="Optimizer step", ylabel="Effective rank", xticks=steps, xlim=(430, 2070))
    ax_rank.legend(loc="upper left", ncol=2)
    ax_rank.text(0.98, 0.06, r"constant dims at 2000: $z=0$, $h=1$", transform=ax_rank.transAxes, ha="right", va="bottom")
    _panel_label(ax_rank, "a", right=True)

    for branch, color, marker, label in (
        ("h_window", colors[5], "o", r"trained $h$"),
        ("z_window", colors[7], "s", r"trained $z$"),
    ):
        intervals = [
            _interval(
                health["trained_encoder"]["representations"][branch]["safe_mask_pair_separation"],
                "paired_minus_unrelated_source_clustered_95_interval",
            )
            for health in healths
        ]
        _plot_interval(ax_paired, steps, intervals, color=color, marker=marker, label=label)
    random_z = _interval(controls["z_window"]["safe_mask_pair_separation"], "paired_minus_unrelated_source_clustered_95_interval")[0]
    scalar = float(final["baselines"]["robust_scalar"]["development_evaluation"]["safe_mask_pair_separation"]["paired_minus_unrelated_mean"])
    pca = float(final["baselines"]["pca"]["development_evaluation"]["safe_mask_pair_separation"]["paired_minus_unrelated_mean"])
    ax_paired.axhline(random_z, color="0.55", ls=":", lw=1.0, label="matched random $z$")
    ax_paired.axhline(scalar, color=colors[2], ls="-.", lw=1.0, label="robust scalar")
    ax_paired.axhline(pca, color="0.3", ls="--", lw=1.0, label="train-fit PCA")
    ax_paired.set(xlabel="Optimizer step", ylabel="Paired − unrelated cosine", xticks=steps, xlim=(430, 2070), yscale="log")
    ax_paired.legend(loc="upper left", ncol=2)
    _panel_label(ax_paired, "b", right=True)

    for branch, color, marker, label in (
        ("h_window", colors[5], "o", r"$h$"),
        ("z_window", colors[7], "s", r"$z$"),
    ):
        intervals = [
            _interval(health["source_paired_trained_minus_random"][branch], "trained_minus_random_source_clustered_95_interval")
            for health in healths
        ]
        _plot_interval(ax_delta, steps, intervals, color=color, marker=marker, label=label)
    ax_delta.axhline(0.0, color="0.25", ls="--", lw=0.9)
    ax_delta.set(xlabel="Optimizer step", ylabel="Trained − random separation", xticks=steps, xlim=(430, 2070))
    ax_delta.legend(loc="upper left")
    _panel_label(ax_delta, "c", right=True)

    excess = np.asarray([100.0 * (float(health["masked_reconstruction"]["model_to_zero_baseline_ratio"]) - 1.0) for health in healths])
    ax_recon.plot(steps, excess, color=colors[6], marker="o", ms=4, lw=1.2)
    ax_recon.axhline(0.0, color="0.25", ls="--", lw=0.9, label="zero predictor / gate")
    ax_recon.set(xlabel="Optimizer step", ylabel="Masked Huber excess vs zero (%)", xticks=steps, xlim=(430, 2070))
    ax_recon.annotate(f"{excess[-1]:.3f}% worse", xy=(steps[-1], excess[-1]), xytext=(-76, 18), textcoords="offset points", arrowprops={"arrowstyle": "-", "color": "0.35", "lw": 0.7})
    ax_recon.legend(loc="upper right")
    _panel_label(ax_recon, "d")

    for branch, color, marker, label in (
        ("h_window", colors[5], "o", r"trained $h$"),
        ("z_window", colors[7], "s", r"trained $z$"),
    ):
        intervals = [
            _interval(health["trained_encoder"]["representations"][branch]["query_sector_excluded_cross_visit_retrieval"], "top1_source_clustered_95_interval")
            for health in healths
        ]
        _plot_interval(ax_retrieval, steps, intervals, color=color, marker=marker, label=label)
    random_retrieval = _interval(controls["z_window"]["query_sector_excluded_cross_visit_retrieval"], "top1_source_clustered_95_interval")[0]
    scalar_retrieval = _interval(final["baselines"]["robust_scalar"]["development_evaluation"]["query_sector_excluded_cross_visit_retrieval"], "top1_source_clustered_95_interval")[0]
    pca_retrieval = _interval(final["baselines"]["pca"]["development_evaluation"]["query_sector_excluded_cross_visit_retrieval"], "top1_source_clustered_95_interval")[0]
    ax_retrieval.axhline(random_retrieval, color="0.55", ls=":", lw=1.0, label="matched random $z$")
    ax_retrieval.axhline(scalar_retrieval, color=colors[2], ls="-.", lw=1.0, label="robust scalar")
    ax_retrieval.axhline(pca_retrieval, color="0.3", ls="--", lw=1.0, label="train-fit PCA")
    ax_retrieval.set(xlabel="Optimizer step", ylabel="Cross-sector top-1 retrieval", xticks=steps, xlim=(430, 2070), ylim=(0.0, 0.10))
    ax_retrieval.legend(loc="upper left", ncol=2)
    _panel_label(ax_retrieval, "e", right=True)

    indices = np.arange(1, 257)
    for health, color, label in zip(healths, (colors[3], colors[5], colors[7]), ("step 500", "step 1000", "step 2000"), strict=True):
        singular = np.asarray(health["trained_encoder"]["projection_spectrum"]["singular_values"], dtype=float)
        ax_spectrum.plot(indices, singular, color=color, lw=1.0, label=label)
    random_singular = np.asarray(final["random_encoder_control"]["projection_spectrum"]["singular_values"], dtype=float)
    ax_spectrum.plot(indices, random_singular, color="0.55", ls=":", lw=1.0, label="matched random")
    ax_spectrum.set(xlabel="Singular-value index", ylabel="Projection singular value", yscale="log", xlim=(1, 256))
    ax_spectrum.legend(loc="upper right")
    _panel_label(ax_spectrum, "f")

    fig.subplots_adjust(left=0.10, right=0.985, bottom=0.075, top=0.985, hspace=0.34, wspace=0.30)
    stem = output_dir / "fm0_2_step2000_representation_trajectory"
    fig.savefig(stem.with_suffix(".png"), dpi=240)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


def _positive_scale(ax: Any, values: list[np.ndarray]) -> None:
    if all(np.all(value > 0.0) for value in values):
        ax.set_yscale("log")
        return
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

    ax_objective.plot(steps, values["total"], color=colors[7], lw=0.85, label="total")
    ax_objective.plot(steps, values["reconstruction_first"], color=colors[4], lw=0.8, label="optimized reconstruction")
    ax_objective.plot(steps, values["vicreg_weighted"], color=colors[1], lw=0.8, label="weighted VICReg")
    for milestone in (500, 1000):
        ax_objective.axvline(milestone, color="0.62", ls=":", lw=0.8)
    ax_objective.set(xlabel="Optimizer step", ylabel="Objective contribution", xlim=(1, 2000))
    ax_objective.legend(loc="upper right")
    _panel_label(ax_objective, "a")

    for name, color, linestyle, label in (
        ("reconstruction_first", colors[2], "-", "first / optimized"),
        ("reconstruction_second", colors[6], "--", "second / diagnostic"),
        ("reconstruction_mean", colors[4], ":", "mean / diagnostic"),
    ):
        ax_reconstruction.plot(steps, values[name], color=color, ls=linestyle, lw=0.8, label=label)
    for milestone in (500, 1000):
        ax_reconstruction.axvline(milestone, color="0.62", ls=":", lw=0.8)
    ax_reconstruction.set(xlabel="Optimizer step", ylabel="Reconstruction loss", xlim=(1, 2000))
    ax_reconstruction.legend(loc="upper right")
    _panel_label(ax_reconstruction, "b")

    vicreg_arrays: list[np.ndarray] = []
    for name, color, linestyle, label in (
        ("vicreg", "0.25", "-", "aggregate"),
        ("vicreg_invariance", colors[2], "-", "invariance"),
        ("vicreg_variance", colors[5], "--", "variance"),
        ("vicreg_covariance", colors[7], ":", "covariance"),
    ):
        array = values[name]
        vicreg_arrays.append(array)
        ax_vicreg.plot(steps, array, color=color, ls=linestyle, lw=0.8, label=label)
    _positive_scale(ax_vicreg, vicreg_arrays)
    for milestone in (500, 1000):
        ax_vicreg.axvline(milestone, color="0.62", ls=":", lw=0.8)
    ax_vicreg.set(xlabel="Optimizer step", ylabel="Raw VICReg term", xlim=(1, 2000))
    ax_vicreg.legend(loc="upper right", ncol=2)
    _panel_label(ax_vicreg, "c")

    gradient = values["embedding_projection_gradient_norm"]
    line_gradient = ax_optimization.plot(steps, gradient, color=colors[6], lw=0.8, label="projection-gradient norm")[0]
    _positive_scale(ax_optimization, [gradient])
    for milestone in (500, 1000):
        ax_optimization.axvline(milestone, color="0.62", ls=":", lw=0.8)
    ax_optimization.set(xlabel="Optimizer step", ylabel="Projection-gradient norm", xlim=(1, 2000))
    ax_lr = ax_optimization.twinx()
    line_lr = ax_lr.plot(steps, values["learning_rate"], color=colors[1], lw=0.9, label="learning rate")[0]
    ax_lr.set_ylabel("Learning rate")
    ax_lr.grid(False)
    ax_optimization.legend([line_gradient, line_lr], [line_gradient.get_label(), line_lr.get_label()], loc="center right")
    _panel_label(ax_optimization, "d")

    fig.subplots_adjust(left=0.10, right=0.86, bottom=0.10, top=0.98, hspace=0.32, wspace=0.34)
    stem = output_dir / "fm0_2_step2000_training_components"
    fig.savefig(stem.with_suffix(".png"), dpi=240)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    args = _parser().parse_args()
    health_by_step = {
        500: _load_json(args.step500_health),
        1000: _load_json(args.step1000_health),
        2000: _load_json(args.step2000_health),
    }
    receipt = _load_json(args.step2000_receipt)
    post_validation = _load_json(args.step2000_post_validation)
    step1000_history = _load_json(args.step1000_loss_history)
    step2000_history = _load_json(args.step2000_loss_history)
    rows = _validate_inputs(
        health_by_step,
        receipt,
        post_validation,
        step1000_history,
        step2000_history,
        args.step2000_health,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(args.output_dir / "training_trace.csv", rows, HISTORY_FIELDS)
    metric_rows = _representation_rows(health_by_step)
    _write_csv(
        args.output_dir / "step500_step1000_step2000_metrics.csv",
        metric_rows,
        list(metric_rows[0]),
    )
    gate_rows = _gate_rows(health_by_step[2000])
    _write_csv(
        args.output_dir / "step2000_representation_criteria.csv",
        gate_rows,
        list(gate_rows[0]),
    )
    _plot_representation(health_by_step, args.output_dir)
    _plot_training(rows, args.output_dir)
    print(args.output_dir / "fm0_2_step2000_representation_trajectory.png")
    print(args.output_dir / "fm0_2_step2000_training_components.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
