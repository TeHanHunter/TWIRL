"""Matched, development-only comparison figures for TWIRL-FM0.1.1/0.1.2."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.plotting.style import apply_twirl_style, get_ordered_palette


COMPARISON_SCHEMA_VERSION = "twirl_fm0_1_development_comparison_v1"
REPRESENTATION_HEALTH_SCHEMA_VERSION = "twirl_fm0_1_representation_health_v1"
VARIANT_ORDER: tuple[str, str] = ("TWIRL-FM0.1.1", "TWIRL-FM0.1.2")
EXPECTED_ARCHITECTURE: Mapping[str, str] = {
    "TWIRL-FM0.1.1": "tcn",
    "TWIRL-FM0.1.2": "conformer",
}
TARGET_STEP = 20_000
EFFECTIVE_BATCH_WINDOWS = 64
PROGRESS_INTERVAL_STEPS = 10
CURVE_BIN_WIDTH_STEPS = 100
FINAL_DIAGNOSTIC_WINDOW_STEPS = 1_000
MINIMUM_PARAMETERS = 8_000_000
MAXIMUM_PARAMETERS = 12_000_000
PARAMETER_MATCH_RELATIVE_TOLERANCE = 0.10
EXPECTED_LOGGED_STEPS: tuple[int, ...] = (
    1,
    *range(PROGRESS_INTERVAL_STEPS, TARGET_STEP + 1, PROGRESS_INTERVAL_STEPS),
)
_FLOAT_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
_TRAINING_LINE = re.compile(
    rf"\[fm0-train\]\s+step=(?P<step>\d+)/(?P<target>\d+)\s+"
    rf"loss=(?P<loss>{_FLOAT_PATTERN})\s+"
    rf"lr=(?P<learning_rate>{_FLOAT_PATTERN})\s+"
    rf"elapsed_s=(?P<elapsed_seconds>{_FLOAT_PATTERN})(?:\s|$)"
)


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, payload: Mapping[str, Any]) -> str:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return _file_sha256(path)


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected a JSON mapping in {path}")
    return payload


def _require_mapping(
    payload: Mapping[str, Any], key: str, *, context: str
) -> Mapping[str, Any]:
    value = payload.get(key)
    if not isinstance(value, Mapping):
        raise ValueError(f"{context}.{key} must be a mapping")
    return value


def _finite_float(value: Any, *, context: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{context} must be a finite number")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{context} must be a finite number") from exc
    if not np.isfinite(number):
        raise ValueError(f"{context} must be a finite number")
    return number


def _positive_int(value: Any, *, context: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{context} must be a positive integer")
    try:
        number = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{context} must be a positive integer") from exc
    if number <= 0 or number != value:
        raise ValueError(f"{context} must be a positive integer")
    return number


def _nonnegative_int(value: Any, *, context: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{context} must be a non-negative integer")
    try:
        number = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{context} must be a non-negative integer") from exc
    if number < 0 or number != value:
        raise ValueError(f"{context} must be a non-negative integer")
    return number


def parse_fm0_training_log(path: str | Path) -> pd.DataFrame:
    """Parse and strictly validate one real-run FM0 progress log.

    Real training reports step 1 and then every tenth step. The complete expected
    schedule therefore has 2,001 records rather than one record per optimizer
    update.
    """

    resolved = Path(path).expanduser().resolve(strict=True)
    records: list[dict[str, float | int]] = []
    malformed_lines: list[int] = []
    for line_number, line in enumerate(
        resolved.read_text(encoding="utf-8", errors="strict").splitlines(), start=1
    ):
        if "[fm0-train]" not in line:
            continue
        match = _TRAINING_LINE.search(line)
        if match is None:
            malformed_lines.append(line_number)
            continue
        records.append(
            {
                "step": int(match.group("step")),
                "target_step": int(match.group("target")),
                "loss": float(match.group("loss")),
                "learning_rate": float(match.group("learning_rate")),
                "elapsed_seconds": float(match.group("elapsed_seconds")),
            }
        )
    if malformed_lines:
        raise ValueError(
            f"training log has malformed [fm0-train] lines at {malformed_lines[:8]}"
        )
    if not records:
        raise ValueError(f"training log has no [fm0-train] records: {resolved}")
    frame = pd.DataFrame.from_records(records)
    if frame["step"].duplicated().any():
        duplicates = frame.loc[frame["step"].duplicated(keep=False), "step"].tolist()
        raise ValueError(
            f"training log has duplicate optimizer steps: {duplicates[:8]}"
        )
    observed_steps = tuple(int(value) for value in frame["step"])
    if observed_steps != EXPECTED_LOGGED_STEPS:
        missing = sorted(set(EXPECTED_LOGGED_STEPS) - set(observed_steps))
        extra = sorted(set(observed_steps) - set(EXPECTED_LOGGED_STEPS))
        raise ValueError(
            "training log does not match the exact step-1 plus 10-step schedule "
            f"through {TARGET_STEP}; missing={missing[:8]}, extra={extra[:8]}"
        )
    if set(frame["target_step"].astype(int)) != {TARGET_STEP}:
        raise ValueError(f"training log target denominator must be {TARGET_STEP}")
    numeric = frame[["loss", "learning_rate", "elapsed_seconds"]].to_numpy(
        dtype=float
    )
    if not np.isfinite(numeric).all():
        raise ValueError("training log contains non-finite numeric values")
    if (frame["loss"] < 0).any() or (frame["learning_rate"] < 0).any():
        raise ValueError("training loss and learning rate must be non-negative")
    elapsed = frame["elapsed_seconds"].to_numpy(dtype=float)
    if elapsed[0] < 0 or np.any(np.diff(elapsed) < 0):
        raise ValueError("training elapsed time must be finite and non-decreasing")
    return frame


def bin_fm0_training_curve(
    frame: pd.DataFrame,
    *,
    variant: str,
) -> pd.DataFrame:
    """Return deterministic, non-overlapping 100-step summaries."""

    required = {"step", "loss", "learning_rate", "elapsed_seconds"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise KeyError(f"training frame lacks {missing}")
    output: list[dict[str, Any]] = []
    for bin_start in range(1, TARGET_STEP + 1, CURVE_BIN_WIDTH_STEPS):
        bin_end = min(TARGET_STEP, bin_start + CURVE_BIN_WIDTH_STEPS - 1)
        selected = frame.loc[frame["step"].between(bin_start, bin_end)]
        if selected.empty:
            raise ValueError(f"training curve bin {bin_start}-{bin_end} is empty")
        loss = selected["loss"].to_numpy(dtype=float)
        output.append(
            {
                "variant": variant,
                "architecture": EXPECTED_ARCHITECTURE[variant],
                "bin_start_step": bin_start,
                "bin_end_step": bin_end,
                "bin_midpoint_step": 0.5 * (bin_start + bin_end),
                "n_logged_steps": int(len(selected)),
                "loss_median": float(np.median(loss)),
                "loss_p10": float(np.quantile(loss, 0.10)),
                "loss_p90": float(np.quantile(loss, 0.90)),
                "learning_rate_median": float(
                    np.median(selected["learning_rate"].to_numpy(dtype=float))
                ),
                "elapsed_seconds_max": float(selected["elapsed_seconds"].iloc[-1]),
            }
        )
    return pd.DataFrame.from_records(output)


def _validate_representation_payload(
    payload: Mapping[str, Any],
    *,
    expected_variant: str,
) -> None:
    if payload.get("schema_version") != REPRESENTATION_HEALTH_SCHEMA_VERSION:
        raise ValueError(f"{expected_variant} representation-health schema mismatch")
    if payload.get("passed") is not True:
        raise ValueError(
            f"{expected_variant} representation-health execution did not pass"
        )
    evaluator_git_sha = payload.get("evaluator_git_sha")
    if not isinstance(evaluator_git_sha, str) or not evaluator_git_sha.strip():
        raise ValueError(f"{expected_variant} lacks an evaluator Git revision")
    run = _require_mapping(payload, "run", context=expected_variant)
    if run.get("variant") != expected_variant:
        raise ValueError(f"representation payload is not for {expected_variant}")
    if run.get("architecture") != EXPECTED_ARCHITECTURE[expected_variant]:
        raise ValueError(f"{expected_variant} architecture binding is incorrect")
    global_step = _positive_int(
        run.get("global_step"), context=f"{expected_variant}.global_step"
    )
    if global_step != TARGET_STEP:
        raise ValueError(
            f"{expected_variant} representation checkpoint is not step {TARGET_STEP}"
        )
    population = _require_mapping(
        payload, "evaluation_population", context=expected_variant
    )
    for key in (
        "source_partition",
        "selected_leakage_components_sha256",
        "selected_observation_keys_sha256",
    ):
        value = population.get(key)
        if not isinstance(value, str) or not value.strip():
            raise ValueError(
                f"{expected_variant}.evaluation_population.{key} is absent"
            )
    _positive_int(
        population.get("selected_leakage_components"),
        context=f"{expected_variant}.selected_leakage_components",
    )
    _positive_int(
        population.get("selected_observation_visits"),
        context=f"{expected_variant}.selected_observation_visits",
    )


def _validate_matched_representation_payloads(
    payloads: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    for variant in VARIANT_ORDER:
        _validate_representation_payload(payloads[variant], expected_variant=variant)
    evaluator_revisions = {str(payloads[v]["evaluator_git_sha"]) for v in VARIANT_ORDER}
    if len(evaluator_revisions) != 1:
        raise ValueError("representation evaluations use different evaluator revisions")
    run_revisions = {
        str(_require_mapping(payloads[v], "run", context=v).get("run_git_sha"))
        for v in VARIANT_ORDER
    }
    if len(run_revisions) != 1 or next(iter(run_revisions)) in {"", "None"}:
        raise ValueError("representation evaluations use different training revisions")
    populations = {
        variant: _require_mapping(
            payloads[variant], "evaluation_population", context=variant
        )
        for variant in VARIANT_ORDER
    }
    matched_keys = (
        "source_partition",
        "selected_leakage_components",
        "selected_observation_visits",
        "selected_leakage_components_sha256",
        "selected_observation_keys_sha256",
    )
    matched: dict[str, Any] = {}
    for key in matched_keys:
        values = [populations[variant].get(key) for variant in VARIANT_ORDER]
        if values[0] != values[1]:
            raise ValueError(f"representation evaluation populations differ for {key}")
        matched[key] = values[0]
    matched["evaluator_git_sha"] = next(iter(evaluator_revisions))
    matched["run_git_sha"] = next(iter(run_revisions))
    if matched["source_partition"] != "poc_development":
        raise ValueError(
            "representation comparison is not restricted to poc_development"
        )
    return matched


def _interval(
    payload: Mapping[str, Any],
    key: str,
    *,
    context: str,
) -> tuple[float, float, float, int]:
    value = _require_mapping(payload, key, context=context)
    mean = _finite_float(value.get("mean"), context=f"{context}.{key}.mean")
    lower = _finite_float(value.get("lower"), context=f"{context}.{key}.lower")
    upper = _finite_float(value.get("upper"), context=f"{context}.{key}.upper")
    clusters = _positive_int(
        value.get("n_source_clusters"),
        context=f"{context}.{key}.n_source_clusters",
    )
    if lower > mean or mean > upper:
        raise ValueError(f"{context}.{key} does not bracket its mean")
    return mean, lower, upper, clusters


def _encoder_metrics(
    encoder: Mapping[str, Any],
    *,
    context: str,
) -> dict[str, Any]:
    health = _require_mapping(encoder, "embedding_health", context=context)
    separation = _require_mapping(
        encoder, "safe_mask_pair_separation", context=context
    )
    retrieval = _require_mapping(
        encoder, "same_component_cross_visit_retrieval", context=context
    )
    separation_raw = _finite_float(
        separation.get("paired_minus_unrelated_mean"),
        context=f"{context}.paired_minus_unrelated_mean",
    )
    separation_interval = _interval(
        separation,
        "paired_minus_unrelated_source_clustered_95_interval",
        context=f"{context}.safe_mask_pair_separation",
    )
    if retrieval.get("status") != "available":
        raise ValueError(f"{context} cross-visit retrieval is unavailable")
    retrieval_raw = _finite_float(
        retrieval.get("top1_same_component_retrieval"),
        context=f"{context}.top1_same_component_retrieval",
    )
    retrieval_interval = _interval(
        retrieval,
        "top1_source_clustered_95_interval",
        context=f"{context}.same_component_cross_visit_retrieval",
    )
    if not 0.0 <= retrieval_raw <= 1.0:
        raise ValueError(f"{context} retrieval is outside [0,1]")
    embedding_dim = _positive_int(
        health.get("embedding_dim"), context=f"{context}.embedding_dim"
    )
    effective_rank = _finite_float(
        health.get("effective_rank"), context=f"{context}.effective_rank"
    )
    if not 0.0 <= effective_rank <= embedding_dim:
        raise ValueError(f"{context} effective rank is outside [0, embedding_dim]")
    leading_share = _finite_float(
        health.get("leading_principal_component_share"),
        context=f"{context}.leading_principal_component_share",
    )
    if not 0.0 <= leading_share <= 1.0:
        raise ValueError(
            f"{context} leading principal-component share is outside [0,1]"
        )
    return {
        "embedding_dim": embedding_dim,
        "effective_rank": effective_rank,
        "zero_or_constant_dimensions": _nonnegative_int(
            health.get("zero_or_constant_dimensions"),
            context=f"{context}.zero_or_constant_dimensions",
        ),
        "near_duplicate_dimension_pairs": _nonnegative_int(
            health.get("near_duplicate_dimension_pairs"),
            context=f"{context}.near_duplicate_dimension_pairs",
        ),
        "leading_principal_component_share": leading_share,
        "separation": separation_interval[0],
        "separation_raw_pair_weighted": separation_raw,
        "separation_ci_lower": separation_interval[1],
        "separation_ci_upper": separation_interval[2],
        "separation_source_clusters": separation_interval[3],
        "separation_pairs": _positive_int(
            separation.get("n_pairs"), context=f"{context}.separation.n_pairs"
        ),
        "retrieval": retrieval_interval[0],
        "retrieval_raw_visit_weighted": retrieval_raw,
        "retrieval_ci_lower": retrieval_interval[1],
        "retrieval_ci_upper": retrieval_interval[2],
        "retrieval_source_clusters": retrieval_interval[3],
        "retrieval_queries": _positive_int(
            retrieval.get("n_repeated_component_queries"),
            context=f"{context}.retrieval.n_repeated_component_queries",
        ),
        "retrieval_visit_embeddings": _positive_int(
            retrieval.get("n_visit_embeddings"),
            context=f"{context}.retrieval.n_visit_embeddings",
        ),
    }


def _load_optional_summary(
    path: Path | None,
    *,
    variant: str,
) -> tuple[dict[str, Any] | None, str | None]:
    if path is None:
        return None, None
    resolved = Path(path).expanduser().resolve(strict=True)
    payload = _load_json(resolved)
    if payload.get("variant") != variant:
        raise ValueError(f"run summary is not for {variant}")
    if payload.get("architecture") != EXPECTED_ARCHITECTURE[variant]:
        raise ValueError(f"run summary architecture is not correct for {variant}")
    global_step = _positive_int(
        payload.get("global_step"), context=f"{variant}.summary.global_step"
    )
    if global_step != TARGET_STEP:
        raise ValueError(f"run summary for {variant} is not at step {TARGET_STEP}")
    parameter_count = _positive_int(
        payload.get("parameter_count"), context=f"{variant}.parameter_count"
    )
    if not MINIMUM_PARAMETERS <= parameter_count <= MAXIMUM_PARAMETERS:
        raise ValueError(
            f"{variant} parameter count is outside the frozen "
            f"{MINIMUM_PARAMETERS:,}--{MAXIMUM_PARAMETERS:,} range"
        )
    _finite_float(
        payload.get("elapsed_seconds_this_invocation"),
        context=f"{variant}.elapsed_seconds_this_invocation",
    )
    return payload, _file_sha256(resolved)


def _comparison_metric_row(
    *,
    variant: str,
    training: pd.DataFrame,
    representation: Mapping[str, Any],
    summary: Mapping[str, Any] | None,
) -> dict[str, Any]:
    run = _require_mapping(representation, "run", context=variant)
    population = _require_mapping(
        representation, "evaluation_population", context=variant
    )
    reconstruction = _require_mapping(
        representation, "masked_reconstruction", context=variant
    )
    masked_huber = _finite_float(
        reconstruction.get("masked_huber_mean"),
        context=f"{variant}.masked_huber_mean",
    )
    zero_huber = _finite_float(
        reconstruction.get("zero_prediction_masked_huber_mean"),
        context=f"{variant}.zero_prediction_masked_huber_mean",
    )
    ratio = _finite_float(
        reconstruction.get("model_to_zero_baseline_ratio"),
        context=f"{variant}.model_to_zero_baseline_ratio",
    )
    if zero_huber <= 0 or not np.isclose(ratio, masked_huber / zero_huber, rtol=2e-6):
        raise ValueError(f"{variant} reconstruction baseline ratio is inconsistent")
    trained = _encoder_metrics(
        _require_mapping(representation, "trained_encoder", context=variant),
        context=f"{variant}.trained_encoder",
    )
    random = _encoder_metrics(
        _require_mapping(representation, "random_encoder_control", context=variant),
        context=f"{variant}.random_encoder_control",
    )
    final_window = training.loc[
        training["step"].ge(TARGET_STEP - FINAL_DIAGNOSTIC_WINDOW_STEPS)
    ]
    post_warmup = training.loc[training["step"].gt(1_000)]
    final_losses = final_window["loss"].to_numpy(dtype=float)
    minimum_row = post_warmup.loc[post_warmup["loss"].idxmin()]
    elapsed_seconds = float(training["elapsed_seconds"].iloc[-1])
    row: dict[str, Any] = {
        "variant": variant,
        "architecture": EXPECTED_ARCHITECTURE[variant],
        "training_target_step": TARGET_STEP,
        "training_logged_records": int(len(training)),
        "training_final_loss": float(training["loss"].iloc[-1]),
        "training_final_learning_rate": float(training["learning_rate"].iloc[-1]),
        "training_elapsed_seconds": elapsed_seconds,
        "training_one_h200_gpu_hours": elapsed_seconds / 3_600.0,
        "training_effective_windows_per_second": (
            TARGET_STEP * EFFECTIVE_BATCH_WINDOWS / elapsed_seconds
        ),
        "training_final_1000_step_median_loss": float(np.median(final_losses)),
        "training_final_1000_step_p10_loss": float(
            np.quantile(final_losses, 0.10)
        ),
        "training_final_1000_step_p90_loss": float(
            np.quantile(final_losses, 0.90)
        ),
        "training_post_warmup_minimum_loss": float(minimum_row["loss"]),
        "training_post_warmup_minimum_loss_step": int(minimum_row["step"]),
        "run_git_sha": run.get("run_git_sha"),
        "checkpoint_sha256": run.get("checkpoint_sha256"),
        "evaluator_git_sha": representation.get("evaluator_git_sha"),
        "source_partition": population.get("source_partition"),
        "selected_leakage_components": int(
            population["selected_leakage_components"]
        ),
        "selected_observation_visits": int(population["selected_observation_visits"]),
        "selected_leakage_components_sha256": population.get(
            "selected_leakage_components_sha256"
        ),
        "selected_observation_keys_sha256": population.get(
            "selected_observation_keys_sha256"
        ),
        "development_masked_valid_target_count": _positive_int(
            reconstruction.get("masked_valid_target_count"),
            context=f"{variant}.masked_valid_target_count",
        ),
        "development_masked_huber_mean": masked_huber,
        "zero_prediction_masked_huber_mean": zero_huber,
        "model_to_zero_baseline_ratio": ratio,
        "masked_reconstruction_improvement_over_zero_percent": 100.0
        * (1.0 - ratio),
    }
    for prefix, metrics in (("trained", trained), ("random_control", random)):
        row.update({f"{prefix}_{key}": value for key, value in metrics.items()})
    if summary is not None:
        row["parameter_count"] = int(summary["parameter_count"])
        row["summary_elapsed_seconds_this_invocation"] = float(
            summary["elapsed_seconds_this_invocation"]
        )
    else:
        row["parameter_count"] = None
        row["summary_elapsed_seconds_this_invocation"] = None
    return row


def _plot_interval_marker(
    axis: Any,
    *,
    x: float,
    value: float,
    lower: float,
    upper: float,
    color: Any,
    trained: bool,
) -> None:
    axis.errorbar(
        [x],
        [value],
        yerr=[[value - lower], [upper - value]],
        fmt="o",
        markersize=5.2,
        markerfacecolor=color if trained else "white",
        markeredgecolor=color,
        markeredgewidth=1.1,
        ecolor=color,
        elinewidth=1.0,
        capsize=2.2,
        zorder=4,
    )


def _plot_comparison(
    *,
    curve_bins: pd.DataFrame,
    metrics: pd.DataFrame,
    output_dir: Path,
    effective_rank_gate: float,
) -> dict[str, Any]:
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    template = apply_twirl_style("full_page")
    figure, axes = plt.subplots(
        2,
        2,
        figsize=(template["figsize"][0], 6.25),
        constrained_layout=True,
    )
    colors = get_ordered_palette(2)
    color_by_variant = dict(zip(VARIANT_ORDER, colors))
    labels = {
        "TWIRL-FM0.1.1": "FM0.1.1 (TCN)",
        "TWIRL-FM0.1.2": "FM0.1.2 (Conformer)",
    }
    development_x = {
        "TWIRL-FM0.1.1": TARGET_STEP - 110,
        "TWIRL-FM0.1.2": TARGET_STEP + 110,
    }

    training_axis = axes[0, 0]
    training_axis.axvspan(0, 1_000, color="0.90", zorder=0)
    training_axis.text(
        500,
        0.98,
        "warm-up",
        ha="center",
        va="top",
        color="0.35",
        transform=training_axis.get_xaxis_transform(),
        fontsize=template["annotation_size"],
    )
    for variant in VARIANT_ORDER:
        subset = curve_bins.loc[curve_bins["variant"].eq(variant)]
        color = color_by_variant[variant]
        x = subset["bin_midpoint_step"].to_numpy(dtype=float)
        median = subset["loss_median"].to_numpy(dtype=float)
        training_axis.fill_between(
            x,
            subset["loss_p10"].to_numpy(dtype=float),
            subset["loss_p90"].to_numpy(dtype=float),
            color=color,
            alpha=0.18,
            linewidth=0,
        )
        training_axis.plot(x, median, color=color, linewidth=1.2, label=labels[variant])
        metric = metrics.loc[metrics["variant"].eq(variant)].iloc[0]
        training_axis.scatter(
            [development_x[variant]],
            [metric["development_masked_huber_mean"]],
            marker="D",
            s=27,
            facecolor=color,
            edgecolor="white",
            linewidth=0.6,
            zorder=5,
        )
        training_axis.scatter(
            [development_x[variant]],
            [metric["zero_prediction_masked_huber_mean"]],
            marker="x",
            s=29,
            color=color,
            linewidth=1.0,
            zorder=5,
        )
    training_axis.set_xlim(0, TARGET_STEP * 1.025)
    post_warmup_bins = curve_bins.loc[curve_bins["bin_start_step"].gt(1_000)]
    diagnostic_low = min(
        float(post_warmup_bins["loss_p10"].min()),
        float(metrics["development_masked_huber_mean"].min()),
    )
    diagnostic_high = max(
        float(post_warmup_bins["loss_p90"].max()),
        float(metrics["zero_prediction_masked_huber_mean"].max()),
    )
    diagnostic_padding = 0.08 * (diagnostic_high - diagnostic_low)
    training_axis.set_ylim(
        diagnostic_low - diagnostic_padding,
        diagnostic_high + diagnostic_padding,
    )
    training_axis.set_xlabel("Optimizer step")
    training_axis.set_ylabel("Masked Huber loss (post-warm-up scale)")
    training_axis.set_title("(a) Reconstruction diagnostics")
    training_axis.text(
        0.02,
        0.03,
        "Curves: 100-step train median (10th–90th percentile)\n"
        "Diamond: dev model; cross: dev zero baseline\n"
        "Early warm-up values are clipped; loss is not model selection.",
        transform=training_axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=template["annotation_size"],
        bbox={"facecolor": "white", "edgecolor": "0.75", "alpha": 0.88, "pad": 2},
    )
    training_axis.legend(loc="upper right")

    x_base = np.arange(len(VARIANT_ORDER), dtype=float)
    offsets = {"trained": -0.09, "random_control": 0.09}
    rank_axis = axes[0, 1]
    rank_axis.axhline(
        effective_rank_gate,
        color="0.35",
        linestyle="--",
        linewidth=0.9,
        label=f"provisional gate = {effective_rank_gate:g}",
    )
    for index, variant in enumerate(VARIANT_ORDER):
        metric = metrics.loc[metrics["variant"].eq(variant)].iloc[0]
        color = color_by_variant[variant]
        trained_rank = float(metric["trained_effective_rank"])
        random_rank = float(metric["random_control_effective_rank"])
        rank_axis.plot(
            [index + offsets["trained"], index + offsets["random_control"]],
            [trained_rank, random_rank],
            color=color,
            alpha=0.35,
            linewidth=0.8,
            zorder=2,
        )
        rank_axis.scatter(
            [index + offsets["trained"]],
            [trained_rank],
            s=34,
            facecolor=color,
            edgecolor=color,
            zorder=4,
        )
        rank_axis.scatter(
            [index + offsets["random_control"]],
            [random_rank],
            s=34,
            facecolor="white",
            edgecolor=color,
            linewidth=1.1,
            zorder=4,
        )
        horizontal_offset = 7 if index == 0 else -7
        horizontal_alignment = "left" if index == 0 else "right"
        rank_axis.annotate(
            f"const={int(metric['trained_zero_or_constant_dimensions'])}; "
            f"PC1={float(metric['trained_leading_principal_component_share']):.2f}",
            (index + offsets["trained"], trained_rank),
            xytext=(horizontal_offset, 5),
            textcoords="offset points",
            ha=horizontal_alignment,
            va="bottom",
            fontsize=template["annotation_size"],
        )
    rank_axis.set_xticks(x_base, [labels[variant] for variant in VARIANT_ORDER])
    rank_axis.set_ylabel("Effective embedding rank")
    rank_axis.set_title("(b) Embedding health")
    rank_axis.set_ylim(bottom=0)
    rank_axis.legend(loc="best")

    separation_axis = axes[1, 0]
    separation_axis.axhline(0.0, color="0.35", linestyle="--", linewidth=0.9)
    retrieval_axis = axes[1, 1]
    for index, variant in enumerate(VARIANT_ORDER):
        metric = metrics.loc[metrics["variant"].eq(variant)].iloc[0]
        color = color_by_variant[variant]
        for prefix, trained in (("trained", True), ("random_control", False)):
            _plot_interval_marker(
                separation_axis,
                x=index + offsets[prefix],
                value=float(metric[f"{prefix}_separation"]),
                lower=float(metric[f"{prefix}_separation_ci_lower"]),
                upper=float(metric[f"{prefix}_separation_ci_upper"]),
                color=color,
                trained=trained,
            )
            _plot_interval_marker(
                retrieval_axis,
                x=index + offsets[prefix],
                value=float(metric[f"{prefix}_retrieval"]),
                lower=float(metric[f"{prefix}_retrieval_ci_lower"]),
                upper=float(metric[f"{prefix}_retrieval_ci_upper"]),
                color=color,
                trained=trained,
            )
    categorical_labels = [labels[variant] for variant in VARIANT_ORDER]
    separation_axis.set_xticks(x_base, categorical_labels)
    separation_axis.set_ylabel("Paired − unrelated cosine similarity")
    separation_axis.set_title("(c) Mask-pair separation (clustered 95% CI)")
    retrieval_axis.set_xticks(x_base, categorical_labels)
    retrieval_axis.set_ylabel("Top-1 same-component retrieval")
    retrieval_axis.set_title("(d) Cross-visit retrieval (clustered 95% CI)")
    retrieval_upper = max(
        float(metrics[f"{prefix}_retrieval_ci_upper"].max())
        for prefix in ("trained", "random_control")
    )
    retrieval_axis.set_ylim(0.0, min(1.03, max(0.10, 1.22 * retrieval_upper)))
    retrieval_axis.text(
        0.02,
        0.78,
        "Cross-visit only\n(no query-sector exclusion)",
        transform=retrieval_axis.transAxes,
        ha="left",
        va="top",
        fontsize=template["annotation_size"],
        bbox={"facecolor": "white", "edgecolor": "0.75", "alpha": 0.88, "pad": 2},
    )
    state_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="0.25",
            markerfacecolor="0.25",
            linestyle="none",
            label="trained",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="0.25",
            markerfacecolor="white",
            linestyle="none",
            label="random initialization",
        ),
    ]
    separation_axis.legend(handles=state_handles, loc="best")
    retrieval_axis.legend(handles=state_handles, loc="best")

    stem = output_dir / "fm0_1_1_vs_fm0_1_2_development"
    png = stem.with_suffix(".png")
    pdf = stem.with_suffix(".pdf")
    figure.savefig(png, dpi=300, bbox_inches="tight")
    figure.savefig(pdf, bbox_inches="tight")
    plt.close(figure)
    return {
        "png": str(png),
        "png_sha256": _file_sha256(png),
        "pdf": str(pdf),
        "pdf_sha256": _file_sha256(pdf),
    }


def make_fm0_development_comparison(
    *,
    fm0_1_1_log_path: str | Path,
    fm0_1_2_log_path: str | Path,
    fm0_1_1_representation_path: str | Path,
    fm0_1_2_representation_path: str | Path,
    output_dir: str | Path,
    fm0_1_1_summary_path: str | Path | None = None,
    fm0_1_2_summary_path: str | Path | None = None,
    effective_rank_gate: float = 26.0,
) -> dict[str, Any]:
    """Build matched numeric artifacts and the four-panel development figure."""

    gate = _finite_float(effective_rank_gate, context="effective_rank_gate")
    if gate <= 0:
        raise ValueError("effective_rank_gate must be positive")
    output = Path(output_dir).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    log_paths = {
        "TWIRL-FM0.1.1": Path(fm0_1_1_log_path).expanduser().resolve(strict=True),
        "TWIRL-FM0.1.2": Path(fm0_1_2_log_path).expanduser().resolve(strict=True),
    }
    representation_paths = {
        "TWIRL-FM0.1.1": Path(fm0_1_1_representation_path)
        .expanduser()
        .resolve(strict=True),
        "TWIRL-FM0.1.2": Path(fm0_1_2_representation_path)
        .expanduser()
        .resolve(strict=True),
    }
    summary_paths: Mapping[str, Path | None] = {
        "TWIRL-FM0.1.1": (
            Path(fm0_1_1_summary_path).expanduser().resolve(strict=True)
            if fm0_1_1_summary_path is not None
            else None
        ),
        "TWIRL-FM0.1.2": (
            Path(fm0_1_2_summary_path).expanduser().resolve(strict=True)
            if fm0_1_2_summary_path is not None
            else None
        ),
    }
    training = {
        variant: parse_fm0_training_log(log_paths[variant])
        for variant in VARIANT_ORDER
    }
    representation = {
        variant: _load_json(representation_paths[variant])
        for variant in VARIANT_ORDER
    }
    matched = _validate_matched_representation_payloads(representation)
    summaries: dict[str, Mapping[str, Any] | None] = {}
    summary_hashes: dict[str, str | None] = {}
    for variant in VARIANT_ORDER:
        summaries[variant], summary_hashes[variant] = _load_optional_summary(
            summary_paths[variant], variant=variant
        )

    curve_bins = pd.concat(
        [
            bin_fm0_training_curve(training[variant], variant=variant)
            for variant in VARIANT_ORDER
        ],
        ignore_index=True,
    )
    metrics = pd.DataFrame.from_records(
        [
            _comparison_metric_row(
                variant=variant,
                training=training[variant],
                representation=representation[variant],
                summary=summaries[variant],
            )
            for variant in VARIANT_ORDER
        ]
    )
    metrics["trained_effective_rank_gate"] = gate
    metrics["trained_effective_rank_gate_pass"] = metrics[
        "trained_effective_rank"
    ].ge(gate)
    metrics["trained_zero_constant_dimensions_gate_pass"] = metrics[
        "trained_zero_or_constant_dimensions"
    ].eq(0)
    metrics["trained_positive_separation_lower_bound_gate_pass"] = metrics[
        "trained_separation_ci_lower"
    ].gt(0.0)
    metrics["all_current_representation_hard_gates_pass"] = (
        metrics["trained_effective_rank_gate_pass"]
        & metrics["trained_zero_constant_dimensions_gate_pass"]
        & metrics["trained_positive_separation_lower_bound_gate_pass"]
    )
    parameter_match: dict[str, Any] | None = None
    if metrics["parameter_count"].notna().all():
        counts = [int(value) for value in metrics["parameter_count"]]
        relative_difference = abs(counts[0] - counts[1]) / max(counts)
        if relative_difference > PARAMETER_MATCH_RELATIVE_TOLERANCE:
            raise ValueError(
                "FM0.1.1 and FM0.1.2 parameter counts exceed the frozen "
                "relative-match tolerance"
            )
        parameter_match = {
            "counts": dict(zip(VARIANT_ORDER, counts)),
            "relative_difference_against_larger_model": relative_difference,
            "relative_tolerance": PARAMETER_MATCH_RELATIVE_TOLERANCE,
            "passed": True,
        }
    curve_path = output / "training_curve_bins.csv"
    metrics_path = output / "comparison_metrics.csv"
    curve_bins.to_csv(
        curve_path, index=False, lineterminator="\n", float_format="%.10g"
    )
    metrics.to_csv(
        metrics_path, index=False, lineterminator="\n", float_format="%.10g"
    )
    figure_artifacts = _plot_comparison(
        curve_bins=curve_bins,
        metrics=metrics,
        output_dir=output,
        effective_rank_gate=gate,
    )
    provenance_path = output / "comparison.provenance.json"
    provenance: dict[str, Any] = {
        "schema_version": COMPARISON_SCHEMA_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "matched development-only FM0.1 architecture comparison",
        "claim_limit": (
            "training loss is diagnostic only; cross-visit retrieval is not a "
            "cross-sector test; no sealed-test, science, or foundation-model claim"
        ),
        "training_curve_contract": {
            "target_step": TARGET_STEP,
            "effective_batch_windows": EFFECTIVE_BATCH_WINDOWS,
            "logged_steps": "step 1 plus every 10th step through 20000",
            "expected_logged_records_per_model": len(EXPECTED_LOGGED_STEPS),
            "nonoverlapping_bin_width_steps": CURVE_BIN_WIDTH_STEPS,
            "statistics": ["median", "10th percentile", "90th percentile"],
            "final_diagnostic_window_steps": FINAL_DIAGNOSTIC_WINDOW_STEPS,
            "minimum_training_loss_is_model_selection_metric": False,
        },
        "effective_rank_gate": gate,
        "architecture_parameter_match": parameter_match,
        "matched_evaluation_binding": matched,
        "sources": {
            variant: {
                "training_log": str(log_paths[variant]),
                "training_log_sha256": _file_sha256(log_paths[variant]),
                "representation_health": str(representation_paths[variant]),
                "representation_health_sha256": _file_sha256(
                    representation_paths[variant]
                ),
                "run_summary": (
                    str(summary_paths[variant])
                    if summary_paths[variant] is not None
                    else None
                ),
                "run_summary_sha256": summary_hashes[variant],
            }
            for variant in VARIANT_ORDER
        },
        "artifacts": {
            "training_curve_bins_csv": str(curve_path),
            "training_curve_bins_csv_sha256": _file_sha256(curve_path),
            "comparison_metrics_csv": str(metrics_path),
            "comparison_metrics_csv_sha256": _file_sha256(metrics_path),
            **figure_artifacts,
        },
    }
    provenance_sha = _write_json(provenance_path, provenance)
    return {
        "schema_version": COMPARISON_SCHEMA_VERSION,
        "comparison_metrics_csv": str(metrics_path),
        "training_curve_bins_csv": str(curve_path),
        **figure_artifacts,
        "provenance_json": str(provenance_path),
        "provenance_json_sha256": provenance_sha,
    }


__all__: Sequence[str] = (
    "COMPARISON_SCHEMA_VERSION",
    "EXPECTED_LOGGED_STEPS",
    "VARIANT_ORDER",
    "bin_fm0_training_curve",
    "make_fm0_development_comparison",
    "parse_fm0_training_log",
)
