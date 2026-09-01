"""Fast oracle-support mechanics control for the FM event-transfer probe.

This module is deliberately separate from :mod:`event_transfer_canary` so the
published v1/v2 evaluator and result paths remain unchanged.  The v3 control
fits and evaluates only the raw and quality-only cadence scorers.  It emits a
logit for every native 200-second cadence, then uses held-out injection support
truth only to select the single center cadence of each exact clean/injected
pair for mechanics metrics.  Support is never appended to the input features.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from .event_transfer_canary import (
    EVENT_CELLS,
    EVENT_DEPTHS,
    EVENT_DURATIONS,
    JITTER_BOUNDS,
    MINIMUM_COMPONENTS,
    SPLIT_SECTORS,
    TARGET_COMPONENTS,
    ProbeFit,
    build_paired_samples_with_token_targets,
    freeze_component_schedule,
    raw_cadence_features,
)
from .temporal_panel import DEVELOPMENT_PARTITION
from .temporal_zero_shot import TEMPORAL_COHORTS, load_temporal_panel

EVENT_TRANSFER_MECHANICS_CONFIG_SCHEMA_VERSION = (
    "twirl_fm0_event_transfer_mechanics_config_v3"
)
EVENT_TRANSFER_MECHANICS_RESULT_SCHEMA_VERSION = (
    "twirl_fm0_event_transfer_mechanics_result_v3"
)
EVENT_TRANSFER_MECHANICS_CAMPAIGN_ID = "twirl_fm0_2_s66_s77_event_transfer_mechanics_v3"
CONTROL_SPECS = (
    "raw_adp_validity_error_cadence_128",
    "quality_only_cadence_128",
)
CONTEXT_CADENCES = 128


@dataclass(frozen=True)
class CenterPairSelection:
    """Exact adjacent clean/injected pair and support-center indices."""

    clean_rows: np.ndarray
    injected_rows: np.ndarray
    center_cadences: np.ndarray
    components: tuple[str, ...]


def _mapping(value: Any, *, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be a mapping")
    return value


def _exact_bool(mapping: Mapping[str, Any], key: str, expected: bool) -> None:
    if mapping.get(key) is not expected:
        raise ValueError(f"frozen field {key!r} must be {expected}")


def _exact_int(value: Any, *, label: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise TypeError(f"{label} must be an integer")
    return int(value)


def _require_sha256(value: Any, *, label: str) -> str:
    text = str(value)
    if len(text) != 64 or any(
        character not in "0123456789abcdef" for character in text
    ):
        raise ValueError(f"{label} must be a lowercase SHA-256")
    return text


def load_event_transfer_mechanics_config(
    config_path: str | Path,
) -> tuple[dict[str, Any], Path, str]:
    """Validate and byte-bind the exact v3 mechanics-only contract."""

    source = Path(config_path).expanduser().resolve(strict=True)
    payload_bytes = source.read_bytes()
    config = dict(
        _mapping(yaml.safe_load(payload_bytes), label="event-transfer mechanics config")
    )
    if (
        config.get("schema_version") != EVENT_TRANSFER_MECHANICS_CONFIG_SCHEMA_VERSION
        or config.get("campaign_id") != EVENT_TRANSFER_MECHANICS_CAMPAIGN_ID
        or config.get("model_family") != "TWIRL-FM0"
        or config.get("status") != "frozen_bls_free_development_mechanics_control"
    ):
        raise ValueError("event-transfer mechanics config identity differs")

    purpose = _mapping(config.get("purpose"), label="purpose")
    if purpose.get("primary") != "verify_raw_and_quality_cadence_probe_mechanics":
        raise ValueError("event-transfer mechanics purpose differs")
    for key, expected in (
        ("mechanics_only", True),
        ("formal_model_gate", False),
        ("fm_metric_interpretation", False),
        ("production_model_claim", False),
    ):
        _exact_bool(purpose, key, expected)

    authorization = _mapping(config.get("authorization"), label="authorization")
    for key, expected in (
        ("development_shard_payload_access", True),
        ("artificial_single_event_injection", True),
        ("low_capacity_probe_fitting", True),
        ("frozen_checkpoint_inference", False),
        ("fm_encoder_optimizer_or_training", False),
        ("bls_period_candidate_or_event_features", False),
        ("sealed_test_access", False),
    ):
        _exact_bool(authorization, key, expected)

    inputs = _mapping(config.get("inputs"), label="inputs")
    _require_sha256(
        inputs.get("temporal_panel_receipt_sha256"),
        label="temporal-panel receipt",
    )
    if (
        tuple(inputs.get("allowed_sectors", ())) != tuple(range(66, 78))
        or inputs.get("excluded_sectors") != [65]
        or inputs.get("source_partition") != DEVELOPMENT_PARTITION
        or inputs.get("variant") != "TWIRL-FM0.2.1"
    ):
        raise ValueError("event-transfer mechanics input identity differs")

    split = _mapping(config.get("split"), label="split")
    for name, sectors in SPLIT_SECTORS.items():
        if tuple(split.get(f"{name}_sectors", ())) != sectors:
            raise ValueError(f"event-transfer mechanics {name} sectors differ")
    if (
        split.get("component_spanning_blocks") != "quarantine"
        or split.get("shortfall_policy") != "fail_closed"
        or split.get("component_and_visit_order") != "deterministic_sha256"
        or tuple(split.get("cohorts", ())) != TEMPORAL_COHORTS
        or _exact_int(
            split.get("target_components_per_cohort_per_block"),
            label="target components",
        )
        != TARGET_COMPONENTS
        or _exact_int(
            split.get("minimum_components_per_cohort_per_block"),
            label="minimum components",
        )
        != MINIMUM_COMPONENTS
    ):
        raise ValueError("event-transfer mechanics split differs")
    _exact_bool(split, "one_visit_per_component", True)

    temporal = _mapping(config.get("temporal_resolution"), label="temporal resolution")
    if (
        _exact_int(temporal.get("nominal_cadence_seconds"), label="nominal cadence")
        != 200
        or _exact_int(temporal.get("context_cadences"), label="context cadences")
        != CONTEXT_CADENCES
        or temporal.get("held_out_reduction")
        != "exact_injected_support_center_vs_exact_paired_clean_cadence"
    ):
        raise ValueError("event-transfer mechanics cadence grid differs")
    for key, expected in (
        ("one_logit_per_native_cadence", True),
        ("direct_centered_crop", True),
        ("temporal_downsampling", False),
        ("patching_or_cadence_averaging", False),
        ("representation_pooling", False),
        ("off_support_used_for_fit_or_evaluation", False),
    ):
        _exact_bool(temporal, key, expected)

    event = _mapping(config.get("artificial_event"), label="artificial event")
    if (
        event.get("profile") != "symmetric_trapezoid_sampled_at_cadence_centers"
        or tuple(event.get("duration_cadences", ())) != EVENT_DURATIONS
        or tuple(float(value) for value in event.get("fractional_depths", ()))
        != EVENT_DEPTHS
        or tuple(event.get("center_jitter_cadences_inclusive", ())) != JITTER_BOUNDS
        or event.get("apply_to_views") != ["adp_1x1", "adp_3x3"]
        or event.get("support_role")
        != "training_target_and_held_out_evaluation_truth_never_input_feature"
    ):
        raise ValueError("event-transfer mechanics event definition differs")
    for key, expected in (
        ("one_clean_and_one_injected_sample_per_component", True),
        ("one_event_per_injected_sample", True),
        ("period_defined", False),
        ("injection_truth_used_as_feature", False),
    ):
        _exact_bool(event, key, expected)

    features = _mapping(config.get("features"), label="features")
    if tuple(features.get("controls", ())) != CONTROL_SPECS:
        raise ValueError("event-transfer mechanics controls differ")
    if tuple(features.get("raw_channels", ())) != (
        "adp_1x1",
        "adp_3x3",
        "flux_valid_1x1",
        "flux_valid_3x3",
        "flux_error_1x1",
        "flux_error_3x3",
        "error_valid_1x1",
        "error_valid_3x3",
        "local_time_cadences",
        "delta_time_cadences",
        "time_valid",
        "segment_boundary",
    ) or tuple(features.get("quality_only_channels", ())) != (
        "flux_valid_1x1",
        "flux_valid_3x3",
        "error_valid_1x1",
        "error_valid_3x3",
        "local_time_cadences",
        "delta_time_cadences",
        "time_valid",
        "segment_boundary",
    ):
        raise ValueError("event-transfer mechanics cadence features differ")
    if tuple(features.get("forbidden_feature_families", ())) != (
        "period",
        "bls",
        "candidate_score",
        "event_duration",
        "event_depth",
        "event_center_numeric_value",
        "injection_support",
    ):
        raise ValueError("event-transfer mechanics forbidden features differ")

    probe = _mapping(config.get("probe"), label="probe")
    if (
        probe.get("family") != "shared_linear_cadence_scorer"
        or probe.get("training_objective")
        != "pair_balanced_center_support_per_cadence_binary_cross_entropy"
        or probe.get("standardization")
        != "training_clean_and_injected_center_cadences_only"
        or _exact_int(probe.get("epochs"), label="probe epochs") != 400
        or _exact_int(probe.get("initialization_seed"), label="probe seed") != 560203
        or not math.isclose(
            float(probe.get("learning_rate", math.nan)),
            0.02,
            rel_tol=0.0,
            abs_tol=0.0,
        )
        or not math.isclose(
            float(probe.get("l2_weight", math.nan)),
            0.001,
            rel_tol=0.0,
            abs_tol=0.0,
        )
    ):
        raise ValueError("event-transfer mechanics probe differs")
    for key, expected in (
        ("one_example_per_clean_or_injected_component_row", True),
        ("support_appended_to_features", False),
        ("off_support_used_for_standardization_or_fit", False),
    ):
        _exact_bool(probe, key, expected)

    evaluation = _mapping(config.get("evaluation"), label="evaluation")
    if (
        evaluation.get("split") != "locked_development_test"
        or evaluation.get("support_use")
        != "select_center_cadence_as_held_out_truth_after_all_logits_are_emitted"
        or evaluation.get("comparison")
        != "injected_center_logit_vs_exact_paired_clean_center_logit"
        or evaluation.get("pair_weighting") != "one_equal_vote_per_component"
        or evaluation.get("duration_weighting")
        != "one_equal_vote_per_component_not_per_support_cadence"
        or evaluation.get("uncertainty")
        != "deterministic_component_pair_bootstrap_95_interval"
        or _exact_int(
            evaluation.get("bootstrap_replicates"), label="bootstrap replicates"
        )
        != 1_000
        or _exact_int(evaluation.get("bootstrap_seed"), label="bootstrap seed")
        != 560205
    ):
        raise ValueError("event-transfer mechanics evaluation differs")

    mechanics = _mapping(config.get("mechanics"), label="mechanics")
    if dict(mechanics) != {
        "raw_minimum_center_paired_ranking_lower_95": 0.95,
        "raw_minimum_each_depth_0_30_center_paired_ranking": 0.90,
        "raw_require_negative_adp_flux_coefficients": True,
        "quality_center_paired_ranking_required_exact": 0.50,
    }:
        raise ValueError("event-transfer mechanics gate differs")
    return config, source, hashlib.sha256(payload_bytes).hexdigest()


def paired_support_center_indices(
    token_valid: np.ndarray,
    token_targets: np.ndarray,
    components: Sequence[str],
) -> CenterPairSelection:
    """Validate adjacent pairs and locate one exact center target per pair."""

    valid = np.asarray(token_valid, dtype=bool)
    targets = np.asarray(token_targets, dtype=bool)
    component_rows = tuple(str(value) for value in components)
    if valid.ndim != 2 or targets.shape != valid.shape:
        raise ValueError("center-pair validity and support must be aligned matrices")
    if valid.shape[0] < 2 or valid.shape[0] % 2:
        raise ValueError("center-pair control requires adjacent clean/injected rows")
    if len(component_rows) != valid.shape[0]:
        raise ValueError("center-pair components do not align with rows")
    if np.any(targets & ~valid):
        raise ValueError("center-pair support includes an invalid injected cadence")

    clean_rows = np.arange(0, valid.shape[0], 2, dtype=np.int64)
    injected_rows = clean_rows + 1
    centers = np.empty(clean_rows.size, dtype=np.int64)
    pair_components: list[str] = []
    seen: set[str] = set()
    for pair_index, (clean_row, injected_row) in enumerate(
        zip(clean_rows, injected_rows, strict=True)
    ):
        component = component_rows[clean_row]
        if (
            not component
            or component_rows[injected_row] != component
            or component in seen
        ):
            raise ValueError("center-pair component identity or adjacency differs")
        seen.add(component)
        if not np.array_equal(valid[clean_row], valid[injected_row]):
            raise ValueError("center-pair clean/injected validity masks differ")
        if np.any(targets[clean_row]):
            raise ValueError("center-pair clean row contains event support")
        support = np.flatnonzero(targets[injected_row])
        if (
            support.size == 0
            or int(support.size) not in EVENT_DURATIONS
            or support.size % 2 == 0
            or not np.all(np.diff(support) == 1)
        ):
            raise ValueError(
                "center-pair injected support must be contiguous with duration 1, 3, or 9"
            )
        center = int(support[support.size // 2])
        if not valid[clean_row, center] or not valid[injected_row, center]:
            raise ValueError("center-pair center cadence is not common-visible")
        centers[pair_index] = center
        pair_components.append(component)
    return CenterPairSelection(
        clean_rows=clean_rows,
        injected_rows=injected_rows,
        center_cadences=centers,
        components=tuple(pair_components),
    )


def center_pair_examples(
    tokens: np.ndarray,
    token_valid: np.ndarray,
    token_targets: np.ndarray,
    components: Sequence[str],
) -> tuple[np.ndarray, np.ndarray, CenterPairSelection]:
    """Return exactly one clean and one injected center token per component."""

    values = np.asarray(tokens, dtype=np.float32)
    valid = np.asarray(token_valid, dtype=bool)
    if values.ndim != 3 or valid.shape != values.shape[:2]:
        raise ValueError("center-pair cadence feature cube is invalid")
    selection = paired_support_center_indices(valid, token_targets, components)
    clean = values[selection.clean_rows, selection.center_cadences]
    injected = values[selection.injected_rows, selection.center_cadences]
    if not np.all(np.isfinite(clean)) or not np.all(np.isfinite(injected)):
        raise ValueError("center-pair selected cadence features are nonfinite")
    examples = np.empty((2 * clean.shape[0], values.shape[2]), dtype=np.float32)
    examples[0::2] = clean
    examples[1::2] = injected
    labels = np.tile(np.asarray([0.0, 1.0], dtype=np.float32), clean.shape[0])
    return examples, labels, selection


def fit_shared_linear_center_probe(
    tokens: np.ndarray,
    token_valid: np.ndarray,
    token_targets: np.ndarray,
    components: Sequence[str],
    *,
    learning_rate: float = 0.02,
    l2_weight: float = 0.001,
    epochs: int = 400,
    seed: int = 560203,
    progress_label: str | None = None,
) -> ProbeFit:
    """Fit on one exact center cadence per clean/injected component row."""

    examples, labels, _selection = center_pair_examples(
        tokens,
        token_valid,
        token_targets,
        components,
    )
    center = np.mean(examples, axis=0)
    scale = np.std(examples, axis=0)
    scale = np.where(np.isfinite(scale) & (scale > 1.0e-8), scale, 1.0)
    standardized = (examples - center[None, :]) / scale[None, :]

    rng = np.random.default_rng(int(seed))
    weight = rng.normal(0.0, 1.0e-3, size=examples.shape[1]).astype(np.float32)
    bias = 0.0
    m_w = np.zeros_like(weight)
    v_w = np.zeros_like(weight)
    m_b = v_b = 0.0
    history: list[float] = []
    beta1, beta2 = 0.9, 0.999
    n_epochs = _exact_int(epochs, label="probe epochs")
    for step in range(1, n_epochs + 1):
        logits = standardized @ weight + bias
        probabilities = np.empty_like(logits, dtype=np.float64)
        nonnegative = logits >= 0.0
        probabilities[nonnegative] = 1.0 / (1.0 + np.exp(-logits[nonnegative]))
        exponential = np.exp(logits[~nonnegative])
        probabilities[~nonnegative] = exponential / (1.0 + exponential)
        residual = probabilities - labels
        grad_w = np.mean(residual[:, None] * standardized, axis=0)
        grad_w += float(l2_weight) * weight
        grad_b = float(np.mean(residual))
        m_w = beta1 * m_w + (1.0 - beta1) * grad_w
        v_w = beta2 * v_w + (1.0 - beta2) * np.square(grad_w)
        m_b = beta1 * m_b + (1.0 - beta1) * grad_b
        v_b = beta2 * v_b + (1.0 - beta2) * grad_b * grad_b
        weight -= (
            float(learning_rate)
            * (m_w / (1.0 - beta1**step))
            / (np.sqrt(v_w / (1.0 - beta2**step)) + 1.0e-8)
        )
        bias -= (
            float(learning_rate)
            * (m_b / (1.0 - beta1**step))
            / (math.sqrt(v_b / (1.0 - beta2**step)) + 1.0e-8)
        )
        if step == 1 or step % 20 == 0 or step == n_epochs:
            history.append(
                float(
                    np.mean(np.logaddexp(0.0, logits) - labels * logits)
                    + 0.5 * float(l2_weight) * np.dot(weight, weight)
                )
            )
        if progress_label is not None and (
            step == 1 or step % 100 == 0 or step == n_epochs
        ):
            print(
                f"FM_EVENT_TRANSFER_MECHANICS control={progress_label} "
                f"phase=fit_center epoch={step}/{n_epochs}",
                flush=True,
            )
    return ProbeFit(
        weight=weight,
        bias=float(bias),
        center=center,
        scale=scale,
        objective_history=tuple(history),
    )


def score_shared_linear_cadence_logits(
    fit: ProbeFit,
    tokens: np.ndarray,
    token_valid: np.ndarray,
) -> np.ndarray:
    """Emit one independently computed logit for every native cadence."""

    values = np.asarray(tokens, dtype=np.float32)
    valid = np.asarray(token_valid, dtype=bool)
    if (
        values.ndim != 3
        or valid.shape != values.shape[:2]
        or not np.all(np.isfinite(values))
        or fit.weight.shape != (values.shape[2],)
        or fit.center.shape != fit.weight.shape
        or fit.scale.shape != fit.weight.shape
    ):
        raise ValueError("cadence-logit feature cube or probe dimensions differ")
    logits = (
        (values - fit.center[None, None, :]) / fit.scale[None, None, :]
    ) @ fit.weight + float(fit.bias)
    if logits.shape != valid.shape or not np.all(np.isfinite(logits)):
        raise ValueError("cadence-logit output lost native cadence alignment")
    return logits


def oracle_support_center_pair_scores(
    cadence_logits: np.ndarray,
    token_valid: np.ndarray,
    token_targets: np.ndarray,
    components: Sequence[str],
) -> tuple[np.ndarray, np.ndarray, CenterPairSelection]:
    """Select held-out center truth after logits exist; never score off-support."""

    logits = np.asarray(cadence_logits, dtype=np.float64)
    valid = np.asarray(token_valid, dtype=bool)
    if logits.shape != valid.shape:
        raise ValueError("center-pair logits and validity differ")
    selection = paired_support_center_indices(valid, token_targets, components)
    clean = logits[selection.clean_rows, selection.center_cadences]
    injected = logits[selection.injected_rows, selection.center_cadences]
    if not np.all(np.isfinite(clean)) or not np.all(np.isfinite(injected)):
        raise ValueError("center-pair selected logits are nonfinite")
    return clean, injected, selection


def paired_center_ranking_bootstrap(
    clean_scores: np.ndarray,
    injected_scores: np.ndarray,
    components: Sequence[str],
    *,
    replicates: int = 1_000,
    seed: int = 560205,
) -> dict[str, Any]:
    """Give each clean/injected component pair exactly one equal vote."""

    clean = np.asarray(clean_scores, dtype=np.float64)
    injected = np.asarray(injected_scores, dtype=np.float64)
    component_values = tuple(str(value) for value in components)
    if (
        clean.ndim != 1
        or injected.shape != clean.shape
        or len(component_values) != clean.size
        or len(set(component_values)) != clean.size
        or clean.size < 2
        or not np.all(np.isfinite(clean))
        or not np.all(np.isfinite(injected))
    ):
        raise ValueError("paired-center score arrays are invalid")
    delta = injected - clean
    votes = (delta > 0.0).astype(np.float64) + 0.5 * (delta == 0.0)
    n_replicates = _exact_int(replicates, label="bootstrap replicates")
    if n_replicates <= 0:
        raise ValueError("bootstrap replicates must be positive")
    rng = np.random.default_rng(int(seed))
    draws = np.empty(n_replicates, dtype=np.float64)
    for index in range(n_replicates):
        sampled = rng.integers(0, votes.size, size=votes.size)
        draws[index] = float(np.mean(votes[sampled]))
    return {
        "estimate": float(np.mean(votes)),
        "component_pair_bootstrap_95_interval": [
            float(np.quantile(draws, 0.025)),
            float(np.quantile(draws, 0.975)),
        ],
        "component_pairs": int(votes.size),
        "injected_higher_pairs": int(np.count_nonzero(delta > 0.0)),
        "tied_pairs": int(np.count_nonzero(delta == 0.0)),
        "injected_lower_pairs": int(np.count_nonzero(delta < 0.0)),
    }


def _control_result(
    fit: ProbeFit,
    tokens: np.ndarray,
    token_valid: np.ndarray,
    token_targets: np.ndarray,
    components: Sequence[str],
    pair_cells: Sequence[str],
    *,
    bootstrap_replicates: int,
    bootstrap_seed: int,
) -> dict[str, Any]:
    logits = score_shared_linear_cadence_logits(fit, tokens, token_valid)
    clean, injected, selection = oracle_support_center_pair_scores(
        logits,
        token_valid,
        token_targets,
        components,
    )
    cells = np.asarray(tuple(str(value) for value in pair_cells), dtype=object)
    if cells.shape != clean.shape:
        raise ValueError("center-pair event cells do not align with components")
    overall = paired_center_ranking_bootstrap(
        clean,
        injected,
        selection.components,
        replicates=bootstrap_replicates,
        seed=bootstrap_seed,
    )
    by_cell: dict[str, Any] = {}
    for duration, depth in EVENT_CELLS:
        name = f"duration_{duration}_depth_{depth:.2f}"
        selected = cells == name
        by_cell[name] = paired_center_ranking_bootstrap(
            clean[selected],
            injected[selected],
            np.asarray(selection.components, dtype=object)[selected],
            replicates=bootstrap_replicates,
            seed=bootstrap_seed + 1 + EVENT_CELLS.index((duration, depth)),
        )
    return {
        "context_cadences": int(tokens.shape[1]),
        "input_native_cadences": int(tokens.shape[1]),
        "emitted_cadence_logits_per_sample": int(logits.shape[1]),
        "one_logit_per_native_cadence": logits.shape == token_valid.shape,
        "temporal_downsampling": False,
        "representation_pooling": False,
        "off_support_used_for_fit_or_evaluation": False,
        "held_out_support_role": "center_cadence_truth_index_only_never_feature",
        "probe_objective_history": list(fit.objective_history),
        "probe_weight_l2": float(np.linalg.norm(fit.weight)),
        "probe_bias": float(fit.bias),
        "locked_development_test_center_pair_ranking": overall,
        "locked_development_test_center_pair_ranking_by_event_cell": by_cell,
    }


def summarize_event_transfer_mechanics(
    control_results: Mapping[str, Any],
    *,
    mechanics_config: Mapping[str, Any],
) -> dict[str, Any]:
    """Apply only the frozen raw/quality center-pair mechanics gates."""

    raw = _mapping(
        control_results["raw_adp_validity_error_cadence_128"],
        label="raw mechanics control",
    )
    quality = _mapping(
        control_results["quality_only_cadence_128"],
        label="quality mechanics control",
    )
    raw_overall = _mapping(
        raw["locked_development_test_center_pair_ranking"],
        label="raw center-pair ranking",
    )
    raw_lower = float(raw_overall["component_pair_bootstrap_95_interval"][0])
    raw_cells = _mapping(
        raw["locked_development_test_center_pair_ranking_by_event_cell"],
        label="raw center-pair cells",
    )
    deep_cells = {
        name: float(value["estimate"])
        for name, value in raw_cells.items()
        if str(name).endswith("_depth_0.30")
    }
    if len(deep_cells) != len(EVENT_DURATIONS):
        raise ValueError("raw center-pair mechanics lacks depth-0.30 cells")
    coefficients = tuple(float(value) for value in raw["raw_flux_coefficients"])
    if len(coefficients) != 2:
        raise ValueError("raw center-pair mechanics lacks ADP coefficients")
    quality_estimate = float(
        quality["locked_development_test_center_pair_ranking"]["estimate"]
    )
    criteria = {
        "raw_center_paired_ranking_lower_95_at_least_floor": raw_lower
        >= float(mechanics_config["raw_minimum_center_paired_ranking_lower_95"]),
        "raw_each_depth_0_30_center_paired_ranking_at_least_floor": all(
            value
            >= float(
                mechanics_config["raw_minimum_each_depth_0_30_center_paired_ranking"]
            )
            for value in deep_cells.values()
        ),
        "raw_adp_flux_coefficients_are_dip_positive": (
            all(value < 0.0 for value in coefficients)
            if mechanics_config["raw_require_negative_adp_flux_coefficients"]
            else True
        ),
        "quality_center_paired_ranking_is_exact_chance": quality_estimate
        == float(mechanics_config["quality_center_paired_ranking_required_exact"]),
    }
    return {
        "passed": all(criteria.values()),
        "criteria": criteria,
        "raw_center_paired_ranking": float(raw_overall["estimate"]),
        "raw_center_paired_ranking_lower_95": raw_lower,
        "raw_depth_0_30_center_paired_ranking": deep_cells,
        "raw_flux_coefficients": list(coefficients),
        "quality_center_paired_ranking": quality_estimate,
        "interpretation": (
            "oracle_support_center_probe_mechanics_valid"
            if all(criteria.values())
            else "oracle_support_center_probe_mechanics_failed"
        ),
    }


def evaluate_event_transfer_mechanics(
    *,
    config_path: str | Path,
    temporal_panel_dir: str | Path,
    temporal_panel_receipt_sha256: str,
) -> dict[str, Any]:
    """Run the fast raw/quality-only v3 mechanics control without an FM."""

    config, config_source, config_sha = load_event_transfer_mechanics_config(
        config_path
    )
    inputs = _mapping(config["inputs"], label="inputs")
    if temporal_panel_receipt_sha256 != inputs["temporal_panel_receipt_sha256"]:
        raise ValueError("event-transfer mechanics temporal-panel hash differs")
    rows, panel_receipt = load_temporal_panel(
        temporal_panel_dir,
        receipt_sha256=temporal_panel_receipt_sha256,
    )
    print("FM_EVENT_TRANSFER_MECHANICS phase=freeze_schedule", flush=True)
    schedule, schedule_audit = freeze_component_schedule(
        rows,
        variant=str(inputs["variant"]),
        target_per_cohort_block=int(
            config["split"]["target_components_per_cohort_per_block"]
        ),
    )

    probe_config = _mapping(config["probe"], label="probe")
    evaluation_config = _mapping(config["evaluation"], label="evaluation")
    train_subset = [item for item in schedule if item["block"] == "train"]
    train_samples, _labels, train_components, _cohorts, train_targets = (
        build_paired_samples_with_token_targets(
            train_subset,
            context_length=CONTEXT_CADENCES,
        )
    )
    fitted_controls: dict[str, ProbeFit] = {}
    for spec, quality_only in (
        ("raw_adp_validity_error_cadence_128", False),
        ("quality_only_cadence_128", True),
    ):
        print(f"FM_EVENT_TRANSFER_MECHANICS control={spec} phase=fit", flush=True)
        train_tokens, train_valid = raw_cadence_features(
            train_samples,
            quality_only=quality_only,
        )
        fitted_controls[spec] = fit_shared_linear_center_probe(
            train_tokens,
            train_valid,
            train_targets,
            train_components,
            learning_rate=float(probe_config["learning_rate"]),
            l2_weight=float(probe_config["l2_weight"]),
            epochs=int(probe_config["epochs"]),
            seed=int(probe_config["initialization_seed"]),
            progress_label=spec,
        )
    del train_samples

    test_subset = [
        item for item in schedule if item["block"] == "locked_development_test"
    ]
    test_samples, _labels, test_components, _cohorts, test_targets = (
        build_paired_samples_with_token_targets(
            test_subset,
            context_length=CONTEXT_CADENCES,
        )
    )
    test_cells = tuple(
        f"duration_{int(item['duration_cadences'])}_depth_{float(item['fractional_depth']):.2f}"
        for item in test_subset
    )
    control_results: dict[str, Any] = {}
    for spec, quality_only in (
        ("raw_adp_validity_error_cadence_128", False),
        ("quality_only_cadence_128", True),
    ):
        print(f"FM_EVENT_TRANSFER_MECHANICS control={spec} phase=test", flush=True)
        test_tokens, test_valid = raw_cadence_features(
            test_samples,
            quality_only=quality_only,
        )
        result = _control_result(
            fitted_controls[spec],
            test_tokens,
            test_valid,
            test_targets,
            test_components,
            test_cells,
            bootstrap_replicates=int(evaluation_config["bootstrap_replicates"]),
            bootstrap_seed=int(evaluation_config["bootstrap_seed"])
            + CONTROL_SPECS.index(spec) * 100,
        )
        if not quality_only:
            result["raw_flux_coefficients"] = [
                float(fitted_controls[spec].weight[0]),
                float(fitted_controls[spec].weight[1]),
            ]
        control_results[spec] = result
        print(f"FM_EVENT_TRANSFER_MECHANICS control={spec} phase=complete", flush=True)

    mechanics = summarize_event_transfer_mechanics(
        control_results,
        mechanics_config=_mapping(config["mechanics"], label="mechanics"),
    )
    schedule_public = []
    for item in schedule:
        public = {key: value for key, value in item.items() if key != "visit"}
        public["window_binding"] = dict(item["visit"]["binding"])
        schedule_public.append(public)
    return {
        "schema_version": EVENT_TRANSFER_MECHANICS_RESULT_SCHEMA_VERSION,
        "campaign_id": EVENT_TRANSFER_MECHANICS_CAMPAIGN_ID,
        "evaluation_completed": True,
        "passed": bool(mechanics["passed"]),
        "config_path": str(config_source),
        "config_sha256": config_sha,
        "temporal_panel_receipt_sha256": temporal_panel_receipt_sha256,
        "temporal_panel_rows": int(panel_receipt["n_panel_rows"]),
        "schedule_sha256": hashlib.sha256(
            json.dumps(schedule_public, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest(),
        "schedule": schedule_public,
        "schedule_audit": schedule_audit,
        "control_results": control_results,
        "probe_mechanics": mechanics,
        "cadence_preservation": {
            "nominal_cadence_seconds": 200,
            "one_logit_per_native_cadence": True,
            "patching_or_cadence_averaging": False,
            "representation_pooling": False,
            "held_out_reduction": (
                "exact_injected_support_center_vs_exact_paired_clean_cadence"
            ),
            "off_support_used_for_fit_or_evaluation": False,
            "support_used_as_input_feature": False,
        },
        "boundaries": {
            "fm_checkpoint_loaded": False,
            "fm_encoder_trained_or_evaluated": False,
            "fm_metrics_interpreted": False,
            "bls_or_period_features_used": False,
            "sealed_test_opened": False,
            "formal_model_gate": False,
            "development_only": True,
        },
    }


__all__ = [
    "CONTROL_SPECS",
    "EVENT_TRANSFER_MECHANICS_CAMPAIGN_ID",
    "EVENT_TRANSFER_MECHANICS_CONFIG_SCHEMA_VERSION",
    "EVENT_TRANSFER_MECHANICS_RESULT_SCHEMA_VERSION",
    "CenterPairSelection",
    "center_pair_examples",
    "evaluate_event_transfer_mechanics",
    "fit_shared_linear_center_probe",
    "load_event_transfer_mechanics_config",
    "oracle_support_center_pair_scores",
    "paired_center_ranking_bootstrap",
    "paired_support_center_indices",
    "score_shared_linear_cadence_logits",
    "summarize_event_transfer_mechanics",
]
