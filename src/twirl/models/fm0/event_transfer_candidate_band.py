"""Fixed-band, cadence-preserving FM event-transfer screen.

The v3 oracle-center artifact verifies only raw/quality probe mechanics.  This
module performs the actual candidate-centered development screen without an
oracle at held-out scoring time.  A shared linear map emits one scalar logit
for each of the 128 native 200-second cadence representations.  Only then is a
sample score formed by taking the maximum over the same predeclared zero-based
cadence band, 36--92 inclusive, for every sample.

Synthetic support may supervise the training-token loss.  It is never appended
to a feature, used to choose a held-out cadence, or passed to the fixed-band
reducer.  No cadence representation is pooled, averaged, patched, resampled,
or merged.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from .centered_event_context_diagnostic import _load_inference_only_trusted_model
from .event_transfer_canary import (
    EVENT_DEPTHS,
    EVENT_DURATIONS,
    JITTER_BOUNDS,
    MINIMUM_COMPONENTS,
    SPLIT_SECTORS,
    TARGET_COMPONENTS,
    ProbeFit,
    build_paired_samples_with_token_targets,
    encode_cadence_tokens,
    fit_shared_linear_token_probe,
    freeze_component_schedule,
    paired_component_bootstrap,
    raw_cadence_features,
    select_balanced_accuracy_threshold,
    select_fpr_threshold,
    summarize_probe_mechanics,
    summarize_readiness,
)
from .event_transfer_mechanics import (
    EVENT_TRANSFER_MECHANICS_CAMPAIGN_ID,
    EVENT_TRANSFER_MECHANICS_RESULT_SCHEMA_VERSION,
    score_shared_linear_cadence_logits,
)
from .temporal_panel import DEVELOPMENT_PARTITION
from .temporal_zero_shot import TEMPORAL_COHORTS, load_temporal_panel

EVENT_TRANSFER_BAND_CONFIG_SCHEMA_VERSION = (
    "twirl_fm0_event_transfer_candidate_band_config_v4"
)
EVENT_TRANSFER_BAND_RESULT_SCHEMA_VERSION = (
    "twirl_fm0_event_transfer_candidate_band_result_v4"
)
EVENT_TRANSFER_BAND_CAMPAIGN_ID = "twirl_fm0_2_s66_s77_event_transfer_candidate_band_v4"

CONTEXT_CADENCES = 128
CANDIDATE_BAND_START = 36
CANDIDATE_BAND_STOP_INCLUSIVE = 92
CANDIDATE_BAND_STOP_EXCLUSIVE = CANDIDATE_BAND_STOP_INCLUSIVE + 1
CANDIDATE_BAND_CADENCES = CANDIDATE_BAND_STOP_EXCLUSIVE - CANDIDATE_BAND_START
EXPECTED_SCHEDULE_SHA256 = (
    "6211057e0a6b8acb4377c11c2ccba92f2fc4a56e4a4f9c2e06f5afc0cef4cc37"
)
FEATURE_SPECS = (
    "raw_adp_validity_error_cadence_128",
    "quality_only_cadence_128",
    "step0_h_cadence_128",
    "step2000_h_cadence_128",
)


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


def fixed_candidate_band_mask(length: int = CONTEXT_CADENCES) -> np.ndarray:
    """Return the one frozen band; the result never depends on event support."""

    cadence_count = _exact_int(length, label="candidate-band context length")
    if cadence_count != CONTEXT_CADENCES:
        raise ValueError("candidate-band screen requires exactly 128 cadences")
    mask = np.zeros(cadence_count, dtype=bool)
    mask[CANDIDATE_BAND_START:CANDIDATE_BAND_STOP_EXCLUSIVE] = True
    if int(np.count_nonzero(mask)) != CANDIDATE_BAND_CADENCES:
        raise RuntimeError("candidate-band geometry is inconsistent")
    return mask


def fixed_band_training_validity(
    token_valid: np.ndarray,
    token_targets: np.ndarray,
) -> np.ndarray:
    """Restrict token-supervised fitting to the fixed band, fail closed."""

    valid = np.asarray(token_valid, dtype=bool)
    targets = np.asarray(token_targets, dtype=bool)
    if valid.ndim != 2 or valid.shape[1] != CONTEXT_CADENCES:
        raise ValueError("candidate-band validity must be [sample, 128]")
    if targets.shape != valid.shape:
        raise ValueError("candidate-band training targets must match validity")
    if valid.shape[0] < 2 or valid.shape[0] % 2:
        raise ValueError("candidate-band training requires adjacent sample pairs")
    band = fixed_candidate_band_mask()
    if np.any(targets[:, ~band]):
        raise ValueError("synthetic support escaped the predeclared fixed band")
    if not np.all(valid[:, band]):
        raise ValueError("every candidate-band cadence must be valid")
    if not np.array_equal(valid[0::2], valid[1::2]):
        raise ValueError("clean/injected validity masks differ")
    fit_valid = np.zeros_like(valid)
    fit_valid[:, band] = True
    return fit_valid


def validate_fixed_band_scoring_validity(token_valid: np.ndarray) -> None:
    """Validate a held-out band without accepting event truth or support."""

    valid = np.asarray(token_valid, dtype=bool)
    if valid.ndim != 2 or valid.shape[1] != CONTEXT_CADENCES:
        raise ValueError("candidate-band validity must be [sample, 128]")
    if valid.shape[0] < 2 or valid.shape[0] % 2:
        raise ValueError("candidate-band scoring requires adjacent sample pairs")
    band = fixed_candidate_band_mask()
    if not np.all(valid[:, band]):
        raise ValueError("every held-out candidate-band cadence must be valid")
    if not np.array_equal(valid[0::2], valid[1::2]):
        raise ValueError("held-out clean/injected validity masks differ")


def reduce_fixed_candidate_band(
    cadence_logits: np.ndarray,
    token_valid: np.ndarray,
) -> np.ndarray:
    """Reduce scalar logits over indices 36--92, with no support argument."""

    logits = np.asarray(cadence_logits, dtype=np.float64)
    valid = np.asarray(token_valid, dtype=bool)
    if (
        logits.ndim != 2
        or logits.shape != valid.shape
        or logits.shape[1] != CONTEXT_CADENCES
        or not np.all(np.isfinite(logits))
    ):
        raise ValueError("candidate-band cadence logits must be finite [sample, 128]")
    band = fixed_candidate_band_mask()
    if not np.all(valid[:, band]):
        raise ValueError("held-out fixed candidate band is not fully valid")
    return np.max(logits[:, band], axis=1)


def score_fixed_candidate_band(
    fit: ProbeFit,
    tokens: np.ndarray,
    token_valid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Score all 128 native tokens, then reduce only the fixed scalar band."""

    logits = score_shared_linear_cadence_logits(fit, tokens, token_valid)
    if logits.shape[1] != CONTEXT_CADENCES:
        raise ValueError("candidate-band scorer did not emit 128 cadence logits")
    scores = reduce_fixed_candidate_band(logits, token_valid)
    return scores, logits


def validate_mechanics_prerequisite(
    result_path: str | Path,
    *,
    expected_sha256: str,
) -> dict[str, Any]:
    """Verify the exact passed v3 mechanics artifact before FM scoring."""

    source = Path(result_path).expanduser().resolve(strict=True)
    payload_bytes = source.read_bytes()
    observed_sha = hashlib.sha256(payload_bytes).hexdigest()
    if observed_sha != _require_sha256(expected_sha256, label="v3 mechanics result"):
        raise ValueError("v3 mechanics result SHA-256 differs")
    payload = _mapping(json.loads(payload_bytes), label="v3 mechanics result")
    if (
        payload.get("schema_version") != EVENT_TRANSFER_MECHANICS_RESULT_SCHEMA_VERSION
        or payload.get("campaign_id") != EVENT_TRANSFER_MECHANICS_CAMPAIGN_ID
        or payload.get("evaluation_completed") is not True
        or payload.get("passed") is not True
        or _mapping(payload.get("probe_mechanics"), label="v3 probe mechanics").get(
            "passed"
        )
        is not True
    ):
        raise ValueError("v3 mechanics prerequisite did not pass its exact contract")
    cadence = _mapping(
        payload.get("cadence_preservation"), label="v3 cadence preservation"
    )
    expected_cadence = {
        "nominal_cadence_seconds": 200,
        "one_logit_per_native_cadence": True,
        "patching_or_cadence_averaging": False,
        "representation_pooling": False,
        "held_out_reduction": (
            "exact_injected_support_center_vs_exact_paired_clean_cadence"
        ),
        "off_support_used_for_fit_or_evaluation": False,
        "support_used_as_input_feature": False,
    }
    for key, expected in expected_cadence.items():
        if cadence.get(key) != expected:
            raise ValueError(f"v3 mechanics cadence field {key!r} differs")
    boundaries = _mapping(payload.get("boundaries"), label="v3 boundaries")
    for key in (
        "fm_checkpoint_loaded",
        "fm_encoder_trained_or_evaluated",
        "fm_metrics_interpreted",
        "bls_or_period_features_used",
        "sealed_test_opened",
        "formal_model_gate",
    ):
        if boundaries.get(key) is not False:
            raise ValueError(f"v3 mechanics boundary {key!r} differs")
    return {
        "path": str(source),
        "sha256": observed_sha,
        "schema_version": str(payload["schema_version"]),
        "campaign_id": str(payload["campaign_id"]),
        "passed": True,
    }


def load_event_transfer_band_config(
    config_path: str | Path,
) -> tuple[dict[str, Any], Path, str]:
    """Validate and byte-bind the exact fixed-band v4 screen contract."""

    source = Path(config_path).expanduser().resolve(strict=True)
    payload_bytes = source.read_bytes()
    config = dict(
        _mapping(yaml.safe_load(payload_bytes), label="candidate-band config")
    )
    if (
        config.get("schema_version") != EVENT_TRANSFER_BAND_CONFIG_SCHEMA_VERSION
        or config.get("campaign_id") != EVENT_TRANSFER_BAND_CAMPAIGN_ID
        or config.get("model_family") != "TWIRL-FM0"
        or config.get("status") != "frozen_bls_free_development_candidate_band_screen"
    ):
        raise ValueError("candidate-band config identity differs")

    purpose = _mapping(config.get("purpose"), label="purpose")
    if (
        purpose.get("primary")
        != "test_fixed_band_cadence_local_information_for_later_classifiers"
        or purpose.get("task_geometry")
        != "candidate_centered_single_event_classification"
    ):
        raise ValueError("candidate-band purpose differs")
    for key, expected in (
        ("probe_only", True),
        ("formal_model_gate", False),
        ("production_model_claim", False),
        ("architecture_selection_claim", False),
    ):
        _exact_bool(purpose, key, expected)

    authorization = _mapping(config.get("authorization"), label="authorization")
    for key, expected in (
        ("development_shard_payload_access", True),
        ("artificial_single_event_injection", True),
        ("frozen_checkpoint_inference", True),
        ("low_capacity_probe_fitting", True),
        ("fm_encoder_optimizer_or_training", False),
        ("bls_period_candidate_or_event_features", False),
        ("sealed_test_access", False),
    ):
        _exact_bool(authorization, key, expected)

    prerequisite = _mapping(config.get("prerequisite"), label="prerequisite")
    _exact_bool(prerequisite, "v3_mechanics_pass_required", True)
    _require_sha256(
        prerequisite.get("v3_mechanics_result_sha256"),
        label="v3 mechanics result",
    )
    _exact_bool(prerequisite, "v1_v2_fm_metrics_interpretable", False)

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
        or tuple(inputs.get("checkpoint_steps", ())) != (0, 2_000)
        or inputs.get("expected_schedule_sha256") != EXPECTED_SCHEDULE_SHA256
    ):
        raise ValueError("candidate-band input identity differs")
    checkpoint_hashes = _mapping(
        inputs.get("checkpoint_sha256"), label="checkpoint hashes"
    )
    _require_sha256(checkpoint_hashes.get("step0"), label="step-0 checkpoint")
    _require_sha256(checkpoint_hashes.get("step2000"), label="step-2000 checkpoint")

    split = _mapping(config.get("split"), label="split")
    for name, sectors in SPLIT_SECTORS.items():
        if tuple(split.get(f"{name}_sectors", ())) != sectors:
            raise ValueError(f"candidate-band {name} sectors differ")
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
        raise ValueError("candidate-band split differs")
    _exact_bool(split, "one_visit_per_component", True)

    temporal = _mapping(config.get("temporal_resolution"), label="temporal resolution")
    if (
        _exact_int(temporal.get("nominal_cadence_seconds"), label="nominal cadence")
        != 200
        or _exact_int(temporal.get("context_cadences"), label="context cadences")
        != CONTEXT_CADENCES
        or _exact_int(
            temporal.get("fixed_band_start_index_zero_based"), label="band start"
        )
        != CANDIDATE_BAND_START
        or _exact_int(
            temporal.get("fixed_band_stop_index_inclusive_zero_based"),
            label="band stop",
        )
        != CANDIDATE_BAND_STOP_INCLUSIVE
        or _exact_int(temporal.get("fixed_band_cadences"), label="band size")
        != CANDIDATE_BAND_CADENCES
        or temporal.get("fixed_band_definition")
        != "union_of_every_allowed_jittered_event_support"
        or temporal.get("held_out_sample_reduction")
        != "maximum_of_scalar_cadence_logits_over_fixed_indices_36_through_92"
    ):
        raise ValueError("candidate-band cadence grid differs")
    for key, expected in (
        ("direct_centered_crop", True),
        ("one_encoder_token_per_native_cadence", True),
        ("one_logit_per_native_cadence", True),
        ("score_all_128_tokens_before_reduction", True),
        ("temporal_downsampling", False),
        ("patching_or_cadence_averaging", False),
        ("representation_pooling", False),
        ("fixed_band_same_for_every_sample", True),
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
        != "per_cadence_training_target_only_never_feature_or_held_out_router"
    ):
        raise ValueError("candidate-band event definition differs")
    for key, expected in (
        ("one_clean_and_one_injected_sample_per_component", True),
        ("one_event_per_injected_sample", True),
        ("period_defined", False),
        ("injection_truth_used_as_feature", False),
        ("held_out_support_used_to_select_score", False),
    ):
        _exact_bool(event, key, expected)

    features = _mapping(config.get("features"), label="features")
    if tuple(features.get("arms", ())) != FEATURE_SPECS:
        raise ValueError("candidate-band feature arms differ")
    if features.get("primary") != "step2000_h_cadence_128":
        raise ValueError("candidate-band primary feature differs")
    if tuple(features.get("forbidden_arms", ())) != (
        "step0_h_cadence_2048",
        "step2000_h_cadence_2048",
        "step2000_z_cadence_128_diagnostic",
    ):
        raise ValueError("candidate-band forbidden arms differ")
    if tuple(features.get("forbidden_feature_families", ())) != (
        "period",
        "bls",
        "candidate_score",
        "event_duration",
        "event_depth",
        "event_center_numeric_value",
        "injection_support",
    ):
        raise ValueError("candidate-band forbidden features differ")

    probe = _mapping(config.get("probe"), label="probe")
    if (
        probe.get("family") != "shared_linear_cadence_scorer"
        or probe.get("training_objective")
        != "pair_balanced_four_stratum_per_cadence_binary_cross_entropy"
        or probe.get("training_target")
        != "synthetic_event_support_boolean_target_not_input_feature"
        or probe.get("standardization")
        != "training_fixed_band_valid_tokens_dimensionwise_center_and_scale"
        or probe.get("held_out_sample_score")
        != "fixed_band_max_after_all_128_scalar_cadence_logits_exist"
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
        raise ValueError("candidate-band probe differs")
    for key, expected in (
        ("support_appended_to_features", False),
        ("held_out_support_used_as_router", False),
        ("representation_pooling_before_score", False),
        ("locked_test_not_used_for_fit_or_threshold", True),
    ):
        _exact_bool(probe, key, expected)

    metrics = _mapping(config.get("metrics"), label="metrics")
    if (
        metrics.get("primary") != "paired_component_ranking_accuracy"
        or tuple(metrics.get("secondary", ()))
        != (
            "roc_auc",
            "average_precision",
            "balanced_accuracy",
            "macro_f1",
            "brier_score",
            "expected_calibration_error",
            "tpr_at_validation_frozen_5_percent_fpr",
            "fpr_at_validation_frozen_5_percent_fpr",
        )
        or not math.isclose(
            float(metrics.get("fpr_operating_point", math.nan)),
            0.05,
            rel_tol=0.0,
            abs_tol=0.0,
        )
        or _exact_int(metrics.get("calibration_bins"), label="calibration bins") != 10
        or _exact_int(metrics.get("bootstrap_replicates"), label="bootstrap replicates")
        != 1_000
        or _exact_int(metrics.get("bootstrap_seed"), label="bootstrap seed") != 560204
    ):
        raise ValueError("candidate-band metrics differ")
    for key in ("report_by_cohort", "report_by_event_cell"):
        _exact_bool(metrics, key, True)

    mechanics = _mapping(config.get("probe_mechanics"), label="probe mechanics")
    if dict(mechanics) != {
        "raw_control_minimum_overall_roc_auc_lower_95": 0.75,
        "raw_control_minimum_each_depth_0_30_roc_auc": 0.90,
        "raw_control_require_negative_adp_flux_coefficients": True,
        "raw_control_paired_ranking_lower_95_strictly_above_chance": True,
        "quality_only_maximum_absolute_roc_auc_from_chance": 0.05,
        "quality_only_paired_ranking_required_exact": 0.50,
    }:
        raise ValueError("candidate-band mechanics gate differs")

    readiness = _mapping(config.get("readiness"), label="readiness")
    if (
        readiness.get("primary_feature") != "step2000_h_cadence_128"
        or tuple(readiness.get("paired_delta_controls", ())) != ("step0_h_cadence_128",)
        or readiness.get("raw_noninferiority_control")
        != "raw_adp_validity_error_cadence_128"
    ):
        raise ValueError("candidate-band readiness identity differs")
    expected_readiness = {
        "minimum_roc_auc_lower_95": 0.75,
        "minimum_paired_roc_auc_delta": 0.02,
        "paired_delta_lower_95_strictly_positive": True,
        "raw_noninferiority_margin": 0.01,
        "minimum_each_cohort_roc_auc": 0.70,
        "minimum_tpr_at_5_percent_fpr": 0.30,
        "minimum_tpr_at_5_percent_fpr_lower_95": 0.20,
        "maximum_locked_test_fpr_at_validation_frozen_threshold": 0.10,
        "require_paired_ranking_lower_95_above_chance": True,
        "raw_noninferiority_supports_readiness_not_superiority": True,
    }
    for key, expected in expected_readiness.items():
        if readiness.get(key) != expected:
            raise ValueError(f"candidate-band readiness field {key!r} differs")
    return config, source, hashlib.sha256(payload_bytes).hexdigest()


def _feature_tokens(
    spec: str,
    *,
    samples: Sequence[Mapping[str, np.ndarray]],
    model0: Any,
    model2000: Any,
    device: Any,
    batch_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    if spec == "raw_adp_validity_error_cadence_128":
        return raw_cadence_features(samples, quality_only=False)
    if spec == "quality_only_cadence_128":
        return raw_cadence_features(samples, quality_only=True)
    if spec == "step0_h_cadence_128":
        model = model0
    elif spec == "step2000_h_cadence_128":
        model = model2000
    else:  # pragma: no cover - guarded by the frozen feature tuple
        raise ValueError(f"unknown candidate-band feature arm: {spec}")
    return encode_cadence_tokens(
        model,
        samples,
        context_length=CONTEXT_CADENCES,
        device=device,
        batch_size=batch_size,
        project_each_cadence=False,
    )


def _band_mechanics(
    feature_results: Mapping[str, Any],
    *,
    mechanics_config: Mapping[str, Any],
) -> dict[str, Any]:
    base_config = dict(mechanics_config)
    required_quality_ranking = float(
        base_config.pop("quality_only_paired_ranking_required_exact")
    )
    result = summarize_probe_mechanics(
        feature_results,
        mechanics_config=base_config,
    )
    quality_ranking = float(
        feature_results["quality_only_cadence_128"]["locked_development_test"][
            "paired_component_ranking_accuracy"
        ]["estimate"]
    )
    exact = quality_ranking == required_quality_ranking
    criteria = dict(result["criteria"])
    criteria["quality_only_paired_ranking_exactly_chance"] = exact
    result.update(
        {
            "passed": all(criteria.values()),
            "criteria": criteria,
            "quality_only_paired_ranking": quality_ranking,
            "fixed_candidate_band_zero_based_inclusive": [
                CANDIDATE_BAND_START,
                CANDIDATE_BAND_STOP_INCLUSIVE,
            ],
            "held_out_support_used_as_router": False,
            "interpretation": (
                "fixed_band_probe_mechanics_valid"
                if all(criteria.values())
                else "fixed_band_probe_mechanics_failed"
            ),
        }
    )
    return result


def evaluate_event_transfer_candidate_band(
    *,
    config_path: str | Path,
    mechanics_result_path: str | Path,
    run_dir: str | Path,
    step0_checkpoint_path: str | Path,
    step2000_checkpoint_path: str | Path,
    temporal_panel_dir: str | Path,
    temporal_panel_receipt_sha256: str,
    batch_size: int = 8,
) -> dict[str, Any]:
    """Run the exact development-only fixed candidate-band screen."""

    try:
        import torch
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("PyTorch is required for candidate-band inference") from exc

    config, config_source, config_sha = load_event_transfer_band_config(config_path)
    prerequisite_config = _mapping(config["prerequisite"], label="prerequisite")
    prerequisite = validate_mechanics_prerequisite(
        mechanics_result_path,
        expected_sha256=str(prerequisite_config["v3_mechanics_result_sha256"]),
    )
    inputs = _mapping(config["inputs"], label="inputs")
    if temporal_panel_receipt_sha256 != inputs["temporal_panel_receipt_sha256"]:
        raise ValueError("candidate-band temporal-panel receipt hash differs")
    rows, panel_receipt = load_temporal_panel(
        temporal_panel_dir,
        receipt_sha256=temporal_panel_receipt_sha256,
    )
    device = torch.device("cpu")
    model0, contract0, validation0 = _load_inference_only_trusted_model(
        Path(run_dir),
        device=device,
        checkpoint_path=Path(step0_checkpoint_path),
    )
    model2000, contract2000, validation2000 = _load_inference_only_trusted_model(
        Path(run_dir),
        device=device,
        checkpoint_path=Path(step2000_checkpoint_path),
    )
    if (
        validation0["global_step"] != 0
        or validation2000["global_step"] != 2_000
        or contract0 != contract2000
    ):
        raise ValueError("candidate-band checkpoint pair differs")
    hashes = _mapping(inputs["checkpoint_sha256"], label="checkpoint hashes")
    if (
        validation0["selected_checkpoint_sha256"] != hashes["step0"]
        or validation2000["selected_checkpoint_sha256"] != hashes["step2000"]
    ):
        raise ValueError("candidate-band checkpoint hash differs")

    print("FM_EVENT_TRANSFER_BAND phase=freeze_schedule", flush=True)
    schedule, schedule_audit = freeze_component_schedule(
        rows,
        variant=str(inputs["variant"]),
    )
    schedule_public = []
    for item in schedule:
        public = {key: value for key, value in item.items() if key != "visit"}
        public["window_binding"] = dict(item["visit"]["binding"])
        schedule_public.append(public)
    schedule_sha = hashlib.sha256(
        json.dumps(schedule_public, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    if schedule_sha != inputs["expected_schedule_sha256"]:
        raise ValueError("candidate-band frozen schedule SHA-256 differs")
    print(
        f"FM_EVENT_TRANSFER_BAND phase=schedule_ready components={len(schedule)}",
        flush=True,
    )

    block_data: dict[
        str,
        tuple[
            list[dict[str, np.ndarray]],
            np.ndarray,
            tuple[str, ...],
            tuple[str, ...],
            np.ndarray,
        ],
    ] = {}
    block_cells: dict[str, tuple[str, ...]] = {}
    for block in SPLIT_SECTORS:
        subset = [item for item in schedule if item["block"] == block]
        block_data[block] = build_paired_samples_with_token_targets(
            subset,
            context_length=CONTEXT_CADENCES,
        )
        block_cells[block] = tuple(
            f"duration_{int(item['duration_cadences'])}_depth_{float(item['fractional_depth']):.2f}"
            for item in subset
            for _ in (0, 1)
        )

    probe_config = _mapping(config["probe"], label="probe")
    metric_config = _mapping(config["metrics"], label="metrics")
    results: dict[str, Any] = {}
    test_scores_by_spec: dict[str, np.ndarray] = {}
    test_reference: tuple[np.ndarray, tuple[str, ...], tuple[str, ...]] | None = None
    for spec in FEATURE_SPECS:
        print(f"FM_EVENT_TRANSFER_BAND feature={spec} phase=encode", flush=True)
        encoded: dict[
            str,
            tuple[
                np.ndarray,
                np.ndarray,
                np.ndarray,
                tuple[str, ...],
                tuple[str, ...],
                np.ndarray,
            ],
        ] = {}
        for block in SPLIT_SECTORS:
            samples, labels, components, cohorts, token_targets = block_data[block]
            tokens, valid = _feature_tokens(
                spec,
                samples=samples,
                model0=model0,
                model2000=model2000,
                device=device,
                batch_size=batch_size,
            )
            if tokens.shape[:2] != (len(samples), CONTEXT_CADENCES):
                raise ValueError("candidate-band feature lost native cadence geometry")
            if block == "train":
                fixed_band_training_validity(valid, token_targets)
            else:
                validate_fixed_band_scoring_validity(valid)
            encoded[block] = (
                tokens,
                valid,
                labels,
                components,
                cohorts,
                token_targets,
            )
            print(
                f"FM_EVENT_TRANSFER_BAND feature={spec} block={block} phase=encoded",
                flush=True,
            )

        train = encoded["train"]
        fit_valid = fixed_band_training_validity(train[1], train[5])
        print(f"FM_EVENT_TRANSFER_BAND feature={spec} phase=fit", flush=True)
        fit = fit_shared_linear_token_probe(
            train[0],
            fit_valid,
            train[5],
            learning_rate=float(probe_config["learning_rate"]),
            l2_weight=float(probe_config["l2_weight"]),
            epochs=int(probe_config["epochs"]),
            seed=int(probe_config["initialization_seed"]),
            progress_label=spec,
        )

        validation = encoded["validation"]
        print(f"FM_EVENT_TRANSFER_BAND feature={spec} phase=threshold", flush=True)
        validation_scores, validation_logits = score_fixed_candidate_band(
            fit,
            validation[0],
            validation[1],
        )
        threshold = select_balanced_accuracy_threshold(validation[2], validation_scores)
        fpr_threshold = select_fpr_threshold(
            validation[2],
            validation_scores,
            maximum_fpr=float(metric_config["fpr_operating_point"]),
        )

        test = encoded["locked_development_test"]
        test_scores, test_logits = score_fixed_candidate_band(
            fit,
            test[0],
            test[1],
        )
        test_scores_by_spec[spec] = test_scores
        current_reference = (test[2], test[3], test[4])
        if test_reference is None:
            test_reference = current_reference
        elif (
            not np.array_equal(test_reference[0], current_reference[0])
            or test_reference[1] != current_reference[1]
            or test_reference[2] != current_reference[2]
        ):
            raise ValueError("candidate-band feature arms changed test identities")

        overall = paired_component_bootstrap(
            test[2],
            test_scores,
            test[3],
            threshold=threshold,
            fpr_threshold=fpr_threshold,
            calibration_bins=int(metric_config["calibration_bins"]),
            replicates=int(metric_config["bootstrap_replicates"]),
            seed=int(metric_config["bootstrap_seed"]) + FEATURE_SPECS.index(spec),
        )
        cohort_values = np.asarray(test[4], dtype=object)
        by_cohort: dict[str, Any] = {}
        for cohort in TEMPORAL_COHORTS:
            selected = cohort_values == cohort
            by_cohort[cohort] = paired_component_bootstrap(
                test[2][selected],
                test_scores[selected],
                np.asarray(test[3], dtype=object)[selected],
                threshold=threshold,
                fpr_threshold=fpr_threshold,
                calibration_bins=int(metric_config["calibration_bins"]),
                replicates=int(metric_config["bootstrap_replicates"]),
                seed=int(metric_config["bootstrap_seed"])
                + 100
                + FEATURE_SPECS.index(spec),
            )

        cell_values = np.asarray(block_cells["locked_development_test"], dtype=object)
        by_cell: dict[str, Any] = {}
        by_cohort_and_cell: dict[str, Any] = {}
        for cell in dict.fromkeys(cell_values.tolist()):
            selected = cell_values == cell
            by_cell[cell] = paired_component_bootstrap(
                test[2][selected],
                test_scores[selected],
                np.asarray(test[3], dtype=object)[selected],
                threshold=threshold,
                fpr_threshold=fpr_threshold,
                calibration_bins=int(metric_config["calibration_bins"]),
                replicates=int(metric_config["bootstrap_replicates"]),
                seed=int(metric_config["bootstrap_seed"])
                + 200
                + FEATURE_SPECS.index(spec),
            )
            for cohort in TEMPORAL_COHORTS:
                joint = selected & (cohort_values == cohort)
                by_cohort_and_cell[f"{cohort}:{cell}"] = paired_component_bootstrap(
                    test[2][joint],
                    test_scores[joint],
                    np.asarray(test[3], dtype=object)[joint],
                    threshold=threshold,
                    fpr_threshold=fpr_threshold,
                    calibration_bins=int(metric_config["calibration_bins"]),
                    replicates=int(metric_config["bootstrap_replicates"]),
                    seed=int(metric_config["bootstrap_seed"])
                    + 300
                    + FEATURE_SPECS.index(spec),
                )

        validation_prediction = validation_scores >= fpr_threshold
        results[spec] = {
            "context_cadences": CONTEXT_CADENCES,
            "input_native_cadences": CONTEXT_CADENCES,
            "emitted_scalar_cadence_logits_per_sample": int(test_logits.shape[1]),
            "fixed_band_zero_based_inclusive": [
                CANDIDATE_BAND_START,
                CANDIDATE_BAND_STOP_INCLUSIVE,
            ],
            "fixed_band_cadences": CANDIDATE_BAND_CADENCES,
            "temporal_downsampling": False,
            "representation_pooling": False,
            "cadence_representations_averaged_or_merged": False,
            "held_out_support_used_as_router": False,
            "held_out_sample_reduction": "max_of_fixed_band_scalar_logits",
            "probe_training_objective": probe_config["training_objective"],
            "validation_threshold": float(threshold),
            "validation_frozen_5_percent_fpr_threshold": float(fpr_threshold),
            "validation_operating_fpr": float(
                np.mean(validation_prediction[validation[2] == 0])
            ),
            "validation_operating_tpr": float(
                np.mean(validation_prediction[validation[2] == 1])
            ),
            "probe_objective_history": list(fit.objective_history),
            "probe_weight_l2": float(np.linalg.norm(fit.weight)),
            "probe_bias": float(fit.bias),
            "locked_development_test": overall,
            "locked_development_test_by_cohort": by_cohort,
            "locked_development_test_by_event_cell": by_cell,
            "locked_development_test_by_cohort_and_event_cell": by_cohort_and_cell,
            "validation_logits_shape": list(validation_logits.shape),
            "locked_development_test_logits_shape": list(test_logits.shape),
        }
        if spec == "raw_adp_validity_error_cadence_128":
            results[spec]["raw_flux_coefficients"] = [
                float(fit.weight[0]),
                float(fit.weight[1]),
            ]
        print(f"FM_EVENT_TRANSFER_BAND feature={spec} phase=complete", flush=True)

    if test_reference is None:
        raise ValueError("candidate-band screen produced no held-out scores")
    mechanics = _band_mechanics(
        results,
        mechanics_config=_mapping(config["probe_mechanics"], label="probe mechanics"),
    )
    readiness = summarize_readiness(
        results,
        test_scores_by_spec,
        test_reference[0],
        test_reference[1],
        readiness_config=_mapping(config["readiness"], label="readiness"),
        bootstrap_replicates=int(metric_config["bootstrap_replicates"]),
        bootstrap_seed=int(metric_config["bootstrap_seed"]),
    )
    metric_criteria_passed = bool(readiness["ready_for_next_real_training"])
    overall_passed = bool(mechanics["passed"] and metric_criteria_passed)
    readiness = dict(readiness)
    readiness.update(
        {
            "fm_metric_criteria_satisfied": metric_criteria_passed,
            "fixed_band_probe_mechanics_gate_required": True,
            "ready_for_next_real_training": overall_passed,
            "useful_representation_claim_supported": False,
            "architecture_selection_supported": False,
            "claim_limit": "synthetic_candidate_band_screen_only",
        }
    )

    return {
        "schema_version": EVENT_TRANSFER_BAND_RESULT_SCHEMA_VERSION,
        "campaign_id": EVENT_TRANSFER_BAND_CAMPAIGN_ID,
        "evaluation_completed": True,
        "passed": overall_passed,
        "probe_mechanics_passed": bool(mechanics["passed"]),
        "fm_metric_criteria_passed": metric_criteria_passed,
        "scientific_readiness_passed": overall_passed,
        "architecture_selection_supported": False,
        "config_path": str(config_source),
        "config_sha256": config_sha,
        "mechanics_prerequisite": prerequisite,
        "temporal_panel_receipt_sha256": temporal_panel_receipt_sha256,
        "temporal_panel_rows": int(panel_receipt["n_panel_rows"]),
        "schedule_sha256": schedule_sha,
        "schedule": schedule_public,
        "schedule_audit": schedule_audit,
        "checkpoint_validation": {"step0": validation0, "step2000": validation2000},
        "feature_results": results,
        "probe_mechanics": mechanics,
        "readiness": readiness,
        "cadence_preservation": {
            "nominal_cadence_seconds": 200,
            "input_native_cadences_per_sample": CONTEXT_CADENCES,
            "emitted_scalar_logits_per_sample": CONTEXT_CADENCES,
            "one_encoder_token_per_input_cadence": True,
            "one_logit_per_native_cadence": True,
            "temporal_downsampling": False,
            "patching_or_cadence_averaging": False,
            "representation_pooling": False,
            "fixed_band_zero_based_inclusive": [
                CANDIDATE_BAND_START,
                CANDIDATE_BAND_STOP_INCLUSIVE,
            ],
            "only_temporal_reduction": "max_of_predeclared_fixed_band_scalar_logits",
            "held_out_support_used_as_router": False,
            "support_used_as_input_feature": False,
        },
        "boundaries": {
            "fm_encoder_trained": False,
            "bls_or_period_features_used": False,
            "sealed_test_opened": False,
            "formal_model_gate": False,
            "development_only": True,
            "candidate_centered_classification_not_blind_detection": True,
            "unknown_location_transfer_established": False,
            "natural_event_performance_established": False,
            "architecture_selected": False,
            "locked_development_test_used_for_fit_or_threshold": False,
        },
    }


__all__ = [
    "CANDIDATE_BAND_CADENCES",
    "CANDIDATE_BAND_START",
    "CANDIDATE_BAND_STOP_INCLUSIVE",
    "CONTEXT_CADENCES",
    "EVENT_TRANSFER_BAND_CAMPAIGN_ID",
    "EVENT_TRANSFER_BAND_CONFIG_SCHEMA_VERSION",
    "EVENT_TRANSFER_BAND_RESULT_SCHEMA_VERSION",
    "EXPECTED_SCHEDULE_SHA256",
    "FEATURE_SPECS",
    "evaluate_event_transfer_candidate_band",
    "fixed_band_training_validity",
    "fixed_candidate_band_mask",
    "load_event_transfer_band_config",
    "reduce_fixed_candidate_band",
    "score_fixed_candidate_band",
    "validate_fixed_band_scoring_validity",
    "validate_mechanics_prerequisite",
]
