"""BLS-free, cadence-preserving event-transfer canary for TWIRL-FM0.

The probe scores every 200-second cadence independently with one shared linear
map and applies a masked maximum only after cadence scoring.  It never patches,
averages, bins, or otherwise merges cadence tokens.  The FM encoder remains
frozen; only the small downstream probe is fitted.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import defaultdict
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from .centered_event_context_diagnostic import (
    BASE_INTERVAL_CADENCES,
    CenteredEventIneligibleError,
    _load_centered_visit,
    _load_inference_only_trusted_model,
    centered_context_bounds,
    centered_trapezoid,
    slice_centered_context,
)
from .dataset import VIEW_NAMES, collate_fm0_samples, move_batch_to_device
from .temporal_panel import DEVELOPMENT_PARTITION
from .temporal_zero_shot import TEMPORAL_COHORTS, load_temporal_panel

EVENT_TRANSFER_CONFIG_SCHEMA_VERSION = "twirl_fm0_event_transfer_canary_config_v1"
EVENT_TRANSFER_CONFIG_SCHEMA_VERSION_V2 = "twirl_fm0_event_transfer_canary_config_v2"
EVENT_TRANSFER_RESULT_SCHEMA_VERSION = "twirl_fm0_event_transfer_canary_result_v1"
EVENT_TRANSFER_RESULT_SCHEMA_VERSION_V2 = "twirl_fm0_event_transfer_canary_result_v2"
EVENT_TRANSFER_CAMPAIGN_ID = "twirl_fm0_2_s66_s77_event_transfer_canary_v1"
EVENT_TRANSFER_CAMPAIGN_ID_V2 = "twirl_fm0_2_s66_s77_event_transfer_canary_v2"
EVENT_TRANSFER_SALT = "twirl_fm0_event_transfer_canary_v1"
SPLIT_SECTORS: Mapping[str, tuple[int, ...]] = {
    "train": tuple(range(66, 72)),
    "validation": (72, 73, 74),
    "locked_development_test": (75, 76, 77),
}
CONTEXT_LENGTHS = (128, 2_048)
EVENT_DURATIONS = (1, 3, 9)
EVENT_DEPTHS = (0.01, 0.03, 0.1, 0.3)
EVENT_CELLS = tuple(
    (duration, depth) for duration in EVENT_DURATIONS for depth in EVENT_DEPTHS
)
TARGET_COMPONENTS = 240
MINIMUM_COMPONENTS = len(EVENT_CELLS)
JITTER_BOUNDS = (-24, 24)
FEATURE_SPECS = (
    "step2000_h_cadence_128",
    "step0_h_cadence_128",
    "step0_h_cadence_2048",
    "step2000_h_cadence_2048",
    "step2000_z_cadence_128_diagnostic",
    "raw_adp_validity_error_cadence_128",
    "quality_only_cadence_128",
)


@dataclass(frozen=True)
class ProbeFit:
    """Fitted shared cadence scorer and training-only standardization."""

    weight: np.ndarray
    bias: float
    center: np.ndarray
    scale: np.ndarray
    objective_history: tuple[float, ...]


def _seed(*parts: object) -> int:
    payload = "\x1f".join(str(part) for part in parts).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


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


def load_event_transfer_config(
    config_path: str | Path,
) -> tuple[dict[str, Any], Path, str]:
    """Validate and byte-bind the exact development canary contract."""

    source = Path(config_path).expanduser().resolve(strict=True)
    payload_bytes = source.read_bytes()
    config = dict(
        _mapping(yaml.safe_load(payload_bytes), label="event-transfer config")
    )
    schema_version = str(config.get("schema_version", ""))
    if schema_version not in {
        EVENT_TRANSFER_CONFIG_SCHEMA_VERSION,
        EVENT_TRANSFER_CONFIG_SCHEMA_VERSION_V2,
    }:
        raise ValueError("event-transfer config schema differs")
    token_supervised = schema_version == EVENT_TRANSFER_CONFIG_SCHEMA_VERSION_V2
    expected_campaign = (
        EVENT_TRANSFER_CAMPAIGN_ID_V2 if token_supervised else EVENT_TRANSFER_CAMPAIGN_ID
    )
    if config.get("campaign_id") != expected_campaign:
        raise ValueError("event-transfer campaign differs")
    if (
        config.get("model_family") != "TWIRL-FM0"
        or config.get("status") != "frozen_bls_free_development_probe"
    ):
        raise ValueError("event-transfer config identity differs")

    purpose = _mapping(config.get("purpose"), label="purpose")
    if (
        purpose.get("primary")
        != "test_cadence_local_single_event_information_for_later_classifiers"
        or purpose.get("task_geometry")
        != "candidate_centered_single_event_classification"
    ):
        raise ValueError("event-transfer purpose differs")
    for key, expected in (
        ("probe_only", True),
        ("formal_model_gate", False),
        ("production_model_claim", False),
    ):
        _exact_bool(purpose, key, expected)
    if token_supervised and purpose.get("evidence_scope") != (
        "injection_support_supervised_synthetic_event_localization_and_linear_decodability"
    ):
        raise ValueError("event-transfer v2 evidence scope differs")

    authorization = _mapping(config.get("authorization"), label="authorization")
    for key in (
        "development_shard_payload_access",
        "artificial_single_event_injection",
        "frozen_checkpoint_inference",
        "low_capacity_probe_fitting",
    ):
        _exact_bool(authorization, key, True)
    for key in (
        "fm_encoder_optimizer_or_training",
        "bls_period_candidate_or_event_features",
        "sealed_test_access",
    ):
        _exact_bool(authorization, key, False)

    inputs = _mapping(config.get("inputs"), label="inputs")
    if tuple(inputs.get("allowed_sectors", ())) != tuple(range(66, 78)) or inputs.get(
        "excluded_sectors"
    ) != [65]:
        raise ValueError("event-transfer sector scope differs")
    if (
        inputs.get("source_partition") != DEVELOPMENT_PARTITION
        or inputs.get("variant") != "TWIRL-FM0.2.1"
    ):
        raise ValueError("event-transfer input identity differs")
    if tuple(inputs.get("checkpoint_steps", ())) != (0, 2_000):
        raise ValueError("event-transfer checkpoint steps differ")
    checkpoint_hashes = _mapping(
        inputs.get("checkpoint_sha256"), label="checkpoint hashes"
    )
    if set(checkpoint_hashes) != {"step0", "step2000"}:
        raise ValueError("event-transfer checkpoint hash keys differ")
    for value in checkpoint_hashes.values():
        if len(str(value)) != 64 or any(
            character not in "0123456789abcdef" for character in str(value)
        ):
            raise ValueError("event-transfer checkpoint SHA-256 differs")

    split = _mapping(config.get("split"), label="split")
    for name, sectors in SPLIT_SECTORS.items():
        key = f"{name}_sectors"
        if tuple(split.get(key, ())) != sectors:
            raise ValueError(f"event-transfer {name} sectors differ")
    if (
        split.get("component_spanning_blocks") != "quarantine"
        or split.get("shortfall_policy") != "fail_closed"
        or split.get("component_and_visit_order") != "deterministic_sha256"
    ):
        raise ValueError("event-transfer split policy differs")
    if tuple(split.get("cohorts", ())) != TEMPORAL_COHORTS:
        raise ValueError("event-transfer cohorts differ")
    for key, expected in (("one_visit_per_component", True),):
        _exact_bool(split, key, expected)
    if (
        _exact_int(
            split.get("target_components_per_cohort_per_block"),
            label="target components",
        )
        != TARGET_COMPONENTS
    ):
        raise ValueError("event-transfer target count differs")
    if (
        _exact_int(
            split.get("minimum_components_per_cohort_per_block"),
            label="minimum components",
        )
        != MINIMUM_COMPONENTS
    ):
        raise ValueError("event-transfer minimum count differs")

    temporal = _mapping(config.get("temporal_resolution"), label="temporal resolution")
    if (
        _exact_int(temporal.get("nominal_cadence_seconds"), label="nominal cadence")
        != 200
        or tuple(temporal.get("context_cadences", ())) != CONTEXT_LENGTHS
    ):
        raise ValueError("event-transfer cadence grid differs")
    for key in (
        "direct_nested_crops",
        "one_encoder_token_per_input_cadence",
        "checkpoint_time_normalization_remains_2048",
    ):
        _exact_bool(temporal, key, True)
    for key in (
        "temporal_downsampling",
        "patching_or_cadence_averaging",
        "pooling_before_token_score",
    ):
        _exact_bool(temporal, key, False)
    if (
        _exact_int(
            temporal.get("conformer_patch_stride_required"),
            label="Conformer patch stride",
        )
        != 1
    ):
        raise ValueError("event-transfer Conformer stride differs")
    if (
        temporal.get("token_aggregation")
        != "masked_max_after_shared_linear_cadence_score"
    ):
        raise ValueError("event-transfer cadence aggregation differs")

    event = _mapping(config.get("artificial_event"), label="artificial event")
    if (
        event.get("profile") != "symmetric_trapezoid_sampled_at_cadence_centers"
        or not math.isclose(
            float(event.get("ingress_egress_fraction_each", math.nan)),
            1.0 / 3.0,
            rel_tol=0.0,
            abs_tol=0.0,
        )
        or tuple(event.get("duration_cadences", ())) != EVENT_DURATIONS
        or tuple(float(value) for value in event.get("fractional_depths", ()))
        != EVENT_DEPTHS
        or _exact_int(event.get("balanced_cells"), label="balanced cells")
        != len(EVENT_CELLS)
        or event.get("apply_to_views") != ["adp_1x1", "adp_3x3"]
    ):
        raise ValueError("event-transfer event grid differs")
    if tuple(event.get("center_jitter_cadences_inclusive", ())) != JITTER_BOUNDS:
        raise ValueError("event-transfer jitter differs")
    for key in (
        "one_clean_and_one_injected_sample_per_component",
        "one_event_per_injected_sample",
        "window_centered_on_synthetic_event",
        "classification_not_blind_detection",
    ):
        _exact_bool(event, key, True)
    for key in ("period_defined", "injection_truth_used_as_feature"):
        _exact_bool(event, key, False)
    if token_supervised:
        _exact_bool(event, "injection_support_used_as_training_target", True)

    features = _mapping(config.get("features"), label="features")
    if (features.get("primary"), *tuple(features.get("controls", ()))) != FEATURE_SPECS:
        raise ValueError("event-transfer feature set differs")
    forbidden = tuple(features.get("forbidden_feature_families", ()))
    if forbidden != (
        "period",
        "bls",
        "candidate_score",
        "event_duration",
        "event_depth",
        "event_center_numeric_value",
        "injection_support",
    ):
        raise ValueError("event-transfer forbidden features differ")
    if features.get("z_cadence_definition") != (
        "apply_frozen_checkpoint_projection_to_each_cadence_independently"
    ) or tuple(features.get("raw_channels", ())) != (
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
    ):
        raise ValueError("event-transfer cadence feature definition differs")
    if tuple(features.get("quality_only_channels", ())) != (
        "flux_valid_1x1",
        "flux_valid_3x3",
        "error_valid_1x1",
        "error_valid_3x3",
        "local_time_cadences",
        "delta_time_cadences",
        "time_valid",
        "segment_boundary",
    ):
        raise ValueError("event-transfer quality-only feature definition differs")

    probe = _mapping(config.get("probe"), label="probe")
    if (
        probe.get("family") != "shared_linear_cadence_scorer"
        or probe.get("sample_logit") != "masked_max_of_per_cadence_logits"
    ):
        raise ValueError("event-transfer probe differs")
    _exact_bool(probe, "temporal_pooling_before_scoring", False)
    _exact_bool(probe, "locked_test_not_used_for_probe_fit_or_threshold", True)
    expected_training_objective = (
        "pair_balanced_four_stratum_per_cadence_binary_cross_entropy"
        if token_supervised
        else "sample_binary_cross_entropy_after_hard_masked_max"
    )
    if probe.get("training_objective", expected_training_objective) != (
        expected_training_objective
    ):
        raise ValueError("event-transfer probe training objective differs")
    if token_supervised and probe.get("training_target") != (
        "synthetic_event_support_boolean_target_not_input_feature"
    ):
        raise ValueError("event-transfer probe training target differs")
    if (
        probe.get("standardization")
        != "training_valid_tokens_dimensionwise_center_and_scale_no_temporal_resampling"
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
        or _exact_int(probe.get("epochs"), label="probe epochs") != 400
        or _exact_int(probe.get("initialization_seed"), label="probe seed") != 560203
        or probe.get("threshold_selection") != "validation_balanced_accuracy"
    ):
        raise ValueError("event-transfer probe epochs differ")

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
        or metrics.get("uncertainty")
        != "deterministic_paired_component_bootstrap_95_interval"
    ):
        raise ValueError("event-transfer metrics differ")
    if (
        _exact_int(metrics.get("bootstrap_replicates"), label="bootstrap replicates")
        != 1_000
    ):
        raise ValueError("event-transfer bootstrap count differs")
    if (
        _exact_int(metrics.get("bootstrap_seed"), label="bootstrap seed") != 560204
        or metrics.get("report_by_cohort") is not True
        or metrics.get("report_by_event_cell") is not True
        or _exact_int(metrics.get("calibration_bins"), label="calibration bins") != 10
        or not math.isclose(
            float(metrics.get("fpr_operating_point", math.nan)),
            0.05,
            rel_tol=0.0,
            abs_tol=0.0,
        )
    ):
        raise ValueError("event-transfer bootstrap policy differs")

    readiness = _mapping(config.get("readiness"), label="readiness")
    expected_readiness = {
        "primary_feature": "step2000_h_cadence_128",
        "minimum_roc_auc_lower_95": 0.75,
        "minimum_paired_roc_auc_delta": 0.02,
        "paired_delta_controls": [
            "step0_h_cadence_128",
            "step2000_h_cadence_2048",
        ],
        "paired_delta_lower_95_strictly_positive": True,
        "raw_noninferiority_control": "raw_adp_validity_error_cadence_128",
        "raw_noninferiority_margin": 0.01,
        "minimum_each_cohort_roc_auc": 0.70,
        "minimum_tpr_at_5_percent_fpr": 0.30,
        "minimum_tpr_at_5_percent_fpr_lower_95": 0.20,
        "maximum_locked_test_fpr_at_validation_frozen_threshold": 0.10,
        "require_paired_ranking_lower_95_above_chance": True,
        "raw_noninferiority_supports_readiness_not_superiority": True,
    }
    if dict(readiness) != expected_readiness:
        raise ValueError("event-transfer readiness rule differs")

    if token_supervised:
        mechanics = _mapping(config.get("probe_mechanics"), label="probe mechanics")
        if dict(mechanics) != {
            "raw_control_minimum_overall_roc_auc_lower_95": 0.75,
            "raw_control_minimum_each_depth_0_30_roc_auc": 0.90,
            "raw_control_require_negative_adp_flux_coefficients": True,
            "raw_control_paired_ranking_lower_95_strictly_above_chance": True,
            "quality_only_maximum_absolute_roc_auc_from_chance": 0.05,
        }:
            raise ValueError("event-transfer probe mechanics rule differs")
    return config, source, hashlib.sha256(payload_bytes).hexdigest()


def sector_block(sector: int) -> str:
    """Return the frozen chronological block for one allowed sector."""

    value = _exact_int(sector, label="sector")
    matches = [name for name, sectors in SPLIT_SECTORS.items() if value in sectors]
    if len(matches) != 1:
        raise ValueError(f"sector S{value} is outside the event-transfer split")
    return matches[0]


def _order_key(*parts: object) -> tuple[str, str]:
    identity = str(parts[-1])
    digest = hashlib.sha256(
        "\x1f".join((EVENT_TRANSFER_SALT, *(str(part) for part in parts))).encode()
    ).hexdigest()
    return digest, identity


def freeze_component_schedule(
    rows: Sequence[Mapping[str, str]],
    *,
    variant: str = "TWIRL-FM0.2.1",
    visit_loader: Callable[[Mapping[str, str]], Mapping[str, Any]] | None = None,
    target_per_cohort_block: int = TARGET_COMPONENTS,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Freeze one eligible visit/component after quarantining cross-block IDs."""

    target = _exact_int(target_per_cohort_block, label="target per cohort/block")
    if target < MINIMUM_COMPONENTS or target % len(EVENT_CELLS):
        raise ValueError("target cannot fit one balanced event-cell cycle")
    loader = visit_loader or (lambda row: _load_centered_visit(row, variant=variant))
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    observations: set[str] = set()
    component_cohort: dict[str, str] = {}
    for source_row in rows:
        row = dict(source_row)
        component = str(row.get("leakage_component_id", ""))
        observation = str(row.get("observation_key", ""))
        cohort = str(row.get("temporal_cohort", ""))
        if (
            not component
            or not observation
            or cohort not in TEMPORAL_COHORTS
            or row.get("source_partition") != DEVELOPMENT_PARTITION
        ):
            raise ValueError("event-transfer schedule received an ineligible row")
        sector_block(int(row.get("sector", 0)))
        if observation in observations:
            raise ValueError("event-transfer schedule received duplicate observations")
        observations.add(observation)
        prior = component_cohort.setdefault(component, cohort)
        if prior != cohort:
            raise ValueError("event-transfer component crosses cohorts")
        grouped[component].append(row)

    quarantined = {
        component
        for component, component_rows in grouped.items()
        if len({sector_block(int(row["sector"])) for row in component_rows}) != 1
    }
    selected: list[dict[str, Any]] = []
    audit: dict[str, Any] = {
        "quarantined_cross_block_components": len(quarantined),
        "quarantined_components_sha256": hashlib.sha256(
            "\n".join(sorted(quarantined)).encode("utf-8")
        ).hexdigest(),
        "blocks": {},
    }
    for block in SPLIT_SECTORS:
        audit["blocks"][block] = {}
        for cohort in TEMPORAL_COHORTS:
            eligible_components = [
                component
                for component, component_rows in grouped.items()
                if component not in quarantined
                and component_cohort[component] == cohort
                and sector_block(int(component_rows[0]["sector"])) == block
            ]
            eligible_components.sort(
                key=lambda value: _order_key("component", block, cohort, value)
            )
            if len(eligible_components) < target:
                raise ValueError(
                    f"event-transfer {block}/{cohort} has only "
                    f"{len(eligible_components)} nonquarantined components; "
                    f"the exact target is {target}"
                )
            desired = target
            retained = 0
            opened = 0
            ineligible = 0
            for component in eligible_components:
                if retained >= desired:
                    break
                chosen: tuple[dict[str, str], Mapping[str, Any]] | None = None
                for row in sorted(
                    grouped[component],
                    key=lambda item: _order_key("visit", item["observation_key"]),
                ):
                    present = json.loads(str(row.get("view_present_json", "")))
                    if (
                        not isinstance(present, list)
                        or len(present) != len(VIEW_NAMES)
                        or not all(
                            type(value) is int and value in {0, 1} for value in present
                        )
                    ):
                        raise ValueError("event-transfer view-present schema differs")
                    if not bool(present[2]) or not bool(present[3]):
                        continue
                    opened += 1
                    try:
                        visit = loader(row)
                        _require_jitter_support(visit["sample"])
                    except CenteredEventIneligibleError:
                        ineligible += 1
                        continue
                    chosen = row, visit
                    break
                if chosen is None:
                    continue
                cell_index = retained % len(EVENT_CELLS)
                duration, depth = EVENT_CELLS[cell_index]
                jitter = (
                    _seed("jitter", block, cohort, component)
                    % (JITTER_BOUNDS[1] - JITTER_BOUNDS[0] + 1)
                    + JITTER_BOUNDS[0]
                )
                row, visit = chosen
                selected.append(
                    {
                        "block": block,
                        "cohort": cohort,
                        "component_id": component,
                        "observation_key": str(row["observation_key"]),
                        "sector": int(row["sector"]),
                        "duration_cadences": duration,
                        "fractional_depth": depth,
                        "jitter_cadences": int(jitter),
                        "visit": visit,
                    }
                )
                retained += 1
            if retained != desired:
                raise ValueError(
                    f"event-transfer {block}/{cohort} yielded {retained} of {desired} frozen components"
                )
            counts = defaultdict(int)
            for item in selected:
                if item["block"] == block and item["cohort"] == cohort:
                    counts[(item["duration_cadences"], item["fractional_depth"])] += 1
            if len(set(counts.values())) != 1 or set(counts) != set(EVENT_CELLS):
                raise ValueError("event-transfer event schedule is not cell-balanced")
            audit["blocks"][block][cohort] = {
                "available_nonquarantined_components": len(eligible_components),
                "selected_components": retained,
                "target_components": target,
                "shortfall_components": target - retained,
                "screening_visits_opened": opened,
                "structurally_ineligible_visits": ineligible,
            }
    return selected, audit


def _require_jitter_support(sample: Mapping[str, np.ndarray]) -> None:
    start128, _ = centered_context_bounds(128)
    half = max(EVENT_DURATIONS) // 2
    start = start128 + 128 // 2 + JITTER_BOUNDS[0] - half
    stop = start128 + 128 // 2 + JITTER_BOUNDS[1] + half + 1
    time_valid = np.asarray(sample["time_valid"], dtype=bool)
    flux_valid = np.asarray(sample["flux_valid"], dtype=bool)
    if (
        time_valid.shape != (BASE_INTERVAL_CADENCES,)
        or flux_valid.shape[-1] != BASE_INTERVAL_CADENCES
    ):
        raise ValueError("event-transfer visit does not retain the 2048-cadence base")
    if not np.all(time_valid[start:stop]) or not np.all(flux_valid[:, start:stop]):
        raise CenteredEventIneligibleError(
            "visit lacks full validity over jittered event support"
        )


def inject_jittered_single_event(
    sample: Mapping[str, np.ndarray],
    *,
    duration_cadences: int,
    fractional_depth: float,
    jitter_cadences: int,
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """Inject one event at a deterministic offset without changing other tensors."""

    length = int(np.asarray(sample["time_valid"]).size)
    duration = _exact_int(duration_cadences, label="event duration")
    jitter = _exact_int(jitter_cadences, label="event jitter")
    if duration not in EVENT_DURATIONS or not (
        JITTER_BOUNDS[0] <= jitter <= JITTER_BOUNDS[1]
    ):
        raise ValueError("event support differs from the frozen schedule")
    center = length // 2 + jitter
    start = center - duration // 2
    stop = start + duration
    if start < 0 or stop > length:
        raise ValueError("jittered event escaped its context")
    profile = centered_trapezoid(duration, float(fractional_depth))
    valid = np.asarray(sample["time_valid"], dtype=bool)[None, :] & np.asarray(
        sample["flux_valid"], dtype=bool
    )
    if not np.all(valid[:, start:stop]):
        raise CenteredEventIneligibleError("jittered event support is not fully valid")
    injected = {key: np.asarray(value).copy() for key, value in sample.items()}
    original = np.asarray(sample["flux"])
    injected["flux"][:, start:stop] = (
        (1.0 + original[:, start:stop].astype(np.float64)) * (1.0 - profile[None, :])
        - 1.0
    ).astype(original.dtype)
    support = np.zeros(length, dtype=bool)
    support[start:stop] = True
    for key in sample:
        if key != "flux" and not np.array_equal(injected[key], np.asarray(sample[key])):
            raise ValueError(f"event injection changed non-flux tensor {key}")
    return injected, support


def build_paired_samples_with_token_targets(
    schedule: Sequence[Mapping[str, Any]], *, context_length: int
) -> tuple[
    list[dict[str, np.ndarray]],
    np.ndarray,
    tuple[str, ...],
    tuple[str, ...],
    np.ndarray,
]:
    """Materialize exact pairs plus training-only event-cadence targets.

    The target mask is never appended to a model or raw-control feature.  It is
    retained separately so the v2 probe can learn a shared cadence scorer
    without using a temporally pooled training objective.
    """

    if context_length not in CONTEXT_LENGTHS:
        raise ValueError("context length differs from the frozen grid")
    samples: list[dict[str, np.ndarray]] = []
    labels: list[int] = []
    components: list[str] = []
    cohorts: list[str] = []
    token_targets: list[np.ndarray] = []
    for item in schedule:
        clean = slice_centered_context(
            item["visit"]["sample"], context_length=context_length
        )
        injected, support = inject_jittered_single_event(
            clean,
            duration_cadences=int(item["duration_cadences"]),
            fractional_depth=float(item["fractional_depth"]),
            jitter_cadences=int(item["jitter_cadences"]),
        )
        for sample, label, target in (
            (clean, 0, np.zeros(context_length, dtype=bool)),
            (injected, 1, support),
        ):
            samples.append(sample)
            labels.append(label)
            components.append(str(item["component_id"]))
            cohorts.append(str(item["cohort"]))
            token_targets.append(np.asarray(target, dtype=bool))
    return (
        samples,
        np.asarray(labels, dtype=np.int8),
        tuple(components),
        tuple(cohorts),
        np.stack(token_targets),
    )


def build_paired_samples(
    schedule: Sequence[Mapping[str, Any]], *, context_length: int
) -> tuple[list[dict[str, np.ndarray]], np.ndarray, tuple[str, ...], tuple[str, ...]]:
    """Materialize clean/injected pairs using direct cadence-preserving crops."""

    samples, labels, components, cohorts, _targets = (
        build_paired_samples_with_token_targets(
            schedule,
            context_length=context_length,
        )
    )
    return samples, labels, components, cohorts


def _validate_token_matrix(
    tokens: np.ndarray, token_valid: np.ndarray, labels: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray]:
    values = np.asarray(tokens, dtype=np.float32)
    valid = np.asarray(token_valid, dtype=bool)
    if (
        values.ndim != 3
        or valid.shape != values.shape[:2]
        or not np.all(np.isfinite(values))
    ):
        raise ValueError("cadence token matrix is invalid")
    if np.any(np.count_nonzero(valid, axis=1) == 0):
        raise ValueError("every sample requires at least one valid cadence token")
    if labels is not None:
        target = np.asarray(labels)
        if target.shape != (values.shape[0],) or not np.all(np.isin(target, (0, 1))):
            raise ValueError("probe labels must be binary and sample-aligned")
    return values, valid


def fit_shared_linear_max_probe(
    tokens: np.ndarray,
    token_valid: np.ndarray,
    labels: np.ndarray,
    *,
    learning_rate: float = 0.02,
    l2_weight: float = 0.001,
    epochs: int = 400,
    seed: int = 560203,
    progress_label: str | None = None,
) -> ProbeFit:
    """Fit a shared linear cadence scorer with hard masked-max aggregation."""

    values, valid = _validate_token_matrix(tokens, token_valid, labels)
    target = np.asarray(labels, dtype=np.float32)
    observed = values[valid]
    center = np.mean(observed, axis=0)
    scale = np.std(observed, axis=0)
    scale = np.where(np.isfinite(scale) & (scale > 1.0e-8), scale, 1.0)
    standardized = (values - center[None, None, :]) / scale[None, None, :]
    rng = np.random.default_rng(int(seed))
    weight = rng.normal(0.0, 1.0e-3, size=values.shape[2]).astype(np.float32)
    bias = 0.0
    m_w = np.zeros_like(weight)
    v_w = np.zeros_like(weight)
    m_b = v_b = 0.0
    history: list[float] = []
    beta1, beta2 = 0.9, 0.999
    for step in range(1, _exact_int(epochs, label="probe epochs") + 1):
        token_logits = standardized @ weight + bias
        token_logits = np.where(valid, token_logits, -np.inf)
        winning = np.argmax(token_logits, axis=1)
        chosen = standardized[np.arange(values.shape[0]), winning]
        logits = token_logits[np.arange(values.shape[0]), winning]
        probabilities = _sigmoid(logits)
        residual = probabilities - target
        grad_w = np.mean(residual[:, None] * chosen, axis=0) + float(l2_weight) * weight
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
        if step == 1 or step % 20 == 0 or step == epochs:
            loss = float(
                np.mean(np.logaddexp(0.0, logits) - target * logits)
                + 0.5 * float(l2_weight) * np.dot(weight, weight)
            )
            history.append(loss)
        if progress_label is not None and (
            step == 1 or step % 100 == 0 or step == epochs
        ):
            print(
                f"FM_EVENT_TRANSFER feature={progress_label} "
                f"phase=fit epoch={step}/{epochs}",
                flush=True,
            )
    return ProbeFit(
        weight=weight,
        bias=float(bias),
        center=center,
        scale=scale,
        objective_history=tuple(history),
    )


def paired_token_probe_examples(
    token_valid: np.ndarray,
    token_targets: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return four-stratum targets and exact per-pair example weights."""

    valid = np.asarray(token_valid, dtype=bool)
    if valid.ndim != 2:
        raise ValueError("cadence probe validity must be two-dimensional")
    targets = np.asarray(token_targets, dtype=bool)
    if targets.shape != valid.shape:
        raise ValueError("cadence probe targets must match token validity")
    if np.any(targets & ~valid):
        raise ValueError("cadence probe has a positive target on an invalid token")
    if valid.shape[0] % 2:
        raise ValueError("cadence probe requires adjacent clean/injected pairs")
    pair_count = valid.shape[0] // 2
    if np.any(targets[0::2]) or np.any(~np.any(targets[1::2], axis=1)):
        raise ValueError("cadence probe target rows do not follow clean/injected pairs")

    example_weights = np.zeros(valid.shape, dtype=np.float64)
    example_targets = np.zeros(valid.shape, dtype=np.float32)
    for pair_index in range(pair_count):
        clean_index = 2 * pair_index
        injected_index = clean_index + 1
        support = targets[injected_index]
        clean_support = support & valid[clean_index]
        injected_support = support & valid[injected_index]
        clean_off = ~support & valid[clean_index]
        injected_off = ~support & valid[injected_index]
        strata = (
            (injected_index, injected_support, 1.0),
            (clean_index, clean_support, 0.0),
            (injected_index, injected_off, 0.0),
            (clean_index, clean_off, 0.0),
        )
        for row_index, mask, target_value in strata:
            count = int(np.count_nonzero(mask))
            if count == 0:
                raise ValueError("cadence probe encountered an empty pair stratum")
            example_weights[row_index, mask] = 1.0 / (4.0 * pair_count * count)
            example_targets[row_index, mask] = target_value
    if not math.isclose(
        float(np.sum(example_weights)),
        1.0,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    ):
        raise ValueError("cadence probe stratum weights do not sum to one")
    return example_targets, example_weights


def fit_shared_linear_token_probe(
    tokens: np.ndarray,
    token_valid: np.ndarray,
    token_targets: np.ndarray,
    *,
    learning_rate: float = 0.02,
    l2_weight: float = 0.001,
    epochs: int = 400,
    seed: int = 560203,
    progress_label: str | None = None,
) -> ProbeFit:
    """Fit one shared scorer from paired cadence labels without pooling.

    Each clean/injected pair contributes equal weight to four strata:
    injected support, paired-clean support, injected off-support, and clean
    off-support.  This prevents event duration or the 128-cadence background
    from diluting the useful gradient.  The scalar loss reduces independently
    scored cadence examples; no representation is averaged across time.
    ``token_targets`` is training truth only and is never supplied to the
    fitted scorer or to :func:`score_shared_linear_max_probe`.
    """

    values, valid = _validate_token_matrix(tokens, token_valid)
    example_targets, example_weights = paired_token_probe_examples(
        valid,
        token_targets,
    )

    observed = values[valid]
    center = np.mean(observed, axis=0)
    scale = np.std(observed, axis=0)
    scale = np.where(np.isfinite(scale) & (scale > 1.0e-8), scale, 1.0)
    standardized = (values - center[None, None, :]) / scale[None, None, :]
    selected = example_weights > 0.0
    training_values = standardized[selected]
    training_targets = example_targets[selected]
    training_weights = example_weights[selected]

    rng = np.random.default_rng(int(seed))
    weight = rng.normal(0.0, 1.0e-3, size=values.shape[2]).astype(np.float32)
    bias = 0.0
    m_w = np.zeros_like(weight)
    v_w = np.zeros_like(weight)
    m_b = v_b = 0.0
    history: list[float] = []
    beta1, beta2 = 0.9, 0.999
    for step in range(1, _exact_int(epochs, label="probe epochs") + 1):
        logits = training_values @ weight + bias
        residual = _sigmoid(logits) - training_targets
        weighted_residual = training_weights * residual
        grad_w = np.sum(weighted_residual[:, None] * training_values, axis=0)
        grad_w += float(l2_weight) * weight
        grad_b = float(np.sum(weighted_residual))
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
        if step == 1 or step % 20 == 0 or step == epochs:
            loss = float(
                np.sum(
                    training_weights
                    * (
                        np.logaddexp(0.0, logits)
                        - training_targets * logits
                    )
                )
                + 0.5 * float(l2_weight) * np.dot(weight, weight)
            )
            history.append(loss)
        if progress_label is not None and (
            step == 1 or step % 100 == 0 or step == epochs
        ):
            print(
                f"FM_EVENT_TRANSFER feature={progress_label} "
                f"phase=fit_token_supervised epoch={step}/{epochs}",
                flush=True,
            )
    return ProbeFit(
        weight=weight,
        bias=float(bias),
        center=center,
        scale=scale,
        objective_history=tuple(history),
    )


def score_shared_linear_max_probe(
    fit: ProbeFit, tokens: np.ndarray, token_valid: np.ndarray
) -> np.ndarray:
    """Score each cadence separately, then take exactly one masked maximum."""

    values, valid = _validate_token_matrix(tokens, token_valid)
    if (
        fit.weight.shape != (values.shape[2],)
        or fit.center.shape != fit.weight.shape
        or fit.scale.shape != fit.weight.shape
    ):
        raise ValueError("probe dimensions differ from cadence tokens")
    standardized = (values - fit.center[None, None, :]) / fit.scale[None, None, :]
    per_cadence = standardized @ fit.weight + float(fit.bias)
    return np.max(np.where(valid, per_cadence, -np.inf), axis=1)


def raw_cadence_features(
    samples: Sequence[Mapping[str, np.ndarray]], *, quality_only: bool
) -> tuple[np.ndarray, np.ndarray]:
    """Return direct per-cadence controls without temporal resampling."""

    rows: list[np.ndarray] = []
    masks: list[np.ndarray] = []
    for sample in samples:
        flux = np.asarray(sample["flux"], dtype=np.float64).T
        flux_valid = np.asarray(sample["flux_valid"], dtype=bool).T
        error = np.asarray(sample["flux_error"], dtype=np.float64).T
        error_valid = np.asarray(sample["error_valid"], dtype=bool).T
        local_time = np.asarray(sample["local_time_cadences"], dtype=np.float64)
        delta_time = np.asarray(sample["delta_time_cadences"], dtype=np.float64)
        time_valid = np.asarray(sample["time_valid"], dtype=bool)
        boundary = np.asarray(sample["segment_boundary"], dtype=bool)
        time_channels = np.column_stack(
            (local_time, delta_time, time_valid.astype(float), boundary.astype(float))
        )
        if quality_only:
            feature = np.column_stack(
                (
                    flux_valid.astype(float),
                    error_valid.astype(float),
                    time_channels,
                )
            )
        else:
            feature = np.column_stack(
                (
                    flux,
                    flux_valid.astype(float),
                    error,
                    error_valid.astype(float),
                    time_channels,
                )
            )
        if feature.shape[0] != time_valid.size:
            raise ValueError("raw cadence feature length differs")
        rows.append(feature)
        masks.append(time_valid)
    return np.stack(rows), np.stack(masks)


def _rank_auc(labels: np.ndarray, scores: np.ndarray) -> float:
    y = np.asarray(labels, dtype=int)
    s = np.asarray(scores, dtype=float)
    n_positive = int(np.count_nonzero(y == 1))
    n_negative = int(np.count_nonzero(y == 0))
    if not n_positive or not n_negative:
        return float("nan")
    order = np.argsort(s, kind="mergesort")
    sorted_scores = s[order]
    ranks = np.empty(s.size, dtype=np.float64)
    start = 0
    while start < s.size:
        stop = start + 1
        while stop < s.size and sorted_scores[stop] == sorted_scores[start]:
            stop += 1
        ranks[order[start:stop]] = 0.5 * ((start + 1) + stop)
        start = stop
    rank_sum = float(np.sum(ranks[y == 1]))
    return (rank_sum - n_positive * (n_positive + 1) / 2.0) / float(
        n_positive * n_negative
    )


def _average_precision(labels: np.ndarray, scores: np.ndarray) -> float:
    order = np.argsort(-np.asarray(scores), kind="mergesort")
    y = np.asarray(labels, dtype=int)[order]
    sorted_scores = np.asarray(scores, dtype=float)[order]
    positives = int(np.count_nonzero(y == 1))
    if positives == 0:
        return float("nan")
    cumulative_positive = 0
    cumulative_rows = 0
    average_precision = 0.0
    start = 0
    while start < y.size:
        stop = start + 1
        while stop < y.size and sorted_scores[stop] == sorted_scores[start]:
            stop += 1
        group_positive = int(np.count_nonzero(y[start:stop] == 1))
        cumulative_positive += group_positive
        cumulative_rows += stop - start
        average_precision += (group_positive / positives) * (
            cumulative_positive / cumulative_rows
        )
        start = stop
    return float(average_precision)


def _sigmoid(values: np.ndarray) -> np.ndarray:
    logits = np.asarray(values, dtype=np.float64)
    result = np.empty_like(logits)
    nonnegative = logits >= 0.0
    result[nonnegative] = 1.0 / (1.0 + np.exp(-logits[nonnegative]))
    exponential = np.exp(logits[~nonnegative])
    result[~nonnegative] = exponential / (1.0 + exponential)
    return result


def select_balanced_accuracy_threshold(labels: np.ndarray, scores: np.ndarray) -> float:
    """Freeze a threshold using validation data only."""

    y = np.asarray(labels, dtype=int)
    s = np.asarray(scores, dtype=float)
    candidates = np.concatenate(
        (
            [np.nextafter(np.min(s), -np.inf)],
            np.unique(s),
            [np.nextafter(np.max(s), np.inf)],
        )
    )
    best = (-np.inf, 0.0)
    for threshold in candidates:
        prediction = s >= threshold
        tpr = np.mean(prediction[y == 1])
        tnr = np.mean(~prediction[y == 0])
        balanced = 0.5 * (tpr + tnr)
        candidate = (float(balanced), -float(threshold))
        best = max(best, candidate)
    return -best[1]


def select_fpr_threshold(
    labels: np.ndarray,
    scores: np.ndarray,
    *,
    maximum_fpr: float = 0.05,
) -> float:
    """Freeze the most sensitive validation threshold at or below an FPR cap."""

    y = np.asarray(labels, dtype=int)
    s = np.asarray(scores, dtype=float)
    if y.shape != s.shape or not np.any(y == 0) or not np.any(y == 1):
        raise ValueError("FPR threshold requires aligned positive and negative rows")
    cap = float(maximum_fpr)
    if not 0.0 <= cap < 1.0:
        raise ValueError("maximum FPR must be in [0, 1)")
    candidates = np.concatenate(
        (
            [np.nextafter(np.max(s), np.inf)],
            np.unique(s)[::-1],
        )
    )
    feasible: list[tuple[float, float]] = []
    for threshold in candidates:
        prediction = s >= threshold
        fpr = float(np.mean(prediction[y == 0]))
        if fpr <= cap:
            tpr = float(np.mean(prediction[y == 1]))
            feasible.append((tpr, -float(threshold)))
    if not feasible:
        raise ValueError("no validation threshold satisfies the frozen FPR cap")
    best = max(feasible)
    return -best[1]


def _macro_f1(labels: np.ndarray, prediction: np.ndarray) -> float:
    y = np.asarray(labels, dtype=int)
    predicted = np.asarray(prediction, dtype=bool)
    values: list[float] = []
    for class_value in (0, 1):
        truth = y == class_value
        assigned = predicted if class_value == 1 else ~predicted
        true_positive = int(np.count_nonzero(truth & assigned))
        false_positive = int(np.count_nonzero(~truth & assigned))
        false_negative = int(np.count_nonzero(truth & ~assigned))
        denominator = 2 * true_positive + false_positive + false_negative
        values.append(0.0 if denominator == 0 else 2.0 * true_positive / denominator)
    return float(np.mean(values))


def _expected_calibration_error(
    labels: np.ndarray,
    probabilities: np.ndarray,
    *,
    bins: int = 10,
) -> float:
    y = np.asarray(labels, dtype=np.float64)
    probability = np.asarray(probabilities, dtype=np.float64)
    n_bins = _exact_int(bins, label="calibration bins")
    if y.shape != probability.shape or np.any(
        (probability < 0.0) | (probability > 1.0)
    ):
        raise ValueError("calibration arrays are invalid")
    indices = np.minimum((probability * n_bins).astype(int), n_bins - 1)
    result = 0.0
    for index in range(n_bins):
        selected = indices == index
        if np.any(selected):
            result += float(np.mean(selected)) * abs(
                float(np.mean(y[selected])) - float(np.mean(probability[selected]))
            )
    return result


def probe_metrics(
    labels: np.ndarray,
    scores: np.ndarray,
    components: Sequence[str],
    *,
    threshold: float,
    fpr_threshold: float | None = None,
    calibration_bins: int = 10,
) -> dict[str, float]:
    y = np.asarray(labels, dtype=int)
    s = np.asarray(scores, dtype=float)
    groups = np.asarray(tuple(components), dtype=object)
    if y.shape != s.shape or groups.shape != y.shape:
        raise ValueError("probe metric arrays do not align")
    pair_hits: list[float] = []
    for component in dict.fromkeys(groups.tolist()):
        mask = groups == component
        if (
            np.count_nonzero(mask & (y == 0)) != 1
            or np.count_nonzero(mask & (y == 1)) != 1
        ):
            raise ValueError("every component must contribute one clean/injected pair")
        delta = float(s[mask & (y == 1)][0] - s[mask & (y == 0)][0])
        pair_hits.append(float(delta > 0.0) + 0.5 * float(delta == 0.0))
    prediction = s >= float(threshold)
    balanced = 0.5 * (np.mean(prediction[y == 1]) + np.mean(~prediction[y == 0]))
    operating_threshold = float(threshold if fpr_threshold is None else fpr_threshold)
    operating_prediction = s >= operating_threshold
    probabilities = _sigmoid(s)
    return {
        "paired_component_ranking_accuracy": float(np.mean(pair_hits)),
        "roc_auc": _rank_auc(y, s),
        "average_precision": _average_precision(y, s),
        "balanced_accuracy": float(balanced),
        "macro_f1": _macro_f1(y, prediction),
        "brier_score": float(np.mean(np.square(probabilities - y))),
        "expected_calibration_error": _expected_calibration_error(
            y, probabilities, bins=calibration_bins
        ),
        "tpr_at_validation_frozen_5_percent_fpr": float(
            np.mean(operating_prediction[y == 1])
        ),
        "fpr_at_validation_frozen_5_percent_fpr": float(
            np.mean(operating_prediction[y == 0])
        ),
    }


def paired_component_bootstrap(
    labels: np.ndarray,
    scores: np.ndarray,
    components: Sequence[str],
    *,
    threshold: float,
    fpr_threshold: float | None = None,
    calibration_bins: int = 10,
    replicates: int = 1_000,
    seed: int = 560204,
) -> dict[str, Any]:
    """Bootstrap whole clean/injected component pairs, never individual rows."""

    y = np.asarray(labels, dtype=int)
    s = np.asarray(scores, dtype=float)
    groups = np.asarray(tuple(components), dtype=object)
    unique = tuple(dict.fromkeys(groups.tolist()))
    if len(unique) < 2:
        raise ValueError("paired bootstrap requires at least two components")
    by_component = {
        component: np.flatnonzero(groups == component) for component in unique
    }
    estimates = probe_metrics(
        y,
        s,
        groups,
        threshold=threshold,
        fpr_threshold=fpr_threshold,
        calibration_bins=calibration_bins,
    )
    rng = np.random.default_rng(int(seed))
    draws = {name: np.empty(int(replicates), dtype=float) for name in estimates}
    for draw_index in range(int(replicates)):
        sampled = rng.choice(unique, size=len(unique), replace=True)
        indices = np.concatenate(
            [by_component[str(component)] for component in sampled]
        )
        resampled_groups = tuple(
            f"{position}:{component}"
            for position, component in enumerate(sampled)
            for _ in by_component[str(component)]
        )
        metrics = probe_metrics(
            y[indices],
            s[indices],
            resampled_groups,
            threshold=threshold,
            fpr_threshold=fpr_threshold,
            calibration_bins=calibration_bins,
        )
        for name, value in metrics.items():
            draws[name][draw_index] = value
    return {
        name: {
            "estimate": value,
            "paired_component_bootstrap_95_interval": [
                float(np.quantile(draws[name], 0.025)),
                float(np.quantile(draws[name], 0.975)),
            ],
        }
        for name, value in estimates.items()
    }


def paired_auc_delta_bootstrap(
    labels: np.ndarray,
    left_scores: np.ndarray,
    right_scores: np.ndarray,
    components: Sequence[str],
    *,
    replicates: int = 1_000,
    seed: int = 560204,
) -> dict[str, Any]:
    """Bootstrap a paired AUROC difference on identical component pairs."""

    y = np.asarray(labels, dtype=int)
    left = np.asarray(left_scores, dtype=float)
    right = np.asarray(right_scores, dtype=float)
    groups = np.asarray(tuple(components), dtype=object)
    if y.shape != left.shape or right.shape != y.shape or groups.shape != y.shape:
        raise ValueError("paired AUROC delta arrays do not align")
    unique = tuple(dict.fromkeys(groups.tolist()))
    if len(unique) < 2:
        raise ValueError("paired AUROC delta requires at least two components")
    by_component = {
        component: np.flatnonzero(groups == component) for component in unique
    }
    estimate = _rank_auc(y, left) - _rank_auc(y, right)
    rng = np.random.default_rng(int(seed))
    draws = np.empty(int(replicates), dtype=np.float64)
    for draw_index in range(int(replicates)):
        sampled = rng.choice(unique, size=len(unique), replace=True)
        indices = np.concatenate(
            [by_component[str(component)] for component in sampled]
        )
        draws[draw_index] = _rank_auc(y[indices], left[indices]) - _rank_auc(
            y[indices], right[indices]
        )
    return {
        "estimate": float(estimate),
        "paired_component_bootstrap_95_interval": [
            float(np.quantile(draws, 0.025)),
            float(np.quantile(draws, 0.975)),
        ],
        "bootstrap_replicates": int(replicates),
    }


def summarize_readiness(
    feature_results: Mapping[str, Any],
    test_scores: Mapping[str, np.ndarray],
    labels: np.ndarray,
    components: Sequence[str],
    *,
    readiness_config: Mapping[str, Any],
    bootstrap_replicates: int,
    bootstrap_seed: int,
) -> dict[str, Any]:
    """Apply the frozen development-canary readiness rule."""

    primary = str(readiness_config["primary_feature"])
    primary_result = _mapping(feature_results[primary], label="primary result")
    overall = _mapping(
        primary_result["locked_development_test"], label="primary test metrics"
    )
    auc = _mapping(overall["roc_auc"], label="primary AUROC")
    paired_ranking = _mapping(
        overall["paired_component_ranking_accuracy"], label="paired ranking"
    )
    tpr = _mapping(
        overall["tpr_at_validation_frozen_5_percent_fpr"], label="primary TPR"
    )
    fpr = _mapping(
        overall["fpr_at_validation_frozen_5_percent_fpr"], label="primary FPR"
    )
    cohort_results = _mapping(
        primary_result["locked_development_test_by_cohort"],
        label="primary cohort metrics",
    )

    deltas: dict[str, Any] = {}
    delta_criteria: dict[str, bool] = {}
    for index, control in enumerate(readiness_config["paired_delta_controls"]):
        control_name = str(control)
        delta = paired_auc_delta_bootstrap(
            labels,
            test_scores[primary],
            test_scores[control_name],
            components,
            replicates=bootstrap_replicates,
            seed=bootstrap_seed + 500 + index,
        )
        deltas[control_name] = delta
        lower = float(delta["paired_component_bootstrap_95_interval"][0])
        delta_criteria[control_name] = (
            float(delta["estimate"])
            >= float(readiness_config["minimum_paired_roc_auc_delta"])
            and lower > 0.0
        )

    raw_control = str(readiness_config["raw_noninferiority_control"])
    raw_delta = paired_auc_delta_bootstrap(
        labels,
        test_scores[primary],
        test_scores[raw_control],
        components,
        replicates=bootstrap_replicates,
        seed=bootstrap_seed + 600,
    )
    raw_lower = float(raw_delta["paired_component_bootstrap_95_interval"][0])
    raw_noninferior = raw_lower >= -float(readiness_config["raw_noninferiority_margin"])
    raw_superior = raw_lower > 0.0

    criteria: dict[str, bool] = {
        "primary_roc_auc_lower_95_at_least_0_75": float(
            auc["paired_component_bootstrap_95_interval"][0]
        )
        >= float(readiness_config["minimum_roc_auc_lower_95"]),
        "paired_roc_auc_delta_requirements_pass": all(delta_criteria.values()),
        "raw_baseline_noninferiority_margin_0_01_pass": raw_noninferior,
        "each_cohort_point_roc_auc_at_least_0_70": all(
            float(cohort_results[cohort]["roc_auc"]["estimate"])
            >= float(readiness_config["minimum_each_cohort_roc_auc"])
            for cohort in TEMPORAL_COHORTS
        ),
        "tpr_at_validation_frozen_5_percent_fpr_point_at_least_0_30": float(
            tpr["estimate"]
        )
        >= float(readiness_config["minimum_tpr_at_5_percent_fpr"]),
        "tpr_at_validation_frozen_5_percent_fpr_lower_95_at_least_0_20": float(
            tpr["paired_component_bootstrap_95_interval"][0]
        )
        >= float(readiness_config["minimum_tpr_at_5_percent_fpr_lower_95"]),
        "locked_test_fpr_at_validation_frozen_threshold_at_most_0_10": float(
            fpr["estimate"]
        )
        <= float(
            readiness_config[
                "maximum_locked_test_fpr_at_validation_frozen_threshold"
            ]
        ),
        "paired_ranking_lower_95_above_chance": float(
            paired_ranking["paired_component_bootstrap_95_interval"][0]
        )
        > 0.5,
    }
    ready = all(criteria.values())
    return {
        "ready_for_next_real_training": ready,
        "useful_representation_claim_supported": ready and raw_superior,
        "raw_noninferiority_is_not_superiority": True,
        "criteria": criteria,
        "paired_delta_criteria": delta_criteria,
        "paired_roc_auc_delta_vs_controls": deltas,
        "paired_roc_auc_delta_vs_raw_baseline": raw_delta,
        "interpretation": (
            "ready_for_controlled_next_training"
            if ready
            else "development_transfer_gate_not_yet_met"
        ),
    }


def summarize_probe_mechanics(
    feature_results: Mapping[str, Any],
    *,
    mechanics_config: Mapping[str, Any],
) -> dict[str, Any]:
    """Verify that flux, but not quality metadata, drives the v2 probe."""

    raw = _mapping(
        feature_results["raw_adp_validity_error_cadence_128"],
        label="raw-control result",
    )
    quality = _mapping(
        feature_results["quality_only_cadence_128"],
        label="quality-only result",
    )
    raw_metrics = _mapping(
        raw["locked_development_test"],
        label="raw-control locked metrics",
    )
    raw_auc = float(raw_metrics["roc_auc"]["estimate"])
    raw_auc_lower = float(
        raw_metrics["roc_auc"]["paired_component_bootstrap_95_interval"][0]
    )
    raw_ranking_lower = float(
        raw_metrics["paired_component_ranking_accuracy"][
            "paired_component_bootstrap_95_interval"
        ][0]
    )
    raw_flux_coefficients = tuple(float(value) for value in raw["raw_flux_coefficients"])
    if len(raw_flux_coefficients) != 2:
        raise ValueError("raw-control result lacks its two ADP flux coefficients")
    quality_auc = float(
        quality["locked_development_test"]["roc_auc"]["estimate"]
    )
    raw_cells = _mapping(
        raw["locked_development_test_by_event_cell"],
        label="raw-control event cells",
    )
    deep_cells = {
        name: float(value["roc_auc"]["estimate"])
        for name, value in raw_cells.items()
        if str(name).endswith("_depth_0.30")
    }
    if len(deep_cells) != len(EVENT_DURATIONS):
        raise ValueError("raw-control mechanics lacks the three depth-0.30 cells")
    criteria = {
        "raw_control_overall_roc_auc_lower_95_at_least_floor": raw_auc_lower
        >= float(mechanics_config["raw_control_minimum_overall_roc_auc_lower_95"]),
        "raw_control_each_depth_0_30_roc_auc_at_least_floor": all(
            value
            >= float(
                mechanics_config["raw_control_minimum_each_depth_0_30_roc_auc"]
            )
            for value in deep_cells.values()
        ),
        "raw_control_adp_flux_coefficients_are_dip_positive": (
            all(value < 0.0 for value in raw_flux_coefficients)
            if mechanics_config["raw_control_require_negative_adp_flux_coefficients"]
            else True
        ),
        "raw_control_paired_ranking_lower_95_above_chance": (
            raw_ranking_lower > 0.5
            if mechanics_config[
                "raw_control_paired_ranking_lower_95_strictly_above_chance"
            ]
            else True
        ),
        "quality_only_roc_auc_near_chance": abs(quality_auc - 0.5)
        <= float(
            mechanics_config["quality_only_maximum_absolute_roc_auc_from_chance"]
        ),
    }
    return {
        "passed": all(criteria.values()),
        "criteria": criteria,
        "raw_control_roc_auc": raw_auc,
        "raw_control_roc_auc_lower_95": raw_auc_lower,
        "raw_control_paired_ranking_lower_95": raw_ranking_lower,
        "raw_control_flux_coefficients": list(raw_flux_coefficients),
        "raw_control_depth_0_30_cell_roc_auc": deep_cells,
        "quality_only_roc_auc": quality_auc,
        "interpretation": (
            "token_supervised_probe_mechanics_valid"
            if all(criteria.values())
            else "token_supervised_probe_mechanics_failed"
        ),
    }


def assert_cadence_preserving_model(
    model: Any, *, context_length: int, output: Mapping[str, Any]
) -> None:
    """Reject any effective Conformer patching or cadence-token count loss."""

    config = model.config
    if config.architecture == "conformer" and int(config.patch_stride) != 1:
        raise ValueError(
            "event-transfer Conformer must use patch_stride=1; cadence averaging is forbidden"
        )
    hidden = output.get("h_cadence")
    if hidden is None or tuple(hidden.shape[1:]) != (
        int(context_length),
        int(config.d_model),
    ):
        raise ValueError(
            "encoder did not emit exactly one hidden token per input cadence"
        )


def encode_cadence_tokens(
    model: Any,
    samples: Sequence[Mapping[str, np.ndarray]],
    *,
    context_length: int,
    device: Any,
    batch_size: int = 8,
    project_each_cadence: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Encode cadence tokens without any temporal pooling or token merging."""

    try:
        import torch
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("PyTorch is required for event-transfer inference") from exc
    tokens: list[np.ndarray] = []
    validity: list[np.ndarray] = []
    model.eval()
    with torch.no_grad():
        for start in range(0, len(samples), batch_size):
            batch = move_batch_to_device(
                collate_fm0_samples(samples[start : start + batch_size]), device
            )
            output = (
                model(batch)
                if context_length == BASE_INTERVAL_CADENCES
                else model.forward_short_context(batch)
            )
            assert_cadence_preserving_model(
                model, context_length=context_length, output=output
            )
            hidden = output["h_cadence"]
            if project_each_cadence:
                hidden = model.embedding_projection(hidden)
            values = hidden.detach().float().cpu().numpy()
            valid = output["token_valid"].detach().cpu().numpy().astype(bool)
            if values.shape[:2] != valid.shape:
                raise ValueError("cadence tokens and validity mask differ")
            tokens.append(values)
            validity.append(valid)
    return np.concatenate(tokens), np.concatenate(validity)


def _feature_tokens(
    spec: str,
    *,
    samples: Sequence[Mapping[str, np.ndarray]],
    model0: Any,
    model2000: Any,
    device: Any,
    batch_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    length = 2_048 if spec.endswith("_2048") else 128
    if spec.startswith("raw_"):
        return raw_cadence_features(samples, quality_only=False)
    if spec.startswith("quality_only"):
        return raw_cadence_features(samples, quality_only=True)
    model = model0 if spec.startswith("step0_") else model2000
    return encode_cadence_tokens(
        model,
        samples,
        context_length=length,
        device=device,
        batch_size=batch_size,
        project_each_cadence="_z_cadence_" in spec,
    )


def evaluate_event_transfer_canary(
    *,
    config_path: str | Path,
    run_dir: str | Path,
    step0_checkpoint_path: str | Path,
    step2000_checkpoint_path: str | Path,
    temporal_panel_dir: str | Path,
    temporal_panel_receipt_sha256: str,
    batch_size: int = 8,
) -> dict[str, Any]:
    """Run the complete development-only cadence-token transfer canary."""

    try:
        import torch
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("PyTorch is required for event-transfer inference") from exc
    config, config_source, config_sha = load_event_transfer_config(config_path)
    token_supervised = (
        config["schema_version"] == EVENT_TRANSFER_CONFIG_SCHEMA_VERSION_V2
    )
    inputs = _mapping(config["inputs"], label="inputs")
    if temporal_panel_receipt_sha256 != inputs["temporal_panel_receipt_sha256"]:
        raise ValueError("event-transfer temporal-panel receipt hash differs")
    rows, panel_receipt = load_temporal_panel(
        temporal_panel_dir, receipt_sha256=temporal_panel_receipt_sha256
    )
    device = torch.device("cpu")
    model0, contract0, validation0 = _load_inference_only_trusted_model(
        Path(run_dir), device=device, checkpoint_path=Path(step0_checkpoint_path)
    )
    model2000, contract2000, validation2000 = _load_inference_only_trusted_model(
        Path(run_dir), device=device, checkpoint_path=Path(step2000_checkpoint_path)
    )
    if (
        validation0["global_step"] != 0
        or validation2000["global_step"] != 2_000
        or contract0 != contract2000
    ):
        raise ValueError("event-transfer checkpoint pair differs")
    hashes = _mapping(inputs["checkpoint_sha256"], label="checkpoint hashes")
    if (
        validation0["selected_checkpoint_sha256"] != hashes["step0"]
        or validation2000["selected_checkpoint_sha256"] != hashes["step2000"]
    ):
        raise ValueError("event-transfer checkpoint hash differs")

    print("FM_EVENT_TRANSFER phase=freeze_schedule", flush=True)
    schedule, schedule_audit = freeze_component_schedule(
        rows, variant=str(inputs["variant"])
    )
    print(
        f"FM_EVENT_TRANSFER phase=schedule_ready components={len(schedule)}",
        flush=True,
    )
    block_data: dict[
        str,
        dict[
            int,
            tuple[
                list[dict[str, np.ndarray]],
                np.ndarray,
                tuple[str, ...],
                tuple[str, ...],
                np.ndarray,
            ],
        ],
    ] = {}
    block_cells: dict[str, tuple[str, ...]] = {}
    for block in SPLIT_SECTORS:
        subset = [item for item in schedule if item["block"] == block]
        block_data[block] = {
            length: build_paired_samples_with_token_targets(
                subset,
                context_length=length,
            )
            for length in CONTEXT_LENGTHS
        }
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
    short_step2000_cache: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    projection_weight = (
        model2000.embedding_projection.weight.detach().float().cpu().numpy()
    )
    projection_bias = (
        model2000.embedding_projection.bias.detach().float().cpu().numpy()
    )
    for spec in FEATURE_SPECS:
        print(f"FM_EVENT_TRANSFER feature={spec} phase=encode", flush=True)
        length = 2_048 if spec.endswith("_2048") else 128
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
            samples, labels, components, cohorts, token_targets = block_data[block][
                length
            ]
            if spec == "step2000_z_cadence_128_diagnostic":
                base_tokens, valid = short_step2000_cache[block]
                tokens = base_tokens @ projection_weight.T + projection_bias
            else:
                tokens, valid = _feature_tokens(
                    spec,
                    samples=samples,
                    model0=model0,
                    model2000=model2000,
                    device=device,
                    batch_size=batch_size,
                )
                if spec == "step2000_h_cadence_128":
                    short_step2000_cache[block] = tokens, valid
            encoded[block] = (
                tokens,
                valid,
                labels,
                components,
                cohorts,
                token_targets,
            )
            print(
                f"FM_EVENT_TRANSFER feature={spec} block={block} phase=encoded",
                flush=True,
            )
        train = encoded["train"]
        print(f"FM_EVENT_TRANSFER feature={spec} phase=fit", flush=True)
        fit_function = (
            fit_shared_linear_token_probe
            if token_supervised
            else fit_shared_linear_max_probe
        )
        fit_target = train[5] if token_supervised else train[2]
        fit = fit_function(
            train[0],
            train[1],
            fit_target,
            learning_rate=float(probe_config["learning_rate"]),
            l2_weight=float(probe_config["l2_weight"]),
            epochs=int(probe_config["epochs"]),
            seed=int(probe_config["initialization_seed"]),
            progress_label=spec,
        )
        validation = encoded["validation"]
        print(f"FM_EVENT_TRANSFER feature={spec} phase=threshold", flush=True)
        validation_scores = score_shared_linear_max_probe(
            fit, validation[0], validation[1]
        )
        threshold = select_balanced_accuracy_threshold(validation[2], validation_scores)
        fpr_threshold = select_fpr_threshold(
            validation[2],
            validation_scores,
            maximum_fpr=float(metric_config["fpr_operating_point"]),
        )
        test = encoded["locked_development_test"]
        test_scores = score_shared_linear_max_probe(fit, test[0], test[1])
        test_scores_by_spec[spec] = test_scores
        current_reference = (test[2], test[3], test[4])
        if test_reference is None:
            test_reference = current_reference
        elif (
            not np.array_equal(test_reference[0], current_reference[0])
            or test_reference[1] != current_reference[1]
            or test_reference[2] != current_reference[2]
        ):
            raise ValueError("event-transfer feature controls changed test identities")
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
        by_cohort: dict[str, Any] = {}
        cohort_values = np.asarray(test[4], dtype=object)
        for cohort in TEMPORAL_COHORTS:
            mask = cohort_values == cohort
            by_cohort[cohort] = paired_component_bootstrap(
                test[2][mask],
                test_scores[mask],
                np.asarray(test[3], dtype=object)[mask],
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
            mask = cell_values == cell
            by_cell[cell] = paired_component_bootstrap(
                test[2][mask],
                test_scores[mask],
                np.asarray(test[3], dtype=object)[mask],
                threshold=threshold,
                fpr_threshold=fpr_threshold,
                calibration_bins=int(metric_config["calibration_bins"]),
                replicates=int(metric_config["bootstrap_replicates"]),
                seed=int(metric_config["bootstrap_seed"])
                + 200
                + FEATURE_SPECS.index(spec),
            )
            for cohort in TEMPORAL_COHORTS:
                joint = mask & (cohort_values == cohort)
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
            "context_cadences": length,
            "input_cadence_tokens": length,
            "temporal_downsampling": False,
            "pooling_before_token_score": False,
            "probe_training_objective": probe_config.get(
                "training_objective",
                "sample_binary_cross_entropy_after_hard_masked_max",
            ),
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
        }
        if spec == "raw_adp_validity_error_cadence_128":
            results[spec]["raw_flux_coefficients"] = [
                float(fit.weight[0]),
                float(fit.weight[1]),
            ]
        print(f"FM_EVENT_TRANSFER feature={spec} phase=complete", flush=True)

    if test_reference is None:
        raise ValueError("event-transfer produced no locked development scores")
    readiness = summarize_readiness(
        results,
        test_scores_by_spec,
        test_reference[0],
        test_reference[1],
        readiness_config=_mapping(config["readiness"], label="readiness"),
        bootstrap_replicates=int(metric_config["bootstrap_replicates"]),
        bootstrap_seed=int(metric_config["bootstrap_seed"]),
    )
    probe_mechanics = (
        summarize_probe_mechanics(
            results,
            mechanics_config=_mapping(
                config["probe_mechanics"],
                label="probe mechanics",
            ),
        )
        if token_supervised
        else {
            "passed": True,
            "interpretation": "legacy_v1_probe_mechanics_not_gated",
        }
    )
    metric_readiness_passed = bool(readiness["ready_for_next_real_training"])
    overall_readiness_passed = bool(
        probe_mechanics["passed"] and metric_readiness_passed
    )
    reported_readiness = dict(readiness)
    if token_supervised:
        reported_readiness["fm_metric_criteria_satisfied"] = metric_readiness_passed
        reported_readiness["probe_mechanics_gate_required"] = True
        reported_readiness["ready_for_next_real_training"] = (
            overall_readiness_passed
        )

    schedule_public = []
    for item in schedule:
        public = {key: value for key, value in item.items() if key != "visit"}
        public["window_binding"] = dict(item["visit"]["binding"])
        schedule_public.append(public)
    return {
        "schema_version": (
            EVENT_TRANSFER_RESULT_SCHEMA_VERSION_V2
            if token_supervised
            else EVENT_TRANSFER_RESULT_SCHEMA_VERSION
        ),
        "campaign_id": (
            EVENT_TRANSFER_CAMPAIGN_ID_V2
            if token_supervised
            else EVENT_TRANSFER_CAMPAIGN_ID
        ),
        "evaluation_completed": True,
        "passed": overall_readiness_passed,
        "probe_mechanics_passed": bool(probe_mechanics["passed"]),
        "fm_metric_criteria_passed": metric_readiness_passed,
        "scientific_readiness_passed": overall_readiness_passed,
        "config_path": str(config_source),
        "config_sha256": config_sha,
        "temporal_panel_receipt_sha256": temporal_panel_receipt_sha256,
        "temporal_panel_rows": int(panel_receipt["n_panel_rows"]),
        "schedule_sha256": hashlib.sha256(
            json.dumps(schedule_public, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest(),
        "schedule": schedule_public,
        "schedule_audit": schedule_audit,
        "checkpoint_validation": {"step0": validation0, "step2000": validation2000},
        "feature_results": results,
        "probe_mechanics": probe_mechanics,
        "readiness": reported_readiness,
        "cadence_preservation": {
            "nominal_cadence_seconds": 200,
            "one_encoder_token_per_input_cadence": True,
            "patching_or_cadence_averaging": False,
            "only_temporal_reduction": "masked_max_after_shared_linear_cadence_score",
            "probe_training_target_used_as_input_feature": False,
            "probe_training_loss_reduction": (
                "class_balanced_scalar_loss_over_independently_scored_cadence_examples"
                if token_supervised
                else "legacy_sample_loss_after_hard_masked_max"
            ),
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
            "locked_development_test_used_for_fit_or_threshold": False,
        },
    }


__all__ = [
    "EVENT_TRANSFER_CONFIG_SCHEMA_VERSION",
    "EVENT_TRANSFER_CONFIG_SCHEMA_VERSION_V2",
    "EVENT_TRANSFER_RESULT_SCHEMA_VERSION",
    "EVENT_TRANSFER_RESULT_SCHEMA_VERSION_V2",
    "FEATURE_SPECS",
    "ProbeFit",
    "assert_cadence_preserving_model",
    "build_paired_samples",
    "build_paired_samples_with_token_targets",
    "encode_cadence_tokens",
    "evaluate_event_transfer_canary",
    "fit_shared_linear_max_probe",
    "fit_shared_linear_token_probe",
    "freeze_component_schedule",
    "inject_jittered_single_event",
    "load_event_transfer_config",
    "paired_auc_delta_bootstrap",
    "paired_component_bootstrap",
    "paired_token_probe_examples",
    "probe_metrics",
    "raw_cadence_features",
    "score_shared_linear_max_probe",
    "sector_block",
    "select_balanced_accuracy_threshold",
    "select_fpr_threshold",
    "summarize_probe_mechanics",
    "summarize_readiness",
]
