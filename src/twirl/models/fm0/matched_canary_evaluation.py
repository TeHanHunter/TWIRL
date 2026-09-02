"""Matched, candidate-centred evaluation for the TWIRL-FM0.3 canary.

The evaluator consumes only the immutable payload-screened schedule.  It
reopens each of the 1,440 named shards exactly once, verifies both the source
and cropped-payload hashes, and then constructs one adjacent clean/injected
pair.  The default controls-only mode needs NumPy but no model checkpoint or
GPU.  Optional checkpoint arms expose only ``h_cadence[:, 64, :]`` to the same
fixed linear-probe implementation.

This module deliberately contains no search, BLS, training, threshold tuning,
or sealed-data path.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from io import StringIO
from pathlib import Path, PurePosixPath
from typing import Any

import numpy as np

from .dataset import collate_fm0_samples, move_batch_to_device, prepare_model_window
from .input_release import load_input_release_bytes
from .matched_canary_payload_plan import (
    BOOTSTRAP_REPLICATES,
    BOOTSTRAP_SEED,
    CONTEXT_CADENCES,
    EVENT_CENTER_INDEX_ZERO_BASED,
    PAYLOAD_SCHEDULE_FIELDS,
    PROBE_EPOCHS,
    PROBE_L2,
    PROBE_LEARNING_RATE,
    PROBE_SEED,
    QUALITY_FEATURE_CHANNELS,
    RAW_FEATURE_CHANNELS,
    crop_arrays,
    crop_payload_sha256,
    inject_centered_event,
    validate_matched_canary_payload_plan,
)
from .matched_canary_plan import (
    EVENT_DEPTHS_TEXT,
    EVENT_DURATIONS,
    SPLIT_ORDER,
    TEMPORAL_COHORTS,
    _frozen_analysis_contract,
)
from .registry import FM0ContractError, sha256_file

EVALUATION_SCHEMA_VERSION = "twirl_fm0_3_matched_canary_evaluation_v1"
EVALUATION_RECEIPT_SCHEMA_VERSION = "twirl_fm0_3_matched_canary_evaluation_receipt_v1"
EVALUATION_READY_STATE = "FM0_3_MATCHED_CANARY_EVALUATION_READY"
TEMPORAL_PANEL_RECEIPT_SHA256 = (
    "78c370e10c556472c5997c20cfe95207a0b334bafe7f024bf7ba4fc7ec4de624"
)

CONTROL_ARMS = (
    "raw_adp_validity_error_exact_center",
    "quality_only_exact_center",
)
MODEL_ARM_CONTRACT: Mapping[str, tuple[str, int, str]] = {
    "TWIRL-FM0.3.1_step0_h_cadence_token64": ("TWIRL-FM0.3.1", 0, "tcn"),
    "TWIRL-FM0.3.1_step2000_h_cadence_token64": (
        "TWIRL-FM0.3.1",
        2_000,
        "tcn",
    ),
    "TWIRL-FM0.3.2_step0_h_cadence_token64": (
        "TWIRL-FM0.3.2",
        0,
        "conformer",
    ),
    "TWIRL-FM0.3.2_step2000_h_cadence_token64": (
        "TWIRL-FM0.3.2",
        2_000,
        "conformer",
    ),
}

SCORE_FIELDS = (
    "evaluation_schema_version",
    "feature_arm",
    "split",
    "cohort",
    "sector",
    "observation_key",
    "leakage_component_id",
    "event_pair_id",
    "duration_cadences",
    "fractional_depth",
    "sample_role",
    "label",
    "score",
)

_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_GIT_SHA = re.compile(r"^[0-9a-f]{40}$")


@dataclass(frozen=True, slots=True)
class ModelArm:
    """One optional immutable FM0.3 checkpoint feature arm."""

    name: str
    run_dir: Path
    checkpoint_path: Path


@dataclass(frozen=True, slots=True)
class FrozenPair:
    """One schedule identity and its exact clean/injected direct-128 tensors."""

    metadata: Mapping[str, str]
    clean: Mapping[str, np.ndarray]
    injected: Mapping[str, np.ndarray]


@dataclass(frozen=True, slots=True)
class LinearProbe:
    """Fixed full-batch Adam probe and train-only standardization."""

    weight: np.ndarray
    bias: float
    center: np.ndarray
    scale: np.ndarray
    objective_history: tuple[float, ...]


@dataclass(frozen=True, slots=True)
class MatchedCanaryEvaluationResult:
    root: Path
    results_path: Path
    scores_path: Path
    receipt_path: Path
    ready_path: Path
    results_sha256: str
    scores_sha256: str
    receipt_sha256: str
    receipt: Mapping[str, Any]


def _digest(value: Any, *, label: str) -> str:
    digest = str(value).strip().lower()
    if _SHA256.fullmatch(digest) is None:
        raise FM0ContractError(f"{label} must be a lowercase SHA-256")
    return digest


def _git_sha(value: Any) -> str:
    sha = str(value).strip().lower()
    if _GIT_SHA.fullmatch(sha) is None:
        raise FM0ContractError("producer_git_sha must be a full lowercase Git SHA")
    return sha


def _canonical_json_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode(
        "utf-8"
    )


def _read_schedule(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != PAYLOAD_SCHEDULE_FIELDS:
            raise FM0ContractError("payload evaluation schedule columns drifted")
        rows = [dict(row) for row in reader]
    if len(rows) != 1_440:
        raise FM0ContractError("payload evaluation schedule must contain 1,440 rows")
    return rows


def evaluation_algorithm_contract() -> dict[str, Any]:
    """Return the exact mechanics implemented by this evaluator."""

    return {
        "source_loading": {
            "schedule_rows": 1_440,
            "each_bound_source_shard_reopened_exactly_once": True,
            "source_shard_sha256_verified_before_npz_load": True,
            "crop_payload_sha256_recomputed_and_verified": True,
            "replacement_or_resampling_after_payload_screen": False,
            "required_temporal_panel_receipt_sha256": (TEMPORAL_PANEL_RECEIPT_SHA256),
        },
        "pair_construction": {
            "one_adjacent_clean_injected_pair_per_schedule_row": True,
            "injection_helper": (
                "twirl.models.fm0.matched_canary_payload_plan.inject_centered_event"
            ),
            "window_local_time_origin": "first_time_valid_cadence_in_crop",
            "cadence_count": CONTEXT_CADENCES,
            "nominal_cadence_seconds": 200,
            "cadence_averaging_merging_resampling_or_downsampling": False,
        },
        "feature_arms": {
            CONTROL_ARMS[0]: list(RAW_FEATURE_CHANNELS),
            CONTROL_ARMS[1]: list(QUALITY_FEATURE_CHANNELS),
            "optional_model_arms": list(MODEL_ARM_CONTRACT),
            "model_feature": "h_cadence_token_64_only",
            "model_emits_all_128_tokens_before_center_selection": True,
            "artificial_temporal_and_reconstruction_masks_all_false": True,
            "duration_depth_support_period_bls_or_search_score_input": False,
        },
        "probe": {
            "family": "full_batch_adam_linear_logit",
            "fit_split": "probe_train",
            "standardization": "dimensionwise_mean_and_std_fit_on_probe_train_only",
            "zero_or_nonfinite_scale_replacement": 1.0,
            "weight_initialization": "normal_mean_0_std_0.001",
            "bias_initialization": 0.0,
            "epochs": PROBE_EPOCHS,
            "learning_rate": PROBE_LEARNING_RATE,
            "l2_weight": PROBE_L2,
            "adam_beta1": 0.9,
            "adam_beta2": 0.999,
            "adam_epsilon": 1.0e-8,
            "seed": PROBE_SEED,
            "validation_or_test_tuning": False,
        },
        "metrics": {
            "primary": "sample_roc_auc",
            "diagnostic": "paired_clean_injected_ranking_accuracy",
            "bootstrap": {
                "method": "whole_component_clean_injected_pair_resampling",
                "replicates": BOOTSTRAP_REPLICATES,
                "seed": BOOTSTRAP_SEED,
                "confidence_interval_quantiles": [0.025, 0.975],
                "same_resample_indices_across_arms_and_paired_deltas": True,
            },
            "breakdowns": [
                "split",
                "cohort_within_split",
                "duration_depth_cell_within_split",
                "depth_aggregated_over_durations_within_split",
            ],
        },
        "frozen_primary_test_gates": _frozen_analysis_contract()["primary_test_gates"],
        "forbidden": {
            "bls_or_periodic_search": True,
            "sealed_test_access": True,
            "fm_encoder_training_or_optimizer": True,
            "hyperparameter_or_threshold_tuning": True,
        },
    }


def _temporal_panel_authority(
    payload_receipt: Mapping[str, Any],
) -> dict[str, Any]:
    """Extract and enforce the sole frozen temporal-panel authority."""

    identity = payload_receipt.get("identity_plan")
    authorities = (
        identity.get("input_authorities") if isinstance(identity, Mapping) else None
    )
    panel = (
        authorities.get("temporal_panel") if isinstance(authorities, Mapping) else None
    )
    if not isinstance(panel, Mapping):
        raise FM0ContractError("payload plan lacks its temporal-panel authority")
    expected_fields = {
        "root",
        "receipt_sha256",
        "panel_sha256",
        "sector_bindings_sha256",
    }
    if set(panel) != expected_fields:
        raise FM0ContractError("payload temporal-panel authority fields drifted")
    receipt_sha = _digest(
        panel.get("receipt_sha256"), label="temporal-panel receipt hash"
    )
    if receipt_sha != TEMPORAL_PANEL_RECEIPT_SHA256:
        raise FM0ContractError("payload binds an unauthorized temporal-panel receipt")
    return {
        "root": str(panel["root"]),
        "receipt_sha256": receipt_sha,
        "panel_sha256": _digest(panel["panel_sha256"], label="temporal-panel hash"),
        "sector_bindings_sha256": _digest(
            panel["sector_bindings_sha256"],
            label="temporal-panel sector-bindings hash",
        ),
    }


def _safe_source_path(row: Mapping[str, str], *, require_read_only: bool) -> Path:
    root_raw = Path(row["source_release_root"]).expanduser()
    if root_raw.is_symlink():
        raise FM0ContractError("evaluation source root must be materialized")
    root = root_raw.resolve(strict=True)
    if not root.is_dir() or (require_read_only and root.stat().st_mode & 0o222):
        raise FM0ContractError("evaluation source root is unavailable or writable")
    relative = PurePosixPath(row["source_relative_path"])
    if relative.is_absolute() or ".." in relative.parts:
        raise FM0ContractError("evaluation source shard path is unsafe")
    raw_path = root.joinpath(*relative.parts)
    if raw_path.is_symlink():
        raise FM0ContractError("evaluation source shard must be materialized")
    path = raw_path.resolve(strict=True)
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise FM0ContractError("evaluation source shard escaped its release") from exc
    if (
        not path.is_file()
        or path.stat().st_size <= 0
        or (require_read_only and path.stat().st_mode & 0o222)
    ):
        raise FM0ContractError("evaluation source shard is unavailable or writable")
    return path


def _normalize_crop_time(arrays: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Apply the ordinary FM window-relative time origin without resampling."""

    normalized = {name: np.asarray(value).copy() for name, value in arrays.items()}
    local_time = normalized["local_time_cadences"]
    time_valid = np.asarray(normalized["time_valid"], dtype=bool)
    valid = time_valid & np.isfinite(local_time)
    if not np.any(valid):
        raise FM0ContractError("frozen crop has no finite time-valid cadence")
    reference = float(local_time[np.flatnonzero(valid)[0]])
    local_time[valid] -= reference
    return normalized


def _load_one_frozen_pair(
    row: Mapping[str, str], *, require_read_only: bool
) -> FrozenPair:
    path = _safe_source_path(row, require_read_only=require_read_only)
    payload = path.read_bytes()
    if hashlib.sha256(payload).hexdigest() != row["source_shard_sha256"]:
        raise FM0ContractError("evaluation source shard SHA-256 drifted")
    release = load_input_release_bytes(payload)
    try:
        source_n_cadences = int(row["source_n_cadences"])
        source_n_segments = int(row["source_n_segments"])
        start = int(row["crop_start_index_zero_based"])
        stop = int(row["crop_stop_index_exclusive"])
        segment = int(row["crop_segment_id"])
        duration = int(row["duration_cadences"])
        depth = float(row["fractional_depth"])
    except (TypeError, ValueError) as exc:
        raise FM0ContractError("evaluation schedule numeric field is invalid") from exc
    if (
        release.n_cadences != source_n_cadences
        or release.n_segments != source_n_segments
        or stop - start != CONTEXT_CADENCES
    ):
        raise FM0ContractError("evaluation source or crop dimensions drifted")
    source_crop = crop_arrays(release, start=start)
    if (
        np.any(np.asarray(source_crop["segment_id"]) != segment)
        or crop_payload_sha256(source_crop) != row["crop_payload_sha256"]
    ):
        raise FM0ContractError("evaluation crop payload binding drifted")
    clean = _normalize_crop_time(source_crop)
    injected = inject_centered_event(
        clean,
        duration_cadences=duration,
        fractional_depth=depth,
    )
    return FrozenPair(metadata=dict(row), clean=clean, injected=injected)


def load_frozen_pairs(
    rows: Sequence[Mapping[str, str]], *, require_read_only: bool = True
) -> tuple[FrozenPair, ...]:
    """Reopen exactly the schedule's unique source shards once apiece."""

    if len(rows) != 1_440:
        raise FM0ContractError("evaluation requires exactly 1,440 frozen crops")
    seen: set[Path] = set()
    pairs: list[FrozenPair] = []
    for index, row in enumerate(rows, start=1):
        path = _safe_source_path(row, require_read_only=require_read_only)
        if path in seen:
            raise FM0ContractError("evaluation schedule repeats a source shard")
        seen.add(path)
        # _load_one_frozen_pair resolves but opens this exact path only once.
        pairs.append(_load_one_frozen_pair(row, require_read_only=require_read_only))
        if index == 1 or index % 50 == 0 or index == len(rows):
            print(
                f"FM0_3_EVALUATION phase=load_pairs reopened={index}/{len(rows)}",
                flush=True,
            )
    if len(seen) != 1_440:
        raise FM0ContractError("evaluation did not reopen exactly 1,440 unique shards")
    return tuple(pairs)


def _paired_metadata(
    pairs: Sequence[FrozenPair],
) -> tuple[np.ndarray, tuple[str, ...], tuple[Mapping[str, str], ...]]:
    labels: list[int] = []
    components: list[str] = []
    metadata: list[Mapping[str, str]] = []
    for pair in pairs:
        for label in (0, 1):
            labels.append(label)
            components.append(pair.metadata["leakage_component_id"])
            metadata.append(pair.metadata)
    return np.asarray(labels, dtype=np.int8), tuple(components), tuple(metadata)


def exact_center_control_features(
    pairs: Sequence[FrozenPair], *, quality_only: bool
) -> np.ndarray:
    """Return the frozen raw or quality-only vector at cadence index 64."""

    rows: list[np.ndarray] = []
    for pair in pairs:
        for sample in (pair.clean, pair.injected):
            center = EVENT_CENTER_INDEX_ZERO_BASED
            flux = np.asarray(sample["flux"], dtype=np.float64)
            flux_valid = np.asarray(sample["flux_valid"], dtype=bool)
            error = np.asarray(sample["flux_error"], dtype=np.float64)
            error_valid = np.asarray(sample["error_valid"], dtype=bool)
            time_valid = np.asarray(sample["time_valid"], dtype=bool)
            local_time = np.asarray(sample["local_time_cadences"], dtype=np.float64)
            delta_time = np.asarray(sample["delta_time_cadences"], dtype=np.float64)
            boundary = np.asarray(sample["segment_boundary"], dtype=bool)
            quality = np.asarray(
                [
                    flux_valid[center, 2],
                    flux_valid[center, 3],
                    error_valid[center, 0],
                    error_valid[center, 1],
                    local_time[center],
                    delta_time[center],
                    time_valid[center],
                    boundary[center],
                ],
                dtype=np.float64,
            )
            if quality_only:
                feature = quality
                expected = len(QUALITY_FEATURE_CHANNELS)
            else:
                feature = np.concatenate(
                    (
                        flux[center, (2, 3)],
                        flux_valid[center, (2, 3)].astype(np.float64),
                        error[center, (0, 1)],
                        error_valid[center, (0, 1)].astype(np.float64),
                        np.asarray(
                            [
                                local_time[center],
                                delta_time[center],
                                time_valid[center],
                                boundary[center],
                            ],
                            dtype=np.float64,
                        ),
                    )
                )
                expected = len(RAW_FEATURE_CHANNELS)
            if feature.shape != (expected,) or not np.all(np.isfinite(feature)):
                raise FM0ContractError("exact-center control feature is invalid")
            rows.append(feature)
    return np.stack(rows)


def _sigmoid(values: np.ndarray) -> np.ndarray:
    logits = np.asarray(values, dtype=np.float64)
    result = np.empty_like(logits)
    nonnegative = logits >= 0.0
    result[nonnegative] = 1.0 / (1.0 + np.exp(-logits[nonnegative]))
    exponential = np.exp(logits[~nonnegative])
    result[~nonnegative] = exponential / (1.0 + exponential)
    return result


def fit_frozen_linear_probe(features: np.ndarray, labels: np.ndarray) -> LinearProbe:
    """Fit the exact no-tuning, full-batch Adam linear-logit probe."""

    values = np.asarray(features, dtype=np.float64)
    target = np.asarray(labels, dtype=np.float64)
    if (
        values.ndim != 2
        or target.shape != (values.shape[0],)
        or not np.all(np.isfinite(values))
        or not np.all(np.isin(target, (0.0, 1.0)))
        or not np.any(target == 0.0)
        or not np.any(target == 1.0)
    ):
        raise FM0ContractError("linear-probe training arrays are invalid")
    center = np.mean(values, axis=0)
    scale = np.std(values, axis=0)
    scale = np.where(np.isfinite(scale) & (scale > 1.0e-8), scale, 1.0)
    standardized = (values - center) / scale

    rng = np.random.default_rng(PROBE_SEED)
    weight = rng.normal(0.0, 1.0e-3, size=values.shape[1]).astype(np.float64)
    bias = 0.0
    m_w = np.zeros_like(weight)
    v_w = np.zeros_like(weight)
    m_b = 0.0
    v_b = 0.0
    beta1 = 0.9
    beta2 = 0.999
    epsilon = 1.0e-8
    history: list[float] = []
    for epoch in range(1, PROBE_EPOCHS + 1):
        logits = standardized @ weight + bias
        residual = _sigmoid(logits) - target
        grad_w = np.mean(residual[:, None] * standardized, axis=0) + PROBE_L2 * weight
        grad_b = float(np.mean(residual))
        m_w = beta1 * m_w + (1.0 - beta1) * grad_w
        v_w = beta2 * v_w + (1.0 - beta2) * np.square(grad_w)
        m_b = beta1 * m_b + (1.0 - beta1) * grad_b
        v_b = beta2 * v_b + (1.0 - beta2) * grad_b * grad_b
        weight -= (
            PROBE_LEARNING_RATE
            * (m_w / (1.0 - beta1**epoch))
            / (np.sqrt(v_w / (1.0 - beta2**epoch)) + epsilon)
        )
        bias -= (
            PROBE_LEARNING_RATE
            * (m_b / (1.0 - beta1**epoch))
            / (math.sqrt(v_b / (1.0 - beta2**epoch)) + epsilon)
        )
        if epoch == 1 or epoch % 20 == 0 or epoch == PROBE_EPOCHS:
            updated_logits = standardized @ weight + bias
            loss = float(
                np.mean(np.logaddexp(0.0, updated_logits) - target * updated_logits)
                + 0.5 * PROBE_L2 * np.dot(weight, weight)
            )
            history.append(loss)
    return LinearProbe(
        weight=weight,
        bias=float(bias),
        center=center,
        scale=scale,
        objective_history=tuple(history),
    )


def score_frozen_linear_probe(fit: LinearProbe, features: np.ndarray) -> np.ndarray:
    values = np.asarray(features, dtype=np.float64)
    if (
        values.ndim != 2
        or values.shape[1:] != fit.weight.shape
        or fit.center.shape != fit.weight.shape
        or fit.scale.shape != fit.weight.shape
        or not np.all(np.isfinite(values))
    ):
        raise FM0ContractError("linear-probe scoring arrays are invalid")
    return (values - fit.center) / fit.scale @ fit.weight + fit.bias


def _rank_auc(labels: np.ndarray, scores: np.ndarray) -> float:
    y = np.asarray(labels, dtype=np.int8)
    s = np.asarray(scores, dtype=np.float64)
    if y.shape != s.shape or not np.all(np.isfinite(s)):
        raise FM0ContractError("AUROC arrays are invalid")
    positive = int(np.count_nonzero(y == 1))
    negative = int(np.count_nonzero(y == 0))
    if positive == 0 or negative == 0:
        raise FM0ContractError("sample AUROC requires both classes")
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
    return (rank_sum - positive * (positive + 1) / 2.0) / (positive * negative)


def _paired_component_scores(
    labels: np.ndarray,
    scores: np.ndarray,
    components: Sequence[str],
) -> tuple[tuple[str, ...], np.ndarray, np.ndarray]:
    y = np.asarray(labels, dtype=np.int8)
    s = np.asarray(scores, dtype=np.float64)
    groups = np.asarray(tuple(components), dtype=object)
    if y.shape != s.shape or groups.shape != y.shape or not np.all(np.isfinite(s)):
        raise FM0ContractError("paired metric arrays do not align")
    unique = tuple(dict.fromkeys(str(value) for value in groups.tolist()))
    if len(unique) < 2:
        raise FM0ContractError("paired bootstrap requires at least two components")
    clean = np.empty(len(unique), dtype=np.float64)
    injected = np.empty(len(unique), dtype=np.float64)
    for index, component in enumerate(unique):
        selected = groups == component
        if (
            np.count_nonzero(selected & (y == 0)) != 1
            or np.count_nonzero(selected & (y == 1)) != 1
        ):
            raise FM0ContractError(
                "every component must contribute one clean/injected pair"
            )
        clean[index] = s[selected & (y == 0)][0]
        injected[index] = s[selected & (y == 1)][0]
    return unique, clean, injected


def _pair_ranking(clean: np.ndarray, injected: np.ndarray) -> float:
    delta = np.asarray(injected) - np.asarray(clean)
    return float(np.mean((delta > 0.0) + 0.5 * (delta == 0.0)))


def paired_component_bootstrap(
    labels: np.ndarray,
    scores: np.ndarray,
    components: Sequence[str],
) -> dict[str, Any]:
    """Compute sample AUROC and resample whole adjacent component pairs."""

    unique, clean, injected = _paired_component_scores(labels, scores, components)
    point_auc = _rank_auc(
        np.tile((0, 1), len(unique)), np.column_stack((clean, injected)).ravel()
    )
    point_ranking = _pair_ranking(clean, injected)
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    draws = rng.integers(
        0,
        len(unique),
        size=(BOOTSTRAP_REPLICATES, len(unique)),
        endpoint=False,
    )
    ranking_draws = np.mean(
        (injected[draws] > clean[draws]) + 0.5 * (injected[draws] == clean[draws]),
        axis=1,
    )
    auc_draws = np.empty(BOOTSTRAP_REPLICATES, dtype=np.float64)
    bootstrap_labels = np.tile((0, 1), len(unique))
    for index, draw in enumerate(draws):
        bootstrap_scores = np.column_stack((clean[draw], injected[draw])).ravel()
        auc_draws[index] = _rank_auc(bootstrap_labels, bootstrap_scores)
    return {
        "n_component_pairs": len(unique),
        "n_samples": 2 * len(unique),
        "sample_roc_auc": {
            "estimate": float(point_auc),
            "component_pair_bootstrap_95_interval": [
                float(np.quantile(auc_draws, 0.025)),
                float(np.quantile(auc_draws, 0.975)),
            ],
        },
        "paired_clean_injected_ranking_accuracy": {
            "estimate": float(point_ranking),
            "component_pair_bootstrap_95_interval": [
                float(np.quantile(ranking_draws, 0.025)),
                float(np.quantile(ranking_draws, 0.975)),
            ],
        },
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "bootstrap_seed": BOOTSTRAP_SEED,
    }


def paired_auc_delta_bootstrap(
    labels: np.ndarray,
    left_scores: np.ndarray,
    right_scores: np.ndarray,
    components: Sequence[str],
) -> dict[str, Any]:
    """Compute a paired AUROC difference with the same frozen resamples."""

    unique, left_clean, left_injected = _paired_component_scores(
        labels, left_scores, components
    )
    right_unique, right_clean, right_injected = _paired_component_scores(
        labels, right_scores, components
    )
    if unique != right_unique:
        raise FM0ContractError("paired AUROC arms have different component order")
    fixed_labels = np.tile((0, 1), len(unique))
    left_point = np.column_stack((left_clean, left_injected)).ravel()
    right_point = np.column_stack((right_clean, right_injected)).ravel()
    estimate = _rank_auc(fixed_labels, left_point) - _rank_auc(
        fixed_labels, right_point
    )
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    draws = rng.integers(
        0,
        len(unique),
        size=(BOOTSTRAP_REPLICATES, len(unique)),
        endpoint=False,
    )
    delta_draws = np.empty(BOOTSTRAP_REPLICATES, dtype=np.float64)
    for index, draw in enumerate(draws):
        left = np.column_stack((left_clean[draw], left_injected[draw])).ravel()
        right = np.column_stack((right_clean[draw], right_injected[draw])).ravel()
        delta_draws[index] = _rank_auc(fixed_labels, left) - _rank_auc(
            fixed_labels, right
        )
    return {
        "estimate": float(estimate),
        "component_pair_bootstrap_95_interval": [
            float(np.quantile(delta_draws, 0.025)),
            float(np.quantile(delta_draws, 0.975)),
        ],
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "bootstrap_seed": BOOTSTRAP_SEED,
    }


def _metric_breakdowns(
    scores: np.ndarray,
    *,
    labels: np.ndarray,
    components: Sequence[str],
    metadata: Sequence[Mapping[str, str]],
) -> dict[str, Any]:
    values = np.asarray(scores, dtype=np.float64)
    if values.shape != labels.shape or len(metadata) != labels.size:
        raise FM0ContractError("feature scores do not align with paired schedule")
    split_values = np.asarray([row["split"] for row in metadata], dtype=object)
    cohort_values = np.asarray([row["cohort"] for row in metadata], dtype=object)
    duration_values = np.asarray(
        [int(row["duration_cadences"]) for row in metadata], dtype=np.int16
    )
    depth_values = np.asarray(
        [f"{float(row['fractional_depth']):.2f}" for row in metadata], dtype=object
    )
    component_array = np.asarray(tuple(components), dtype=object)
    result: dict[str, Any] = {}
    for split in SPLIT_ORDER:
        split_mask = split_values == split
        by_cohort = {}
        for cohort in TEMPORAL_COHORTS:
            selected = split_mask & (cohort_values == cohort)
            by_cohort[cohort] = paired_component_bootstrap(
                labels[selected], values[selected], component_array[selected]
            )
        by_cell = {}
        for duration in EVENT_DURATIONS:
            for depth_text in EVENT_DEPTHS_TEXT:
                depth = f"{float(depth_text):.2f}"
                selected = (
                    split_mask & (duration_values == duration) & (depth_values == depth)
                )
                by_cell[f"duration_{duration}_depth_{depth}"] = (
                    paired_component_bootstrap(
                        labels[selected], values[selected], component_array[selected]
                    )
                )
        by_depth = {}
        for depth_text in EVENT_DEPTHS_TEXT:
            depth = f"{float(depth_text):.2f}"
            selected = split_mask & (depth_values == depth)
            by_depth[f"depth_{depth}"] = paired_component_bootstrap(
                labels[selected], values[selected], component_array[selected]
            )
        result[split] = {
            "overall": paired_component_bootstrap(
                labels[split_mask], values[split_mask], component_array[split_mask]
            ),
            "by_cohort": by_cohort,
            "by_duration_depth_cell": by_cell,
            "by_depth_aggregated_over_durations": by_depth,
        }
    return result


def _model_samples(
    pairs: Sequence[FrozenPair], *, variant: str
) -> list[dict[str, np.ndarray]]:
    samples: list[dict[str, np.ndarray]] = []
    for pair in pairs:
        for source in (pair.clean, pair.injected):
            window = {
                name: np.asarray(value)
                for name, value in source.items()
                if name != "segment_id"
            }
            sample = prepare_model_window(
                window,
                variant=variant,
                mask_seed=0,
                window_length=CONTEXT_CADENCES,
                mask_target_fraction=0.0,
            )
            if (
                sample["flux"].shape != (2, CONTEXT_CADENCES)
                or np.any(sample["temporal_mask"])
                or np.any(sample["reconstruction_mask"])
            ):
                raise FM0ContractError(
                    "model evaluation sample is not unmasked native-128"
                )
            samples.append(sample)
    return samples


def _device(value: str) -> Any:
    try:
        import torch
    except ImportError as exc:  # pragma: no cover - optional locally
        raise RuntimeError("PyTorch is required for FM0.3 checkpoint arms") from exc
    requested = str(value)
    if requested == "auto":
        requested = "cuda" if torch.cuda.is_available() else "cpu"
    if requested == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but is unavailable")
    if requested not in {"cpu", "cuda"}:
        raise ValueError("model evaluation device must be auto, cpu, or cuda")
    return torch.device(requested)


def _encode_model_arm(
    arm: ModelArm,
    pairs: Sequence[FrozenPair],
    *,
    payload_plan_root: Path,
    payload_plan_receipt_sha256: str,
    device: str,
    batch_size: int,
) -> tuple[np.ndarray, Mapping[str, Any]]:
    try:
        import torch
    except ImportError as exc:  # pragma: no cover - optional locally
        raise RuntimeError("PyTorch is required for FM0.3 checkpoint arms") from exc
    from .centered_event_context_diagnostic import _load_inference_only_trusted_model

    try:
        variant, expected_step, architecture = MODEL_ARM_CONTRACT[arm.name]
    except KeyError as exc:
        raise FM0ContractError(
            f"unknown FM0.3 model evaluation arm: {arm.name}"
        ) from exc
    target_device = _device(device)
    model, contract, validation = _load_inference_only_trusted_model(
        Path(arm.run_dir),
        device=target_device,
        checkpoint_path=Path(arm.checkpoint_path),
    )
    evaluation_binding = contract.get("evaluation_plan")
    if (
        contract.get("variant") != variant
        or contract.get("architecture") != architecture
        or validation.get("variant") != variant
        or validation.get("architecture") != architecture
        or validation.get("global_step") != expected_step
        or not isinstance(evaluation_binding, Mapping)
        or Path(str(evaluation_binding.get("root", ""))).resolve(strict=True)
        != payload_plan_root
        or evaluation_binding.get("receipt_sha256") != payload_plan_receipt_sha256
        or int(model.config.window_length) != CONTEXT_CADENCES
        or int(model.config.patch_stride) != 1
    ):
        raise FM0ContractError(
            "model arm differs from the matched payload/checkpoint contract"
        )
    samples = _model_samples(pairs, variant=variant)
    if batch_size <= 0:
        raise ValueError("model batch_size must be positive")
    rows: list[np.ndarray] = []
    model.eval()
    with torch.no_grad():
        for start in range(0, len(samples), batch_size):
            batch = move_batch_to_device(
                collate_fm0_samples(samples[start : start + batch_size]),
                target_device,
            )
            if bool(torch.any(batch["temporal_mask"])) or bool(
                torch.any(batch["reconstruction_mask"])
            ):
                raise FM0ContractError("checkpoint inference received a nonzero mask")
            output = model(batch)
            hidden = output.get("h_cadence")
            token_valid = output.get("token_valid")
            if (
                hidden is None
                or token_valid is None
                or tuple(hidden.shape[1:])
                != (CONTEXT_CADENCES, int(model.config.d_model))
                or tuple(token_valid.shape[1:]) != (CONTEXT_CADENCES,)
                or not bool(torch.all(token_valid[:, EVENT_CENTER_INDEX_ZERO_BASED]))
                or not bool(torch.isfinite(hidden).all())
            ):
                raise FM0ContractError(
                    "checkpoint did not emit valid native-128 center tokens"
                )
            rows.append(
                hidden[:, EVENT_CENTER_INDEX_ZERO_BASED, :]
                .detach()
                .float()
                .cpu()
                .numpy()
            )
    features = np.concatenate(rows).astype(np.float64, copy=False)
    if features.shape != (2 * len(pairs), int(model.config.d_model)):
        raise FM0ContractError("model center-token feature matrix shape drifted")
    return features, {
        "name": arm.name,
        "variant": variant,
        "architecture": architecture,
        "expected_checkpoint_step": expected_step,
        "run_dir": str(Path(arm.run_dir).resolve(strict=True)),
        "checkpoint_path": str(Path(arm.checkpoint_path).resolve(strict=True)),
        "checkpoint_sha256": validation["selected_checkpoint_sha256"],
        "run_contract_sha256": validation["artifact_sha256"]["run_contract.json"],
        "payload_plan_receipt_sha256": payload_plan_receipt_sha256,
        "device": str(target_device),
        "all_128_tokens_emitted": True,
        "selected_token_index_zero_based": EVENT_CENTER_INDEX_ZERO_BASED,
        "temporal_or_reconstruction_mask_nonzero": False,
    }


def _control_gates(feature_results: Mapping[str, Any]) -> dict[str, Any]:
    gates = _frozen_analysis_contract()["primary_test_gates"]
    raw = feature_results[CONTROL_ARMS[0]]["metrics"]["fresh_s77_test"]["overall"]
    quality = feature_results[CONTROL_ARMS[1]]["metrics"]["fresh_s77_test"]["overall"]
    raw_auc = raw["sample_roc_auc"]
    quality_auc = quality["sample_roc_auc"]
    quality_interval = quality_auc["component_pair_bootstrap_95_interval"]
    criteria = {
        "raw_overall_sample_roc_auc_lower_95_at_least_0_90": (
            float(raw_auc["component_pair_bootstrap_95_interval"][0])
            >= float(gates["raw_control"]["minimum_overall_sample_roc_auc_lower_95"])
        ),
        "quality_overall_sample_roc_auc_within_0_03_of_chance": (
            abs(float(quality_auc["estimate"]) - 0.5)
            <= float(
                gates["quality_only_control"][
                    "maximum_absolute_overall_sample_roc_auc_from_chance"
                ]
            )
        ),
        "quality_overall_sample_roc_auc_interval_contains_chance": (
            float(quality_interval[0]) <= 0.5 <= float(quality_interval[1])
        ),
    }
    return {
        "passed": all(criteria.values()),
        "criteria": criteria,
        "raw_control_overall_sample_roc_auc": raw_auc,
        "quality_only_overall_sample_roc_auc": quality_auc,
    }


def _model_gates(
    variant: str,
    *,
    feature_results: Mapping[str, Any],
    scores: Mapping[str, np.ndarray],
    labels: np.ndarray,
    components: Sequence[str],
    metadata: Sequence[Mapping[str, str]],
) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
    step0 = f"{variant}_step0_h_cadence_token64"
    step2000 = f"{variant}_step2000_h_cadence_token64"
    if step0 not in feature_results or step2000 not in feature_results:
        return None, None
    test_mask = np.asarray([row["split"] == "fresh_s77_test" for row in metadata])
    test_components = np.asarray(tuple(components), dtype=object)[test_mask]
    delta = paired_auc_delta_bootstrap(
        labels[test_mask],
        scores[step2000][test_mask],
        scores[step0][test_mask],
        test_components,
    )
    gate = _frozen_analysis_contract()["primary_test_gates"][
        "each_architecture_step2000"
    ]
    metrics = feature_results[step2000]["metrics"]["fresh_s77_test"]
    overall = metrics["overall"]["sample_roc_auc"]
    by_cohort = metrics["by_cohort"]
    by_depth = metrics["by_depth_aggregated_over_durations"]
    criteria = {
        "overall_sample_roc_auc_lower_95_at_least_0_75": (
            float(overall["component_pair_bootstrap_95_interval"][0])
            >= float(gate["minimum_overall_sample_roc_auc_lower_95"])
        ),
        "each_cohort_sample_roc_auc_estimate_at_least_0_70": all(
            float(by_cohort[name]["sample_roc_auc"]["estimate"])
            >= float(gate["minimum_each_cohort_sample_roc_auc_estimate"])
            for name in TEMPORAL_COHORTS
        ),
        "each_cohort_sample_roc_auc_lower_95_strictly_above_chance": all(
            float(
                by_cohort[name]["sample_roc_auc"][
                    "component_pair_bootstrap_95_interval"
                ][0]
            )
            > float(gate["minimum_each_cohort_sample_roc_auc_lower_95_strictly_above"])
            for name in TEMPORAL_COHORTS
        ),
        "step2000_minus_own_step0_auc_estimate_at_least_0_02": (
            float(delta["estimate"])
            >= float(gate["paired_step2000_minus_own_step0_auc"]["minimum_estimate"])
        ),
        "step2000_minus_own_step0_auc_lower_95_strictly_above_zero": (
            float(delta["component_pair_bootstrap_95_interval"][0])
            > float(
                gate["paired_step2000_minus_own_step0_auc"]["lower_95_strictly_above"]
            )
        ),
        "depth_0_10_and_0_30_auc_lower_95_at_least_0_80": all(
            float(
                by_depth[f"depth_{depth:.2f}"]["sample_roc_auc"][
                    "component_pair_bootstrap_95_interval"
                ][0]
            )
            >= float(
                gate["blocking_depth_aggregates"][
                    "minimum_each_depth_sample_roc_auc_lower_95"
                ]
            )
            for depth in gate["blocking_depth_aggregates"]["fractional_depths"]
        ),
    }
    return {"passed": all(criteria.values()), "criteria": criteria}, delta


def _scores_csv(
    arm_order: Sequence[str],
    scores: Mapping[str, np.ndarray],
    metadata: Sequence[Mapping[str, str]],
) -> bytes:
    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(SCORE_FIELDS), lineterminator="\n")
    writer.writeheader()
    for arm in arm_order:
        arm_scores = np.asarray(scores[arm], dtype=np.float64)
        if arm_scores.shape != (len(metadata),):
            raise FM0ContractError("score vector does not align with schedule pairs")
        for index, (row, score) in enumerate(zip(metadata, arm_scores, strict=True)):
            label = index % 2
            writer.writerow(
                {
                    "evaluation_schema_version": EVALUATION_SCHEMA_VERSION,
                    "feature_arm": arm,
                    "split": row["split"],
                    "cohort": row["cohort"],
                    "sector": row["sector"],
                    "observation_key": row["observation_key"],
                    "leakage_component_id": row["leakage_component_id"],
                    "event_pair_id": row["event_pair_id"],
                    "duration_cadences": row["duration_cadences"],
                    "fractional_depth": row["fractional_depth"],
                    "sample_role": "clean" if label == 0 else "injected",
                    "label": str(label),
                    "score": f"{float(score):.17g}",
                }
            )
    return stream.getvalue().encode("utf-8")


def _make_tree_read_only(root: Path) -> None:
    for path in root.rglob("*"):
        path.chmod(0o555 if path.is_dir() else 0o444)
    root.chmod(0o555)


def validate_matched_canary_evaluation(
    root: str | Path,
    *,
    expected_receipt_sha256: str,
    require_read_only: bool = True,
) -> MatchedCanaryEvaluationResult:
    """Validate the compact immutable result without reopening source shards."""

    output_raw = Path(root).expanduser()
    if output_raw.is_symlink():
        raise FM0ContractError("evaluation result root must be materialized")
    output = output_raw.resolve(strict=True)
    expected_files = {"results.json", "scores.csv", "receipt.json", "READY"}
    if (
        not output.is_dir()
        or {path.name for path in output.iterdir()} != expected_files
    ):
        raise FM0ContractError("evaluation result closure drifted")
    results_path = output / "results.json"
    scores_path = output / "scores.csv"
    receipt_path = output / "receipt.json"
    ready_path = output / "READY"
    for path in (results_path, scores_path, receipt_path, ready_path):
        if path.is_symlink() or not path.is_file():
            raise FM0ContractError("evaluation result artifact must be materialized")
        if require_read_only and path.stat().st_mode & 0o222:
            raise FM0ContractError("evaluation result artifact must be read-only")
    if require_read_only and output.stat().st_mode & 0o222:
        raise FM0ContractError("evaluation result root must be read-only")
    receipt_hash = sha256_file(receipt_path)
    if receipt_hash != _digest(
        expected_receipt_sha256, label="evaluation receipt hash"
    ):
        raise FM0ContractError("evaluation receipt SHA-256 drifted")
    if ready_path.read_text(encoding="utf-8").strip() != receipt_hash:
        raise FM0ContractError("evaluation READY/receipt binding drifted")
    try:
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        results = json.loads(results_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FM0ContractError("evaluation JSON artifact is invalid") from exc
    if (
        not isinstance(receipt, Mapping)
        or receipt.get("schema_version") != EVALUATION_RECEIPT_SCHEMA_VERSION
        or receipt.get("ready_state") != EVALUATION_READY_STATE
        or receipt.get("artifact_validation_passed") is not True
        or _GIT_SHA.fullmatch(str(receipt.get("producer_git_sha", ""))) is None
        or receipt.get("algorithm_contract") != evaluation_algorithm_contract()
    ):
        raise FM0ContractError("evaluation receipt contract drifted")
    artifacts = receipt.get("artifacts")
    results_hash = sha256_file(results_path)
    scores_hash = sha256_file(scores_path)
    if artifacts != {
        "results": {
            "path": "results.json",
            "sha256": results_hash,
        },
        "scores": {
            "path": "scores.csv",
            "sha256": scores_hash,
            "n_rows": receipt.get("n_score_rows"),
        },
    }:
        raise FM0ContractError("evaluation artifact hashes drifted")
    payload = receipt.get("payload_plan")
    if not isinstance(payload, Mapping):
        raise FM0ContractError("evaluation lacks its payload-plan binding")
    try:
        validated_payload = validate_matched_canary_payload_plan(
            payload["root"],
            expected_receipt_sha256=payload["receipt_sha256"],
            require_read_only=require_read_only,
        )
    except (KeyError, OSError, TypeError, ValueError) as exc:
        raise FM0ContractError("evaluation payload-plan binding failed") from exc
    temporal_panel = _temporal_panel_authority(validated_payload.receipt)
    expected_payload = {
        "root": str(validated_payload.root),
        "receipt_sha256": validated_payload.receipt_sha256,
        "schedule_sha256": validated_payload.schedule_sha256,
        "producer_git_sha": validated_payload.receipt["producer_git_sha"],
        "source_shard_bindings_sha256": validated_payload.receipt["payload_bindings"][
            "source_shard_bindings_sha256"
        ],
        "crop_payload_bindings_sha256": validated_payload.receipt["payload_bindings"][
            "crop_payload_bindings_sha256"
        ],
        "temporal_panel_authority": temporal_panel,
        "n_crops": 1_440,
    }
    if dict(payload) != expected_payload:
        raise FM0ContractError("evaluation payload-plan aggregate binding drifted")
    if (
        not isinstance(results, Mapping)
        or results.get("schema_version") != EVALUATION_SCHEMA_VERSION
        or results.get("producer_git_sha") != receipt["producer_git_sha"]
        or results.get("payload_plan") != expected_payload
        or results.get("algorithm_contract") != receipt["algorithm_contract"]
        or results.get("gate_summary") != receipt.get("gate_summary")
    ):
        raise FM0ContractError("evaluation result/receipt contract drifted")
    arms = receipt.get("feature_arms")
    if (
        not isinstance(arms, list)
        or tuple(arms[:2]) != CONTROL_ARMS
        or len(set(arms)) != len(arms)
        or any(name not in MODEL_ARM_CONTRACT for name in arms[2:])
        or set(results.get("feature_results", {})) != set(arms)
        or receipt.get("n_score_rows") != 2 * 1_440 * len(arms)
    ):
        raise FM0ContractError("evaluation feature-arm closure drifted")
    with scores_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != SCORE_FIELDS:
            raise FM0ContractError("evaluation score columns drifted")
        score_rows = [dict(row) for row in reader]
    if len(score_rows) != receipt["n_score_rows"]:
        raise FM0ContractError("evaluation score row count drifted")
    counts = {arm: 0 for arm in arms}
    for row in score_rows:
        arm = row["feature_arm"]
        if (
            row["evaluation_schema_version"] != EVALUATION_SCHEMA_VERSION
            or arm not in counts
            or row["sample_role"] not in {"clean", "injected"}
            or row["label"] != ("0" if row["sample_role"] == "clean" else "1")
            or row["split"] not in SPLIT_ORDER
            or row["cohort"] not in TEMPORAL_COHORTS
            or not math.isfinite(float(row["score"]))
        ):
            raise FM0ContractError("evaluation score row contract drifted")
        counts[arm] += 1
    if any(count != 2 * 1_440 for count in counts.values()):
        raise FM0ContractError("evaluation score arm count drifted")
    limits = receipt.get("limits")
    if limits != {
        "source_shards_reopened": 1_440,
        "events_injected": 1_440,
        "clean_injected_samples_scored_per_arm": 2_880,
        "bls_or_search_executed": False,
        "sealed_test_opened": False,
        "fm_encoder_training_executed": False,
        "hyperparameter_or_threshold_tuning_executed": False,
    }:
        raise FM0ContractError("evaluation execution boundary drifted")
    return MatchedCanaryEvaluationResult(
        root=output,
        results_path=results_path,
        scores_path=scores_path,
        receipt_path=receipt_path,
        ready_path=ready_path,
        results_sha256=results_hash,
        scores_sha256=scores_hash,
        receipt_sha256=receipt_hash,
        receipt=receipt,
    )


def evaluate_matched_canary(
    output_dir: str | Path,
    *,
    payload_plan_root: str | Path,
    payload_plan_receipt_sha256: str,
    producer_git_sha: str,
    model_arms: Sequence[ModelArm] = (),
    device: str = "auto",
    batch_size: int = 32,
    require_read_only: bool = True,
) -> MatchedCanaryEvaluationResult:
    """Run controls and any explicitly supplied immutable checkpoint arms."""

    producer = _git_sha(producer_git_sha)
    payload_hash = _digest(
        payload_plan_receipt_sha256, label="payload-plan receipt hash"
    )
    payload = validate_matched_canary_payload_plan(
        payload_plan_root,
        expected_receipt_sha256=payload_hash,
        require_read_only=require_read_only,
    )
    temporal_panel = _temporal_panel_authority(payload.receipt)
    if payload.receipt["producer_git_sha"] != producer:
        raise FM0ContractError("payload plan was not frozen by this exact Git revision")
    requested = tuple(model_arms)
    names = tuple(arm.name for arm in requested)
    if len(set(names)) != len(names) or any(
        name not in MODEL_ARM_CONTRACT for name in names
    ):
        raise FM0ContractError("model evaluation arms are duplicated or unknown")
    output = Path(output_dir).expanduser().absolute()
    if output.exists():
        ready = output / "READY"
        if not ready.is_file():
            raise FM0ContractError("existing evaluation result has no READY binding")
        result = validate_matched_canary_evaluation(
            output,
            expected_receipt_sha256=ready.read_text(encoding="utf-8").strip(),
            require_read_only=require_read_only,
        )
        if (
            result.receipt["producer_git_sha"] != producer
            or result.receipt["payload_plan"]["receipt_sha256"] != payload_hash
            or tuple(result.receipt["feature_arms"][2:]) != names
        ):
            raise FM0ContractError(
                "existing evaluation result binds another invocation"
            )
        return result
    if output.is_symlink() or output.name in {"", ".", ".."}:
        raise FM0ContractError("evaluation output must be a named directory")

    schedule = _read_schedule(payload.schedule_path)
    pairs = load_frozen_pairs(schedule, require_read_only=require_read_only)
    labels, components, metadata = _paired_metadata(pairs)
    train_mask = np.asarray([row["split"] == "probe_train" for row in metadata])
    features: dict[str, np.ndarray] = {
        CONTROL_ARMS[0]: exact_center_control_features(pairs, quality_only=False),
        CONTROL_ARMS[1]: exact_center_control_features(pairs, quality_only=True),
    }
    checkpoint_bindings: dict[str, Mapping[str, Any]] = {}
    for arm in requested:
        print(f"FM0_3_EVALUATION phase=encode feature={arm.name}", flush=True)
        encoded, binding = _encode_model_arm(
            arm,
            pairs,
            payload_plan_root=payload.root,
            payload_plan_receipt_sha256=payload.receipt_sha256,
            device=device,
            batch_size=batch_size,
        )
        features[arm.name] = encoded
        checkpoint_bindings[arm.name] = binding

    feature_results: dict[str, Any] = {}
    all_scores: dict[str, np.ndarray] = {}
    arm_order = (*CONTROL_ARMS, *names)
    for arm_name in arm_order:
        print(f"FM0_3_EVALUATION phase=fit feature={arm_name}", flush=True)
        fit = fit_frozen_linear_probe(
            features[arm_name][train_mask], labels[train_mask]
        )
        arm_scores = score_frozen_linear_probe(fit, features[arm_name])
        all_scores[arm_name] = arm_scores
        feature_results[arm_name] = {
            "feature_dimension": int(features[arm_name].shape[1]),
            "probe": {
                "weight": [float(value) for value in fit.weight],
                "bias": float(fit.bias),
                "train_center": [float(value) for value in fit.center],
                "train_scale": [float(value) for value in fit.scale],
                "objective_history": [float(value) for value in fit.objective_history],
            },
            "metrics": _metric_breakdowns(
                arm_scores,
                labels=labels,
                components=components,
                metadata=metadata,
            ),
        }
        if arm_name in checkpoint_bindings:
            feature_results[arm_name]["checkpoint"] = checkpoint_bindings[arm_name]
        print(f"FM0_3_EVALUATION phase=complete feature={arm_name}", flush=True)

    controls = _control_gates(feature_results)
    model_gates: dict[str, Any] = {}
    paired_deltas: dict[str, Any] = {}
    for variant in ("TWIRL-FM0.3.1", "TWIRL-FM0.3.2"):
        gate, delta = _model_gates(
            variant,
            feature_results=feature_results,
            scores=all_scores,
            labels=labels,
            components=components,
            metadata=metadata,
        )
        model_gates[variant] = gate
        if delta is not None:
            paired_deltas[f"{variant}_step2000_minus_own_step0"] = delta
    left = "TWIRL-FM0.3.2_step2000_h_cadence_token64"
    right = "TWIRL-FM0.3.1_step2000_h_cadence_token64"
    if left in all_scores and right in all_scores:
        test_mask = np.asarray([row["split"] == "fresh_s77_test" for row in metadata])
        paired_deltas["TWIRL-FM0.3.2_minus_TWIRL-FM0.3.1_step2000"] = (
            paired_auc_delta_bootstrap(
                labels[test_mask],
                all_scores[left][test_mask],
                all_scores[right][test_mask],
                np.asarray(tuple(components), dtype=object)[test_mask],
            )
        )
    complete_model_set = set(names) == set(MODEL_ARM_CONTRACT)
    architecture_interpretable = bool(
        complete_model_set
        and controls["passed"]
        and all(
            isinstance(model_gates[variant], Mapping) and model_gates[variant]["passed"]
            for variant in ("TWIRL-FM0.3.1", "TWIRL-FM0.3.2")
        )
    )
    gate_summary = {
        "evaluation_scope": "controls_only" if not names else "controls_and_model_arms",
        "controls_preflight": controls,
        "all_four_model_arms_present": complete_model_set,
        "model_step2000_gates": model_gates,
        "paired_test_auc_deltas": paired_deltas,
        "architecture_comparison_interpretable": architecture_interpretable,
        "one_seed_pass_authorizes_only": (
            "longer_matched_two_seed_architecture_comparison"
            if architecture_interpretable
            else None
        ),
        "architecture_promotion_authorized": False,
        "foundation_model_claim_authorized": False,
    }
    payload_binding = {
        "root": str(payload.root),
        "receipt_sha256": payload.receipt_sha256,
        "schedule_sha256": payload.schedule_sha256,
        "producer_git_sha": payload.receipt["producer_git_sha"],
        "source_shard_bindings_sha256": payload.receipt["payload_bindings"][
            "source_shard_bindings_sha256"
        ],
        "crop_payload_bindings_sha256": payload.receipt["payload_bindings"][
            "crop_payload_bindings_sha256"
        ],
        "temporal_panel_authority": temporal_panel,
        "n_crops": 1_440,
    }
    algorithm = evaluation_algorithm_contract()
    results = {
        "schema_version": EVALUATION_SCHEMA_VERSION,
        "evaluation_completed": True,
        "producer_git_sha": producer,
        "payload_plan": payload_binding,
        "algorithm_contract": algorithm,
        "feature_results": feature_results,
        "gate_summary": gate_summary,
        "checkpoint_bindings": checkpoint_bindings,
        "limits": {
            "source_shards_reopened": 1_440,
            "events_injected": 1_440,
            "clean_injected_samples_scored_per_arm": 2_880,
            "bls_or_search_executed": False,
            "sealed_test_opened": False,
            "fm_encoder_training_executed": False,
            "hyperparameter_or_threshold_tuning_executed": False,
        },
    }
    scores_payload = _scores_csv(arm_order, all_scores, metadata)
    results_payload = _canonical_json_bytes(results)
    results_hash = hashlib.sha256(results_payload).hexdigest()
    scores_hash = hashlib.sha256(scores_payload).hexdigest()
    limits = dict(results["limits"])
    receipt = {
        "schema_version": EVALUATION_RECEIPT_SCHEMA_VERSION,
        "ready_state": EVALUATION_READY_STATE,
        "artifact_validation_passed": True,
        "producer_git_sha": producer,
        "payload_plan": payload_binding,
        "algorithm_contract": algorithm,
        "feature_arms": list(arm_order),
        "model_checkpoint_bindings": checkpoint_bindings,
        "gate_summary": gate_summary,
        "n_score_rows": 2 * 1_440 * len(arm_order),
        "artifacts": {
            "results": {"path": "results.json", "sha256": results_hash},
            "scores": {
                "path": "scores.csv",
                "sha256": scores_hash,
                "n_rows": 2 * 1_440 * len(arm_order),
            },
        },
        "limits": limits,
        "claim_limit": (
            "Development-only candidate-centred synthetic-event evaluation. "
            "A controls-only pass validates probe mechanics but does not score an FM. "
            "A one-seed full pass can authorize only a longer matched two-seed "
            "comparison; it cannot promote an architecture or support a foundation-"
            "model, discovery, or completeness claim."
        ),
    }
    receipt_payload = _canonical_json_bytes(receipt)
    receipt_hash = hashlib.sha256(receipt_payload).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = Path(
        tempfile.mkdtemp(prefix=f".{output.name}.partial.", dir=output.parent)
    )
    try:
        (partial / "results.json").write_bytes(results_payload)
        (partial / "scores.csv").write_bytes(scores_payload)
        (partial / "receipt.json").write_bytes(receipt_payload)
        (partial / "READY").write_text(receipt_hash + "\n", encoding="utf-8")
        _make_tree_read_only(partial)
        os.replace(partial, output)
    except BaseException:
        if partial.exists():
            partial.chmod(0o755)
            for path in partial.rglob("*"):
                path.chmod(0o755 if path.is_dir() else 0o644)
            shutil.rmtree(partial)
        raise
    return validate_matched_canary_evaluation(
        output,
        expected_receipt_sha256=receipt_hash,
        require_read_only=require_read_only,
    )


__all__ = [
    "CONTROL_ARMS",
    "EVALUATION_RECEIPT_SCHEMA_VERSION",
    "EVALUATION_SCHEMA_VERSION",
    "MODEL_ARM_CONTRACT",
    "FrozenPair",
    "LinearProbe",
    "MatchedCanaryEvaluationResult",
    "ModelArm",
    "evaluate_matched_canary",
    "evaluation_algorithm_contract",
    "exact_center_control_features",
    "fit_frozen_linear_probe",
    "load_frozen_pairs",
    "paired_auc_delta_bootstrap",
    "paired_component_bootstrap",
    "score_frozen_linear_probe",
    "validate_matched_canary_evaluation",
]
