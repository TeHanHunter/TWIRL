"""Training-free centered-event context diagnostic for TWIRL-FM0.2.1.

The diagnostic measures how one visible, artificial, short event changes the
frozen FM representation and decoder response when the same real development
interval is presented at 128, 256, 512, and 2,048 cadences.  It deliberately
contains no period, BLS, label, candidate, probe, optimizer, or training path.

Every injected sample contains exactly one event at the designated center
cadence.  The event truth is retained only by this evaluator; model batches
continue to use the ordinary FM0 tensor allowlist.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import fields
from pathlib import Path
from typing import Any

import numpy as np

from .dataset import (
    MODEL_WINDOW_KEYS,
    VIEW_NAMES,
    collate_fm0_samples,
    move_batch_to_device,
    prepare_model_window,
    variant_view_indices,
)
from .input_release import (
    WindowSpec,
    deterministic_training_window,
    extract_window,
    load_input_release_bytes,
)
from .representation_health import (
    _identity_sha256,
    _model_projection_spectrum,
    cosine_similarity_rows,
    source_clustered_mean_interval,
)
from .temporal_panel import DEVELOPMENT_PARTITION
from .temporal_zero_shot import (
    EXPECTED_CHECKPOINT_STEPS,
    TEMPORAL_COHORTS,
    _sample_set_sha256,
    load_temporal_panel,
    select_temporal_rows,
)
from .validation import read_json, sha256_file, validate_real_run_release

CENTERED_EVENT_CONTEXT_CONFIG_SCHEMA_VERSION = (
    "twirl_fm0_centered_event_context_config_v1"
)
CENTERED_EVENT_CONTEXT_RESULT_SCHEMA_VERSION = (
    "twirl_fm0_centered_event_context_diagnostic_v1"
)
CENTERED_EVENT_SELECTION_SALT = "twirl_fm0_centered_event_context_v1"
BASE_INTERVAL_CADENCES = 2_048
CONTEXT_LENGTHS = (128, 256, 512, 2_048)
EVENT_DURATIONS_CADENCES = (1, 3, 9)
EVENT_FRACTIONAL_DEPTHS = (0.01, 0.03, 0.1, 0.3)
EVENT_INGRESS_EGRESS_FRACTION = 1.0 / 3.0
COMPONENTS_PER_COHORT = 48
MINIMUM_VALID_FRACTION = 0.8
_REPRESENTATIONS = ("h_window", "z_window")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
# Exact schema carried by the frozen FM0.2.1 checkpoint artifacts.  Keep this
# local so the inference-only evaluator never imports the training module.
_TRUSTED_CHECKPOINT_SCHEMA_VERSION = "twirl_fm0_1_checkpoint_v1"

_CADENCE_KEYS = (
    "flux",
    "flux_valid",
    "flux_error",
    "error_valid",
    "local_time_cadences",
    "delta_time_cadences",
    "time_valid",
    "segment_boundary",
    "temporal_mask",
    "reconstruction_mask",
)
_MODEL_SAMPLE_KEYS = frozenset(
    (*MODEL_WINDOW_KEYS, "temporal_mask", "reconstruction_mask")
)
_AUTHORIZED_TRUE = (
    "development_shard_payload_access",
    "artificial_single_event_injection",
    "checkpoint_inference",
)
_AUTHORIZED_FALSE = (
    "labels_or_candidates",
    "bls_or_period_features",
    "repeated_synthetic_events",
    "classifier_or_probe_fitting",
    "optimizer_or_model_training",
    "sealed_test_access",
    "formal_model_gate",
    "production_model_claim",
)


class CenteredEventIneligibleError(ValueError):
    """A development visit cannot support the frozen centered experiment."""


def _seed(*parts: object) -> int:
    payload = "\x1f".join(str(part) for part in parts).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


def _mapping(value: Any, *, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be a mapping")
    return value


def _exact_sha256(value: Any, *, label: str) -> str:
    digest = str(value).strip()
    if _SHA256.fullmatch(digest) is None:
        raise ValueError(f"{label} must be a lowercase SHA-256 digest")
    return digest


def _exact_bool(mapping: Mapping[str, Any], key: str, expected: bool) -> None:
    if mapping.get(key) is not expected:
        raise ValueError(f"frozen config field {key!r} must be {expected}")


def _exact_int(value: Any, *, label: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise TypeError(f"{label} must be an integer")
    return int(value)


def _load_frozen_config(
    config_path: str | Path,
) -> tuple[dict[str, Any], Path, str]:
    """Load, validate, and byte-bind the frozen diagnostic YAML."""

    try:
        import yaml
    except ImportError as exc:  # pragma: no cover - PyYAML is a project dependency
        raise RuntimeError("PyYAML is required for the centered-event config") from exc

    source = Path(config_path).expanduser().resolve(strict=True)
    payload_bytes = source.read_bytes()
    payload = yaml.safe_load(payload_bytes)
    config = dict(_mapping(payload, label="centered-event config"))
    if config.get("schema_version") != CENTERED_EVENT_CONTEXT_CONFIG_SCHEMA_VERSION:
        raise ValueError("centered-event config schema version differs")
    if config.get("model_family") != "TWIRL-FM0" or config.get("status") != (
        "frozen_training_free_development_diagnostic"
    ):
        raise ValueError("centered-event config identity differs")
    if config.get("campaign_id") != "twirl_fm0_2_s66_s77_centered_event_context_v1":
        raise ValueError("centered-event campaign identifier differs")

    purpose = _mapping(config.get("purpose"), label="centered-event purpose")
    if purpose.get("primary") != (
        "measure_single_short_event_information_retained_by_frozen_fm_representations"
    ) or purpose.get("downstream_goal") != (
        "support_later_bls_free_supervised_classification_and_triage"
    ):
        raise ValueError("centered-event purpose differs")
    _exact_bool(purpose, "diagnostic_only", True)
    authorization = _mapping(
        config.get("authorization"), label="centered-event authorization"
    )
    for key in _AUTHORIZED_TRUE:
        _exact_bool(authorization, key, True)
    for key in _AUTHORIZED_FALSE:
        _exact_bool(authorization, key, False)

    inputs = _mapping(config.get("inputs"), label="centered-event inputs")
    _exact_sha256(
        inputs.get("temporal_panel_receipt_sha256"),
        label="temporal-panel receipt SHA-256",
    )
    if tuple(inputs.get("allowed_sectors", ())) != tuple(range(66, 78)):
        raise ValueError("centered-event allowed sectors differ")
    if inputs.get("excluded_sectors") != [65]:
        raise ValueError("centered-event excluded sectors differ")
    if inputs.get("source_partition") != DEVELOPMENT_PARTITION:
        raise ValueError("centered-event source partition differs")
    if tuple(inputs.get("checkpoint_steps", ())) != EXPECTED_CHECKPOINT_STEPS:
        raise ValueError("centered-event checkpoint steps differ")
    checkpoint_hashes = _mapping(
        inputs.get("checkpoint_sha256"), label="centered-event checkpoint hashes"
    )
    for key in ("step0", "step2000"):
        _exact_sha256(checkpoint_hashes.get(key), label=f"{key} checkpoint SHA-256")

    selection = _mapping(config.get("selection"), label="centered-event selection")
    if tuple(selection.get("cohorts", ())) != TEMPORAL_COHORTS:
        raise ValueError("centered-event cohorts differ")
    if (
        _exact_int(
            selection.get("max_components_per_cohort"),
            label="max components per cohort",
        )
        != COMPONENTS_PER_COHORT
    ):
        raise ValueError("centered-event component count differs")
    if (
        _exact_int(
            selection.get("base_interval_cadences"), label="base interval cadences"
        )
        != BASE_INTERVAL_CADENCES
    ):
        raise ValueError("centered-event base interval differs")
    if not math.isclose(
        float(selection.get("minimum_valid_fraction_each_required_view", math.nan)),
        MINIMUM_VALID_FRACTION,
        rel_tol=0.0,
        abs_tol=0.0,
    ):
        raise ValueError("centered-event minimum valid fraction differs")
    for key in (
        "one_visit_per_component",
        "require_full_validity_over_maximum_event_support",
        "only_structural_context_ineligibility_is_skippable",
        "path_checksum_schema_and_partition_failures_are_fatal",
    ):
        _exact_bool(selection, key, True)
    if selection.get("component_and_visit_order") != "deterministic_sha256":
        raise ValueError("centered-event selection order differs")

    contexts = _mapping(config.get("contexts"), label="centered-event contexts")
    if tuple(contexts.get("cadence_lengths", ())) != CONTEXT_LENGTHS:
        raise ValueError("centered-event context grid differs")
    if (
        _exact_int(
            contexts.get("nominal_cadence_seconds"), label="nominal cadence seconds"
        )
        != 200
    ):
        raise ValueError("centered-event nominal cadence differs")
    for key in (
        "nested_and_centered_on_identical_base_interval",
        "checkpoint_time_normalization_remains_2048",
        "effective_time_valid_counts_reported",
    ):
        _exact_bool(contexts, key, True)
    _exact_bool(contexts, "tile_or_cross_crop_averaging", False)

    event = _mapping(config.get("artificial_event"), label="artificial-event config")
    if event.get("profile") != "symmetric_trapezoid_sampled_at_cadence_centers":
        raise ValueError("centered-event profile differs")
    if not math.isclose(
        float(event.get("ingress_egress_fraction_each", math.nan)),
        EVENT_INGRESS_EGRESS_FRACTION,
        rel_tol=0.0,
        abs_tol=0.0,
    ):
        raise ValueError("centered-event ingress fraction differs")
    if tuple(event.get("duration_cadences", ())) != EVENT_DURATIONS_CADENCES:
        raise ValueError("centered-event duration grid differs")
    if tuple(float(value) for value in event.get("fractional_depths", ())) != (
        EVENT_FRACTIONAL_DEPTHS
    ):
        raise ValueError("centered-event depth grid differs")
    if event.get("event_center") != "designated_center_cadence_base_index_1024":
        raise ValueError("centered-event center convention differs")
    if event.get("apply_to_views") != ["adp_1x1", "adp_3x3"]:
        raise ValueError("centered-event injected views differ")
    if event.get("relative_flux_equation") != (
        "injected=(1+original)*(1-depth_profile)-1"
    ):
        raise ValueError("centered-event flux equation differs")
    for key in (
        "all_depth_duration_conditions_applied_to_each_source",
        "one_event_in_each_injected_sample",
        "time_errors_validity_and_masks_unchanged",
    ):
        _exact_bool(event, key, True)
    _exact_bool(event, "period_defined", False)

    metrics = _mapping(config.get("metrics"), label="centered-event metrics")
    if tuple(metrics.get("representations", ())) != _REPRESENTATIONS:
        raise ValueError("centered-event representation stages differ")
    if tuple(metrics.get("primary", ())) != (
        "original_population_robust_scaled_rms_displacement",
        "raw_l2_displacement",
        "one_minus_cosine_similarity",
        "coherent_shift_and_direction",
        "paired_context_ratio_to_2048",
        "paired_step2000_minus_step0",
    ):
        raise ValueError("centered-event primary metrics differ")
    if tuple(metrics.get("secondary_reconstruction", ())) != (
        "visible_event_projection_gain",
        "event_support_normalized_residual",
    ) or metrics.get("uncertainty") != (
        "deterministic_source_clustered_bootstrap_95_interval"
    ):
        raise ValueError("centered-event secondary metrics differ")

    return config, source, hashlib.sha256(payload_bytes).hexdigest()


def _load_inference_only_trusted_model(
    run_dir: Path,
    *,
    device: Any,
    checkpoint_path: Path,
) -> tuple[Any, dict[str, Any], dict[str, Any]]:
    """Load one checksum-bound model without constructing training state.

    The ordinary full release validator proves optimizer resumability by
    constructing an optimizer and scheduler.  That action is outside this
    inference-only diagnostic.  Here the release, selected checkpoint hash,
    run contract, global step, exact model-config schema, and strict finite
    model state are validated while optimizer and scheduler payloads remain
    untouched.
    """

    try:
        import torch
    except ImportError as exc:  # pragma: no cover - Torch is optional locally
        raise RuntimeError("PyTorch is required for centered-event inference") from exc
    from .model import TWIRLFM0, FM0ModelConfig, count_trainable_parameters

    root = Path(run_dir).resolve(strict=True)
    validation = validate_real_run_release(root, inspect_checkpoint=False)
    contract = read_json(root / "run_contract.json")
    summary = read_json(root / "summary.json")
    selected = Path(checkpoint_path).resolve(strict=True)
    if selected.parent != root:
        raise ValueError("centered-event checkpoint must be directly inside run_dir")
    milestone_match = re.fullmatch(r"checkpoint_step_([0-9]{8})\.pt", selected.name)
    if selected.name == "checkpoint.pt":
        digest = str(validation["artifact_sha256"]["checkpoint.pt"])
    elif milestone_match is not None:
        sidecar = selected.with_name(selected.name + ".sha256")
        if not sidecar.is_file():
            raise ValueError("centered-event milestone lacks its SHA-256 sidecar")
        sidecar_fields = sidecar.read_text(encoding="utf-8").strip().split()
        digest = sha256_file(selected)
        if (
            len(sidecar_fields) != 2
            or sidecar_fields[0] != digest
            or sidecar_fields[1] != selected.name
        ):
            raise ValueError("centered-event milestone SHA-256 sidecar differs")
    else:
        raise ValueError("centered-event checkpoint is not a trusted release artifact")

    try:
        checkpoint = torch.load(selected, map_location=device, weights_only=False)
    except Exception as exc:
        raise ValueError(
            "centered-event checkpoint is not structurally loadable"
        ) from exc
    required = {
        "schema_version",
        "model_state",
        "model_config",
        "progress",
        "run_contract",
    }
    if not isinstance(checkpoint, Mapping) or not required.issubset(checkpoint):
        raise ValueError("centered-event checkpoint model schema is incomplete")
    if checkpoint.get("schema_version") != _TRUSTED_CHECKPOINT_SCHEMA_VERSION:
        raise ValueError("centered-event checkpoint schema differs")
    if checkpoint.get("run_contract") != contract:
        raise ValueError("centered-event checkpoint run contract differs")

    progress = checkpoint.get("progress")
    if not isinstance(progress, Mapping):
        raise TypeError("centered-event checkpoint progress is malformed")
    global_step = progress.get("global_step")
    if (
        isinstance(global_step, bool)
        or not isinstance(global_step, int)
        or global_step < 0
    ):
        raise ValueError("centered-event checkpoint step is invalid")
    if milestone_match is not None:
        filename_step = int(milestone_match.group(1))
        allowed_steps = contract.get("immutable_milestone_steps")
        if (
            global_step != filename_step
            or not isinstance(allowed_steps, list)
            or global_step not in allowed_steps
        ):
            raise ValueError("centered-event milestone step was not predeclared")
    elif global_step != validation.get("global_step"):
        raise ValueError("centered-event release checkpoint step differs")

    model_payload = checkpoint.get("model_config")
    expected_fields = {field.name for field in fields(FM0ModelConfig)}
    if not isinstance(model_payload, Mapping) or set(model_payload) != expected_fields:
        raise ValueError("centered-event model configuration schema differs")
    try:
        model_config = FM0ModelConfig(**dict(model_payload))
    except (TypeError, ValueError) as exc:
        raise ValueError("centered-event model configuration is invalid") from exc
    variant = contract.get("variant")
    if (
        not isinstance(variant, str)
        or model_config.n_flux_views != len(variant_view_indices(variant))
        or model_config.architecture != contract.get("architecture")
        or model_config.architecture != validation.get("architecture")
        or variant != validation.get("variant")
    ):
        raise ValueError("centered-event checkpoint model identity differs")

    model_state = checkpoint.get("model_state")
    if (
        not isinstance(model_state, Mapping)
        or not model_state
        or not all(isinstance(value, torch.Tensor) for value in model_state.values())
    ):
        raise ValueError("centered-event checkpoint model state is malformed")
    try:
        model = TWIRLFM0(model_config)
        model.load_state_dict(model_state, strict=True)
    except (RuntimeError, TypeError, ValueError) as exc:
        raise ValueError("centered-event model state differs from its config") from exc
    for value in model.state_dict().values():
        if (value.is_floating_point() or value.is_complex()) and not bool(
            torch.isfinite(value).all()
        ):
            raise ValueError("centered-event model state contains non-finite values")
    parameter_count = summary.get("parameter_count")
    if (
        isinstance(parameter_count, bool)
        or not isinstance(parameter_count, int)
        or parameter_count != count_trainable_parameters(model)
    ):
        raise ValueError("centered-event model parameter count differs")

    model = model.to(device)
    model.eval()
    selected_validation = dict(validation)
    selected_validation.update(
        {
            "global_step": int(global_step),
            "selected_checkpoint_path": str(selected),
            "selected_checkpoint_sha256": digest,
            "selected_checkpoint_is_immutable_milestone": (milestone_match is not None),
            "checkpoint_model_state_inspected": True,
            "checkpoint_inspection_mode": "model_only_inference",
            "optimizer_or_scheduler_constructed": False,
            "model_state_tensor_count": len(model_state),
        }
    )
    return model, contract, selected_validation


def centered_context_bounds(
    context_length: int,
    *,
    base_length: int = BASE_INTERVAL_CADENCES,
) -> tuple[int, int]:
    """Return the frozen nested crop around the designated center cadence."""

    length = _exact_int(context_length, label="context length")
    base = _exact_int(base_length, label="base length")
    if base != BASE_INTERVAL_CADENCES or length not in CONTEXT_LENGTHS:
        raise ValueError("context bounds are outside the frozen grid")
    center = base // 2
    return center - length // 2, center + length // 2


def centered_event_bounds(length: int, duration_cadences: int) -> tuple[int, int]:
    """Return one odd-duration support centered on cadence ``length // 2``."""

    size = _exact_int(length, label="sample length")
    duration = _exact_int(duration_cadences, label="event duration")
    if size not in CONTEXT_LENGTHS and size != BASE_INTERVAL_CADENCES:
        raise ValueError("sample length is outside the frozen context grid")
    if duration not in EVENT_DURATIONS_CADENCES:
        raise ValueError("event duration is outside the frozen grid")
    center = size // 2
    start = center - duration // 2
    stop = start + duration
    if start < 0 or stop > size:
        raise ValueError("event support escaped the sample")
    return start, stop


def centered_trapezoid(
    duration_cadences: int,
    fractional_depth: float,
    *,
    ingress_egress_fraction: float = EVENT_INGRESS_EGRESS_FRACTION,
) -> np.ndarray:
    """Return the frozen cadence-center-sampled symmetric depth profile."""

    duration = _exact_int(duration_cadences, label="event duration")
    depth = float(fractional_depth)
    fraction = float(ingress_egress_fraction)
    if duration not in EVENT_DURATIONS_CADENCES:
        raise ValueError("event duration is outside the frozen grid")
    if depth not in EVENT_FRACTIONAL_DEPTHS:
        raise ValueError("event depth is outside the frozen grid")
    if not math.isclose(
        fraction, EVENT_INGRESS_EGRESS_FRACTION, rel_tol=0.0, abs_tol=0.0
    ):
        raise ValueError("event ingress fraction differs from the frozen contract")
    if duration == 1:
        return np.asarray([depth], dtype=np.float64)
    phase = (np.arange(duration, dtype=np.float64) + 0.5) / float(duration)
    amplitude = np.minimum.reduce(
        (
            phase / fraction,
            (1.0 - phase) / fraction,
            np.ones(duration, dtype=np.float64),
        )
    )
    return depth * amplitude


def _validate_model_sample(sample: Mapping[str, np.ndarray], *, length: int) -> None:
    if frozenset(sample) != _MODEL_SAMPLE_KEYS:
        raise ValueError("centered-event sample differs from the FM0 tensor allowlist")
    for key in _CADENCE_KEYS:
        array = np.asarray(sample[key])
        if array.shape[-1] != int(length):
            raise ValueError(f"centered-event sample field {key} has wrong length")
    flux = np.asarray(sample["flux"])
    if flux.ndim != 2 or np.asarray(sample["flux_valid"]).shape != flux.shape:
        raise ValueError("centered-event flux tensors differ")
    if np.asarray(sample["reconstruction_mask"]).shape != flux.shape:
        raise ValueError("centered-event reconstruction mask differs")
    if np.asarray(sample["time_valid"]).shape != (int(length),):
        raise ValueError("centered-event time-valid tensor differs")


def slice_centered_context(
    sample: Mapping[str, np.ndarray], *, context_length: int
) -> dict[str, np.ndarray]:
    """Slice one direct nested context without tiling or embedding averaging."""

    _validate_model_sample(sample, length=BASE_INTERVAL_CADENCES)
    start, stop = centered_context_bounds(context_length)
    cropped: dict[str, np.ndarray] = {}
    for key, value in sample.items():
        array = np.asarray(value)
        cropped[key] = (
            array[..., start:stop].copy() if key in _CADENCE_KEYS else array.copy()
        )
    time_valid = np.asarray(cropped["time_valid"], dtype=bool)
    local = np.asarray(cropped["local_time_cadences"])
    rebased = np.zeros_like(local)
    valid_indices = np.flatnonzero(time_valid)
    if valid_indices.size:
        rebased[time_valid] = local[time_valid] - local[valid_indices[0]]
    cropped["local_time_cadences"] = rebased
    _validate_model_sample(cropped, length=int(context_length))
    return cropped


def inject_centered_event(
    sample: Mapping[str, np.ndarray],
    *,
    duration_cadences: int,
    fractional_depth: float,
    ingress_egress_fraction: float = EVENT_INGRESS_EGRESS_FRACTION,
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """Inject exactly one visible event and return its one-dimensional support."""

    length = int(np.asarray(sample.get("time_valid", ())).size)
    _validate_model_sample(sample, length=length)
    start, stop = centered_event_bounds(length, duration_cadences)
    profile = centered_trapezoid(
        duration_cadences,
        fractional_depth,
        ingress_egress_fraction=ingress_egress_fraction,
    )
    flux_valid = np.asarray(sample["flux_valid"], dtype=bool)
    time_valid = np.asarray(sample["time_valid"], dtype=bool)
    if not np.all(flux_valid[:, start:stop] & time_valid[None, start:stop]):
        raise CenteredEventIneligibleError(
            "centered event support is not fully valid in every required view"
        )
    injected = {key: np.asarray(value).copy() for key, value in sample.items()}
    original = np.asarray(sample["flux"])
    values = (1.0 + original[:, start:stop].astype(np.float64)) * (
        1.0 - profile[None, :]
    ) - 1.0
    injected["flux"][:, start:stop] = values.astype(original.dtype, copy=False)
    support = np.zeros(length, dtype=bool)
    support[start:stop] = True
    for key in sample:
        if key != "flux" and not np.array_equal(injected[key], np.asarray(sample[key])):
            raise ValueError(f"event injection changed non-flux tensor {key}")
    return injected, support


def _window_valid_counts(
    release: Any,
    *,
    segment_id: int,
    start_offset: int,
    required_views: Sequence[int],
) -> tuple[bool, dict[int, tuple[int, ...]]]:
    indices = np.flatnonzero(release.segment_id == int(segment_id))
    selected = indices[int(start_offset) : int(start_offset) + BASE_INTERVAL_CADENCES]
    if selected.size != BASE_INTERVAL_CADENCES:
        return False, {}
    valid = (
        release.time_valid[selected, None]
        & release.flux_valid[
            selected[:, None], np.asarray(required_views, dtype=int)[None, :]
        ]
    )
    counts: dict[int, tuple[int, ...]] = {}
    for length in CONTEXT_LENGTHS:
        start, stop = centered_context_bounds(length)
        observed = tuple(
            int(value) for value in np.count_nonzero(valid[start:stop], axis=0)
        )
        counts[length] = observed
        if any(value / float(length) < MINIMUM_VALID_FRACTION for value in observed):
            return False, counts
    event_start, event_stop = centered_event_bounds(
        BASE_INTERVAL_CADENCES, max(EVENT_DURATIONS_CADENCES)
    )
    if not np.all(valid[event_start:event_stop]):
        return False, counts
    return True, counts


def _centered_window_spec(
    release: Any,
    *,
    observation_key: str,
    variant: str,
) -> tuple[WindowSpec, dict[int, tuple[int, ...]]]:
    """Choose a deterministic unpadded base supporting every nested context."""

    required = tuple(int(value) for value in variant_view_indices(variant))
    initial = deterministic_training_window(
        release,
        observation_key=observation_key,
        epoch=0,
        draw_index=0,
    )
    if initial.n_padded == 0:
        eligible, counts = _window_valid_counts(
            release,
            segment_id=initial.segment_id,
            start_offset=initial.start_offset,
            required_views=required,
        )
        if eligible:
            return initial, counts

    candidates: list[tuple[int, int]] = []
    for segment_id in np.unique(release.segment_id):
        indices = np.flatnonzero(release.segment_id == int(segment_id))
        if indices.size and np.any(np.diff(indices) != 1):
            raise ValueError("FM0 release segment is non-contiguous")
        n_starts = int(indices.size) - BASE_INTERVAL_CADENCES + 1
        if n_starts <= 0:
            continue
        valid = (
            release.time_valid[indices, None]
            & release.flux_valid[
                indices[:, None], np.asarray(required, dtype=int)[None, :]
            ]
        )
        cumulative = np.vstack(
            [np.zeros((1, len(required)), dtype=np.int64), np.cumsum(valid, axis=0)]
        )
        starts = np.arange(n_starts, dtype=np.int64)
        eligible = np.ones(n_starts, dtype=bool)
        for length in CONTEXT_LENGTHS:
            offset, _ = centered_context_bounds(length)
            counts = cumulative[starts + offset + length] - cumulative[starts + offset]
            eligible &= np.all(counts / float(length) >= MINIMUM_VALID_FRACTION, axis=1)
        event_start, event_stop = centered_event_bounds(
            BASE_INTERVAL_CADENCES, max(EVENT_DURATIONS_CADENCES)
        )
        event_counts = (
            cumulative[starts + event_stop] - cumulative[starts + event_start]
        )
        eligible &= np.all(event_counts == event_stop - event_start, axis=1)
        candidates.extend((int(segment_id), int(start)) for start in starts[eligible])
    if not candidates:
        raise CenteredEventIneligibleError(
            "development visit has no unpadded 2048-cadence interval satisfying "
            "all centered-context and event-support validity requirements: "
            f"{observation_key}"
        )
    selected_segment, selected_start = candidates[
        _seed(CENTERED_EVENT_SELECTION_SALT, "base", observation_key) % len(candidates)
    ]
    eligible, selected_counts = _window_valid_counts(
        release,
        segment_id=selected_segment,
        start_offset=selected_start,
        required_views=required,
    )
    if not eligible:
        raise ValueError("vectorized centered-window eligibility did not reproduce")
    return (
        WindowSpec(selected_segment, selected_start, BASE_INTERVAL_CADENCES, 0),
        selected_counts,
    )


def _load_centered_visit(row: Mapping[str, str], *, variant: str) -> dict[str, Any]:
    if row.get("source_partition") != DEVELOPMENT_PARTITION:
        raise ValueError("centered-event diagnostic may open only development rows")
    root = Path(str(row["sector_release_root"])).expanduser().resolve(strict=True)
    relative = Path(str(row["relative_path"]))
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("centered-event diagnostic received an unsafe shard path")
    shard = (root / relative).resolve(strict=True)
    if root not in shard.parents:
        raise ValueError("centered-event diagnostic shard escaped its release root")
    payload = shard.read_bytes()
    observed_hash = hashlib.sha256(payload).hexdigest()
    if observed_hash != str(row["sha256"]):
        raise ValueError("centered-event shard differs from its panel binding")
    release = load_input_release_bytes(payload)
    spec, counts = _centered_window_spec(
        release,
        observation_key=str(row["observation_key"]),
        variant=variant,
    )
    window = extract_window(
        release,
        segment_id=spec.segment_id,
        start_offset=spec.start_offset,
    )
    sample = prepare_model_window(
        window,
        variant=variant,
        mask_seed=_seed(
            CENTERED_EVENT_SELECTION_SALT, "unmasked", row["observation_key"]
        ),
        window_length=BASE_INTERVAL_CADENCES,
    )
    sample["temporal_mask"] = np.zeros(BASE_INTERVAL_CADENCES, dtype=bool)
    sample["reconstruction_mask"] = np.zeros_like(
        sample["reconstruction_mask"], dtype=bool
    )
    _validate_model_sample(sample, length=BASE_INTERVAL_CADENCES)
    time_counts = {
        str(length): int(
            np.count_nonzero(
                sample["time_valid"][slice(*centered_context_bounds(length))]
            )
        )
        for length in CONTEXT_LENGTHS
    }
    return {
        "sample": sample,
        "binding": {
            "observation_key": str(row["observation_key"]),
            "leakage_component_id": str(row["leakage_component_id"]),
            "sector": int(row["sector"]),
            "shard_sha256": observed_hash,
            "segment_id": int(spec.segment_id),
            "start_offset": int(spec.start_offset),
            "n_observed": int(spec.n_observed),
            "n_padded": int(spec.n_padded),
            "valid_counts_by_context_and_view": {
                str(length): list(counts[length]) for length in CONTEXT_LENGTHS
            },
            "effective_time_valid_count_by_context": time_counts,
        },
    }


def _component_order_key(cohort: str, component: str) -> tuple[str, str]:
    digest = hashlib.sha256(
        f"{CENTERED_EVENT_SELECTION_SALT}\x1fcomponent\x1f{cohort}\x1f{component}".encode()
    ).hexdigest()
    return digest, component


def _visit_order_key(row: Mapping[str, str]) -> tuple[str, str]:
    observation = str(row["observation_key"])
    digest = hashlib.sha256(
        f"{CENTERED_EVENT_SELECTION_SALT}\x1fvisit\x1f{observation}".encode()
    ).hexdigest()
    return digest, observation


def select_centered_event_visits(
    rows: Sequence[Mapping[str, str]],
    *,
    variant: str,
    max_components_per_cohort: int = COMPONENTS_PER_COHORT,
) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    """Select one structurally eligible development visit per component."""

    limit = _exact_int(max_components_per_cohort, label="component limit")
    if limit <= 1:
        raise ValueError("each centered-event cohort requires at least two components")
    validated = select_temporal_rows(
        rows,
        max_repeated_components=len(rows),
        max_new_components=len(rows),
        required_view_indices=(),
    )
    required_views = variant_view_indices(variant)
    selected: dict[str, dict[str, Any]] = {}
    audit: dict[str, Any] = {}
    opened_observations: set[str] = set()
    for cohort in TEMPORAL_COHORTS:
        grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in validated[cohort]:
            try:
                present = json.loads(str(row.get("view_present_json", "")))
            except json.JSONDecodeError as exc:
                raise ValueError(
                    "centered-event row has invalid view_present_json"
                ) from exc
            if (
                not isinstance(present, list)
                or len(present) != len(VIEW_NAMES)
                or any(
                    type(value) is not int or value not in {0, 1} for value in present
                )
            ):
                raise ValueError(
                    "centered-event row has invalid view_present_json schema"
                )
            grouped[str(row["leakage_component_id"])].append(dict(row))

        chosen_rows: list[dict[str, str]] = []
        chosen_visits: list[dict[str, Any]] = []
        skipped_components: list[str] = []
        skipped_visits = 0
        missing_view_visits = 0
        opened_for_cohort = 0
        for component in sorted(
            grouped, key=lambda value: _component_order_key(cohort, value)
        ):
            if len(chosen_rows) >= limit:
                break
            selected_pair: tuple[dict[str, str], dict[str, Any]] | None = None
            for row in sorted(grouped[component], key=_visit_order_key):
                present = json.loads(str(row["view_present_json"]))
                if any(not bool(present[index]) for index in required_views):
                    missing_view_visits += 1
                    continue
                observation = str(row["observation_key"])
                if observation in opened_observations:
                    raise ValueError("centered-event selection reopened a shard")
                opened_observations.add(observation)
                opened_for_cohort += 1
                try:
                    visit = _load_centered_visit(row, variant=variant)
                except CenteredEventIneligibleError:
                    skipped_visits += 1
                    continue
                selected_pair = row, visit
                break
            if selected_pair is None:
                skipped_components.append(component)
                continue
            chosen_rows.append(selected_pair[0])
            chosen_visits.append(selected_pair[1])
        if len(chosen_rows) != limit:
            raise ValueError(
                f"centered-event cohort {cohort} yielded {len(chosen_rows)} of "
                f"{limit} required components"
            )
        selected[cohort] = {"rows": chosen_rows, "visits": chosen_visits}
        audit[cohort] = {
            "selected_components": len(chosen_rows),
            "selected_visits": len(chosen_visits),
            "one_visit_per_component": True,
            "screening_shards_opened": opened_for_cohort,
            "structurally_ineligible_visits_skipped": skipped_visits,
            "missing_required_view_visits_skipped": missing_view_visits,
            "components_without_eligible_visit_skipped": len(skipped_components),
            "skipped_components_sha256": _identity_sha256(
                tuple(sorted(skipped_components))
            ),
        }
    return selected, audit


def _development_shard_access_summary(
    selection_audit: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    if set(selection_audit) != set(TEMPORAL_COHORTS):
        raise ValueError("centered-event shard audit cohorts differ")
    screening = 0
    selected = 0
    for cohort in TEMPORAL_COHORTS:
        audit = selection_audit[cohort]
        opened = audit.get("screening_shards_opened")
        retained = audit.get("selected_visits")
        if (
            isinstance(opened, bool)
            or not isinstance(opened, int)
            or isinstance(retained, bool)
            or not isinstance(retained, int)
            or opened < retained
            or retained < 0
        ):
            raise ValueError("centered-event shard audit counts are invalid")
        screening += opened
        selected += retained
    return {
        "screening_development_shards_opened": screening,
        "selected_development_shards": selected,
        "screening_shards_opened_but_not_selected": screening - selected,
        "each_screened_development_shard_opened_at_most_once": True,
    }


def _encode_samples(
    *,
    model: Any,
    samples: Sequence[Mapping[str, np.ndarray]],
    context_length: int,
    device: Any,
    batch_size: int,
) -> dict[str, np.ndarray]:
    """Encode direct contexts while retaining full visible-input decoder output."""

    try:
        import torch
    except ImportError as exc:  # pragma: no cover - Torch is optional locally
        raise RuntimeError(
            "PyTorch is required for the centered-event diagnostic"
        ) from exc
    if not samples or batch_size <= 0:
        raise ValueError("centered-event samples and batch size must be positive")
    model.eval()
    collected: dict[str, list[np.ndarray]] = {
        "h_window": [],
        "z_window": [],
        "reconstruction": [],
    }
    with torch.no_grad():
        for start in range(0, len(samples), batch_size):
            batch_samples = samples[start : start + batch_size]
            batch = move_batch_to_device(collate_fm0_samples(batch_samples), device)
            output = (
                model(batch)
                if int(context_length) == BASE_INTERVAL_CADENCES
                else model.forward_short_context(batch)
            )
            for name, destination in collected.items():
                if name not in output:
                    raise ValueError(f"centered-event encoder output lacks {name}")
                values = output[name].detach().float().cpu().numpy()
                expected_shape = (
                    (
                        len(batch_samples),
                        np.asarray(batch_samples[0]["flux"]).shape[0],
                        int(context_length),
                    )
                    if name == "reconstruction"
                    else (len(batch_samples), values.shape[-1])
                )
                if values.shape != expected_shape or not np.all(np.isfinite(values)):
                    raise ValueError(f"centered-event encoder output {name} is invalid")
                destination.append(values)
    return {name: np.concatenate(values, axis=0) for name, values in collected.items()}


def _fit_original_robust_scale(
    embeddings: np.ndarray,
) -> tuple[np.ndarray, dict[str, Any]]:
    array = np.asarray(embeddings, dtype=np.float64)
    if array.ndim != 2 or array.shape[0] < 2 or not np.all(np.isfinite(array)):
        raise ValueError("original embedding population is invalid")
    center = np.median(array, axis=0)
    scale = 1.4826 * np.median(np.abs(array - center), axis=0)
    usable = np.isfinite(scale) & (scale > 1.0e-12)
    if not np.any(usable):
        raise ValueError("original embedding population has no positive robust scale")
    result = np.full(scale.shape, np.nan, dtype=np.float64)
    result[usable] = scale[usable]
    return result, {
        "definition": "1.4826_times_dimensionwise_original_population_mad",
        "n_original_embeddings": int(array.shape[0]),
        "embedding_dimensions": int(array.shape[1]),
        "positive_scale_dimensions": int(np.count_nonzero(usable)),
        "zero_or_nonfinite_scale_dimensions_excluded": int(
            usable.size - np.count_nonzero(usable)
        ),
        "median_positive_scale": float(np.median(scale[usable])),
    }


def _scalar_summary(
    values: np.ndarray,
    *,
    component_ids: Sequence[str],
    condition_ids: Sequence[str],
    seed_label: str,
) -> dict[str, Any]:
    array = np.asarray(values, dtype=np.float64)
    components = np.asarray(tuple(str(value) for value in component_ids), dtype=object)
    conditions = np.asarray(tuple(str(value) for value in condition_ids), dtype=object)
    if (
        array.ndim != 1
        or components.shape != array.shape
        or conditions.shape != array.shape
    ):
        raise ValueError("centered-event scalar metric identities do not align")

    def one(mask: np.ndarray, label: str) -> dict[str, Any]:
        finite = mask & np.isfinite(array)
        if not np.any(finite):
            return {
                "status": "not_available",
                "reason": "no_finite_pairs",
                "n_pairs": int(np.count_nonzero(mask)),
                "n_finite_pairs": 0,
            }
        selected = array[finite]
        selected_components = tuple(str(value) for value in components[finite])
        return {
            "status": "available",
            "n_pairs": int(np.count_nonzero(mask)),
            "n_finite_pairs": int(selected.size),
            "mean": float(np.mean(selected)),
            "median": float(np.median(selected)),
            "source_clustered_95_interval": source_clustered_mean_interval(
                selected,
                selected_components,
                seed=_seed(CENTERED_EVENT_SELECTION_SALT, seed_label, label),
            ),
        }

    return {
        "overall": one(np.ones(array.shape, dtype=bool), "overall"),
        "by_condition": {
            condition: one(conditions == condition, condition)
            for condition in dict.fromkeys(str(value) for value in conditions)
        },
    }


def _coherent_shift_summary(delta: np.ndarray, scale: np.ndarray) -> dict[str, Any]:
    values = np.asarray(delta, dtype=np.float64)
    robust_scale = np.asarray(scale, dtype=np.float64)
    if (
        values.ndim != 2
        or robust_scale.shape != (values.shape[1],)
        or not np.all(np.isfinite(values))
    ):
        raise ValueError("coherent-shift inputs differ")
    mean_delta = np.mean(values, axis=0)
    usable = np.isfinite(robust_scale) & (robust_scale > 0.0)
    norms = np.linalg.norm(values, axis=1)
    nonzero = norms > 0.0
    direction = (
        float(np.linalg.norm(np.mean(values[nonzero] / norms[nonzero, None], axis=0)))
        if np.any(nonzero)
        else 0.0
    )
    return {
        "n_pairs": int(values.shape[0]),
        "zero_displacement_pairs": int(np.count_nonzero(~nonzero)),
        "raw_coherent_centroid_shift_l2": float(np.linalg.norm(mean_delta)),
        "original_population_robust_scaled_coherent_shift_rms": float(
            np.sqrt(np.mean(np.square(mean_delta[usable] / robust_scale[usable])))
        ),
        "unit_delta_direction_coherence": direction,
    }


def _embedding_response(
    original: np.ndarray,
    injected: np.ndarray,
    *,
    robust_scale: np.ndarray,
    component_ids: Sequence[str],
    condition_ids: Sequence[str],
    seed_label: str,
    pooling_normalizer: np.ndarray | None = None,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    left = np.asarray(original, dtype=np.float64)
    right = np.asarray(injected, dtype=np.float64)
    scale = np.asarray(robust_scale, dtype=np.float64)
    if left.ndim != 2 or right.shape != left.shape or scale.shape != (left.shape[1],):
        raise ValueError("centered-event paired embeddings differ")
    if not np.all(np.isfinite(left)) or not np.all(np.isfinite(right)):
        raise ValueError("centered-event paired embeddings are non-finite")
    delta = right - left
    usable = np.isfinite(scale) & (scale > 0.0)
    raw_l2 = np.linalg.norm(delta, axis=1)
    scaled_rms = np.sqrt(np.mean(np.square(delta[:, usable] / scale[usable]), axis=1))
    cosine_distance = 1.0 - cosine_similarity_rows(left, right)
    metrics = {
        "original_population_robust_scaled_rms_displacement": scaled_rms,
        "raw_l2_displacement": raw_l2,
        "one_minus_cosine_similarity": cosine_distance,
    }
    if pooling_normalizer is not None:
        normalizer = np.asarray(pooling_normalizer, dtype=np.float64)
        if (
            normalizer.shape != raw_l2.shape
            or not np.all(np.isfinite(normalizer))
            or np.any(normalizer <= 0.0)
        ):
            raise ValueError("centered-event pooling normalizer is invalid")
        metrics[
            "effective_time_valid_pooling_corrected_l2_per_integrated_event_signal"
        ] = raw_l2 * normalizer
    conditions = np.asarray(tuple(str(value) for value in condition_ids), dtype=object)
    summary = {
        name: _scalar_summary(
            values,
            component_ids=component_ids,
            condition_ids=condition_ids,
            seed_label=f"{seed_label}:{name}",
        )
        for name, values in metrics.items()
    }
    summary["coherent_shift_and_direction"] = {
        "overall": _coherent_shift_summary(delta, scale),
        "by_condition": {
            condition: _coherent_shift_summary(delta[conditions == condition], scale)
            for condition in dict.fromkeys(str(value) for value in conditions)
        },
    }
    return summary, metrics


def summarize_embedding_response(
    original: np.ndarray,
    injected: np.ndarray,
    *,
    robust_scale: np.ndarray,
    component_ids: Sequence[str],
    condition_ids: Sequence[str],
    seed_label: str = "embedding-response",
    pooling_normalizer: np.ndarray | None = None,
) -> dict[str, Any]:
    """Summarize paired event information in one embedding stage."""

    return _embedding_response(
        original,
        injected,
        robust_scale=robust_scale,
        component_ids=component_ids,
        condition_ids=condition_ids,
        seed_label=seed_label,
        pooling_normalizer=pooling_normalizer,
    )[0]


def visible_event_decoder_response(
    original_flux: np.ndarray,
    injected_flux: np.ndarray,
    original_reconstruction: np.ndarray,
    injected_reconstruction: np.ndarray,
    *,
    flux_valid: np.ndarray,
    event_support: np.ndarray,
) -> dict[str, np.ndarray]:
    """Measure visible-input decoder sensitivity on the event support.

    This is not reconstruction fidelity: FM0 trains the decoder on hidden
    targets.  The pairwise difference only asks whether the visible event
    changes decoder output in the same direction and amplitude.
    """

    original = np.asarray(original_flux, dtype=np.float64)
    injected = np.asarray(injected_flux, dtype=np.float64)
    original_prediction = np.asarray(original_reconstruction, dtype=np.float64)
    injected_prediction = np.asarray(injected_reconstruction, dtype=np.float64)
    valid = np.asarray(flux_valid, dtype=bool)
    support = np.asarray(event_support, dtype=bool)
    if (
        original.ndim != 3
        or injected.shape != original.shape
        or original_prediction.shape != original.shape
        or injected_prediction.shape != original.shape
        or valid.shape != original.shape
        or support.shape != (original.shape[0], original.shape[2])
    ):
        raise ValueError("visible-event decoder tensors differ")
    truth = injected - original
    predicted = injected_prediction - original_prediction
    gains = np.empty((original.shape[0], original.shape[1]), dtype=np.float64)
    residuals = np.empty_like(gains)
    for pair in range(original.shape[0]):
        for view in range(original.shape[1]):
            mask = support[pair] & valid[pair, view]
            truth_values = truth[pair, view, mask]
            predicted_values = predicted[pair, view, mask]
            denominator = float(np.sum(np.square(truth_values)))
            if denominator <= 0.0 or not np.all(np.isfinite(predicted_values)):
                raise ValueError(
                    "visible-event decoder pair has no finite truth support"
                )
            gains[pair, view] = float(
                np.sum(predicted_values * truth_values) / denominator
            )
            residuals[pair, view] = float(
                np.linalg.norm(predicted_values - truth_values) / math.sqrt(denominator)
            )
    return {
        "visible_event_projection_gain": np.mean(gains, axis=1),
        "event_support_normalized_residual": np.mean(residuals, axis=1),
    }


def _summarize_decoder_values(
    values: Mapping[str, np.ndarray],
    *,
    component_ids: Sequence[str],
    condition_ids: Sequence[str],
    seed_label: str,
) -> dict[str, Any]:
    return {
        name: _scalar_summary(
            array,
            component_ids=component_ids,
            condition_ids=condition_ids,
            seed_label=f"{seed_label}:{name}",
        )
        for name, array in values.items()
    }


def _condition_grid() -> tuple[dict[str, Any], ...]:
    return tuple(
        {
            "condition_id": f"duration_{duration:02d}_depth_{depth:.2f}",
            "duration_cadences": duration,
            "fractional_depth": depth,
        }
        for duration in EVENT_DURATIONS_CADENCES
        for depth in EVENT_FRACTIONAL_DEPTHS
    )


def _paired_metric_comparisons(
    current: Mapping[str, np.ndarray],
    reference: Mapping[str, np.ndarray],
    *,
    component_ids: Sequence[str],
    condition_ids: Sequence[str],
    seed_label: str,
    include_ratio: bool,
) -> dict[str, Any]:
    if set(current) != set(reference):
        raise ValueError("paired metric stages differ")
    result: dict[str, Any] = {}
    for name in current:
        left = np.asarray(current[name], dtype=np.float64)
        right = np.asarray(reference[name], dtype=np.float64)
        if left.shape != right.shape:
            raise ValueError("paired metric arrays differ")
        payload: dict[str, Any] = {
            "paired_difference": _scalar_summary(
                left - right,
                component_ids=component_ids,
                condition_ids=condition_ids,
                seed_label=f"{seed_label}:{name}:difference",
            )
        }
        if include_ratio:
            ratio = np.divide(
                left,
                right,
                out=np.full(left.shape, np.nan, dtype=np.float64),
                where=np.isfinite(right) & (np.abs(right) > 1.0e-12),
            )
            payload["paired_ratio"] = _scalar_summary(
                ratio,
                component_ids=component_ids,
                condition_ids=condition_ids,
                seed_label=f"{seed_label}:{name}:ratio",
            )
        result[name] = payload
    return result


def evaluate_centered_event_context(
    *,
    config_path: str | Path,
    run_dir: str | Path,
    step0_checkpoint_path: str | Path,
    step2000_checkpoint_path: str | Path,
    temporal_panel_dir: str | Path,
    temporal_panel_receipt_sha256: str,
    batch_size: int = 32,
) -> dict[str, Any]:
    """Evaluate the complete frozen centered-event context experiment."""

    try:
        import torch
    except ImportError as exc:  # pragma: no cover - Torch is optional locally
        raise RuntimeError(
            "PyTorch is required for the centered-event diagnostic"
        ) from exc
    if batch_size <= 0:
        raise ValueError("centered-event batch size must be positive")
    config, config_source, config_hash = _load_frozen_config(config_path)
    inputs = _mapping(config["inputs"], label="centered-event inputs")
    root = Path(run_dir).resolve(strict=True)
    target_device = torch.device("cpu")
    configured_panel_hash = _exact_sha256(
        inputs["temporal_panel_receipt_sha256"],
        label="configured temporal-panel receipt SHA-256",
    )
    supplied_panel_hash = _exact_sha256(
        temporal_panel_receipt_sha256,
        label="supplied temporal-panel receipt SHA-256",
    )
    if supplied_panel_hash != configured_panel_hash:
        raise ValueError("supplied temporal-panel receipt differs from frozen config")
    panel_rows, panel_summary = load_temporal_panel(
        temporal_panel_dir,
        receipt_sha256=supplied_panel_hash,
    )
    model0, contract0, validation0 = _load_inference_only_trusted_model(
        root,
        device=target_device,
        checkpoint_path=Path(step0_checkpoint_path),
    )
    model2000, contract2000, validation2000 = _load_inference_only_trusted_model(
        root,
        device=target_device,
        checkpoint_path=Path(step2000_checkpoint_path),
    )
    if validation0.get("global_step") != EXPECTED_CHECKPOINT_STEPS[0]:
        raise ValueError("centered-event initial checkpoint is not exact step 0")
    if validation2000.get("global_step") != EXPECTED_CHECKPOINT_STEPS[1]:
        raise ValueError("centered-event trained checkpoint is not exact step 2000")
    checkpoint_hashes = _mapping(
        inputs["checkpoint_sha256"], label="centered-event checkpoint hashes"
    )
    if validation0.get("selected_checkpoint_sha256") != checkpoint_hashes["step0"]:
        raise ValueError("centered-event step-0 checkpoint hash differs")
    if (
        validation2000.get("selected_checkpoint_sha256")
        != checkpoint_hashes["step2000"]
    ):
        raise ValueError("centered-event step-2000 checkpoint hash differs")
    if contract0 != contract2000:
        raise ValueError("centered-event checkpoints have different run contracts")
    variant = str(contract0.get("variant", ""))
    if variant != "TWIRL-FM0.2.1":
        raise ValueError("centered-event diagnostic requires TWIRL-FM0.2.1")
    for model in (model0, model2000):
        if int(model.config.window_length) != BASE_INTERVAL_CADENCES or not callable(
            getattr(model, "forward_short_context", None)
        ):
            raise ValueError("centered-event checkpoint lacks short-context support")

    cohort_data, selection_audit = select_centered_event_visits(
        panel_rows,
        variant=variant,
        max_components_per_cohort=COMPONENTS_PER_COHORT,
    )
    for cohort in TEMPORAL_COHORTS:
        data = cohort_data[cohort]
        data["component_ids"] = [
            str(row["leakage_component_id"]) for row in data["rows"]
        ]
        data["base_sample_sha256"] = _sample_set_sha256(
            data["rows"], [visit["sample"] for visit in data["visits"]]
        )

    conditions = _condition_grid()
    models = {0: model0, 2_000: model2000}
    checkpoint_results: dict[int, dict[str, Any]] = {
        step: {"contexts": {}} for step in EXPECTED_CHECKPOINT_STEPS
    }
    pair_values: dict[int, dict[int, dict[str, Any]]] = {
        step: {} for step in EXPECTED_CHECKPOINT_STEPS
    }
    context_bindings: dict[str, Any] = {}

    for context_length in CONTEXT_LENGTHS:
        prepared: dict[str, dict[str, Any]] = {}
        for cohort in TEMPORAL_COHORTS:
            data = cohort_data[cohort]
            originals = [
                slice_centered_context(visit["sample"], context_length=context_length)
                for visit in data["visits"]
            ]
            injected_samples: list[dict[str, np.ndarray]] = []
            source_indices: list[int] = []
            pair_components: list[str] = []
            condition_ids: list[str] = []
            event_support: list[np.ndarray] = []
            effective_time_valid_counts: list[int] = []
            integrated_event_signals: list[float] = []
            pseudo_rows: list[dict[str, str]] = []
            for source_index, (row, original) in enumerate(
                zip(data["rows"], originals, strict=True)
            ):
                for condition in conditions:
                    injected, support = inject_centered_event(
                        original,
                        duration_cadences=int(condition["duration_cadences"]),
                        fractional_depth=float(condition["fractional_depth"]),
                    )
                    injected_samples.append(injected)
                    source_indices.append(source_index)
                    pair_components.append(str(row["leakage_component_id"]))
                    condition_ids.append(str(condition["condition_id"]))
                    event_support.append(support)
                    effective_time_valid_counts.append(
                        int(np.count_nonzero(original["time_valid"]))
                    )
                    truth_delta = injected["flux"] - original["flux"]
                    integrated_signal = float(
                        np.mean(np.sum(np.abs(truth_delta[:, support]), axis=1))
                    )
                    if not math.isfinite(integrated_signal) or integrated_signal <= 0.0:
                        raise ValueError(
                            "centered event has no finite integrated signal"
                        )
                    integrated_event_signals.append(integrated_signal)
                    pseudo_rows.append(
                        {
                            "observation_key": (
                                f"{row['observation_key']}:{condition['condition_id']}"
                            )
                        }
                    )
            if len(injected_samples) != COMPONENTS_PER_COHORT * len(conditions):
                raise ValueError("centered-event condition grid did not close")
            prepared[cohort] = {
                "originals": originals,
                "injected": injected_samples,
                "source_indices": np.asarray(source_indices, dtype=np.int64),
                "pair_components": tuple(pair_components),
                "condition_ids": tuple(condition_ids),
                "event_support": np.stack(event_support),
                "pooling_normalizer": np.asarray(
                    effective_time_valid_counts, dtype=np.float64
                )
                / np.asarray(integrated_event_signals, dtype=np.float64),
            }
            exact_time_counts = [
                int(np.count_nonzero(sample["time_valid"])) for sample in originals
            ]
            exact_flux_counts = [
                [
                    int(np.count_nonzero(sample["flux_valid"][view]))
                    for view in range(len(variant_view_indices(variant)))
                ]
                for sample in originals
            ]
            context_bindings.setdefault(str(context_length), {})[cohort] = {
                "bounds_in_base_interval": list(
                    centered_context_bounds(context_length)
                ),
                "n_original_samples": len(originals),
                "n_injected_samples": len(injected_samples),
                "original_sample_set_sha256": _sample_set_sha256(
                    data["rows"], originals
                ),
                "injected_sample_set_sha256": _sample_set_sha256(
                    pseudo_rows, injected_samples
                ),
                "minimum_time_valid_cadences": int(min(exact_time_counts)),
                "effective_time_valid_cadences_by_sample": exact_time_counts,
                "minimum_flux_valid_cadences_by_view": [
                    int(
                        min(
                            np.count_nonzero(sample["flux_valid"][view])
                            for sample in originals
                        )
                    )
                    for view in range(len(variant_view_indices(variant)))
                ],
                "effective_flux_valid_cadences_by_sample_and_view": exact_flux_counts,
            }

        encoded: dict[int, dict[str, dict[str, np.ndarray]]] = {}
        for step, model in models.items():
            encoded[step] = {}
            for cohort in TEMPORAL_COHORTS:
                data = prepared[cohort]
                original_output = _encode_samples(
                    model=model,
                    samples=data["originals"],
                    context_length=context_length,
                    device=target_device,
                    batch_size=batch_size,
                )
                injected_output = _encode_samples(
                    model=model,
                    samples=data["injected"],
                    context_length=context_length,
                    device=target_device,
                    batch_size=batch_size,
                )
                encoded[step][cohort] = {
                    "original": original_output,
                    "injected": injected_output,
                }

            scales: dict[str, np.ndarray] = {}
            scale_audit: dict[str, Any] = {}
            for representation in _REPRESENTATIONS:
                scales[representation], scale_audit[representation] = (
                    _fit_original_robust_scale(
                        np.concatenate(
                            [
                                encoded[step][cohort]["original"][representation]
                                for cohort in TEMPORAL_COHORTS
                            ],
                            axis=0,
                        )
                    )
                )
            context_result: dict[str, Any] = {
                "context_length": context_length,
                "original_population_robust_scale": scale_audit,
                "cohorts": {},
            }
            pair_values[step][context_length] = {}
            for cohort in TEMPORAL_COHORTS:
                data = prepared[cohort]
                source_indices = data["source_indices"]
                original_output = encoded[step][cohort]["original"]
                injected_output = encoded[step][cohort]["injected"]
                cohort_result: dict[str, Any] = {"representations": {}}
                stored: dict[str, Any] = {"representations": {}}
                for representation in _REPRESENTATIONS:
                    summary, values = _embedding_response(
                        original_output[representation][source_indices],
                        injected_output[representation],
                        robust_scale=scales[representation],
                        component_ids=data["pair_components"],
                        condition_ids=data["condition_ids"],
                        seed_label=(
                            f"step{step}:context{context_length}:{cohort}:{representation}"
                        ),
                        pooling_normalizer=data["pooling_normalizer"],
                    )
                    cohort_result["representations"][representation] = summary
                    stored["representations"][representation] = values

                original_flux = np.stack(
                    [data["originals"][index]["flux"] for index in source_indices]
                )
                injected_flux = np.stack(
                    [sample["flux"] for sample in data["injected"]]
                )
                flux_valid = np.stack(
                    [data["originals"][index]["flux_valid"] for index in source_indices]
                )
                decoder_values = visible_event_decoder_response(
                    original_flux,
                    injected_flux,
                    original_output["reconstruction"][source_indices],
                    injected_output["reconstruction"],
                    flux_valid=flux_valid,
                    event_support=data["event_support"],
                )
                cohort_result["visible_input_decoder_response"] = (
                    _summarize_decoder_values(
                        decoder_values,
                        component_ids=data["pair_components"],
                        condition_ids=data["condition_ids"],
                        seed_label=f"step{step}:context{context_length}:{cohort}:decoder",
                    )
                )
                stored["decoder"] = decoder_values
                stored["component_ids"] = data["pair_components"]
                stored["condition_ids"] = data["condition_ids"]
                context_result["cohorts"][cohort] = cohort_result
                pair_values[step][context_length][cohort] = stored
            checkpoint_results[step]["contexts"][str(context_length)] = context_result

    context_comparisons: dict[str, Any] = {}
    for step in EXPECTED_CHECKPOINT_STEPS:
        context_comparisons[f"step{step}"] = {}
        for context_length in CONTEXT_LENGTHS:
            context_comparisons[f"step{step}"][str(context_length)] = {}
            for cohort in TEMPORAL_COHORTS:
                current = pair_values[step][context_length][cohort]
                reference = pair_values[step][BASE_INTERVAL_CADENCES][cohort]
                payload: dict[str, Any] = {"representations": {}}
                for representation in _REPRESENTATIONS:
                    payload["representations"][representation] = (
                        _paired_metric_comparisons(
                            current["representations"][representation],
                            reference["representations"][representation],
                            component_ids=current["component_ids"],
                            condition_ids=current["condition_ids"],
                            seed_label=(
                                f"step{step}:context{context_length}:vs2048:"
                                f"{cohort}:{representation}"
                            ),
                            include_ratio=True,
                        )
                    )
                payload["visible_input_decoder_response"] = _paired_metric_comparisons(
                    current["decoder"],
                    reference["decoder"],
                    component_ids=current["component_ids"],
                    condition_ids=current["condition_ids"],
                    seed_label=(
                        f"step{step}:context{context_length}:vs2048:{cohort}:decoder"
                    ),
                    include_ratio=False,
                )
                context_comparisons[f"step{step}"][str(context_length)][cohort] = (
                    payload
                )

    checkpoint_deltas: dict[str, Any] = {}
    for context_length in CONTEXT_LENGTHS:
        checkpoint_deltas[str(context_length)] = {}
        for cohort in TEMPORAL_COHORTS:
            initial = pair_values[0][context_length][cohort]
            trained = pair_values[2_000][context_length][cohort]
            payload = {"representations": {}}
            for representation in _REPRESENTATIONS:
                payload["representations"][representation] = _paired_metric_comparisons(
                    trained["representations"][representation],
                    initial["representations"][representation],
                    component_ids=initial["component_ids"],
                    condition_ids=initial["condition_ids"],
                    seed_label=(
                        f"step2000-minus-step0:context{context_length}:"
                        f"{cohort}:{representation}"
                    ),
                    include_ratio=False,
                )
            payload["visible_input_decoder_response"] = _paired_metric_comparisons(
                trained["decoder"],
                initial["decoder"],
                component_ids=initial["component_ids"],
                condition_ids=initial["condition_ids"],
                seed_label=(
                    f"step2000-minus-step0:context{context_length}:{cohort}:decoder"
                ),
                include_ratio=False,
            )
            checkpoint_deltas[str(context_length)][cohort] = payload

    run_contract_path = root / "run_contract.json"
    selection_population = {
        cohort: {
            "selected_components": COMPONENTS_PER_COHORT,
            "selected_visits": COMPONENTS_PER_COHORT,
            "selected_component_ids_sha256": _identity_sha256(
                tuple(sorted(cohort_data[cohort]["component_ids"]))
            ),
            "selected_observation_keys_sha256": _identity_sha256(
                tuple(row["observation_key"] for row in cohort_data[cohort]["rows"])
            ),
            "base_sample_set_sha256": cohort_data[cohort]["base_sample_sha256"],
            "visits": [visit["binding"] for visit in cohort_data[cohort]["visits"]],
        }
        for cohort in TEMPORAL_COHORTS
    }
    shard_access = _development_shard_access_summary(selection_audit)
    return {
        "schema_version": CENTERED_EVENT_CONTEXT_RESULT_SCHEMA_VERSION,
        "status": "descriptive_training_free_development_diagnostic_complete",
        "claim_limit": (
            "Centered, artificial, visible single-event sensitivity of a frozen "
            "FM0.2.1 backbone; not blind detection, real classification, training, "
            "sealed testing, reconstruction fidelity, or model promotion"
        ),
        "config": {
            "path": str(config_source),
            "sha256": config_hash,
            "schema_version": config["schema_version"],
            "campaign_id": config.get("campaign_id"),
        },
        "run": {
            "run_dir": str(root),
            "variant": variant,
            "architecture": contract0.get("architecture"),
            "run_git_sha": contract0.get("expected_git_sha"),
            "run_contract_path": str(run_contract_path),
            "run_contract_sha256": hashlib.sha256(
                run_contract_path.read_bytes()
            ).hexdigest(),
            "device": "cpu",
            "exact_checkpoint_steps": list(EXPECTED_CHECKPOINT_STEPS),
        },
        "temporal_panel": panel_summary,
        "selection": {
            "one_visit_per_component": True,
            "minimum_valid_fraction_each_required_view": MINIMUM_VALID_FRACTION,
            "full_validity_over_maximum_event_support": True,
            "only_structural_ineligibility_is_skippable": True,
            "cohorts": selection_audit,
            "population": selection_population,
        },
        "experiment": {
            "base_interval_cadences": BASE_INTERVAL_CADENCES,
            "designated_center_cadence_base_index": BASE_INTERVAL_CADENCES // 2,
            "context_lengths": list(CONTEXT_LENGTHS),
            "event_conditions": list(conditions),
            "conditions_per_source": len(conditions),
            "one_event_per_injected_sample": True,
            "all_conditions_applied_to_each_source": True,
            "period_defined": False,
            "bls_or_period_features_used": False,
            "cross_crop_embedding_averaging": False,
            "pooling_law_metric_definition": (
                "raw_l2_displacement_times_effective_time_valid_count_divided_by_"
                "equal_view_mean_integrated_absolute_event_flux_change"
            ),
            "context_bindings": context_bindings,
        },
        "checkpoints": {
            "step0": {
                "global_step": 0,
                "checkpoint_path": validation0["selected_checkpoint_path"],
                "checkpoint_sha256": validation0["selected_checkpoint_sha256"],
                "projection_spectrum": _model_projection_spectrum(model0),
                **checkpoint_results[0],
            },
            "step2000": {
                "global_step": 2_000,
                "checkpoint_path": validation2000["selected_checkpoint_path"],
                "checkpoint_sha256": validation2000["selected_checkpoint_sha256"],
                "projection_spectrum": _model_projection_spectrum(model2000),
                **checkpoint_results[2_000],
            },
        },
        "paired_context_comparisons_to_2048": context_comparisons,
        "paired_step2000_minus_step0": checkpoint_deltas,
        "metric_interpretation": {
            "h_and_z_are_not_independent": True,
            "z_is_a_linear_projection_of_h": True,
            "h_stage": "pre_projection_masked_mean_pooled_hidden_state",
            "z_stage": "post_linear_projection_window_embedding",
            "visible_input_decoder_response_is_not_reconstruction_fidelity": True,
            "short_contexts_keep_checkpoint_2048_time_normalization": True,
            "fm0_2_1_tcn_receptive_field_cadences": 1_543,
            "contexts_below_trained_receptive_field": [128, 256, 512],
            "short_contexts_cross_the_trained_context_and_receptive_field_regime": True,
            "short_context_result_does_not_predict_a_retrained_short_model": True,
            "centered_result_is_an_aligned_event_upper_bound": True,
        },
        "data_access": {
            **shard_access,
            "sealed_shards_opened": 0,
            "labels_candidates_or_events_opened": False,
            "artificial_events_injected": True,
            "bls_or_period_features_opened": False,
            "classifier_or_probe_fitted": False,
            "optimizer_or_training_created": False,
        },
    }


__all__ = [
    "BASE_INTERVAL_CADENCES",
    "CENTERED_EVENT_CONTEXT_CONFIG_SCHEMA_VERSION",
    "CENTERED_EVENT_CONTEXT_RESULT_SCHEMA_VERSION",
    "COMPONENTS_PER_COHORT",
    "CONTEXT_LENGTHS",
    "EVENT_DURATIONS_CADENCES",
    "EVENT_FRACTIONAL_DEPTHS",
    "CenteredEventIneligibleError",
    "centered_context_bounds",
    "centered_event_bounds",
    "centered_trapezoid",
    "evaluate_centered_event_context",
    "inject_centered_event",
    "select_centered_event_visits",
    "slice_centered_context",
    "summarize_embedding_response",
    "visible_event_decoder_response",
]
