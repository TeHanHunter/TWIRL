"""Artifact-backed FM0.3 cadence-collapse preflight.

This module is deliberately inference-only.  It loads one hash-bound panel and
the exact step-0 and step-2,000 checkpoints, runs both models on the same
native-cadence tensors, and passes only their cadence-indexed hidden states to
the frozen collapse math.  No cadence representation is binned, merged,
resampled, or temporally averaged.  Until a frozen panel receipt and validated
real-run/input-release contract exist, this module emits a non-authoritative
preflight only.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from dataclasses import fields
from pathlib import Path
from typing import Any

import numpy as np

from .cadence_collapse_gate import _evaluate_cadence_collapse_math

PANEL_SCHEMA_VERSION = "twirl_fm0_3_cadence_collapse_panel_v1"
RESULT_SCHEMA_VERSION = "twirl_fm0_3_cadence_collapse_artifact_preflight_v1"
CHECKPOINT_SCHEMA_VERSION = "twirl_fm0_1_checkpoint_v1"
OBJECTIVE_STATE_SCHEMA_VERSION = "twirl_fm0_objective_state_v2"
CADENCE_OBJECTIVE_IDENTITY = (
    "masked_reconstruction_plus_position_centered_cadence_vicreg_v1"
)
FM0_3_VARIANTS = frozenset({"TWIRL-FM0.3.1", "TWIRL-FM0.3.2"})
WINDOW_LENGTH = 128
PATCH_STRIDE = 1
NATIVE_CADENCE_SECONDS = 200
MASK_TARGET_FRACTION = 0.15
MASK_SPAN_RANGE = (1, 4)

_FLOAT_FIELDS = frozenset(
    {
        "flux",
        "flux_error",
        "local_time_cadences",
        "delta_time_cadences",
    }
)
_BOOL_FIELDS = frozenset(
    {
        "flux_valid",
        "error_valid",
        "time_valid",
        "segment_boundary",
        "view_present",
        "reconstruction_mask",
    }
)
_MODEL_BATCH_FIELDS = tuple(sorted(_FLOAT_FIELDS | _BOOL_FIELDS))
_PANEL_FIELDS = frozenset(
    {
        "schema_version",
        "variant",
        "native_cadence_seconds",
        "cadence_axis_operation",
        "sample_identity",
        *_MODEL_BATCH_FIELDS,
    }
)
_OBJECTIVE_STATE = {
    "schema_version": OBJECTIVE_STATE_SCHEMA_VERSION,
    "identity": CADENCE_OBJECTIVE_IDENTITY,
    "use_vicreg": True,
    "reconstruct_second_view": False,
}
_SHA256_HEX = frozenset("0123456789abcdef")


def sha256_file(path: str | Path) -> str:
    """Return the SHA-256 computed from the actual file bytes."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _exact_sha256(value: Any, *, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or value != value.lower()
        or any(character not in _SHA256_HEX for character in value)
    ):
        raise ValueError(f"{label} must be a lowercase SHA-256 digest")
    return value


def _resolve_hash_bound_file(
    path: str | Path,
    *,
    expected_sha256: str,
    label: str,
) -> tuple[Path, str]:
    source = Path(path).expanduser().resolve(strict=True)
    if not source.is_file():
        raise ValueError(f"{label} is not a regular file")
    expected = _exact_sha256(expected_sha256, label=f"expected {label} SHA-256")
    observed = sha256_file(source)
    if observed != expected:
        raise ValueError(f"{label} SHA-256 differs from the actual artifact bytes")
    return source, observed


def _unicode_scalar(value: np.ndarray, *, label: str) -> str:
    if value.shape != () or value.dtype.kind != "U":
        raise ValueError(f"panel {label} must be one Unicode scalar")
    result = str(value.item())
    if not result:
        raise ValueError(f"panel {label} cannot be empty")
    return result


def _load_panel(path: Path) -> dict[str, Any]:
    """Load the exact pickle-free native-cadence NPZ panel schema."""

    try:
        with np.load(path, allow_pickle=False) as archive:
            if frozenset(archive.files) != _PANEL_FIELDS:
                raise ValueError("cadence-collapse panel fields differ from the schema")
            arrays = {name: np.asarray(archive[name]).copy() for name in archive.files}
    except (OSError, TypeError, ValueError) as exc:
        if isinstance(exc, ValueError) and "panel fields" in str(exc):
            raise
        raise ValueError("cadence-collapse panel is not a pickle-free NPZ") from exc

    if _unicode_scalar(arrays["schema_version"], label="schema_version") != (
        PANEL_SCHEMA_VERSION
    ):
        raise ValueError("cadence-collapse panel schema version differs")
    variant = _unicode_scalar(arrays["variant"], label="variant")
    if variant not in FM0_3_VARIANTS:
        raise ValueError("cadence-collapse panel variant is not FM0.3")
    operation = _unicode_scalar(
        arrays["cadence_axis_operation"], label="cadence_axis_operation"
    )
    if operation != "none":
        raise ValueError("cadence-collapse panel must not transform the cadence axis")
    cadence_seconds = arrays["native_cadence_seconds"]
    if cadence_seconds.shape != () or cadence_seconds.dtype != np.dtype("int64"):
        raise ValueError("panel native_cadence_seconds must be one int64 scalar")
    if int(cadence_seconds.item()) != NATIVE_CADENCE_SECONDS:
        raise ValueError("cadence-collapse panel is not native 200-second cadence data")

    identities = arrays["sample_identity"]
    if identities.ndim != 1 or identities.dtype.kind != "U" or identities.size < 2:
        raise ValueError(
            "panel sample_identity must be a Unicode vector of length >= 2"
        )
    sample_identity = tuple(str(value) for value in identities.tolist())
    if any(not value for value in sample_identity):
        raise ValueError("panel sample identities must be nonempty")
    if len(set(sample_identity)) != len(sample_identity):
        raise ValueError("panel sample identities must be unique")

    flux = arrays["flux"]
    if flux.dtype != np.dtype("float32") or flux.ndim != 3:
        raise ValueError("panel flux must have float32 shape [B,V,L]")
    batch_size, n_views, cadence_length = flux.shape
    if batch_size != len(sample_identity) or cadence_length != WINDOW_LENGTH:
        raise ValueError("panel flux must align to B identities and 128 cadences")
    expected_shapes = {
        "flux": (batch_size, n_views, cadence_length),
        "flux_valid": (batch_size, n_views, cadence_length),
        "flux_error": (batch_size, 2, cadence_length),
        "error_valid": (batch_size, 2, cadence_length),
        "local_time_cadences": (batch_size, cadence_length),
        "delta_time_cadences": (batch_size, cadence_length),
        "time_valid": (batch_size, cadence_length),
        "segment_boundary": (batch_size, cadence_length),
        "view_present": (batch_size, n_views),
        "reconstruction_mask": (batch_size, n_views, cadence_length),
    }
    for name in _FLOAT_FIELDS:
        value = arrays[name]
        if value.dtype != np.dtype("float32") or value.shape != expected_shapes[name]:
            raise ValueError(f"panel {name} has the wrong float32 shape")
        if not np.all(np.isfinite(value)):
            raise ValueError(f"panel {name} contains non-finite values")
    for name in _BOOL_FIELDS:
        value = arrays[name]
        if value.dtype != np.dtype("bool") or value.shape != expected_shapes[name]:
            raise ValueError(f"panel {name} has the wrong boolean shape")
    if np.any(arrays["reconstruction_mask"]):
        raise ValueError("cadence-collapse inference panel must be unmasked")
    if np.any(arrays["flux_valid"] & ~arrays["view_present"][:, :, None]):
        raise ValueError("panel marks flux valid in an absent view")
    if np.any(arrays["error_valid"] & ~arrays["time_valid"][:, None, :]):
        raise ValueError("panel marks an error valid at an invalid cadence")

    batch = {name: arrays[name] for name in _MODEL_BATCH_FIELDS}
    return {
        "variant": variant,
        "sample_identity": sample_identity,
        "batch": batch,
        "batch_size": batch_size,
        "n_views": n_views,
        "cadence_length": cadence_length,
    }


def _dataset_geometry(dataset_contract: Mapping[str, Any]) -> Mapping[str, Any]:
    kind = dataset_contract.get("kind")
    if kind == "synthetic":
        config = dataset_contract.get("config")
        if not isinstance(config, Mapping):
            raise ValueError("FM0.3 synthetic dataset contract is malformed")
        return config
    if kind in {"fm0_input_release", "fm0_3_composite_release"}:
        return dataset_contract
    raise ValueError("FM0.3 checkpoint dataset contract kind is unsupported")


def _validate_checkpoint_contract(
    payload: Mapping[str, Any],
    *,
    expected_step: int,
) -> tuple[Any, dict[str, Any]]:
    """Validate one checkpoint and build only its inference model."""

    try:
        import torch
    except ImportError as exc:  # pragma: no cover - Torch is platform optional
        raise RuntimeError(
            "PyTorch is required for cadence-collapse inference"
        ) from exc
    from .dataset import variant_view_indices
    from .model import FM0ModelConfig, architecture_for_variant, build_fm0_model

    required = {
        "schema_version",
        "model_state",
        "model_config",
        "progress",
        "run_contract",
        "dataset_contract",
        "objective_state",
    }
    if not isinstance(payload, Mapping) or not required.issubset(payload):
        raise ValueError("FM0.3 cadence-collapse checkpoint schema is incomplete")
    if payload.get("schema_version") != CHECKPOINT_SCHEMA_VERSION:
        raise ValueError("FM0.3 cadence-collapse checkpoint schema version differs")

    progress = payload.get("progress")
    if not isinstance(progress, Mapping):
        raise TypeError("FM0.3 checkpoint progress is malformed")
    global_step = progress.get("global_step")
    if isinstance(global_step, bool) or not isinstance(global_step, int):
        raise TypeError("FM0.3 checkpoint global step is malformed")
    if global_step != expected_step:
        raise ValueError(
            f"FM0.3 cadence-collapse checkpoint must be exact step {expected_step}"
        )

    run_contract = payload.get("run_contract")
    dataset_contract = payload.get("dataset_contract")
    objective_state = payload.get("objective_state")
    if not isinstance(run_contract, Mapping) or not isinstance(
        dataset_contract, Mapping
    ):
        raise TypeError("FM0.3 checkpoint run or dataset contract is malformed")
    if not isinstance(objective_state, Mapping):
        raise TypeError("FM0.3 checkpoint cadence objective contract is malformed")
    if dict(objective_state) != _OBJECTIVE_STATE:
        raise ValueError("FM0.3 checkpoint cadence objective contract differs")

    variant = run_contract.get("variant")
    if variant not in FM0_3_VARIANTS:
        raise ValueError("cadence-collapse checkpoint variant is not FM0.3")
    architecture = architecture_for_variant(str(variant))
    if (
        run_contract.get("architecture") != architecture
        or run_contract.get("objective") != CADENCE_OBJECTIVE_IDENTITY
        or run_contract.get("use_vicreg") is not True
        or run_contract.get("reconstruct_second_view") is not False
        or run_contract.get("mask_target_fraction") != MASK_TARGET_FRACTION
        or run_contract.get("mask_span_range") != list(MASK_SPAN_RANGE)
    ):
        raise ValueError("FM0.3 checkpoint run contract differs from the frozen design")

    model_payload = payload.get("model_config")
    expected_model_fields = {field.name for field in fields(FM0ModelConfig)}
    if not isinstance(model_payload, Mapping) or set(model_payload) != (
        expected_model_fields
    ):
        raise ValueError("FM0.3 checkpoint model configuration schema differs")
    try:
        model_config = FM0ModelConfig(**dict(model_payload))
    except (TypeError, ValueError) as exc:
        raise ValueError("FM0.3 checkpoint model configuration is invalid") from exc
    if (
        model_config.window_length != WINDOW_LENGTH
        or model_config.patch_stride != PATCH_STRIDE
        or model_config.architecture != architecture
        or model_config.n_flux_views != len(variant_view_indices(str(variant)))
    ):
        raise ValueError(
            "FM0.3 checkpoint must preserve one token for each of 128 cadences"
        )

    dataset_geometry = _dataset_geometry(dataset_contract)
    if (
        dataset_geometry.get("variant") != variant
        or dataset_geometry.get("window_length") != WINDOW_LENGTH
        or dataset_geometry.get("mask_target_fraction") != MASK_TARGET_FRACTION
        or tuple(dataset_geometry.get("mask_span_range", ())) != MASK_SPAN_RANGE
    ):
        raise ValueError("FM0.3 checkpoint dataset contract geometry differs")

    model_state = payload.get("model_state")
    if (
        not isinstance(model_state, Mapping)
        or not model_state
        or not all(isinstance(value, torch.Tensor) for value in model_state.values())
    ):
        raise ValueError("FM0.3 checkpoint model state is malformed")
    try:
        model = build_fm0_model(str(variant), enforce_parameter_budget=True)
    except (RuntimeError, TypeError, ValueError) as exc:
        raise ValueError("canonical FM0.3 model construction failed") from exc
    if model_config != model.config:
        raise ValueError(
            "FM0.3 checkpoint model configuration differs from the canonical build"
        )
    try:
        model.load_state_dict(model_state, strict=True)
    except (RuntimeError, TypeError, ValueError) as exc:
        raise ValueError(
            "FM0.3 checkpoint model state differs from its config"
        ) from exc
    if any(
        (value.is_floating_point() or value.is_complex())
        and not bool(torch.isfinite(value).all())
        for value in model.state_dict().values()
    ):
        raise ValueError("FM0.3 checkpoint model state contains non-finite values")
    model.eval()
    contracts = {
        "run_contract": dict(run_contract),
        "model_config": dict(model_payload),
        "dataset_contract": dict(dataset_contract),
        "objective_state": dict(objective_state),
        "global_step": int(global_step),
        "variant": str(variant),
        "architecture": architecture,
    }
    return model, contracts


def _load_checkpoint(path: Path, *, expected_step: int) -> tuple[Any, dict[str, Any]]:
    try:
        import torch
    except ImportError as exc:  # pragma: no cover - Torch is platform optional
        raise RuntimeError(
            "PyTorch is required for cadence-collapse inference"
        ) from exc
    try:
        payload = torch.load(path, map_location="cpu", weights_only=False)
    except Exception as exc:
        raise ValueError("FM0.3 cadence-collapse checkpoint is not loadable") from exc
    return _validate_checkpoint_contract(payload, expected_step=expected_step)


def _canonical_contract_sha256(value: Mapping[str, Any]) -> str:
    try:
        encoded = json.dumps(
            dict(value),
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    except (TypeError, ValueError) as exc:
        raise ValueError("FM0.3 checkpoint contract is not canonical JSON") from exc
    return hashlib.sha256(encoded).hexdigest()


def _sample_order_sha256(identities: tuple[str, ...]) -> str:
    encoded = json.dumps(
        identities,
        ensure_ascii=False,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _infer_h_cadence(model: Any, batch: Mapping[str, np.ndarray]) -> tuple[Any, Any]:
    try:
        import torch
    except ImportError as exc:  # pragma: no cover - Torch is platform optional
        raise RuntimeError(
            "PyTorch is required for cadence-collapse inference"
        ) from exc
    tensor_batch = {
        name: torch.from_numpy(np.ascontiguousarray(value))
        for name, value in batch.items()
    }
    with torch.inference_mode():
        output = model(tensor_batch)
    hidden = output.get("h_cadence")
    token_valid = output.get("token_valid")
    if (
        not isinstance(hidden, torch.Tensor)
        or not isinstance(token_valid, torch.Tensor)
        or hidden.ndim != 3
        or token_valid.shape != hidden.shape[:2]
    ):
        raise ValueError("FM0.3 checkpoint did not emit cadence-indexed tokens")
    if hidden.shape[1] != WINDOW_LENGTH:
        raise ValueError("FM0.3 checkpoint changed the native cadence-token lattice")
    if not bool(torch.isfinite(hidden).all()):
        raise ValueError("FM0.3 checkpoint emitted non-finite cadence tokens")
    return hidden.detach().cpu(), token_valid.detach().cpu().bool()


def evaluate_cadence_collapse_artifacts(
    *,
    panel_path: str | Path,
    panel_sha256: str,
    step0_checkpoint_path: str | Path,
    step0_checkpoint_sha256: str,
    step2000_checkpoint_path: str | Path,
    step2000_checkpoint_sha256: str,
) -> dict[str, Any]:
    """Run non-authoritative collapse preflight from hash-verified artifacts."""

    panel_file, panel_digest = _resolve_hash_bound_file(
        panel_path,
        expected_sha256=panel_sha256,
        label="cadence-collapse panel",
    )
    step0_file, step0_digest = _resolve_hash_bound_file(
        step0_checkpoint_path,
        expected_sha256=step0_checkpoint_sha256,
        label="step-0 checkpoint",
    )
    step2000_file, step2000_digest = _resolve_hash_bound_file(
        step2000_checkpoint_path,
        expected_sha256=step2000_checkpoint_sha256,
        label="step-2000 checkpoint",
    )
    if step0_file == step2000_file or step0_digest == step2000_digest:
        raise ValueError("step-0 and step-2000 checkpoints must be distinct artifacts")

    panel = _load_panel(panel_file)
    step0_model, step0_contracts = _load_checkpoint(step0_file, expected_step=0)
    step2000_model, step2000_contracts = _load_checkpoint(
        step2000_file, expected_step=2_000
    )
    for name in (
        "run_contract",
        "model_config",
        "dataset_contract",
        "objective_state",
    ):
        if step0_contracts[name] != step2000_contracts[name]:
            raise ValueError(f"step-0 and step-2000 {name.replace('_', ' ')} differ")
    if panel["variant"] != step0_contracts["variant"]:
        raise ValueError("cadence-collapse panel and checkpoint variants differ")
    if panel["n_views"] != step0_model.config.n_flux_views:
        raise ValueError("cadence-collapse panel and checkpoint view counts differ")

    step0_hidden, step0_valid = _infer_h_cadence(step0_model, panel["batch"])
    step2000_hidden, step2000_valid = _infer_h_cadence(step2000_model, panel["batch"])
    if not np.array_equal(step0_valid.numpy(), step2000_valid.numpy()):
        raise ValueError("step-0 and step-2000 token_valid differ")
    panel_time_valid = panel["batch"]["time_valid"]
    if not np.array_equal(step0_valid.numpy(), panel_time_valid):
        raise ValueError("checkpoint token_valid differs from the exact panel")
    if step0_hidden.shape != step2000_hidden.shape:
        raise ValueError("step-0 and step-2000 cadence-token shapes differ")

    identities = panel["sample_identity"]
    natural_available = panel_time_valid & np.any(
        panel["batch"]["flux_valid"] & panel["batch"]["view_present"][:, :, None],
        axis=1,
    )
    gate = _evaluate_cadence_collapse_math(
        step0_hidden.numpy(),
        step2000_hidden.numpy(),
        natural_available,
        panel["batch"]["segment_boundary"],
        panel["batch"]["delta_time_cadences"],
    )
    contract_hashes = {
        name: _canonical_contract_sha256(step0_contracts[name])
        for name in (
            "run_contract",
            "model_config",
            "dataset_contract",
            "objective_state",
        )
    }
    return {
        "schema_version": RESULT_SCHEMA_VERSION,
        "status": "non_authoritative_preflight",
        "authoritative_artifact_gate": False,
        "artifact_bindings": {
            "panel_path": str(panel_file),
            "panel_sha256": panel_digest,
            "step0_checkpoint_path": str(step0_file),
            "step0_checkpoint_sha256": step0_digest,
            "step2000_checkpoint_path": str(step2000_file),
            "step2000_checkpoint_sha256": step2000_digest,
            "sample_order_sha256": _sample_order_sha256(identities),
            "sample_count": len(identities),
            "contract_sha256": contract_hashes,
        },
        "checkpoint_contract": {
            "variant": step0_contracts["variant"],
            "architecture": step0_contracts["architecture"],
            "objective": CADENCE_OBJECTIVE_IDENTITY,
            "exact_steps": [0, 2_000],
            "window_length": WINDOW_LENGTH,
            "patch_stride": PATCH_STRIDE,
        },
        "panel": {
            "schema_version": PANEL_SCHEMA_VERSION,
            "sample_count": panel["batch_size"],
            "cadence_length": panel["cadence_length"],
            "native_cadence_seconds": NATIVE_CADENCE_SECONDS,
            "cadence_axis_operation": "none",
        },
        "inference": {
            "representation": "h_cadence",
            "step0_token_shape": list(step0_hidden.shape),
            "step2000_token_shape": list(step2000_hidden.shape),
            "identical_token_valid": True,
            "gate_eligibility": ("time_valid_and_any_flux_valid_in_present_model_view"),
            "naturally_available_token_count": int(np.count_nonzero(natural_available)),
            "one_native_cadence_per_token": True,
            "temporal_averaging_or_pooling": False,
            "optimizer_or_scheduler_constructed": False,
        },
        "gate": gate,
        "criteria_satisfied": bool(gate["criteria_satisfied"]),
        "formal_gate_prerequisites": [
            "separately_frozen_panel_receipt",
            "validated_real_run_and_input_release_contract",
            "immutable_step0_and_step2000_milestone_checkpoints",
        ],
    }


__all__ = [
    "CADENCE_OBJECTIVE_IDENTITY",
    "CHECKPOINT_SCHEMA_VERSION",
    "PANEL_SCHEMA_VERSION",
    "RESULT_SCHEMA_VERSION",
    "evaluate_cadence_collapse_artifacts",
    "sha256_file",
]
