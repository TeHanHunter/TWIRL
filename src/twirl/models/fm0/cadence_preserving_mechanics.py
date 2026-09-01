"""Training-free mechanics for cadence-preserving short-context FM0 models."""

from __future__ import annotations

import hashlib
import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from .centered_event_context_diagnostic import inject_centered_event
from .dataset import VIEW_NAMES, collate_fm0_samples, prepare_model_window
from .model import TWIRLFM0, FM0ModelConfig, count_trainable_parameters, require_torch

CONFIG_SCHEMA_VERSION = "twirl_fm0_3_stride1_mechanics_config_v1"
RESULT_SCHEMA_VERSION = "twirl_fm0_3_stride1_mechanics_result_v1"
CAMPAIGN_ID = "twirl_fm0_3_stride1_mechanics_v1"
CONTEXT_LENGTH = 128
INITIALIZATION_SEEDS = (560067, 560068)
EVENT_DURATIONS = (1, 3, 9)
EVENT_DEPTH = 0.03
SYNTHETIC_SOURCES = 8
CANDIDATES = {
    "TWIRL-FM0.3.1": ("tcn", "TWIRL-FM0.1.1", 8_825_602),
    "TWIRL-FM0.3.2": ("conformer", "TWIRL-FM0.1.2", 9_345_282),
}


def _mapping(value: Any, *, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be a mapping")
    return value


def load_mechanics_config(path: str | Path) -> tuple[dict[str, Any], Path, str]:
    """Load and strictly validate the small, training-free mechanics contract."""

    try:
        import yaml
    except ImportError as exc:  # pragma: no cover - project dependency
        raise RuntimeError("PyYAML is required for FM0 mechanics") from exc

    source = Path(path).expanduser().resolve(strict=True)
    raw = source.read_bytes()
    config = dict(_mapping(yaml.safe_load(raw), label="mechanics config"))
    if (
        config.get("schema_version") != CONFIG_SCHEMA_VERSION
        or config.get("campaign_id") != CAMPAIGN_ID
        or config.get("model_family") != "TWIRL-FM0"
        or config.get("status") != "frozen_training_free_mechanics"
    ):
        raise ValueError("FM0 mechanics config identity differs")

    purpose = _mapping(config.get("purpose"), label="mechanics purpose")
    if (
        purpose.get("primary")
        != "verify_cadence_preserving_short_context_model_mechanics"
        or purpose.get("learned_utility_claim") is not False
        or purpose.get("architecture_selection_claim") is not False
    ):
        raise ValueError("FM0 mechanics purpose differs")

    authorization = _mapping(
        config.get("authorization"), label="mechanics authorization"
    )
    expected_authorization = {
        "deterministic_synthetic_inputs": True,
        "random_initialization_only": True,
        "checkpoint_access": False,
        "development_shard_access": False,
        "probe_fitting": False,
        "optimizer_or_model_training": False,
        "sealed_test_access": False,
        "gpu_use": False,
        "formal_model_gate": False,
    }
    if dict(authorization) != expected_authorization:
        raise ValueError("FM0 mechanics authorization differs")

    cadence = _mapping(config.get("cadence_contract"), label="cadence contract")
    expected_cadence = {
        "nominal_cadence_seconds": 200,
        "context_length_cadences": CONTEXT_LENGTH,
        "context_duration_hours": CONTEXT_LENGTH * 200.0 / 3600.0,
        "one_encoder_token_per_cadence": True,
        "patch_stride": 1,
        "cadence_averaging": False,
        "temporal_downsampling": False,
        "cadence_sequence_output": "h_cadence",
        "window_mean_outputs_diagnostic_only": ["h_window", "z_window"],
    }
    if dict(cadence) != expected_cadence:
        raise ValueError("FM0 cadence-preservation contract differs")

    candidates = _mapping(config.get("candidates"), label="mechanics candidates")
    if set(candidates) != set(CANDIDATES):
        raise ValueError("FM0 mechanics candidate set differs")
    for name, (architecture, template, parameters) in CANDIDATES.items():
        if dict(_mapping(candidates[name], label=name)) != {
            "architecture": architecture,
            "template_variant": template,
            "expected_parameter_count": parameters,
        }:
            raise ValueError(f"FM0 mechanics candidate {name} differs")

    evaluation = _mapping(config.get("evaluation"), label="mechanics evaluation")
    if dict(evaluation) != {
        "initialization_seeds": list(INITIALIZATION_SEEDS),
        "synthetic_sources": SYNTHETIC_SOURCES,
        "event_duration_cadences": list(EVENT_DURATIONS),
        "event_fractional_depth": EVENT_DEPTH,
        "identical_inputs_across_architectures": True,
        "invalid_fill_invariance": True,
        "require_finite_outputs": True,
    }:
        raise ValueError("FM0 mechanics evaluation contract differs")
    return config, source, hashlib.sha256(raw).hexdigest()


def _stable_seed(*parts: object) -> int:
    payload = "\x1f".join(str(part) for part in parts).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


def deterministic_base_samples(
    *, n_sources: int = SYNTHETIC_SOURCES
) -> list[dict[str, np.ndarray]]:
    """Build identical unmasked 128-cadence inputs for both architectures."""

    if n_sources <= 1:
        raise ValueError("mechanics screen requires at least two synthetic sources")
    samples: list[dict[str, np.ndarray]] = []
    cadence = np.arange(CONTEXT_LENGTH, dtype=np.float32)
    for source_index in range(n_sources):
        rng = np.random.default_rng(_stable_seed(CAMPAIGN_ID, source_index))
        common = (
            0.0015 * np.sin(2.0 * np.pi * cadence / (41.0 + source_index))
            + rng.normal(0.0, 0.001, CONTEXT_LENGTH)
        )
        flux = np.empty((CONTEXT_LENGTH, len(VIEW_NAMES)), dtype=np.float32)
        for view in range(len(VIEW_NAMES)):
            flux[:, view] = common + rng.normal(
                0.0, 0.0005 * (1.0 + 0.05 * view), CONTEXT_LENGTH
            )
        flux_valid = np.ones_like(flux, dtype=bool)
        # Exercise invalid-fill masking away from the centered event support.
        flux_valid[source_index % 8, :] = False
        flux_error = np.full((CONTEXT_LENGTH, 2), 0.001, dtype=np.float32)
        error_valid = np.ones_like(flux_error, dtype=bool)
        release_window = {
            "flux": flux,
            "flux_valid": flux_valid,
            "flux_error": flux_error,
            "error_valid": error_valid,
            "local_time_cadences": cadence.copy(),
            "delta_time_cadences": np.r_[0.0, np.ones(CONTEXT_LENGTH - 1)].astype(
                np.float32
            ),
            "time_valid": np.ones(CONTEXT_LENGTH, dtype=bool),
            "segment_boundary": np.r_[True, np.zeros(CONTEXT_LENGTH - 1, dtype=bool)],
            "view_present": np.ones(len(VIEW_NAMES), dtype=bool),
        }
        sample = prepare_model_window(
            release_window,
            variant="TWIRL-FM0.1.1",
            mask_seed=_stable_seed(CAMPAIGN_ID, "mask", source_index),
            window_length=CONTEXT_LENGTH,
        )
        sample["temporal_mask"] = np.zeros(CONTEXT_LENGTH, dtype=bool)
        sample["reconstruction_mask"] = np.zeros_like(
            sample["reconstruction_mask"], dtype=bool
        )
        samples.append(sample)
    return samples


def _model_config(architecture: str) -> FM0ModelConfig:
    return FM0ModelConfig(
        architecture=architecture,  # type: ignore[arg-type]
        n_flux_views=2,
        window_length=CONTEXT_LENGTH,
        patch_stride=1,
    )


def _finite_tensor(value: Any) -> bool:
    import torch

    return bool(torch.isfinite(value).all())


def _response_summary(
    original: np.ndarray,
    injected: np.ndarray,
    supports: Sequence[np.ndarray],
) -> dict[str, float]:
    delta = np.asarray(injected, dtype=np.float64) - np.asarray(
        original, dtype=np.float64
    )
    if delta.ndim != 3 or delta.shape[:2] != (len(supports), CONTEXT_LENGTH):
        raise ValueError("cadence response tensors differ from the mechanics contract")
    cadence_l2 = np.linalg.norm(delta, axis=2)
    event_values: list[float] = []
    off_values: list[float] = []
    for values, support in zip(cadence_l2, supports, strict=True):
        event_values.extend(values[np.asarray(support, dtype=bool)].tolist())
        off_values.extend(values[~np.asarray(support, dtype=bool)].tolist())
    event_rms = float(np.sqrt(np.mean(np.square(event_values))))
    off_rms = float(np.sqrt(np.mean(np.square(off_values))))
    return {
        "event_support_rms_displacement": event_rms,
        "off_support_rms_displacement": off_rms,
        "event_to_off_support_ratio": (
            event_rms / off_rms if off_rms > 0.0 else math.inf
        ),
        "maximum_cadence_displacement": float(np.max(cadence_l2)),
    }


def evaluate_cadence_preserving_mechanics(
    *, config_path: str | Path
) -> dict[str, Any]:
    """Run the exact no-optimizer TCN/stride-1-Conformer mechanics screen."""

    require_torch()
    import torch

    config, config_source, config_sha256 = load_mechanics_config(config_path)
    originals = deterministic_base_samples()
    injected: list[dict[str, np.ndarray]] = []
    supports: list[np.ndarray] = []
    source_indices: list[int] = []
    conditions: list[dict[str, int | float]] = []
    for source_index, sample in enumerate(originals):
        for duration in EVENT_DURATIONS:
            changed, support = inject_centered_event(
                sample,
                duration_cadences=duration,
                fractional_depth=EVENT_DEPTH,
            )
            injected.append(changed)
            supports.append(support)
            source_indices.append(source_index)
            conditions.append(
                {"duration_cadences": duration, "fractional_depth": EVENT_DEPTH}
            )
    original_batch = collate_fm0_samples(originals)
    injected_batch = collate_fm0_samples(injected)

    results: dict[str, Any] = {}
    for candidate, (architecture, template, expected_parameters) in CANDIDATES.items():
        seed_results: dict[str, Any] = {}
        for seed in INITIALIZATION_SEEDS:
            torch.manual_seed(seed)
            model = TWIRLFM0(_model_config(architecture)).eval()
            parameters = count_trainable_parameters(model)
            if parameters != expected_parameters:
                raise ValueError(f"{candidate} parameter count differs")
            with torch.no_grad():
                original_output = model(original_batch)
                injected_output = model(injected_batch)
                repeated_original = original_output["h_cadence"][source_indices]

                invalid_batch = {key: value.clone() for key, value in original_batch.items()}
                invalid_batch["flux"][~invalid_batch["flux_valid"]] = 1.0e6
                invalid_output = model(invalid_batch)

                packed, valid, _local, _delta = model._pack_input(original_batch)
                stem = model.stem(packed, valid)
                patched, patched_valid = model._patch_mean(stem, valid)
                patched_time = model._patch_time(
                    original_batch["local_time_cadences"].float(), valid
                )

            required_outputs = ("reconstruction", "h_cadence", "h_window", "z_window")
            if not all(_finite_tensor(original_output[name]) for name in required_outputs):
                raise ValueError(f"{candidate} emitted non-finite tensors")
            if original_output["token_valid"].shape != (
                SYNTHETIC_SOURCES,
                CONTEXT_LENGTH,
            ):
                raise ValueError(f"{candidate} reduced the cadence-token axis")
            if original_output["h_cadence"].shape != (
                SYNTHETIC_SOURCES,
                CONTEXT_LENGTH,
                256,
            ):
                raise ValueError(f"{candidate} cadence sequence shape differs")
            if not torch.equal(patched, stem) or not torch.equal(patched_valid, valid):
                raise ValueError(f"{candidate} patch operation is not identity")
            if not torch.equal(
                patched_time, original_batch["local_time_cadences"].float()
            ):
                raise ValueError(f"{candidate} time patch operation is not identity")
            for name in required_outputs:
                if not torch.equal(original_output[name], invalid_output[name]):
                    raise ValueError(f"{candidate} is sensitive to invalid fill in {name}")

            h_response = _response_summary(
                repeated_original.detach().cpu().numpy(),
                injected_output["h_cadence"].detach().cpu().numpy(),
                supports,
            )
            if not math.isfinite(h_response["event_support_rms_displacement"]) or (
                h_response["event_support_rms_displacement"] <= 0.0
            ):
                raise ValueError(f"{candidate} has no cadence event response")
            seed_results[str(seed)] = {
                "parameter_count": parameters,
                "encoder_token_count": CONTEXT_LENGTH,
                "input_cadence_count": CONTEXT_LENGTH,
                "patch_stride": model.config.patch_stride,
                "h_cadence_shape": list(original_output["h_cadence"].shape),
                "h_window_shape": list(original_output["h_window"].shape),
                "z_window_shape": list(original_output["z_window"].shape),
                "patch_mean_is_exact_identity": True,
                "patch_time_is_exact_identity": True,
                "invalid_fill_invariant": True,
                "finite_outputs": True,
                "cadence_event_response": h_response,
                "h_window_event_l2_mean": float(
                    torch.linalg.vector_norm(
                        injected_output["h_window"]
                        - original_output["h_window"][source_indices],
                        dim=1,
                    ).mean()
                ),
                "z_window_event_l2_mean": float(
                    torch.linalg.vector_norm(
                        injected_output["z_window"]
                        - original_output["z_window"][source_indices],
                        dim=1,
                    ).mean()
                ),
            }
        results[candidate] = {
            "architecture": architecture,
            "template_variant": template,
            "seeds": seed_results,
        }

    return {
        "schema_version": RESULT_SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "status": "passed",
        "config": {
            "path": str(config_source),
            "sha256": config_sha256,
            "schema_version": config["schema_version"],
        },
        "cadence_contract": dict(config["cadence_contract"]),
        "experiment": {
            "synthetic_sources": SYNTHETIC_SOURCES,
            "event_conditions": conditions,
            "optimizer_or_training_steps": 0,
            "checkpoint_access": False,
            "development_or_sealed_data_access": False,
        },
        "candidates": results,
        "claim_limit": config["claim_limit"],
    }
