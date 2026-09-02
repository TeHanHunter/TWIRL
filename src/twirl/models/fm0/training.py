"""Losses and deterministic synthetic training for TWIRL-FM0.1."""
from __future__ import annotations

import hashlib
import math
import os
import random
import time
from collections.abc import Mapping
from contextlib import nullcontext
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .cadence_objective import position_centered_cadence_vicreg_loss
from .dataset import (
    FM0ReleaseDataset,
    SyntheticFM0Dataset,
    collate_fm0_samples,
    move_batch_to_device,
)
from .model import FM0ModelConfig, require_torch

try:  # Optional platform dependency.
    import torch
    from torch.nn import functional as F
except ImportError:  # pragma: no cover - local base environment omits Torch
    torch = None  # type: ignore[assignment]
    F = None  # type: ignore[assignment]


CHECKPOINT_SCHEMA_VERSION = "twirl_fm0_1_checkpoint_v1"
OBJECTIVE_STATE_SCHEMA_V2 = "twirl_fm0_objective_state_v2"
CADENCE_VICREG_OBJECTIVE_IDENTITY = (
    "masked_reconstruction_plus_position_centered_cadence_vicreg_v1"
)
FM0_3_VARIANTS = frozenset({"TWIRL-FM0.3.1", "TWIRL-FM0.3.2"})
FM0_3_WINDOW_LENGTH = 128
FM0_3_MASK_TARGET_FRACTION = 0.15
FM0_3_MASK_SPAN_RANGE = (1, 4)


def _objective_state(
    *,
    use_vicreg: bool,
    reconstruct_second_view: bool,
    objective_identity: str | None = None,
) -> dict[str, Any]:
    """Return the checkpoint-bound objective switches for one invocation."""

    if type(use_vicreg) is not bool or type(reconstruct_second_view) is not bool:
        raise ValueError("FM0 objective switches must be boolean")
    if reconstruct_second_view and not use_vicreg:
        raise ValueError("a second reconstruction view requires paired VICReg input")
    state: dict[str, Any] = {
        "use_vicreg": use_vicreg,
        "reconstruct_second_view": reconstruct_second_view,
    }
    if objective_identity is None:
        return state
    if objective_identity != CADENCE_VICREG_OBJECTIVE_IDENTITY:
        raise ValueError(f"unsupported explicit FM0 objective: {objective_identity!r}")
    if not use_vicreg or reconstruct_second_view:
        raise ValueError(
            "FM0.3 cadence VICReg requires paired views and first-mask-only "
            "reconstruction"
        )
    return {
        "schema_version": OBJECTIVE_STATE_SCHEMA_V2,
        "identity": objective_identity,
        **state,
    }


def _validated_objective_state(
    value: Any,
    *,
    label: str,
) -> dict[str, Any]:
    """Validate the exact, intentionally small objective-state schema."""

    legacy_keys = {"use_vicreg", "reconstruct_second_view"}
    explicit_keys = legacy_keys | {"schema_version", "identity"}
    observed_keys = frozenset(value) if isinstance(value, Mapping) else frozenset()
    if not isinstance(value, Mapping) or observed_keys not in {
        frozenset(legacy_keys),
        frozenset(explicit_keys),
    }:
        raise ValueError(
            f"{label} must contain exactly {sorted(legacy_keys)} or "
            f"{sorted(explicit_keys)}"
        )
    if any(type(value[name]) is not bool for name in legacy_keys):
        raise ValueError(f"{label} switches must be boolean")
    if observed_keys == frozenset(explicit_keys) and value.get(
        "schema_version"
    ) != OBJECTIVE_STATE_SCHEMA_V2:
        raise ValueError(f"{label} explicit schema version is invalid")
    return _objective_state(
        use_vicreg=value["use_vicreg"],
        reconstruct_second_view=value["reconstruct_second_view"],
        objective_identity=value.get("identity"),
    )


def _legacy_objective_state(payload: Mapping[str, Any]) -> dict[str, bool] | None:
    """Infer objective switches for trusted checkpoints predating the field.

    The frozen FM0.1 ladder makes this inference unambiguous.  Unknown ad-hoc
    checkpoints fail closed instead of being resumed under a possibly different
    loss.
    """

    contract = payload.get("run_contract")
    if not isinstance(contract, Mapping):
        return None
    if type(contract.get("use_vicreg")) is bool and type(
        contract.get("reconstruct_second_view")
    ) is bool:
        return _objective_state(
            use_vicreg=contract["use_vicreg"],
            reconstruct_second_view=contract["reconstruct_second_view"],
        )
    variant = contract.get("variant")
    if variant == "TWIRL-FM0.1.5":
        return _objective_state(use_vicreg=True, reconstruct_second_view=True)
    if variant in {
        "TWIRL-FM0.1.1",
        "TWIRL-FM0.1.2",
        "TWIRL-FM0.1.3",
        "TWIRL-FM0.1.4",
    }:
        return _objective_state(use_vicreg=False, reconstruct_second_view=False)
    return None


def _optional_contract_bool(
    contract: Mapping[str, Any], name: str
) -> bool | None:
    if name not in contract:
        return None
    value = contract[name]
    if type(value) is not bool:
        raise ValueError(f"FM0 run contract {name} must be boolean")
    return value


def _immutable_milestone_steps(contract: Mapping[str, Any]) -> tuple[int, ...]:
    raw = contract.get("immutable_milestone_steps", ())
    if not isinstance(raw, (list, tuple)):
        raise ValueError("FM0 immutable milestone steps must be a sequence")
    steps: list[int] = []
    for value in raw:
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(
                "FM0 immutable milestone steps must be nonnegative integers"
            )
        steps.append(value)
    if len(set(steps)) != len(steps) or steps != sorted(steps):
        raise ValueError("FM0 immutable milestone steps must be unique and sorted")
    return tuple(steps)


def _synthetic_dataset_config_contract(config: Any) -> dict[str, Any]:
    """Serialize synthetic config without drifting legacy checkpoint keys."""

    serialized = asdict(config)
    for name in ("mask_target_fraction", "mask_span_range"):
        if serialized.get(name) is None:
            serialized.pop(name, None)
    return serialized


def _dataset_contract(dataset: Any) -> dict[str, Any]:
    if isinstance(dataset, SyntheticFM0Dataset):
        return {
            "kind": "synthetic",
            "config": _synthetic_dataset_config_contract(dataset.config),
        }
    if isinstance(dataset, FM0ReleaseDataset):
        return dict(dataset.contract)
    raise TypeError("unsupported FM0 training dataset")


def _validate_cadence_objective_contract(
    *,
    model: Any,
    dataset: Any,
    run_contract: Mapping[str, Any],
) -> None:
    """Fail closed unless the complete FM0.3 cadence contract is explicit."""

    variant = run_contract.get("variant")
    if variant not in FM0_3_VARIANTS:
        raise ValueError("cadence VICReg requires an explicit FM0.3 variant")
    if run_contract.get("objective") != CADENCE_VICREG_OBJECTIVE_IDENTITY:
        raise ValueError("FM0.3 run contract lacks the cadence objective identity")
    if run_contract.get("mask_target_fraction") != FM0_3_MASK_TARGET_FRACTION:
        raise ValueError("FM0.3 run contract mask target fraction is not frozen")
    if run_contract.get("mask_span_range") != list(FM0_3_MASK_SPAN_RANGE):
        raise ValueError("FM0.3 run contract mask span range is not frozen")

    dataset_config = getattr(dataset, "config", None)
    if dataset_config is None or dataset_config.variant != variant:
        raise ValueError("FM0.3 dataset variant differs from the run contract")
    if dataset_config.window_length != FM0_3_WINDOW_LENGTH:
        raise ValueError("FM0.3 dataset context length must be 128 cadences")
    if (
        dataset_config.mask_target_fraction != FM0_3_MASK_TARGET_FRACTION
        or dataset_config.mask_span_range != FM0_3_MASK_SPAN_RANGE
    ):
        raise ValueError("FM0.3 dataset mask geometry differs from the run contract")

    model_config = getattr(model, "config", None)
    if (
        model_config is None
        or model_config.window_length != FM0_3_WINDOW_LENGTH
        or model_config.patch_stride != 1
    ):
        raise ValueError(
            "FM0.3 cadence VICReg requires 128 inputs and one token per cadence"
        )


def _cadence_visibility_mask(batch: Mapping[str, Any]) -> Any:
    """Return natural, present, and unmasked cadence eligibility for one view."""

    required = ("flux_valid", "time_valid", "view_present", "reconstruction_mask")
    missing = [name for name in required if name not in batch]
    if missing:
        raise ValueError(f"FM0 cadence objective batch lacks fields: {missing}")
    flux_valid = batch["flux_valid"]
    reconstruction_mask = batch["reconstruction_mask"]
    time_valid = batch["time_valid"]
    view_present = batch["view_present"]
    if flux_valid.ndim != 3 or reconstruction_mask.shape != flux_valid.shape:
        raise ValueError("FM0 cadence objective flux masks must have shape [B,V,L]")
    if time_valid.shape != (flux_valid.shape[0], flux_valid.shape[2]):
        raise ValueError("FM0 cadence objective time_valid must have shape [B,L]")
    if view_present.shape != flux_valid.shape[:2]:
        raise ValueError("FM0 cadence objective view_present must have shape [B,V]")

    available = flux_valid.bool() & view_present.bool()[:, :, None]
    natural = time_valid.bool() & torch.any(available, dim=1)
    reconstruction = reconstruction_mask.bool()
    return natural & ~torch.any(reconstruction, dim=1)


def _paired_cadence_visibility_masks(
    first_batch: Mapping[str, Any],
    second_batch: Mapping[str, Any],
) -> tuple[Any, Any]:
    """Validate paired natural masks and return each independently visible set."""

    source_fields = (
        "flux",
        "flux_valid",
        "flux_error",
        "error_valid",
        "local_time_cadences",
        "delta_time_cadences",
        "time_valid",
        "segment_boundary",
        "view_present",
    )
    for name in source_fields:
        if name not in first_batch or name not in second_batch:
            raise ValueError(f"FM0 cadence objective paired batch lacks {name}")
        if not torch.equal(first_batch[name], second_batch[name]):
            raise ValueError(
                f"FM0 cadence objective paired source tensor differs: {name}"
            )
    return (
        _cadence_visibility_mask(first_batch),
        _cadence_visibility_mask(second_batch),
    )


def _project_cadence_tokens(model: Any, output: Mapping[str, Any]) -> Any:
    """Apply the shared projection tokenwise without constructing a window mean."""

    hidden = output.get("h_cadence")
    token_valid = output.get("token_valid")
    if hidden is None or token_valid is None or hidden.ndim != 3:
        raise ValueError("FM0.3 model output lacks cadence-indexed hidden states")
    if token_valid.shape != hidden.shape[:2]:
        raise ValueError("FM0.3 cadence hidden states and validity mask differ")
    projected = model.embedding_projection(hidden)
    if projected.shape[:2] != hidden.shape[:2] or projected.ndim != 3:
        raise ValueError("FM0.3 cadence projection changed the token lattice")
    return projected


@dataclass(frozen=True)
class FM0OptimizationConfig:
    learning_rate: float = 3.0e-4
    weight_decay: float = 0.01
    warmup_steps: int = 1000
    max_optimizer_steps: int = 20_000
    effective_batch_windows: int = 64
    huber_delta: float = 0.01
    vicreg_total_weight: float = 0.1
    vicreg_invariance_weight: float = 25.0
    vicreg_variance_weight: float = 25.0
    vicreg_covariance_weight: float = 1.0

    def __post_init__(self) -> None:
        if self.learning_rate <= 0 or self.weight_decay < 0:
            raise ValueError("invalid optimizer hyperparameters")
        if self.warmup_steps < 0 or self.max_optimizer_steps <= 0:
            raise ValueError("invalid scheduler step counts")
        if self.effective_batch_windows <= 0 or self.huber_delta <= 0:
            raise ValueError("batch size and Huber delta must be positive")
        if any(
            not math.isfinite(value) or value < 0
            for value in (
                self.vicreg_total_weight,
                self.vicreg_invariance_weight,
                self.vicreg_variance_weight,
                self.vicreg_covariance_weight,
            )
        ):
            raise ValueError("VICReg weights must be nonnegative")


def seed_everything(seed: int, *, deterministic: bool = True) -> None:
    """Seed Python, NumPy, Torch, and every visible CUDA device."""

    require_torch()
    random.seed(int(seed))
    np.random.seed(int(seed) % (2**32))
    torch.manual_seed(int(seed))
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(int(seed))
    if deterministic:
        torch.use_deterministic_algorithms(True)
        if hasattr(torch.backends, "cudnn"):
            torch.backends.cudnn.benchmark = False
            torch.backends.cudnn.deterministic = True


def masked_huber_reconstruction_loss(
    prediction: Any,
    target: Any,
    target_mask: Any,
    *,
    delta: float = 0.01,
) -> tuple[Any, dict[str, Any]]:
    """Average per-view masked Huber losses with equal view weight."""

    require_torch()
    if prediction.shape != target.shape or prediction.shape != target_mask.shape:
        raise ValueError("prediction, target, and target_mask shapes must match")
    if prediction.ndim != 3:
        raise ValueError("reconstruction tensors must have shape [batch, view, cadence]")
    mask = target_mask.bool()
    point_loss = F.huber_loss(
        prediction,
        target,
        reduction="none",
        delta=float(delta),
    )
    counts = mask.sum(dim=(0, 2))
    sums = (point_loss * mask.to(dtype=point_loss.dtype)).sum(dim=(0, 2))
    per_view = sums / counts.clamp_min(1).to(dtype=sums.dtype)
    active = counts > 0
    if not bool(torch.any(active)):
        raise ValueError("masked reconstruction batch has no valid targets")
    loss = per_view[active].mean()
    return loss, {
        "loss": loss,
        "per_view": per_view,
        "target_counts": counts,
        "active_views": active,
    }


def _off_diagonal(matrix: Any) -> Any:
    size = matrix.shape[0]
    if matrix.ndim != 2 or matrix.shape[1] != size:
        raise ValueError("covariance matrix must be square")
    return matrix.flatten()[:-1].view(size - 1, size + 1)[:, 1:].flatten()


def vicreg_loss(
    first: Any,
    second: Any,
    *,
    invariance_weight: float = 25.0,
    variance_weight: float = 25.0,
    covariance_weight: float = 1.0,
) -> tuple[Any, dict[str, Any]]:
    """Standard same-window VICReg objective with frozen 25:25:1 weights."""

    require_torch()
    if first.shape != second.shape or first.ndim != 2:
        raise ValueError("VICReg embeddings must share shape [batch, feature]")
    invariance = F.mse_loss(first, second)
    centered_first = first - first.mean(dim=0)
    centered_second = second - second.mean(dim=0)
    std_first = torch.sqrt(centered_first.var(dim=0, unbiased=False) + 1.0e-4)
    std_second = torch.sqrt(centered_second.var(dim=0, unbiased=False) + 1.0e-4)
    variance = 0.5 * (
        F.relu(1.0 - std_first).mean() + F.relu(1.0 - std_second).mean()
    )
    denominator = max(1, first.shape[0] - 1)
    covariance_first = centered_first.T @ centered_first / denominator
    covariance_second = centered_second.T @ centered_second / denominator
    feature_count = max(1, first.shape[1])
    covariance = (
        _off_diagonal(covariance_first).pow(2).sum()
        + _off_diagonal(covariance_second).pow(2).sum()
    ) / feature_count
    total = (
        float(invariance_weight) * invariance
        + float(variance_weight) * variance
        + float(covariance_weight) * covariance
    )
    return total, {
        "loss": total,
        "invariance": invariance,
        "variance": variance,
        "covariance": covariance,
    }


def fm0_objective(
    model: Any,
    batch: dict[str, Any],
    *,
    second_batch: dict[str, Any] | None,
    config: FM0OptimizationConfig,
    reconstruct_second_view: bool = True,
) -> tuple[Any, dict[str, Any]]:
    """Compute masked reconstruction and optional same-window VICReg.

    This helper treats the supplied tensors as one complete batch.  The
    synthetic training driver does not call it once per microbatch for
    ``TWIRL-FM0.1.5``: VICReg variance and covariance must see the complete
    effective batch, so that path is implemented explicitly below.
    """

    first_output = model(batch)
    reconstruction, reconstruction_stats = masked_huber_reconstruction_loss(
        first_output["reconstruction"],
        batch["flux"].float(),
        batch["reconstruction_mask"],
        delta=config.huber_delta,
    )
    diagnostics: dict[str, Any] = {
        "reconstruction": reconstruction,
        "reconstruction_first": reconstruction,
    }
    if second_batch is None:
        diagnostics["vicreg"] = reconstruction.new_zeros(())
        diagnostics["total"] = reconstruction
        diagnostics["target_counts"] = reconstruction_stats["target_counts"]
        return reconstruction, diagnostics

    second_output = model(second_batch)
    second_reconstruction, _ = masked_huber_reconstruction_loss(
        second_output["reconstruction"],
        second_batch["flux"].float(),
        second_batch["reconstruction_mask"],
        delta=config.huber_delta,
    )
    reconstruction_mean = 0.5 * (reconstruction + second_reconstruction)
    optimized_reconstruction = (
        reconstruction_mean if reconstruct_second_view else reconstruction
    )
    if first_output["z_window"].shape[0] != config.effective_batch_windows:
        raise ValueError(
            "VICReg variance and covariance require the complete effective batch"
        )
    vicreg, vicreg_stats = vicreg_loss(
        first_output["z_window"],
        second_output["z_window"],
        invariance_weight=config.vicreg_invariance_weight,
        variance_weight=config.vicreg_variance_weight,
        covariance_weight=config.vicreg_covariance_weight,
    )
    total = optimized_reconstruction + config.vicreg_total_weight * vicreg
    diagnostics.update(
        {
            "reconstruction": optimized_reconstruction,
            "reconstruction_second": second_reconstruction,
            "reconstruction_mean": reconstruction_mean,
            "vicreg": vicreg,
            "vicreg_invariance": vicreg_stats["invariance"],
            "vicreg_variance": vicreg_stats["variance"],
            "vicreg_covariance": vicreg_stats["covariance"],
            "total": total,
            "target_counts": reconstruction_stats["target_counts"],
        }
    )
    return total, diagnostics


def _scheduler_multiplier(step: int, config: FM0OptimizationConfig) -> float:
    if config.warmup_steps and step < config.warmup_steps:
        return float(step + 1) / float(config.warmup_steps)
    remaining = max(1, config.max_optimizer_steps - config.warmup_steps)
    progress = min(1.0, max(0.0, (step - config.warmup_steps) / remaining))
    return 0.5 * (1.0 + math.cos(math.pi * progress))


def make_optimizer_and_scheduler(
    model: Any,
    config: FM0OptimizationConfig,
) -> tuple[Any, Any]:
    require_torch()
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=config.learning_rate,
        weight_decay=config.weight_decay,
    )
    scheduler = torch.optim.lr_scheduler.LambdaLR(
        optimizer,
        lr_lambda=lambda step: _scheduler_multiplier(step, config),
    )
    return optimizer, scheduler


def _capture_rng_state() -> dict[str, Any]:
    require_torch()
    return {
        "python": random.getstate(),
        "numpy": np.random.get_state(),
        "torch_cpu": torch.get_rng_state(),
        "torch_cuda": torch.cuda.get_rng_state_all() if torch.cuda.is_available() else [],
    }


def _restore_rng_state(state: Mapping[str, Any]) -> None:
    require_torch()
    random.setstate(state["python"])
    np.random.set_state(state["numpy"])
    torch.set_rng_state(state["torch_cpu"])
    if torch.cuda.is_available() and state.get("torch_cuda"):
        torch.cuda.set_rng_state_all(state["torch_cuda"])


def _atomic_torch_save(payload: Mapping[str, Any], path: Path) -> None:
    require_torch()
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    torch.save(dict(payload), temporary)
    os.replace(temporary, path)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_write_sha256_sidecar(path: Path) -> str:
    digest = _sha256_file(path)
    sidecar = path.with_name(path.name + ".sha256")
    temporary = sidecar.with_name(f".{sidecar.name}.tmp-{os.getpid()}")
    temporary.write_text(f"{digest}  {path.name}\n", encoding="utf-8")
    os.replace(temporary, sidecar)
    return digest


def _verify_sha256_sidecar(path: Path, *, required: bool) -> str | None:
    sidecar = path.with_name(path.name + ".sha256")
    if not sidecar.is_file():
        if required:
            raise ValueError(f"missing FM0 checkpoint SHA256 sidecar: {sidecar}")
        return None
    fields = sidecar.read_text(encoding="utf-8").strip().split()
    if len(fields) != 2 or fields[1] != path.name:
        raise ValueError(f"malformed FM0 checkpoint SHA256 sidecar: {sidecar}")
    observed = _sha256_file(path)
    if fields[0] != observed:
        raise ValueError(f"FM0 checkpoint SHA256 mismatch: {path}")
    return observed


def checkpoint_payload(
    *,
    model: Any,
    optimizer: Any,
    scheduler: Any,
    global_step: int,
    run_contract: Mapping[str, Any],
    optimization: FM0OptimizationConfig,
    dataset: Any,
    loss_history: list[dict[str, float]],
    objective_state: Mapping[str, Any],
) -> dict[str, Any]:
    """Capture every state needed for exact optimizer-step continuation."""

    model_config: FM0ModelConfig = model.config
    return {
        "schema_version": CHECKPOINT_SCHEMA_VERSION,
        "model_state": model.state_dict(),
        "optimizer_state": optimizer.state_dict(),
        "scheduler_state": scheduler.state_dict(),
        "scaler_state": None,
        "rng_state": _capture_rng_state(),
        "progress": {"global_step": int(global_step), "epoch": int(dataset.epoch)},
        "sampler_state": {
            "next_absolute_sample": int(
                global_step * optimization.effective_batch_windows
            ),
            "dataset_epoch": int(dataset.epoch),
        },
        "model_config": asdict(model_config),
        "optimization_config": asdict(optimization),
        "dataset_contract": _dataset_contract(dataset),
        "synthetic_dataset_config": (
            _synthetic_dataset_config_contract(dataset.config)
            if isinstance(dataset, SyntheticFM0Dataset)
            else None
        ),
        "run_contract": dict(run_contract),
        "objective_state": _validated_objective_state(
            objective_state,
            label="FM0 checkpoint objective state",
        ),
        "loss_history": list(loss_history),
    }


def save_rotating_checkpoint(
    output_dir: str | Path,
    payload: Mapping[str, Any],
) -> Path:
    output = Path(output_dir)
    latest = output / "checkpoint_latest.pt"
    previous = output / "checkpoint_previous.pt"
    if latest.is_file():
        os.replace(latest, previous)
    _atomic_torch_save(payload, latest)
    return latest


def save_immutable_milestone_checkpoint(
    output_dir: str | Path,
    payload: Mapping[str, Any],
) -> Path:
    """Write one step-addressed checkpoint without replacing an existing file."""

    progress = payload.get("progress")
    raw_step = progress.get("global_step") if isinstance(progress, Mapping) else None
    if isinstance(raw_step, bool) or not isinstance(raw_step, int) or raw_step < 0:
        raise ValueError("FM0 milestone checkpoint has no valid global step")
    step = raw_step
    destination = Path(output_dir) / f"checkpoint_step_{step:08d}.pt"
    sidecar = destination.with_name(destination.name + ".sha256")
    if destination.exists() or sidecar.exists():
        raise FileExistsError(
            "refusing to replace immutable FM0 milestone or sidecar: "
            f"{destination}"
        )
    _atomic_torch_save(payload, destination)
    _atomic_write_sha256_sidecar(destination)
    return destination


def _immutable_milestone_artifacts(
    output_dir: str | Path,
    milestone_steps: tuple[int, ...],
    *,
    global_step: int,
) -> list[dict[str, Any]]:
    """Return verified bindings for every declared milestone already reached."""

    output = Path(output_dir)
    artifacts: list[dict[str, Any]] = []
    for step in milestone_steps:
        if step > global_step:
            break
        checkpoint = output / f"checkpoint_step_{step:08d}.pt"
        if not checkpoint.is_file():
            raise ValueError(f"missing reached FM0 milestone checkpoint: {checkpoint}")
        digest = _verify_sha256_sidecar(checkpoint, required=True)
        artifacts.append(
            {
                "step": step,
                "checkpoint": str(checkpoint.resolve()),
                "sha256_sidecar": str(
                    checkpoint.with_name(checkpoint.name + ".sha256").resolve()
                ),
                "sha256": digest,
            }
        )
    return artifacts


def load_checkpoint(
    path: str | Path,
    *,
    model: Any,
    optimizer: Any,
    scheduler: Any,
    expected_run_contract: Mapping[str, Any],
    expected_optimization: FM0OptimizationConfig,
    expected_objective_state: Mapping[str, Any],
    dataset: Any,
    expected_global_step: int | None = None,
) -> tuple[int, list[dict[str, float]]]:
    """Load and validate a full FM0 checkpoint, including RNG and sampler state."""

    require_torch()
    checkpoint_path = Path(path)
    _verify_sha256_sidecar(
        checkpoint_path,
        required=checkpoint_path.name.startswith("checkpoint_step_"),
    )
    payload = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    if not isinstance(payload, Mapping):
        raise ValueError("FM0 checkpoint is not a mapping")
    if payload.get("schema_version") != CHECKPOINT_SCHEMA_VERSION:
        raise ValueError("FM0 checkpoint schema mismatch")
    if payload.get("run_contract") != dict(expected_run_contract):
        raise ValueError("FM0 checkpoint run contract mismatch")
    if payload.get("optimization_config") != asdict(expected_optimization):
        raise ValueError("FM0 checkpoint optimization contract mismatch")
    expected_objective = _validated_objective_state(
        expected_objective_state,
        label="expected FM0 objective state",
    )
    raw_objective = payload.get("objective_state")
    if raw_objective is None:
        observed_objective = _legacy_objective_state(payload)
    else:
        observed_objective = _validated_objective_state(
            raw_objective,
            label="FM0 checkpoint objective state",
        )
    if observed_objective is None:
        raise ValueError("FM0 checkpoint lacks a provable objective state")
    if observed_objective != expected_objective:
        raise ValueError("FM0 checkpoint objective state mismatch")
    if payload.get("model_config") != asdict(model.config):
        raise ValueError("FM0 checkpoint model configuration mismatch")
    observed_dataset_contract = payload.get("dataset_contract")
    if observed_dataset_contract is None and isinstance(dataset, SyntheticFM0Dataset):
        observed_dataset_contract = {
            "kind": "synthetic",
            "config": payload.get("synthetic_dataset_config"),
        }
    if observed_dataset_contract != _dataset_contract(dataset):
        raise ValueError("FM0 checkpoint dataset configuration mismatch")
    model.load_state_dict(payload["model_state"], strict=True)
    optimizer.load_state_dict(payload["optimizer_state"])
    scheduler.load_state_dict(payload["scheduler_state"])
    _restore_rng_state(payload["rng_state"])
    progress = payload.get("progress", {})
    global_step = int(progress.get("global_step", -1))
    if global_step < 0:
        raise ValueError("FM0 checkpoint has invalid global step")
    if expected_global_step is not None and global_step != expected_global_step:
        raise ValueError(
            "FM0 checkpoint global step differs from the required resume step"
        )
    progress_epoch = int(progress.get("epoch", 0))
    sampler_epoch = int(payload.get("sampler_state", {}).get("dataset_epoch", -1))
    if sampler_epoch != progress_epoch:
        raise ValueError("FM0 checkpoint dataset epoch is inconsistent")
    dataset.set_epoch(progress_epoch)
    expected_sample = global_step * expected_optimization.effective_batch_windows
    if int(payload.get("sampler_state", {}).get("next_absolute_sample", -1)) != expected_sample:
        raise ValueError("FM0 checkpoint sampler state is inconsistent")
    history = [dict(item) for item in payload.get("loss_history", [])]
    return global_step, history


def _sample_for_absolute_index(
    dataset: Any,
    absolute_index: int,
    *,
    mask_view: int,
) -> dict[str, np.ndarray]:
    epoch, local_index = divmod(int(absolute_index), len(dataset))
    dataset.set_epoch(epoch)
    return dataset.sample(local_index, mask_view=mask_view)


def _autocast_context(device: Any, precision: str) -> Any:
    if precision == "fp32":
        return nullcontext()
    if precision != "bf16":
        raise ValueError("precision must be 'fp32' or 'bf16'")
    return torch.autocast(device_type=device.type, dtype=torch.bfloat16)


def _capture_forward_rng_state(device: Any) -> dict[str, Any]:
    """Capture only the Torch RNG state consumed by a model forward pass."""

    require_torch()
    state: dict[str, Any] = {"cpu": torch.get_rng_state().clone()}
    if device.type == "cuda":
        state["cuda"] = torch.cuda.get_rng_state(device).clone()
    return state


def _restore_forward_rng_state(state: Mapping[str, Any], device: Any) -> None:
    """Restore a forward-pass RNG state for exact dropout replay."""

    require_torch()
    torch.set_rng_state(state["cpu"])
    if device.type == "cuda":
        torch.cuda.set_rng_state(state["cuda"], device=device)


def _embedding_projection_gradient_norm(model: Any) -> float:
    """Return the L2 norm reaching the shared window/cadence projection."""

    require_torch()
    squared = None
    for parameter in model.embedding_projection.parameters():
        if parameter.grad is None:
            continue
        value = parameter.grad.detach().float().pow(2).sum()
        squared = value if squared is None else squared + value
    if squared is None:
        return 0.0
    return float(torch.sqrt(squared).cpu())


def _backprop_effective_batch_vicreg(
    *,
    model: Any,
    replay_records: list[dict[str, Any]],
    first_embeddings: list[Any],
    second_embeddings: list[Any],
    config: FM0OptimizationConfig,
    device: Any,
    precision: str,
) -> tuple[Any, dict[str, Any]]:
    """Backpropagate one VICReg objective over the full effective batch.

    Reconstruction gradients are already released one microbatch at a time.
    Holding every encoder activation until the effective batch is complete
    would defeat gradient accumulation, so this routine keeps only detached
    window embeddings.  It computes the exact full-batch VICReg gradient with
    respect to those embeddings, then replays one microbatch at a time and
    supplies the corresponding embedding gradient to the encoder.  Replaying
    the saved Torch RNG states makes dropout identical to the first forward
    pass while keeping peak activation memory at one microbatch.
    """

    require_torch()
    if not replay_records or len(first_embeddings) != len(replay_records):
        raise ValueError("VICReg replay records are incomplete")
    if len(second_embeddings) != len(replay_records):
        raise ValueError("VICReg second-view replay records are incomplete")

    full_first = torch.cat(first_embeddings, dim=0).float().requires_grad_(True)
    full_second = torch.cat(second_embeddings, dim=0).float().requires_grad_(True)
    if full_first.shape[0] != config.effective_batch_windows:
        raise ValueError("VICReg did not receive the complete effective batch")
    vicreg, diagnostics = vicreg_loss(
        full_first,
        full_second,
        invariance_weight=config.vicreg_invariance_weight,
        variance_weight=config.vicreg_variance_weight,
        covariance_weight=config.vicreg_covariance_weight,
    )
    first_gradient, second_gradient = torch.autograd.grad(
        vicreg, (full_first, full_second)
    )

    # Computing the detached full-batch statistics consumes no model RNG.  The
    # final state after the original forward passes is nevertheless restored
    # after replay so the next optimizer step follows the uninterrupted stream.
    post_forward_rng = _capture_forward_rng_state(device)
    offset = 0
    try:
        for record in replay_records:
            count = int(record["batch_windows"])
            _restore_forward_rng_state(record["first_rng"], device)
            scale = float(config.vicreg_total_weight)
            with _autocast_context(device, precision):
                first_z = model(record["first_batch"])["z_window"]
            first_z.backward(
                (first_gradient[offset : offset + count] * scale).to(
                    first_z.dtype
                )
            )
            del first_z
            _restore_forward_rng_state(record["second_rng"], device)
            with _autocast_context(device, precision):
                second_z = model(record["second_batch"])["z_window"]
            second_z.backward(
                (second_gradient[offset : offset + count] * scale).to(
                    second_z.dtype
                )
            )
            del second_z
            offset += count
    finally:
        _restore_forward_rng_state(post_forward_rng, device)
    if offset != config.effective_batch_windows:
        raise ValueError("VICReg replay did not cover the effective batch")
    return vicreg.detach(), {
        name: value.detach() for name, value in diagnostics.items()
    }


def _backprop_effective_batch_cadence_vicreg(
    *,
    model: Any,
    replay_records: list[dict[str, Any]],
    first_embeddings: list[Any],
    second_embeddings: list[Any],
    first_visible_masks: list[Any],
    second_visible_masks: list[Any],
    config: FM0OptimizationConfig,
    device: Any,
    precision: str,
) -> tuple[Any, dict[str, Any]]:
    """Replay exact forwards for full-batch cadence VICReg gradients."""

    require_torch()
    record_count = len(replay_records)
    if not replay_records or any(
        len(values) != record_count
        for values in (
            first_embeddings,
            second_embeddings,
            first_visible_masks,
            second_visible_masks,
        )
    ):
        raise ValueError("cadence VICReg replay records are incomplete")

    full_first = torch.cat(first_embeddings, dim=0).float().requires_grad_(True)
    full_second = torch.cat(second_embeddings, dim=0).float().requires_grad_(True)
    full_first_visible = torch.cat(first_visible_masks, dim=0).bool()
    full_second_visible = torch.cat(second_visible_masks, dim=0).bool()
    if full_first.shape[0] != config.effective_batch_windows:
        raise ValueError("cadence VICReg did not receive the complete effective batch")
    cadence_vicreg, diagnostics = position_centered_cadence_vicreg_loss(
        full_first,
        full_second,
        full_first_visible,
        full_second_visible,
        invariance_weight=config.vicreg_invariance_weight,
        variance_weight=config.vicreg_variance_weight,
        covariance_weight=config.vicreg_covariance_weight,
    )
    first_gradient, second_gradient = torch.autograd.grad(
        cadence_vicreg, (full_first, full_second)
    )

    post_forward_rng = _capture_forward_rng_state(device)
    offset = 0
    try:
        for record in replay_records:
            count = int(record["batch_windows"])
            _restore_forward_rng_state(record["first_rng"], device)
            with _autocast_context(device, precision):
                first_output = model(record["first_batch"])
                first_z = _project_cadence_tokens(model, first_output)
            first_z.backward(
                (
                    first_gradient[offset : offset + count]
                    * float(config.vicreg_total_weight)
                ).to(first_z.dtype)
            )
            del first_output, first_z

            _restore_forward_rng_state(record["second_rng"], device)
            with _autocast_context(device, precision):
                second_output = model(record["second_batch"])
                second_z = _project_cadence_tokens(model, second_output)
            second_z.backward(
                (
                    second_gradient[offset : offset + count]
                    * float(config.vicreg_total_weight)
                ).to(second_z.dtype)
            )
            del second_output, second_z
            offset += count
    finally:
        _restore_forward_rng_state(post_forward_rng, device)
    if offset != config.effective_batch_windows:
        raise ValueError("cadence VICReg replay did not cover the effective batch")
    return cadence_vicreg.detach(), {
        name: value.detach() for name, value in diagnostics.items()
    }


def _run_training(
    *,
    model: Any,
    dataset: Any,
    output_dir: str | Path,
    run_contract: Mapping[str, Any],
    optimization: FM0OptimizationConfig,
    target_step: int,
    micro_batch_windows: int,
    device: str,
    precision: str,
    use_vicreg: bool,
    reconstruct_second_view: bool,
    objective_identity: str | None,
    synthetic_only: bool,
    checkpoint_interval_seconds: float,
    progress_interval_steps: int,
    resume_checkpoint: str | Path | None = None,
    expected_resume_step: int | None = None,
) -> dict[str, Any]:
    """Run deterministic FM0 training and write full rotating checkpoints."""

    require_torch()
    invocation_objective = _objective_state(
        use_vicreg=use_vicreg,
        reconstruct_second_view=reconstruct_second_view,
        objective_identity=objective_identity,
    )
    cadence_vicreg = objective_identity == CADENCE_VICREG_OBJECTIVE_IDENTITY
    if run_contract.get("objective") == CADENCE_VICREG_OBJECTIVE_IDENTITY and not (
        cadence_vicreg
    ):
        raise ValueError("FM0.3 run contract and explicit objective path disagree")
    if cadence_vicreg:
        _validate_cadence_objective_contract(
            model=model,
            dataset=dataset,
            run_contract=run_contract,
        )
        if "target_step" in run_contract:
            raise ValueError("FM0.3 invocation stop must remain runtime-only")
        if run_contract.get("training_horizon_step") != (
            optimization.max_optimizer_steps
        ):
            raise ValueError("FM0.3 run contract optimizer horizon differs")
        if (
            run_contract.get(
                "stop_after_step_is_execution_state_not_scientific_contract"
            )
            is not True
        ):
            raise ValueError("FM0.3 runtime-only stop declaration is missing")
    contract_use_vicreg = _optional_contract_bool(run_contract, "use_vicreg")
    if contract_use_vicreg is not None and contract_use_vicreg != use_vicreg:
        raise ValueError("FM0 run contract and objective path disagree")
    contract_reconstruct_second = _optional_contract_bool(
        run_contract, "reconstruct_second_view"
    )
    if (
        contract_reconstruct_second is not None
        and contract_reconstruct_second != reconstruct_second_view
    ):
        raise ValueError("FM0 run contract and reconstruction path disagree")
    milestone_steps = _immutable_milestone_steps(run_contract)
    if target_step <= 0 or target_step > optimization.max_optimizer_steps:
        raise ValueError("target_step is outside the frozen optimizer budget")
    if milestone_steps and milestone_steps[-1] > optimization.max_optimizer_steps:
        raise ValueError("immutable milestone schedule exceeds optimizer horizon")
    if micro_batch_windows <= 0:
        raise ValueError("micro_batch_windows must be positive")
    if optimization.effective_batch_windows % micro_batch_windows:
        raise ValueError("effective batch must be divisible by micro batch")
    if checkpoint_interval_seconds < 0 or progress_interval_steps <= 0:
        raise ValueError("checkpoint/progress intervals are invalid")
    resolved_device = torch.device(device)
    if resolved_device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but torch.cuda.is_available() is false")
    model.to(resolved_device)
    optimizer, scheduler = make_optimizer_and_scheduler(model, optimization)
    global_step = 0
    loss_history: list[dict[str, float]] = []
    if expected_resume_step is not None and resume_checkpoint is None:
        raise ValueError("required FM0 resume checkpoint is absent")
    if resume_checkpoint is not None:
        global_step, loss_history = load_checkpoint(
            resume_checkpoint,
            model=model,
            optimizer=optimizer,
            scheduler=scheduler,
            expected_run_contract=run_contract,
            expected_optimization=optimization,
            expected_objective_state=invocation_objective,
            dataset=dataset,
            expected_global_step=expected_resume_step,
        )
    if global_step > target_step:
        raise ValueError("resume checkpoint is beyond requested target step")

    accumulation = optimization.effective_batch_windows // micro_batch_windows
    started = time.monotonic()
    last_checkpoint_at = started
    model.train()
    if global_step == 0 and resume_checkpoint is None and 0 in milestone_steps:
        save_immutable_milestone_checkpoint(
            output_dir,
            checkpoint_payload(
                model=model,
                optimizer=optimizer,
                scheduler=scheduler,
                global_step=global_step,
                run_contract=run_contract,
                optimization=optimization,
                dataset=dataset,
                loss_history=loss_history,
                objective_state=invocation_objective,
            ),
        )
    while global_step < target_step:
        optimizer.zero_grad(set_to_none=True)
        metric_sums: dict[str, float] = {}
        replay_records: list[dict[str, Any]] = []
        first_embeddings: list[Any] = []
        second_embeddings: list[Any] = []
        first_visible_masks: list[Any] = []
        second_visible_masks: list[Any] = []
        for micro_index in range(accumulation):
            batch_start = (
                global_step * optimization.effective_batch_windows
                + micro_index * micro_batch_windows
            )
            first_samples = [
                _sample_for_absolute_index(dataset, batch_start + offset, mask_view=0)
                for offset in range(micro_batch_windows)
            ]
            first_batch = move_batch_to_device(
                collate_fm0_samples(first_samples), resolved_device
            )
            second_batch = None
            if use_vicreg:
                second_samples = [
                    _sample_for_absolute_index(
                        dataset, batch_start + offset, mask_view=1
                    )
                    for offset in range(micro_batch_windows)
                ]
                second_batch = move_batch_to_device(
                    collate_fm0_samples(second_samples), resolved_device
                )
            if use_vicreg:
                if second_batch is None:  # pragma: no cover - defensive
                    raise RuntimeError("VICReg requires two masked views")
                if cadence_vicreg:
                    first_visible, second_visible = (
                        _paired_cadence_visibility_masks(first_batch, second_batch)
                    )
                first_rng = _capture_forward_rng_state(resolved_device)
                with _autocast_context(resolved_device, precision):
                    first_output = model(first_batch)
                    first_embedding = (
                        _project_cadence_tokens(model, first_output)
                        if cadence_vicreg
                        else first_output["z_window"]
                    )
                    first_reconstruction, _ = masked_huber_reconstruction_loss(
                        first_output["reconstruction"],
                        first_batch["flux"].float(),
                        first_batch["reconstruction_mask"],
                        delta=optimization.huber_delta,
                    )
                    second_rng = _capture_forward_rng_state(resolved_device)
                    if reconstruct_second_view:
                        second_output = model(second_batch)
                        second_embedding = (
                            _project_cadence_tokens(model, second_output)
                            if cadence_vicreg
                            else second_output["z_window"]
                        )
                        second_reconstruction, _ = masked_huber_reconstruction_loss(
                            second_output["reconstruction"],
                            second_batch["flux"].float(),
                            second_batch["reconstruction_mask"],
                            delta=optimization.huber_delta,
                        )
                    else:
                        # FM0.2 and FM0.3 use the second mask only to define the
                        # paired representation view. Its reconstruction is a
                        # no-gradient diagnostic; only the first mask is
                        # optimized in these objective contracts.
                        with torch.no_grad():
                            second_output = model(second_batch)
                            second_embedding = (
                                _project_cadence_tokens(model, second_output)
                                if cadence_vicreg
                                else second_output["z_window"]
                            )
                            second_reconstruction, _ = (
                                masked_huber_reconstruction_loss(
                                    second_output["reconstruction"],
                                    second_batch["flux"].float(),
                                    second_batch["reconstruction_mask"],
                                    delta=optimization.huber_delta,
                                )
                            )
                    reconstruction_mean = 0.5 * (
                        first_reconstruction + second_reconstruction
                    )
                    reconstruction = (
                        reconstruction_mean
                        if reconstruct_second_view
                        else first_reconstruction
                    )
                (reconstruction / accumulation).backward()
                first_embeddings.append(first_embedding.detach())
                second_embeddings.append(second_embedding.detach())
                if cadence_vicreg:
                    first_visible_masks.append(first_visible.detach())
                    second_visible_masks.append(second_visible.detach())
                replay_records.append(
                    {
                        "first_batch": first_batch,
                        "second_batch": second_batch,
                        "first_rng": first_rng,
                        "second_rng": second_rng,
                        "batch_windows": micro_batch_windows,
                    }
                )
                del first_output, second_output
                metric_sums["reconstruction"] = metric_sums.get(
                    "reconstruction", 0.0
                ) + float(reconstruction.detach().cpu()) / accumulation
                metric_sums["reconstruction_first"] = metric_sums.get(
                    "reconstruction_first", 0.0
                ) + float(first_reconstruction.detach().cpu()) / accumulation
                metric_sums["reconstruction_second"] = metric_sums.get(
                    "reconstruction_second", 0.0
                ) + float(second_reconstruction.detach().cpu()) / accumulation
                metric_sums["reconstruction_mean"] = metric_sums.get(
                    "reconstruction_mean", 0.0
                ) + float(reconstruction_mean.detach().cpu()) / accumulation
            else:
                with _autocast_context(resolved_device, precision):
                    loss, diagnostics = fm0_objective(
                        model,
                        first_batch,
                        second_batch=None,
                        config=optimization,
                    )
                    scaled_loss = loss / accumulation
                scaled_loss.backward()
                for name in ("total", "reconstruction", "vicreg"):
                    metric_sums[name] = metric_sums.get(name, 0.0) + float(
                        diagnostics[name].detach().cpu()
                    ) / accumulation
        if use_vicreg:
            if cadence_vicreg:
                vicreg, vicreg_diagnostics = (
                    _backprop_effective_batch_cadence_vicreg(
                        model=model,
                        replay_records=replay_records,
                        first_embeddings=first_embeddings,
                        second_embeddings=second_embeddings,
                        first_visible_masks=first_visible_masks,
                        second_visible_masks=second_visible_masks,
                        config=optimization,
                        device=resolved_device,
                        precision=precision,
                    )
                )
            else:
                vicreg, vicreg_diagnostics = _backprop_effective_batch_vicreg(
                    model=model,
                    replay_records=replay_records,
                    first_embeddings=first_embeddings,
                    second_embeddings=second_embeddings,
                    config=optimization,
                    device=resolved_device,
                    precision=precision,
                )
            metric_sums["vicreg"] = float(vicreg.cpu())
            metric_sums["vicreg_weighted"] = (
                optimization.vicreg_total_weight * metric_sums["vicreg"]
            )
            for name in ("invariance", "variance", "covariance"):
                metric_sums[f"vicreg_{name}"] = float(
                    vicreg_diagnostics[name].cpu()
                )
            if cadence_vicreg:
                for name in (
                    "common_visible_tokens",
                    "statistical_cadence_positions",
                    "statistical_degrees_of_freedom",
                ):
                    metric_sums[f"vicreg_{name}"] = float(
                        vicreg_diagnostics[name].cpu()
                    )
            metric_sums["total"] = (
                metric_sums["reconstruction"]
                + metric_sums["vicreg_weighted"]
            )
        metric_sums["embedding_projection_gradient_norm"] = (
            _embedding_projection_gradient_norm(model)
        )
        optimizer.step()
        scheduler.step()
        global_step += 1
        metrics = {
            "step": float(global_step),
            "learning_rate": float(optimizer.param_groups[0]["lr"]),
            "window_draws_seen": float(
                global_step * optimization.effective_batch_windows
            ),
            "masked_views_seen": float(
                global_step
                * optimization.effective_batch_windows
                * (2 if use_vicreg else 1)
            ),
            **metric_sums,
        }
        if not all(math.isfinite(value) for value in metrics.values()):
            raise RuntimeError("non-finite FM0 training metric")
        loss_history.append(metrics)
        now = time.monotonic()
        if global_step in milestone_steps:
            save_immutable_milestone_checkpoint(
                output_dir,
                checkpoint_payload(
                    model=model,
                    optimizer=optimizer,
                    scheduler=scheduler,
                    global_step=global_step,
                    run_contract=run_contract,
                    optimization=optimization,
                    dataset=dataset,
                    loss_history=loss_history,
                    objective_state=invocation_objective,
                ),
            )
        if global_step == 1 or global_step % progress_interval_steps == 0:
            print(
                "[fm0-train] "
                f"step={global_step}/{target_step} "
                f"loss={metrics['total']:.8g} "
                f"lr={metrics['learning_rate']:.8g} "
                f"elapsed_s={now-started:.1f}",
                flush=True,
            )
        if (
            synthetic_only
            or checkpoint_interval_seconds == 0
            or now - last_checkpoint_at >= checkpoint_interval_seconds
        ):
            payload = checkpoint_payload(
                model=model,
                optimizer=optimizer,
                scheduler=scheduler,
                global_step=global_step,
                run_contract=run_contract,
                optimization=optimization,
                dataset=dataset,
                loss_history=loss_history,
                objective_state=invocation_objective,
            )
            save_rotating_checkpoint(output_dir, payload)
            last_checkpoint_at = now

    final_payload = checkpoint_payload(
        model=model,
        optimizer=optimizer,
        scheduler=scheduler,
        global_step=global_step,
        run_contract=run_contract,
        optimization=optimization,
        dataset=dataset,
        loss_history=loss_history,
        objective_state=invocation_objective,
    )
    final_checkpoint = Path(output_dir) / "checkpoint.pt"
    _atomic_torch_save(final_payload, final_checkpoint)
    milestone_artifacts = _immutable_milestone_artifacts(
        output_dir,
        milestone_steps,
        global_step=global_step,
    )
    elapsed = time.monotonic() - started
    return {
        "global_step": global_step,
        "loss_history": loss_history,
        "final_metrics": loss_history[-1] if loss_history else {},
        "checkpoint": str(final_checkpoint.resolve()),
        "immutable_milestone_checkpoints": milestone_artifacts,
        "elapsed_seconds_this_invocation": elapsed,
        "precision": precision,
        "device": str(resolved_device),
        "synthetic_only": bool(synthetic_only),
    }


def run_synthetic_training(
    *,
    model: Any,
    dataset: SyntheticFM0Dataset,
    output_dir: str | Path,
    run_contract: Mapping[str, Any],
    optimization: FM0OptimizationConfig,
    target_step: int,
    micro_batch_windows: int,
    device: str,
    precision: str,
    use_vicreg: bool,
    reconstruct_second_view: bool | None = None,
    objective_identity: str | None = None,
    resume_checkpoint: str | Path | None = None,
    expected_resume_step: int | None = None,
) -> dict[str, Any]:
    """Run an explicitly synthetic numerical smoke."""

    resolved_reconstruct_second = (
        (
            False
            if objective_identity == CADENCE_VICREG_OBJECTIVE_IDENTITY
            else bool(use_vicreg)
        )
        if reconstruct_second_view is None
        else bool(reconstruct_second_view)
    )
    return _run_training(
        model=model,
        dataset=dataset,
        output_dir=output_dir,
        run_contract=run_contract,
        optimization=optimization,
        target_step=target_step,
        micro_batch_windows=micro_batch_windows,
        device=device,
        precision=precision,
        use_vicreg=use_vicreg,
        reconstruct_second_view=resolved_reconstruct_second,
        objective_identity=objective_identity,
        synthetic_only=True,
        checkpoint_interval_seconds=0,
        progress_interval_steps=1,
        resume_checkpoint=resume_checkpoint,
        expected_resume_step=expected_resume_step,
    )


def run_real_training(
    *,
    model: Any,
    dataset: FM0ReleaseDataset,
    output_dir: str | Path,
    run_contract: Mapping[str, Any],
    optimization: FM0OptimizationConfig,
    target_step: int,
    micro_batch_windows: int,
    device: str,
    precision: str,
    use_vicreg: bool = False,
    reconstruct_second_view: bool | None = None,
    objective_identity: str | None = None,
    resume_checkpoint: str | Path | None = None,
    expected_resume_step: int | None = None,
) -> dict[str, Any]:
    """Train one declared FM0 objective on a checksum-bound input release."""
    resolved_reconstruct_second = (
        (
            False
            if objective_identity == CADENCE_VICREG_OBJECTIVE_IDENTITY
            else bool(use_vicreg)
        )
        if reconstruct_second_view is None
        else bool(reconstruct_second_view)
    )
    return _run_training(
        model=model,
        dataset=dataset,
        output_dir=output_dir,
        run_contract=run_contract,
        optimization=optimization,
        target_step=target_step,
        micro_batch_windows=micro_batch_windows,
        device=device,
        precision=precision,
        use_vicreg=use_vicreg,
        reconstruct_second_view=resolved_reconstruct_second,
        objective_identity=objective_identity,
        synthetic_only=False,
        checkpoint_interval_seconds=1800,
        progress_interval_steps=10,
        resume_checkpoint=resume_checkpoint,
        expected_resume_step=expected_resume_step,
    )
