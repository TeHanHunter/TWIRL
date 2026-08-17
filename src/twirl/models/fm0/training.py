"""Losses and deterministic synthetic training for TWIRL-FM0.1."""
from __future__ import annotations

from contextlib import nullcontext
from dataclasses import asdict, dataclass
import math
import os
from pathlib import Path
import random
import time
from typing import Any, Mapping

import numpy as np

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


def _dataset_contract(dataset: Any) -> dict[str, Any]:
    if isinstance(dataset, SyntheticFM0Dataset):
        return {"kind": "synthetic", "config": asdict(dataset.config)}
    if isinstance(dataset, FM0ReleaseDataset):
        return dict(dataset.contract)
    raise TypeError("unsupported FM0 training dataset")


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
    reconstruction = 0.5 * (reconstruction + second_reconstruction)
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
    total = reconstruction + config.vicreg_total_weight * vicreg
    diagnostics.update(
        {
            "reconstruction": reconstruction,
            "reconstruction_second": second_reconstruction,
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
            asdict(dataset.config) if isinstance(dataset, SyntheticFM0Dataset) else None
        ),
        "run_contract": dict(run_contract),
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


def load_checkpoint(
    path: str | Path,
    *,
    model: Any,
    optimizer: Any,
    scheduler: Any,
    expected_run_contract: Mapping[str, Any],
    expected_optimization: FM0OptimizationConfig,
    dataset: Any,
) -> tuple[int, list[dict[str, float]]]:
    """Load and validate a full FM0 checkpoint, including RNG and sampler state."""

    require_torch()
    payload = torch.load(Path(path), map_location="cpu", weights_only=False)
    if not isinstance(payload, Mapping):
        raise ValueError("FM0 checkpoint is not a mapping")
    if payload.get("schema_version") != CHECKPOINT_SCHEMA_VERSION:
        raise ValueError("FM0 checkpoint schema mismatch")
    if payload.get("run_contract") != dict(expected_run_contract):
        raise ValueError("FM0 checkpoint run contract mismatch")
    if payload.get("optimization_config") != asdict(expected_optimization):
        raise ValueError("FM0 checkpoint optimization contract mismatch")
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
    dataset.set_epoch(int(progress.get("epoch", 0)))
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
            with _autocast_context(device, precision):
                first_z = model(record["first_batch"])["z_window"]
            _restore_forward_rng_state(record["second_rng"], device)
            with _autocast_context(device, precision):
                second_z = model(record["second_batch"])["z_window"]
            scale = float(config.vicreg_total_weight)
            torch.autograd.backward(
                (first_z, second_z),
                (
                    first_gradient[offset : offset + count].to(first_z.dtype)
                    * scale,
                    second_gradient[offset : offset + count].to(second_z.dtype)
                    * scale,
                ),
            )
            offset += count
    finally:
        _restore_forward_rng_state(post_forward_rng, device)
    if offset != config.effective_batch_windows:
        raise ValueError("VICReg replay did not cover the effective batch")
    return vicreg.detach(), {
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
    synthetic_only: bool,
    checkpoint_interval_seconds: float,
    progress_interval_steps: int,
    resume_checkpoint: str | Path | None = None,
) -> dict[str, Any]:
    """Run deterministic FM0 training and write full rotating checkpoints."""

    require_torch()
    if target_step <= 0 or target_step > optimization.max_optimizer_steps:
        raise ValueError("target_step is outside the frozen optimizer budget")
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
    if resume_checkpoint is not None:
        global_step, loss_history = load_checkpoint(
            resume_checkpoint,
            model=model,
            optimizer=optimizer,
            scheduler=scheduler,
            expected_run_contract=run_contract,
            expected_optimization=optimization,
            dataset=dataset,
        )
    if global_step > target_step:
        raise ValueError("resume checkpoint is beyond requested target step")

    accumulation = optimization.effective_batch_windows // micro_batch_windows
    started = time.monotonic()
    last_checkpoint_at = started
    model.train()
    while global_step < target_step:
        optimizer.zero_grad(set_to_none=True)
        metric_sums: dict[str, float] = {}
        replay_records: list[dict[str, Any]] = []
        first_embeddings: list[Any] = []
        second_embeddings: list[Any] = []
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
                first_rng = _capture_forward_rng_state(resolved_device)
                with _autocast_context(resolved_device, precision):
                    first_output = model(first_batch)
                    first_reconstruction, _ = masked_huber_reconstruction_loss(
                        first_output["reconstruction"],
                        first_batch["flux"].float(),
                        first_batch["reconstruction_mask"],
                        delta=optimization.huber_delta,
                    )
                    second_rng = _capture_forward_rng_state(resolved_device)
                    second_output = model(second_batch)
                    second_reconstruction, _ = masked_huber_reconstruction_loss(
                        second_output["reconstruction"],
                        second_batch["flux"].float(),
                        second_batch["reconstruction_mask"],
                        delta=optimization.huber_delta,
                    )
                    reconstruction = 0.5 * (
                        first_reconstruction + second_reconstruction
                    )
                (reconstruction / accumulation).backward()
                first_embeddings.append(first_output["z_window"].detach())
                second_embeddings.append(second_output["z_window"].detach())
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
            vicreg, _ = _backprop_effective_batch_vicreg(
                model=model,
                replay_records=replay_records,
                first_embeddings=first_embeddings,
                second_embeddings=second_embeddings,
                config=optimization,
                device=resolved_device,
                precision=precision,
            )
            metric_sums["vicreg"] = float(vicreg.cpu())
            metric_sums["total"] = (
                metric_sums["reconstruction"]
                + optimization.vicreg_total_weight * metric_sums["vicreg"]
            )
        optimizer.step()
        scheduler.step()
        global_step += 1
        metrics = {
            "step": float(global_step),
            "learning_rate": float(optimizer.param_groups[0]["lr"]),
            **metric_sums,
        }
        if not all(math.isfinite(value) for value in metrics.values()):
            raise RuntimeError("non-finite FM0 training metric")
        loss_history.append(metrics)
        now = time.monotonic()
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
    )
    final_checkpoint = Path(output_dir) / "checkpoint.pt"
    _atomic_torch_save(final_payload, final_checkpoint)
    elapsed = time.monotonic() - started
    return {
        "global_step": global_step,
        "loss_history": loss_history,
        "final_metrics": loss_history[-1] if loss_history else {},
        "checkpoint": str(final_checkpoint.resolve()),
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
    resume_checkpoint: str | Path | None = None,
) -> dict[str, Any]:
    """Run an explicitly synthetic numerical smoke."""

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
        synthetic_only=True,
        checkpoint_interval_seconds=0,
        progress_interval_steps=1,
        resume_checkpoint=resume_checkpoint,
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
    resume_checkpoint: str | Path | None = None,
) -> dict[str, Any]:
    """Train FM0.1.1--0.1.4 on a checksum-bound real input release."""

    if dataset.config.variant == "TWIRL-FM0.1.5":
        raise ValueError("FM0.1.5 real training must explicitly enable VICReg")
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
        use_vicreg=False,
        synthetic_only=False,
        checkpoint_interval_seconds=1800,
        progress_interval_steps=10,
        resume_checkpoint=resume_checkpoint,
    )
