"""Cadence-preserving noncollapse losses for TWIRL-FM0.3.

This module intentionally operates on ``[batch, cadence, feature]`` tensors.
It never constructs a window embedding or reduces a representation across the
cadence axis.  Variance and covariance pool only scalar sufficient statistics
after every feature has been batch-centered at each fixed cadence position.
"""
from __future__ import annotations

import math
from contextlib import nullcontext
from typing import Any

from .model import require_torch

try:  # Optional platform dependency.
    import torch
    from torch.nn import functional as F
except ImportError:  # pragma: no cover - local base environment omits Torch
    torch = None  # type: ignore[assignment]
    F = None  # type: ignore[assignment]


def _validate_weight(value: float, *, name: str) -> float:
    weight = float(value)
    if not math.isfinite(weight) or weight < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return weight


def _fp32_context(device: Any) -> Any:
    if device.type in {"cpu", "cuda"}:
        return torch.autocast(device_type=device.type, enabled=False)
    return nullcontext()


def _pooled_covariance_penalty(
    centered: Any,
    statistical_df: Any,
) -> Any:
    """Return off-diagonal penalty from position-centered sufficient stats."""

    covariance = torch.einsum("bld,ble->de", centered, centered)
    covariance = covariance / statistical_df.to(dtype=centered.dtype)
    squared = covariance.square()
    off_diagonal = squared.sum() - squared.diagonal().sum()
    return off_diagonal / float(centered.shape[-1])


def position_centered_cadence_vicreg_loss(
    first: Any,
    second: Any,
    first_visible: Any,
    second_visible: Any,
    *,
    invariance_weight: float = 25.0,
    variance_weight: float = 25.0,
    covariance_weight: float = 1.0,
    variance_target: float = 1.0,
    variance_epsilon: float = 1.0e-4,
) -> tuple[Any, dict[str, Any]]:
    """Compute position-centered cadence VICReg without a window embedding.

    ``first`` and ``second`` must contain one token per input cadence with
    shape ``[batch, cadence, feature]``.  Each visibility mask has shape
    ``[batch, cadence]`` and should already encode natural validity and whether
    the token was left visible by its view's independent reconstruction mask.
    Only the intersection of the two masks enters any term.

    Batch means are computed separately for each fixed cadence index.  The
    centered residuals then contribute pooled variance and covariance
    sufficient statistics with
    ``df = sum_t max(common_count[t] - 1, 0)``.  Positions with fewer than two
    common-visible samples are ignored by those statistics.  The objective
    fails closed when there are no common tokens or when ``df`` is smaller than
    the embedding dimension.  All loss arithmetic is promoted to FP32 for
    stable BF16/FP16 training.  The cadence axis is never pooled into a feature
    representation.
    """

    require_torch()
    if first.shape != second.shape or first.ndim != 3:
        raise ValueError(
            "cadence embeddings must share shape [batch, cadence, feature]"
        )
    if first.shape[-1] <= 0:
        raise ValueError("cadence embeddings must contain at least one feature")
    expected_mask_shape = first.shape[:2]
    if first_visible.shape != expected_mask_shape or second_visible.shape != (
        expected_mask_shape
    ):
        raise ValueError("visibility masks must share shape [batch, cadence]")
    if first_visible.dtype != torch.bool or second_visible.dtype != torch.bool:
        raise ValueError("visibility masks must be boolean")
    tensors = (first, second, first_visible, second_visible)
    if len({value.device for value in tensors}) != 1:
        raise ValueError("cadence embeddings and visibility masks must share a device")
    if not first.is_floating_point() or not second.is_floating_point():
        raise ValueError("cadence embeddings must be floating point")

    invariance_weight = _validate_weight(
        invariance_weight, name="invariance_weight"
    )
    variance_weight = _validate_weight(variance_weight, name="variance_weight")
    covariance_weight = _validate_weight(
        covariance_weight, name="covariance_weight"
    )
    variance_target = float(variance_target)
    variance_epsilon = float(variance_epsilon)
    if not math.isfinite(variance_target) or variance_target <= 0.0:
        raise ValueError("variance_target must be finite and positive")
    if not math.isfinite(variance_epsilon) or variance_epsilon <= 0.0:
        raise ValueError("variance_epsilon must be finite and positive")

    common_visible = first_visible & second_visible
    position_counts = common_visible.sum(dim=0)
    common_visible_tokens = position_counts.sum()
    if int(common_visible_tokens) == 0:
        raise ValueError("cadence VICReg has no common-visible tokens")
    statistical_df = (position_counts - 1).clamp_min(0).sum()
    if int(statistical_df) < first.shape[-1]:
        raise ValueError(
            "cadence VICReg has insufficient degrees of freedom for its embedding"
        )

    with _fp32_context(first.device):
        first_fp32 = first.float()
        second_fp32 = second.float()
        expanded_visible = common_visible[:, :, None]
        eligible_finite = (
            (torch.isfinite(first_fp32) & torch.isfinite(second_fp32))
            | ~expanded_visible
        )
        if not bool(torch.all(eligible_finite)):
            raise ValueError(
                "cadence VICReg received a nonfinite common-visible token"
            )

        zeros = torch.zeros((), dtype=torch.float32, device=first.device)
        safe_first = torch.where(expanded_visible, first_fp32, zeros)
        safe_second = torch.where(expanded_visible, second_fp32, zeros)
        count_fp32 = position_counts.to(dtype=torch.float32).clamp_min(1.0)

        first_mean = safe_first.sum(dim=0) / count_fp32[:, None]
        second_mean = safe_second.sum(dim=0) / count_fp32[:, None]
        statistical_positions = position_counts >= 2
        statistical_visible = expanded_visible & statistical_positions[
            None, :, None
        ]
        centered_first = torch.where(
            statistical_visible,
            first_fp32 - first_mean[None, :, :],
            zeros,
        )
        centered_second = torch.where(
            statistical_visible,
            second_fp32 - second_mean[None, :, :],
            zeros,
        )

        difference = torch.where(
            expanded_visible, first_fp32 - second_fp32, zeros
        )
        invariance = difference.square().sum() / (
            common_visible_tokens.to(dtype=torch.float32)
            * float(first.shape[-1])
        )

        first_std = torch.sqrt(
            centered_first.square().sum(dim=(0, 1))
            / statistical_df.to(dtype=torch.float32)
            + variance_epsilon
        )
        second_std = torch.sqrt(
            centered_second.square().sum(dim=(0, 1))
            / statistical_df.to(dtype=torch.float32)
            + variance_epsilon
        )
        variance = 0.5 * (
            F.relu(variance_target - first_std).mean()
            + F.relu(variance_target - second_std).mean()
        )

        covariance = _pooled_covariance_penalty(
            centered_first, statistical_df
        ) + _pooled_covariance_penalty(centered_second, statistical_df)
        total = (
            invariance_weight * invariance
            + variance_weight * variance
            + covariance_weight * covariance
        )
    return total, {
        "loss": total,
        "invariance": invariance,
        "variance": variance,
        "covariance": covariance,
        "common_visible_tokens": common_visible_tokens,
        "statistical_cadence_positions": (position_counts >= 2).sum(),
        "statistical_degrees_of_freedom": statistical_df,
        "position_batch_counts": position_counts,
    }


__all__ = ["position_centered_cadence_vicreg_loss"]
