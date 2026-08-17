"""Deterministic window and masking helpers for the TWIRL-FM0.1 smoke.

This module deliberately contains no survey input builder.  The authoritative
FM input release is built and validated separately; this layer accepts only its
model-facing tensor allowlist, or generates explicitly synthetic windows for a
numerical smoke.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import hashlib
import math
from typing import Any

import numpy as np


WINDOW_LENGTH = 2048
EVALUATION_STRIDE = 1024
MASK_TARGET_FRACTION = 0.15
MASK_SPAN_RANGE = (1, 64)

VIEW_NAMES = (
    "raw_relative_1x1",
    "raw_relative_3x3",
    "adp_1x1",
    "adp_3x3",
    "adp015_1x1",
    "adp015_3x3",
)
VARIANT_VIEW_INDICES = {
    "TWIRL-FM0.1.1": (2, 3),
    "TWIRL-FM0.1.2": (2, 3),
    "TWIRL-FM0.1.3": (2, 3, 4, 5),
    "TWIRL-FM0.1.4": (0, 1, 2, 3, 4, 5),
    "TWIRL-FM0.1.5": (0, 1, 2, 3, 4, 5),
}

MODEL_WINDOW_KEYS = (
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


def _stable_seed(*parts: object) -> int:
    digest = hashlib.sha256(
        "\x1f".join(str(part) for part in parts).encode("utf-8")
    ).digest()
    return int.from_bytes(digest[:8], "big", signed=False)


def variant_view_indices(variant: str) -> tuple[int, ...]:
    """Return the frozen six-view indices exposed by one FM0.1 variant."""

    try:
        return VARIANT_VIEW_INDICES[str(variant)]
    except KeyError as exc:
        raise ValueError(f"unknown TWIRL-FM0.1 variant: {variant!r}") from exc


def deterministic_window_start(
    segment_length: int,
    *,
    sample_index: int,
    epoch: int,
    seed: int,
    window_length: int = WINDOW_LENGTH,
) -> int:
    """Choose a reproducible training-window start within one segment."""

    if segment_length < 0:
        raise ValueError("segment_length must be nonnegative")
    if window_length <= 0:
        raise ValueError("window_length must be positive")
    maximum = max(0, int(segment_length) - int(window_length))
    if maximum == 0:
        return 0
    rng = np.random.default_rng(
        _stable_seed("fm0-window", seed, epoch, sample_index)
    )
    return int(rng.integers(0, maximum + 1))


def evaluation_window_starts(
    segment_length: int,
    *,
    window_length: int = WINDOW_LENGTH,
    stride: int = EVALUATION_STRIDE,
) -> tuple[int, ...]:
    """Tile one segment without crossing its boundary, including the last edge."""

    if segment_length < 0:
        raise ValueError("segment_length must be nonnegative")
    if window_length <= 0 or stride <= 0:
        raise ValueError("window_length and stride must be positive")
    if segment_length == 0:
        return ()
    # Evaluation is a strict stride lattice.  The final start may yield a
    # padded window; do not shift it backward to manufacture an exact edge.
    return tuple(range(0, segment_length, stride))


def synchronized_temporal_mask(
    flux_valid: np.ndarray,
    time_valid: np.ndarray,
    *,
    seed: int,
    target_fraction: float = MASK_TARGET_FRACTION,
    span_range: tuple[int, int] = MASK_SPAN_RANGE,
) -> tuple[np.ndarray, np.ndarray]:
    """Make one temporal mask shared by every flux view.

    Returns ``(temporal_mask[L], reconstruction_mask[V,L])``.  The latter is
    intersected with per-view validity, so an invalid cadence is never a target.
    """

    valid = np.asarray(flux_valid, dtype=bool)
    time_ok = np.asarray(time_valid, dtype=bool)
    if valid.ndim != 2:
        raise ValueError("flux_valid must have shape [view, cadence]")
    if time_ok.shape != (valid.shape[1],):
        raise ValueError("time_valid must have shape [cadence]")
    if not 0.0 <= float(target_fraction) <= 1.0:
        raise ValueError("target_fraction must be in [0, 1]")
    minimum_span, maximum_span = (int(span_range[0]), int(span_range[1]))
    if minimum_span <= 0 or maximum_span < minimum_span:
        raise ValueError("span_range must contain positive ordered bounds")

    eligible = time_ok & np.any(valid, axis=0)
    eligible_count = int(np.count_nonzero(eligible))
    temporal = np.zeros(valid.shape[1], dtype=bool)
    if eligible_count == 0 or target_fraction == 0:
        return temporal, valid & temporal[None, :]

    target_count = max(1, int(math.ceil(target_fraction * eligible_count)))
    eligible_indices = np.flatnonzero(eligible)
    rng = np.random.default_rng(int(seed))
    attempts = 0
    maximum_attempts = max(64, 8 * target_count)
    while int(np.count_nonzero(temporal & eligible)) < target_count:
        anchor = int(rng.choice(eligible_indices))
        span = int(rng.integers(minimum_span, maximum_span + 1))
        left = int(rng.integers(0, span))
        start = max(0, anchor - left)
        stop = min(valid.shape[1], start + span)
        temporal[start:stop] = True
        attempts += 1
        if attempts >= maximum_attempts:
            remaining = eligible_indices[~temporal[eligible_indices]]
            need = target_count - int(np.count_nonzero(temporal & eligible))
            temporal[rng.choice(remaining, size=need, replace=False)] = True
            break

    temporal &= eligible
    return temporal, valid & temporal[None, :]


def _as_cadence_view(
    value: Any,
    *,
    name: str,
    expected_views: int,
) -> np.ndarray:
    array = np.asarray(value)
    if array.ndim != 2:
        raise ValueError(f"{name} must be two-dimensional")
    if array.shape[1] == expected_views:
        return array
    if array.shape[0] == expected_views:
        return array.T
    raise ValueError(
        f"{name} must have one axis of length {expected_views}; got {array.shape}"
    )


def prepare_model_window(
    window: Mapping[str, Any],
    *,
    variant: str,
    mask_seed: int,
    window_length: int = WINDOW_LENGTH,
) -> dict[str, np.ndarray]:
    """Validate, orient, select views, pad, and mask one release window.

    The input follows the model-tensor allowlist emitted by
    :mod:`twirl.models.fm0.input_release`.  Identifier, detector, BLS, label,
    and injection-truth fields are never copied into the returned sample.
    """

    missing = [name for name in MODEL_WINDOW_KEYS if name not in window]
    if missing:
        raise ValueError(f"model window lacks required fields: {missing}")

    indices = variant_view_indices(variant)
    flux_cv = _as_cadence_view(
        window["flux"], name="flux", expected_views=len(VIEW_NAMES)
    ).astype(np.float32, copy=False)
    valid_cv = _as_cadence_view(
        window["flux_valid"],
        name="flux_valid",
        expected_views=len(VIEW_NAMES),
    ).astype(bool, copy=False)
    error_cv = _as_cadence_view(
        window["flux_error"], name="flux_error", expected_views=2
    ).astype(np.float32, copy=False)
    error_valid_cv = _as_cadence_view(
        window["error_valid"], name="error_valid", expected_views=2
    ).astype(bool, copy=False)
    n_cadences = flux_cv.shape[0]
    if valid_cv.shape != flux_cv.shape:
        raise ValueError("flux and flux_valid shapes differ")
    if error_cv.shape != (n_cadences, 2) or error_valid_cv.shape != error_cv.shape:
        raise ValueError("error arrays must have shape [cadence, 2]")
    if n_cadences > window_length:
        raise ValueError("prepare_model_window accepts one already bounded window")

    one_dimensional: dict[str, np.ndarray] = {}
    for name, dtype in (
        ("local_time_cadences", np.float32),
        ("delta_time_cadences", np.float32),
        ("time_valid", bool),
        ("segment_boundary", bool),
    ):
        value = np.asarray(window[name], dtype=dtype)
        if value.shape != (n_cadences,):
            raise ValueError(f"{name} must have shape [cadence]")
        one_dimensional[name] = value

    present = np.asarray(window["view_present"], dtype=bool)
    if present.shape != (len(VIEW_NAMES),):
        raise ValueError("view_present must have shape [6]")
    selected_present = present[np.asarray(indices, dtype=int)]
    if not np.all(selected_present):
        missing_views = [
            VIEW_NAMES[index]
            for index in indices
            if not bool(present[index])
        ]
        raise ValueError(
            f"variant {variant} requires missing flux views: {missing_views}"
        )
    selected_has_valid_cadence = np.any(valid_cv[:, indices], axis=0)
    if not np.all(selected_has_valid_cadence):
        missing_views = [
            VIEW_NAMES[index]
            for index, has_valid in zip(indices, selected_has_valid_cadence)
            if not bool(has_valid)
        ]
        raise ValueError(
            f"variant {variant} has no valid cadences for required flux views: "
            f"{missing_views}"
        )

    flux = np.zeros((len(indices), window_length), dtype=np.float32)
    flux_valid = np.zeros((len(indices), window_length), dtype=bool)
    flux[:, :n_cadences] = flux_cv[:, indices].T
    flux_valid[:, :n_cadences] = valid_cv[:, indices].T
    flux_valid &= np.isfinite(flux)
    flux[~flux_valid] = 0.0

    errors = np.zeros((2, window_length), dtype=np.float32)
    errors_valid = np.zeros((2, window_length), dtype=bool)
    errors[:, :n_cadences] = error_cv.T
    errors_valid[:, :n_cadences] = error_valid_cv.T
    errors_valid &= np.isfinite(errors)
    errors[~errors_valid] = 0.0

    sample: dict[str, np.ndarray] = {
        "flux": flux,
        "flux_valid": flux_valid,
        "flux_error": errors,
        "error_valid": errors_valid,
        "local_time_cadences": np.zeros(window_length, dtype=np.float32),
        "delta_time_cadences": np.zeros(window_length, dtype=np.float32),
        "time_valid": np.zeros(window_length, dtype=bool),
        "segment_boundary": np.zeros(window_length, dtype=bool),
        "view_present": selected_present,
    }
    for name in (
        "local_time_cadences",
        "delta_time_cadences",
        "time_valid",
        "segment_boundary",
    ):
        sample[name][:n_cadences] = one_dimensional[name]

    sample["flux_valid"] &= sample["view_present"][:, None]
    sample["flux"][~sample["flux_valid"]] = 0.0
    temporal, reconstruction = synchronized_temporal_mask(
        sample["flux_valid"],
        sample["time_valid"],
        seed=mask_seed,
    )
    sample["temporal_mask"] = temporal
    sample["reconstruction_mask"] = reconstruction
    return sample


@dataclass(frozen=True)
class SyntheticFM0Config:
    """Explicitly non-scientific synthetic data for numerical smokes."""

    variant: str
    seed: int = 560067
    n_sources: int = 64
    visits_per_source: int = 2
    windows_per_epoch: int = 1024
    window_length: int = WINDOW_LENGTH
    noise_scale: float = 0.002
    event_depth: float = 0.01

    def __post_init__(self) -> None:
        variant_view_indices(self.variant)
        if self.n_sources <= 0 or self.visits_per_source <= 0:
            raise ValueError("synthetic source and visit counts must be positive")
        if self.windows_per_epoch <= 0 or self.window_length <= 0:
            raise ValueError("synthetic window counts and length must be positive")


class SyntheticFM0Dataset:
    """Hash-indexed synthetic windows with source-first sampling semantics."""

    def __init__(self, config: SyntheticFM0Config) -> None:
        self.config = config
        self.epoch = 0

    def __len__(self) -> int:
        return self.config.windows_per_epoch

    def set_epoch(self, epoch: int) -> None:
        if epoch < 0:
            raise ValueError("epoch must be nonnegative")
        self.epoch = int(epoch)

    def _source_visit(self, index: int) -> tuple[int, int]:
        rng = np.random.default_rng(
            _stable_seed("fm0-host-first", self.config.seed, self.epoch, index)
        )
        source = int(rng.integers(0, self.config.n_sources))
        visit = int(rng.integers(0, self.config.visits_per_source))
        return source, visit

    def sample(self, index: int, *, mask_view: int = 0) -> dict[str, np.ndarray]:
        if index < 0:
            index += len(self)
        if not 0 <= index < len(self):
            raise IndexError(index)
        source, visit = self._source_visit(index)
        cadence_count = self.config.window_length
        rng = np.random.default_rng(
            _stable_seed(
                "fm0-synthetic-signal",
                self.config.seed,
                self.epoch,
                source,
                visit,
                index,
            )
        )
        time = np.arange(cadence_count, dtype=np.float32)
        delta = np.ones(cadence_count, dtype=np.float32)
        delta[0] = 0.0
        common = rng.normal(0.0, self.config.noise_scale, cadence_count)
        event_center = int(rng.integers(128, max(129, cadence_count - 128)))
        event_width = int(rng.integers(2, 13))
        event = np.zeros(cadence_count, dtype=np.float32)
        event[event_center : min(cadence_count, event_center + event_width)] = (
            -self.config.event_depth
        )
        flux = np.empty((cadence_count, len(VIEW_NAMES)), dtype=np.float32)
        for view in range(len(VIEW_NAMES)):
            view_noise = rng.normal(
                0.0,
                self.config.noise_scale * (1.0 + 0.08 * view),
                cadence_count,
            )
            flux[:, view] = common + view_noise + event * (1.0 + 0.02 * view)

        time_valid = np.ones(cadence_count, dtype=bool)
        flux_valid = rng.random(flux.shape) >= 0.02
        flux_error = np.full((cadence_count, 2), self.config.noise_scale, np.float32)
        error_valid = rng.random(flux_error.shape) >= 0.01
        segment_boundary = np.zeros(cadence_count, dtype=bool)
        segment_boundary[0] = True
        release_window = {
            "flux": flux,
            "flux_valid": flux_valid,
            "flux_error": flux_error,
            "error_valid": error_valid,
            "local_time_cadences": time,
            "delta_time_cadences": delta,
            "time_valid": time_valid,
            "segment_boundary": segment_boundary,
            "view_present": np.ones(len(VIEW_NAMES), dtype=bool),
        }
        return prepare_model_window(
            release_window,
            variant=self.config.variant,
            mask_seed=_stable_seed(
                "fm0-mask",
                self.config.seed,
                self.epoch,
                index,
                mask_view,
            ),
            window_length=self.config.window_length,
        )

    def __getitem__(self, index: int) -> dict[str, np.ndarray]:
        return self.sample(index, mask_view=0)


def collate_fm0_samples(samples: Sequence[Mapping[str, np.ndarray]]) -> dict[str, Any]:
    """Stack samples into Torch tensors without importing Torch at module load."""

    if not samples:
        raise ValueError("cannot collate an empty FM0 batch")
    try:
        import torch
    except ImportError as exc:  # pragma: no cover - exercised without Torch
        raise RuntimeError("PyTorch is required to collate FM0 samples") from exc

    keys = tuple(samples[0])
    if any(tuple(sample) != keys for sample in samples[1:]):
        raise ValueError("all FM0 samples must expose the same ordered fields")
    batch: dict[str, Any] = {}
    for key in keys:
        array = np.stack([np.asarray(sample[key]) for sample in samples], axis=0)
        batch[key] = torch.from_numpy(array)
    return batch


def move_batch_to_device(batch: Mapping[str, Any], device: Any) -> dict[str, Any]:
    """Move only tensor values, preserving the exact model-facing key set."""

    try:
        import torch
    except ImportError as exc:  # pragma: no cover - exercised without Torch
        raise RuntimeError("PyTorch is required for FM0 training") from exc
    return {
        key: value.to(device, non_blocking=True) if torch.is_tensor(value) else value
        for key, value in batch.items()
    }
