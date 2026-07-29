"""Leakage-safe self-supervision helpers for the harmonic teacher.

The primary SSL pool is intentionally real-only and fold-local.  This module
does not inspect morphology labels: it selects development observations by the
frozen TIC split, creates two event-preserving views of the current harmonic
CNN tensors, and provides the VICReg objective used to pretrain the encoder.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.vetting.harmonic_inputs import (
    HARMONIC_FACTORS,
    MODEL_INPUT_NUMERIC_ABS_MAX,
)


HARMONIC_SSL_CONTRACT_VERSION = "twirl_harmonic_ssl_v1"
_CACHE_IDENTITY_SCHEMA = "twirl_harmonic_ssl_cache_identity_v1"
_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_MODEL_PHOTOMETRY_ERROR_CHANNELS = slice(0, 5)
_MODEL_ERROR_CHANNELS = (3, 4)
_MODEL_PHASE_CHANNEL = 5
_MODEL_QUALITY_CHANNEL = 6
_MODEL_PHASE_MIN = -0.5
_MODEL_PHASE_MAX = 0.5
_SSL_SELECTION_COLUMNS = (
    "review_id",
    "tic",
    "sector",
    "fixed_split",
    "cv_fold",
    "is_injected_row",
)


@dataclass(frozen=True)
class EventPreservingAugmentationConfig:
    """Configuration for shape-preserving native-tensor augmentations.

    Cadence dropout applies only to full harmonic folds and never to samples
    inside the protected candidate-event window.  Local folds retain their
    masks exactly.  Independent noise is applied to the two aperture-flux
    channels, scaled by their uncertainty channels, and the
    aperture-difference channel is then recomputed.  Phase, quality, order,
    length, and sign are never transformed.  Metadata is always zeroed so SSL
    cannot solve the view-agreement task through scalar BLS shortcuts.
    """

    harmonic_cadence_dropout_probability: float = 0.02
    harmonic_flux_noise_scale: float = 0.05
    local_flux_noise_scale: float = 0.03
    periodogram_bin_dropout_probability: float = 0.02
    event_protection_duration_multiplier: float = 2.0
    max_event_phase_half_width: float = 0.25
    phase_channel: int = 5
    flux_channels: tuple[int, int, int] = (0, 1, 2)
    error_channels: tuple[int, int] = (3, 4)
    harmonic_factors: tuple[float, ...] = HARMONIC_FACTORS
    metadata_policy: str = "zero"
    contract_version: str = HARMONIC_SSL_CONTRACT_VERSION

    def __post_init__(self) -> None:
        probability = float(self.harmonic_cadence_dropout_probability)
        if not 0.0 <= probability <= 1.0:
            raise ValueError("harmonic cadence dropout probability must be in [0, 1]")
        for name in ("harmonic_flux_noise_scale", "local_flux_noise_scale"):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be finite and nonnegative")
        periodogram_probability = float(self.periodogram_bin_dropout_probability)
        if not 0.0 <= periodogram_probability <= 1.0:
            raise ValueError("periodogram bin dropout probability must be in [0, 1]")
        if (
            not np.isfinite(self.event_protection_duration_multiplier)
            or self.event_protection_duration_multiplier <= 0
        ):
            raise ValueError(
                "event_protection_duration_multiplier must be finite and positive"
            )
        if (
            not np.isfinite(self.max_event_phase_half_width)
            or not 0 < self.max_event_phase_half_width <= 0.5
        ):
            raise ValueError("max_event_phase_half_width must be in (0, 0.5]")
        if self.phase_channel < 0:
            raise ValueError("phase_channel must be nonnegative")
        if len(self.flux_channels) != 3 or len(set(self.flux_channels)) != 3:
            raise ValueError("flux_channels must contain three distinct channel indices")
        if len(self.error_channels) != 2 or len(set(self.error_channels)) != 2:
            raise ValueError("error_channels must contain two distinct channel indices")
        all_channels = self.flux_channels + self.error_channels + (self.phase_channel,)
        if any(int(value) < 0 for value in all_channels):
            raise ValueError("augmentation channel indices must be nonnegative")
        if set(self.flux_channels) & set(self.error_channels):
            raise ValueError("flux and uncertainty channels must be distinct")
        factors = tuple(float(value) for value in self.harmonic_factors)
        if not factors or any(not np.isfinite(value) or value <= 0 for value in factors):
            raise ValueError("harmonic_factors must be finite and positive")
        if self.metadata_policy != "zero":
            raise ValueError("SSL metadata_policy must be exactly 'zero'")
        if self.contract_version != HARMONIC_SSL_CONTRACT_VERSION:
            raise ValueError(
                f"augmentation contract_version must be {HARMONIC_SSL_CONTRACT_VERSION!r}"
            )


@dataclass(frozen=True)
class VICRegConfig:
    """Weights and numerical constants for the VICReg representation loss."""

    invariance_weight: float = 25.0
    variance_weight: float = 25.0
    covariance_weight: float = 1.0
    target_std: float = 1.0
    eps: float = 1.0e-4
    contract_version: str = HARMONIC_SSL_CONTRACT_VERSION

    def __post_init__(self) -> None:
        for name in (
            "invariance_weight",
            "variance_weight",
            "covariance_weight",
        ):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be finite and nonnegative")
        if not np.isfinite(self.target_std) or self.target_std <= 0:
            raise ValueError("target_std must be finite and positive")
        if not np.isfinite(self.eps) or self.eps <= 0:
            raise ValueError("eps must be finite and positive")
        if self.contract_version != HARMONIC_SSL_CONTRACT_VERSION:
            raise ValueError(
                f"VICReg contract_version must be {HARMONIC_SSL_CONTRACT_VERSION!r}"
            )


@dataclass(frozen=True)
class SSLCacheIdentity:
    """Complete hashable identity for one fold-local encoder cache."""

    training_table_sha256: str
    native_registry_sha256: str
    split_registry_sha256: str
    selected_rows_sha256: str
    selected_tics_sha256: str
    model_config_sha256: str
    augmentation_config_sha256: str
    vicreg_config_sha256: str
    profile: str
    held_out_fold: int
    seed: int
    code_revision: str
    n_selected_rows: int
    n_selected_tics: int
    schema_version: str = _CACHE_IDENTITY_SCHEMA

    def __post_init__(self) -> None:
        for name in (
            "training_table_sha256",
            "native_registry_sha256",
            "split_registry_sha256",
            "selected_rows_sha256",
            "selected_tics_sha256",
            "model_config_sha256",
            "augmentation_config_sha256",
            "vicreg_config_sha256",
        ):
            value = str(getattr(self, name))
            if _SHA256_PATTERN.fullmatch(value) is None:
                raise ValueError(f"{name} must be a lowercase 64-character SHA-256")
        if self.schema_version != _CACHE_IDENTITY_SCHEMA:
            raise ValueError(f"schema_version must be {_CACHE_IDENTITY_SCHEMA!r}")
        if not str(self.profile).strip():
            raise ValueError("profile must be nonempty")
        if int(self.held_out_fold) not in range(5):
            raise ValueError("held_out_fold must be in [0, 4]")
        if int(self.seed) < 0:
            raise ValueError("seed must be nonnegative")
        if not str(self.code_revision).strip():
            raise ValueError("code_revision must be nonempty")
        if int(self.n_selected_rows) <= 0:
            raise ValueError("n_selected_rows must be positive")
        if int(self.n_selected_tics) <= 0:
            raise ValueError("n_selected_tics must be positive")
        if int(self.n_selected_tics) > int(self.n_selected_rows):
            raise ValueError("n_selected_tics cannot exceed n_selected_rows")

    def as_dict(self) -> dict[str, Any]:
        """Return the canonical identity fields without the outer digest."""

        return asdict(self)

    def digest(self) -> str:
        """Return the canonical SHA-256 of this identity."""

        return _canonical_sha256(self.as_dict())

    def to_manifest(self) -> dict[str, Any]:
        """Return strict cache-manifest content including its own digest."""

        return {**self.as_dict(), "identity_sha256": self.digest()}


def _canonical_value(value: Any) -> Any:
    if hasattr(value, "__dataclass_fields__"):
        return _canonical_value(asdict(value))
    if isinstance(value, Mapping):
        if any(not isinstance(key, str) for key in value):
            raise TypeError("canonical mappings require string keys")
        return {
            key: _canonical_value(item)
            for key, item in sorted(value.items(), key=lambda pair: pair[0])
        }
    if isinstance(value, (list, tuple)):
        return [_canonical_value(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return _canonical_value(value.item())
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("canonical payloads may not contain nonfinite floats")
        return value
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    raise TypeError(f"unsupported canonical payload value: {type(value).__name__}")


def _canonical_sha256(payload: Any) -> str:
    encoded = json.dumps(
        _canonical_value(payload),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _integer_series(series: pd.Series, *, name: str, positive: bool) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce").to_numpy(dtype=np.float64)
    if not np.isfinite(numeric).all() or not np.equal(numeric, np.rint(numeric)).all():
        raise ValueError(f"{name} must contain finite integers")
    if positive and np.any(numeric <= 0):
        raise ValueError(f"{name} must contain positive integers")
    return pd.Series(np.rint(numeric).astype(np.int64), index=series.index, name=name)


def _strict_bool_series(series: pd.Series, *, name: str) -> pd.Series:
    true_values = {"1", "true", "t", "yes", "y"}
    false_values = {"0", "false", "f", "no", "n"}
    parsed: list[bool] = []
    for value in series.tolist():
        if isinstance(value, (bool, np.bool_)):
            parsed.append(bool(value))
            continue
        if isinstance(value, (int, np.integer)) and int(value) in (0, 1):
            parsed.append(bool(value))
            continue
        if isinstance(value, str):
            text = value.strip().lower()
            if text in true_values:
                parsed.append(True)
                continue
            if text in false_values:
                parsed.append(False)
                continue
        raise ValueError(f"{name} contains an invalid boolean value: {value!r}")
    return pd.Series(parsed, index=series.index, dtype=bool, name=name)


def _selection_hashes(rows: pd.DataFrame) -> tuple[str, str]:
    missing = sorted(set(_SSL_SELECTION_COLUMNS) - set(rows.columns))
    if missing:
        raise KeyError(f"SSL selection rows are missing columns: {missing}")
    if rows.empty:
        raise ValueError("SSL selection is empty")
    review_id = rows["review_id"].fillna("").astype(str).str.strip()
    if review_id.eq("").any() or review_id.duplicated().any():
        raise ValueError("selected SSL review_id values must be nonempty and unique")
    tic = _integer_series(rows["tic"], name="tic", positive=True)
    sector = _integer_series(rows["sector"], name="sector", positive=True)
    fold = _integer_series(rows["cv_fold"], name="cv_fold", positive=False)
    fixed_split = rows["fixed_split"].fillna("").astype(str).str.strip().str.lower()
    injected = _strict_bool_series(rows["is_injected_row"], name="is_injected_row")
    if not fixed_split.eq("development").all():
        raise ValueError("selected SSL rows must all be development observations")
    if injected.any():
        raise ValueError("selected SSL rows must all be real observations")
    records = [
        {
            "review_id": str(review),
            "tic": int(host),
            "sector": int(observed_sector),
            "fixed_split": str(split),
            "cv_fold": int(observed_fold),
        }
        for review, host, observed_sector, split, observed_fold in zip(
            review_id, tic, sector, fixed_split, fold
        )
    ]
    records.sort(
        key=lambda row: (
            row["tic"],
            row["sector"],
            row["cv_fold"],
            row["review_id"],
        )
    )
    unique_tics = sorted({int(value) for value in tic})
    return _canonical_sha256(records), _canonical_sha256(unique_tics)


def select_ssl_fold_rows(
    rows: pd.DataFrame,
    *,
    held_out_fold: int,
    reserved_sectors: Sequence[int] = (63,),
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Select a real-only, TIC-disjoint SSL pool for one development fold.

    The fixed test, the held validation fold, every row belonging to a TIC
    observed in a reserved sector, and every injected row are excluded.
    Reserved-sector exclusion is host-wide so an earlier-sector observation
    cannot leak a future S63 host into pretraining.
    """

    if int(held_out_fold) not in range(5):
        raise ValueError("held_out_fold must be in [0, 4]")
    missing = sorted(
        {"tic", "sector", "fixed_split", "cv_fold", "is_injected_row", "review_id"}
        - set(rows.columns)
    )
    if missing:
        raise KeyError(f"SSL fold selection is missing required columns: {missing}")
    if rows.empty:
        raise ValueError("cannot select SSL rows from an empty table")

    reserved = tuple(sorted({int(value) for value in reserved_sectors}))
    if not reserved or any(value <= 0 for value in reserved):
        raise ValueError("reserved_sectors must contain positive sector integers")

    work = rows.copy()
    work["_ssl_tic"] = _integer_series(work["tic"], name="tic", positive=True)
    work["_ssl_sector"] = _integer_series(work["sector"], name="sector", positive=True)
    work["_ssl_fold"] = _integer_series(work["cv_fold"], name="cv_fold", positive=False)
    work["_ssl_injected"] = _strict_bool_series(
        work["is_injected_row"], name="is_injected_row"
    )
    work["_ssl_split"] = (
        work["fixed_split"].fillna("").astype(str).str.strip().str.lower()
    )
    if not work["_ssl_split"].isin({"development", "test"}).all():
        bad = sorted(
            work.loc[
                ~work["_ssl_split"].isin({"development", "test"}), "_ssl_split"
            ].unique()
        )
        raise ValueError(f"fixed_split contains unsupported values: {bad}")
    development = work["_ssl_split"].eq("development")
    test = work["_ssl_split"].eq("test")
    if not work.loc[development, "_ssl_fold"].isin(range(5)).all():
        raise ValueError("development rows require cv_fold in [0, 4]")
    if not work.loc[test, "_ssl_fold"].eq(-1).all():
        raise ValueError("fixed-test rows require cv_fold=-1")
    if not (development & work["_ssl_fold"].eq(int(held_out_fold))).any():
        raise ValueError(f"held validation fold {held_out_fold} has no development rows")

    assignment_counts = work.groupby("_ssl_tic", sort=False).agg(
        n_fixed_splits=("_ssl_split", "nunique"),
        n_cv_folds=("_ssl_fold", "nunique"),
    )
    inconsistent = assignment_counts[
        assignment_counts["n_fixed_splits"].ne(1)
        | assignment_counts["n_cv_folds"].ne(1)
    ]
    if not inconsistent.empty:
        raise ValueError(
            "TIC split leakage across fixed_split/cv_fold assignments; "
            f"first={inconsistent.index.astype(int).tolist()[:10]}"
        )

    test_tics = set(work.loc[test, "_ssl_tic"].astype(int))
    validation_tics = set(
        work.loc[
            development & work["_ssl_fold"].eq(int(held_out_fold)), "_ssl_tic"
        ].astype(int)
    )
    reserved_tics = set(
        work.loc[work["_ssl_sector"].isin(reserved), "_ssl_tic"].astype(int)
    )
    forbidden_tics = test_tics | validation_tics | reserved_tics
    selected_mask = (
        development
        & work["_ssl_fold"].ne(int(held_out_fold))
        & ~work["_ssl_injected"]
        & ~work["_ssl_tic"].isin(forbidden_tics)
    )
    selected = rows.loc[selected_mask].copy()
    if selected.empty:
        raise ValueError("fold-local real-only SSL selection is empty")

    selected_tics = set(
        _integer_series(selected["tic"], name="tic", positive=True).astype(int)
    )
    disjoint = {
        "fixed_test_tics": not bool(selected_tics & test_tics),
        "held_validation_tics": not bool(selected_tics & validation_tics),
        "reserved_sector_tics": not bool(selected_tics & reserved_tics),
    }
    if not all(disjoint.values()):
        raise RuntimeError(f"SSL TIC-disjointness check failed: {disjoint}")
    selected_rows_sha256, selected_tics_sha256 = _selection_hashes(selected)
    audit = {
        "contract_version": HARMONIC_SSL_CONTRACT_VERSION,
        "held_out_fold": int(held_out_fold),
        "reserved_sectors": list(reserved),
        "primary_ssl_real_only": True,
        "reserved_sector_exclusion_scope": "whole_tic",
        "n_input_rows": int(len(rows)),
        "n_input_tics": int(work["_ssl_tic"].nunique()),
        "n_selected_rows": int(len(selected)),
        "n_selected_tics": int(len(selected_tics)),
        "n_fixed_test_rows_excluded": int(test.sum()),
        "n_held_validation_rows_excluded": int(
            (development & work["_ssl_fold"].eq(int(held_out_fold))).sum()
        ),
        "n_injected_rows_excluded": int(work["_ssl_injected"].sum()),
        "n_reserved_sector_rows": int(work["_ssl_sector"].isin(reserved).sum()),
        "n_reserved_sector_tics": int(len(reserved_tics)),
        "selected_rows_sha256": selected_rows_sha256,
        "selected_tics_sha256": selected_tics_sha256,
        "tic_disjoint": disjoint,
    }
    return selected, audit


def _first_true_index(mask: Any) -> tuple[int, ...]:
    import torch

    indices = torch.nonzero(mask, as_tuple=False)
    if int(indices.shape[0]) == 0:
        raise ValueError("cannot identify an index in an empty failure mask")
    return tuple(int(value) for value in indices[0].detach().cpu().tolist())


def _validate_active_numeric_values(
    *,
    name: str,
    values: Any,
    mask: Any,
    stage: str,
) -> None:
    import torch

    payload_nonfinite = ~torch.isfinite(values)
    if bool(payload_nonfinite.any()):
        raise ValueError(
            f"{stage} {name} contains nonfinite tensor payload; "
            f"first_index={_first_true_index(payload_nonfinite)}"
        )
    nonfinite = mask & ~torch.isfinite(values)
    if bool(nonfinite.any()):
        raise ValueError(
            f"{stage} {name} contains nonfinite model-active values; "
            f"first_index={_first_true_index(nonfinite)}"
        )
    envelope_exceeded = (
        mask
        & torch.isfinite(values)
        & values.abs().gt(float(MODEL_INPUT_NUMERIC_ABS_MAX))
    )
    if bool(envelope_exceeded.any()):
        raise ValueError(
            f"{stage} {name} exceeds the model-active absolute bound "
            f"{MODEL_INPUT_NUMERIC_ABS_MAX:g}; "
            f"first_index={_first_true_index(envelope_exceeded)}"
        )


def _validate_sequence_numeric_domains(
    *,
    name: str,
    values: Any,
    mask: Any,
    stage: str,
) -> None:
    import torch

    if int(values.shape[2]) <= _MODEL_QUALITY_CHANNEL:
        raise ValueError(
            f"{stage} {name} must provide the seven locked harmonic channels"
        )
    _validate_active_numeric_values(
        name=name,
        values=values,
        mask=mask,
        stage=stage,
    )

    phase_mask = mask[:, :, _MODEL_PHASE_CHANNEL, :]
    quality_mask = mask[:, :, _MODEL_QUALITY_CHANNEL, :]
    coordinate_mask_mismatch = phase_mask.ne(quality_mask)
    if bool(coordinate_mask_mismatch.any()):
        raise ValueError(
            f"{stage} {name} phase/quality masks must match; "
            f"first_index={_first_true_index(coordinate_mask_mismatch)}"
        )
    photometry_without_coordinates = (
        mask[:, :, _MODEL_PHOTOMETRY_ERROR_CHANNELS, :].any(dim=2)
        & ~(phase_mask & quality_mask)
    )
    if bool(photometry_without_coordinates.any()):
        raise ValueError(
            f"{stage} {name} active photometry requires phase and quality; "
            f"first_index={_first_true_index(photometry_without_coordinates)}"
        )

    for channel in _MODEL_ERROR_CHANNELS:
        negative = mask[:, :, channel, :] & values[:, :, channel, :].lt(0.0)
        if bool(negative.any()):
            raise ValueError(
                f"{stage} {name} error channels must be nonnegative; "
                f"channel={channel}, first_index={_first_true_index(negative)}"
            )

    phase = values[:, :, _MODEL_PHASE_CHANNEL, :]
    phase_active = mask[:, :, _MODEL_PHASE_CHANNEL, :] & torch.isfinite(phase)
    phase_invalid = phase_active & (
        phase.lt(float(_MODEL_PHASE_MIN))
        | phase.gt(float(_MODEL_PHASE_MAX))
    )
    if bool(phase_invalid.any()):
        raise ValueError(
            f"{stage} {name} phase channel must be in "
            f"[{_MODEL_PHASE_MIN}, {_MODEL_PHASE_MAX}]; "
            f"first_index={_first_true_index(phase_invalid)}"
        )

    quality = values[:, :, _MODEL_QUALITY_CHANNEL, :]
    quality_active = (
        quality_mask & torch.isfinite(quality)
    )
    quality_invalid = quality_active & quality.ne(0.0) & quality.ne(1.0)
    if bool(quality_invalid.any()):
        raise ValueError(
            f"{stage} {name} quality channel must be binary; "
            f"first_index={_first_true_index(quality_invalid)}"
        )

    quality_bad = quality_active & quality.eq(1.0)
    bad_photometry_active = (
        mask[:, :, _MODEL_PHOTOMETRY_ERROR_CHANNELS, :]
        & quality_bad.unsqueeze(2)
    )
    if bool(bad_photometry_active.any()):
        raise ValueError(
            f"{stage} {name} quality-bad cadences must mask photometry/error "
            "channels 0:5; "
            f"first_index={_first_true_index(bad_photometry_active)}"
        )


def _validate_model_input_numeric_contract(
    batch: Mapping[str, Any],
    *,
    stage: str,
) -> None:
    import torch

    for name, mask_name in (
        ("harmonic_values", "harmonic_mask"),
        ("local_values", "local_mask"),
    ):
        _validate_sequence_numeric_domains(
            name=name,
            values=batch[name],
            mask=batch[mask_name],
            stage=stage,
        )
    _validate_active_numeric_values(
        name="periodogram_values",
        values=batch["periodogram_values"],
        mask=batch["periodogram_mask"],
        stage=stage,
    )
    metadata = batch["metadata"]
    _validate_active_numeric_values(
        name="metadata",
        values=metadata,
        mask=torch.ones_like(metadata, dtype=torch.bool),
        stage=stage,
    )


def _validate_batch_tensor_contract(batch: Mapping[str, Any], duration_min: Any) -> tuple[Any, Any]:
    import torch

    required = {
        "period_d",
        "harmonic_values",
        "harmonic_mask",
        "local_values",
        "local_mask",
        "periodogram_values",
        "periodogram_mask",
        "metadata",
    }
    missing = sorted(required - set(batch))
    if missing:
        raise KeyError(f"SSL batch is missing tensors: {missing}")
    for name in required:
        if not torch.is_tensor(batch[name]):
            raise TypeError(f"batch[{name!r}] must be a torch.Tensor")
    harmonic = batch["harmonic_values"]
    harmonic_mask = batch["harmonic_mask"]
    local = batch["local_values"]
    local_mask = batch["local_mask"]
    periodogram = batch["periodogram_values"]
    periodogram_mask = batch["periodogram_mask"]
    metadata = batch["metadata"]
    period = batch["period_d"]
    if harmonic.ndim != 4 or local.ndim != 4:
        raise ValueError("harmonic_values and local_values must have shape (B,V,C,T)")
    if periodogram.ndim != 3:
        raise ValueError("periodogram_values must have shape (B,C,T)")
    if metadata.ndim != 2 or period.ndim != 1:
        raise ValueError("metadata and period_d must have shapes (B,M) and (B,)")
    if harmonic_mask.shape != harmonic.shape or local_mask.shape != local.shape:
        raise ValueError("harmonic/local masks must match their value tensors")
    if periodogram_mask.shape != periodogram.shape:
        raise ValueError("periodogram_mask must match periodogram_values")
    if (
        harmonic_mask.dtype != torch.bool
        or local_mask.dtype != torch.bool
        or periodogram_mask.dtype != torch.bool
    ):
        raise TypeError("SSL sequence masks must have torch.bool dtype")
    for name, value in (
        ("harmonic_values", harmonic),
        ("local_values", local),
        ("periodogram_values", periodogram),
        ("metadata", metadata),
        ("period_d", period),
    ):
        if not value.is_floating_point():
            raise TypeError(f"{name} must be floating point")
    batch_size = int(harmonic.shape[0])
    if batch_size == 0:
        raise ValueError("SSL batch must be nonempty")
    if any(
        int(value.shape[0]) != batch_size
        for value in (local, periodogram, metadata, period)
    ):
        raise ValueError("all SSL batch tensors must have the same batch dimension")
    devices = {
        value.device
        for value in (
            harmonic,
            harmonic_mask,
            local,
            local_mask,
            periodogram,
            periodogram_mask,
            metadata,
            period,
        )
    }
    if len(devices) != 1:
        raise ValueError("all SSL batch tensors must be on the same device")
    duration = torch.as_tensor(
        duration_min,
        dtype=period.dtype,
        device=period.device,
    )
    if duration.ndim == 0:
        duration = duration.repeat(batch_size)
    if duration.shape != period.shape:
        raise ValueError("duration_min must be scalar or have shape (B,)")
    if not torch.isfinite(period).all() or not (period > 0).all():
        raise ValueError("period_d must be finite and positive")
    if not torch.isfinite(duration).all() or not (duration > 0).all():
        raise ValueError("duration_min must be finite and positive")
    _validate_model_input_numeric_contract(batch, stage="pre-augmentation")
    return period, duration


def _generator(device: Any, *, seed: int, view_index: int) -> Any:
    import torch

    if int(seed) < 0:
        raise ValueError("seed must be nonnegative")
    if int(view_index) < 0:
        raise ValueError("view_index must be nonnegative")
    modulus = 2**63 - 1
    effective_seed = (
        int(seed) + 0x1E35A7BD * (int(view_index) + 1)
    ) % modulus
    generator = torch.Generator(device=device)
    generator.manual_seed(effective_seed)
    return generator


def _event_phase_centers(factor: float) -> tuple[float, ...]:
    """Return primary and secondary event aliases for one harmonic view."""

    value = float(factor)
    if not np.isfinite(value) or value <= 0:
        raise ValueError("harmonic factor must be finite and positive")
    if value >= 1.0:
        multiplicity = int(round(value))
        if not np.isclose(value, multiplicity, rtol=0.0, atol=1.0e-10):
            raise ValueError(
                "event protection requires integer or reciprocal-integer "
                f"harmonic factors, observed {value}"
            )
        raw = [
            (event_number + offset) / value
            for event_number in range(multiplicity)
            for offset in (0.0, 0.5)
        ]
    else:
        reciprocal = int(round(1.0 / value))
        if not np.isclose(
            value,
            1.0 / reciprocal,
            rtol=0.0,
            atol=1.0e-10,
        ):
            raise ValueError(
                "event protection requires integer or reciprocal-integer "
                f"harmonic factors, observed {value}"
            )
        raw = [offset / value for offset in (0.0, 0.5)]
    wrapped = {
        round(((center + 0.5) % 1.0) - 0.5, 12)
        for center in raw
    }
    return tuple(sorted(float(center) for center in wrapped))


def _harmonic_event_protection_mask(
    *,
    phase: Any,
    phase_valid: Any,
    period: Any,
    duration: Any,
    config: EventPreservingAugmentationConfig,
) -> Any:
    """Protect every repeated primary/secondary alias in each harmonic fold."""

    import torch

    if phase.ndim != 3 or phase_valid.shape != phase.shape:
        raise ValueError("phase and phase_valid must have shape (B,V,T)")
    if int(phase.shape[1]) != len(config.harmonic_factors):
        raise ValueError(
            "phase view count does not match the configured harmonic factors"
        )
    factors = phase.new_tensor(config.harmonic_factors).view(1, -1, 1)
    event_half_width = (
        float(config.event_protection_duration_multiplier)
        * duration.view(-1, 1, 1)
        / (2.0 * 1440.0 * period.view(-1, 1, 1) * factors)
    ).clamp(max=float(config.max_event_phase_half_width))
    protected = torch.zeros_like(phase_valid)
    for view_index, factor in enumerate(config.harmonic_factors):
        view_phase = phase[:, view_index, :]
        view_valid = phase_valid[:, view_index, :]
        view_width = event_half_width[:, view_index, :]
        for center in _event_phase_centers(float(factor)):
            distance = torch.abs(
                torch.remainder(view_phase - float(center) + 0.5, 1.0) - 0.5
            )
            protected[:, view_index, :] |= view_valid & distance.le(view_width)
    return protected


def _flux_noise_scale(
    values: Any,
    mask: Any,
    config: EventPreservingAugmentationConfig,
) -> tuple[Any, Any]:
    import torch

    small_error_channel, primary_error_channel = config.error_channels
    small_error = values[:, :, small_error_channel, :].abs()
    primary_error = values[:, :, primary_error_channel, :].abs()
    error_values = torch.stack([small_error, primary_error], dim=2)
    small_valid = mask[:, :, small_error_channel, :] & torch.isfinite(small_error)
    primary_valid = mask[:, :, primary_error_channel, :] & torch.isfinite(primary_error)
    error_valid = torch.stack([small_valid, primary_valid], dim=2)
    return error_values, error_valid


def _add_flux_noise(
    values: Any,
    mask: Any,
    *,
    scale: float,
    generator: Any,
    config: EventPreservingAugmentationConfig,
    protected: Any | None = None,
) -> Any:
    import torch

    if float(scale) == 0.0:
        return values
    small_channel, primary_channel, difference_channel = config.flux_channels
    independently_perturbed = [small_channel, primary_channel]
    error_values, error_valid = _flux_noise_scale(values, mask, config)
    active = mask[:, :, independently_perturbed, :] & error_valid
    if protected is not None:
        active = active & ~protected.unsqueeze(2)
    noise = torch.randn(
        error_values.shape,
        generator=generator,
        device=values.device,
        dtype=values.dtype,
    )
    perturbed = (
        values[:, :, independently_perturbed, :]
        + float(scale) * noise * error_values
    )
    values[:, :, independently_perturbed, :] = torch.where(
        active,
        perturbed,
        values[:, :, independently_perturbed, :],
    )
    both_flux_valid = (
        mask[:, :, small_channel, :]
        & mask[:, :, primary_channel, :]
        & torch.isfinite(values[:, :, small_channel, :])
        & torch.isfinite(values[:, :, primary_channel, :])
    )
    if protected is not None:
        both_flux_valid = both_flux_valid & ~protected
    coherent_difference = (
        values[:, :, primary_channel, :] - values[:, :, small_channel, :]
    )
    values[:, :, difference_channel, :] = torch.where(
        both_flux_valid,
        coherent_difference,
        values[:, :, difference_channel, :],
    )
    return values


def augment_ssl_batch(
    batch: Mapping[str, Any],
    *,
    duration_min: Any,
    config: EventPreservingAugmentationConfig = EventPreservingAugmentationConfig(),
    seed: int,
    view_index: int,
) -> dict[str, Any]:
    """Create one deterministic, event-preserving SSL view of a native batch.

    The input mapping and tensors are never mutated.  The caller must pass the
    BLS duration explicitly because the current supervised collate contract
    carries ``period_d`` but not ``duration_min``.
    """

    import torch

    period, duration = _validate_batch_tensor_contract(batch, duration_min)
    harmonic = batch["harmonic_values"]
    harmonic_mask = batch["harmonic_mask"]
    maximum_channel = max(
        config.flux_channels + config.error_channels + (config.phase_channel,)
    )
    if (
        int(harmonic.shape[2]) <= maximum_channel
        or int(batch["local_values"].shape[2]) <= maximum_channel
    ):
        raise ValueError("harmonic/local tensors do not satisfy augmentation channel indices")
    if int(harmonic.shape[1]) != len(config.harmonic_factors):
        raise ValueError(
            "harmonic view count does not match the configured harmonic factors"
        )

    generator = _generator(harmonic.device, seed=int(seed), view_index=int(view_index))
    out = dict(batch)
    augmented_harmonic = harmonic.clone()
    augmented_harmonic_mask = harmonic_mask.clone()
    augmented_local = batch["local_values"].clone()
    augmented_local_mask = batch["local_mask"].clone()
    augmented_periodogram = batch["periodogram_values"].clone()
    augmented_periodogram_mask = batch["periodogram_mask"].clone()

    phase = harmonic[:, :, config.phase_channel, :]
    phase_valid = harmonic_mask[:, :, config.phase_channel, :] & torch.isfinite(phase)
    protected = _harmonic_event_protection_mask(
        phase=phase,
        phase_valid=phase_valid,
        period=period,
        duration=duration,
        config=config,
    )

    cadence_valid = harmonic_mask.any(dim=2)
    dropout_draw = torch.rand(
        cadence_valid.shape,
        generator=generator,
        device=harmonic.device,
        dtype=harmonic.dtype,
    )
    dropped = (
        dropout_draw.lt(float(config.harmonic_cadence_dropout_probability))
        & cadence_valid
        & ~protected
    )
    augmented_harmonic_mask = augmented_harmonic_mask & ~dropped.unsqueeze(2)
    augmented_harmonic = augmented_harmonic.masked_fill(
        dropped.unsqueeze(2),
        0.0,
    )
    augmented_harmonic = _add_flux_noise(
        augmented_harmonic,
        augmented_harmonic_mask,
        scale=float(config.harmonic_flux_noise_scale),
        generator=generator,
        config=config,
        protected=protected,
    )
    augmented_local = _add_flux_noise(
        augmented_local,
        augmented_local_mask,
        scale=float(config.local_flux_noise_scale),
        generator=generator,
        config=config,
    )
    if float(config.periodogram_bin_dropout_probability) > 0:
        periodogram_dropout = torch.rand(
            augmented_periodogram_mask.shape,
            generator=generator,
            device=augmented_periodogram.device,
            dtype=augmented_periodogram.dtype,
        )
        periodogram_dropped = (
            periodogram_dropout.lt(float(config.periodogram_bin_dropout_probability))
            & augmented_periodogram_mask
        )
        augmented_periodogram_mask = (
            augmented_periodogram_mask & ~periodogram_dropped
        )
        augmented_periodogram = augmented_periodogram.masked_fill(
            periodogram_dropped, 0.0
        )

    out["harmonic_values"] = augmented_harmonic
    out["harmonic_mask"] = augmented_harmonic_mask
    out["local_values"] = augmented_local
    out["local_mask"] = augmented_local_mask
    out["periodogram_values"] = augmented_periodogram
    out["periodogram_mask"] = augmented_periodogram_mask
    out["metadata"] = torch.zeros_like(batch["metadata"])
    _validate_model_input_numeric_contract(out, stage="post-augmentation")
    return out


def _off_diagonal(matrix: Any) -> Any:
    import torch

    size = int(matrix.shape[0])
    if matrix.ndim != 2 or int(matrix.shape[1]) != size:
        raise ValueError("covariance matrix must be square")
    if size <= 1:
        return matrix.new_empty((0,))
    mask = ~torch.eye(size, dtype=torch.bool, device=matrix.device)
    return matrix[mask]


def vicreg_loss(
    embedding_a: Any,
    embedding_b: Any,
    *,
    config: VICRegConfig = VICRegConfig(),
) -> tuple[Any, dict[str, Any]]:
    """Return VICReg loss and tensor-valued optimization diagnostics."""

    import torch
    from torch.nn import functional

    if not torch.is_tensor(embedding_a) or not torch.is_tensor(embedding_b):
        raise TypeError("VICReg embeddings must be torch tensors")
    if embedding_a.shape != embedding_b.shape:
        raise ValueError("VICReg embedding views must have identical shapes")
    if embedding_a.ndim != 2:
        raise ValueError("VICReg embeddings must have shape (batch, features)")
    if int(embedding_a.shape[0]) < 2:
        raise ValueError("VICReg requires at least two observations per batch")
    if int(embedding_a.shape[1]) < 1:
        raise ValueError("VICReg embeddings require at least one feature")
    if not embedding_a.is_floating_point() or not embedding_b.is_floating_point():
        raise TypeError("VICReg embeddings must be floating point")
    if embedding_a.device != embedding_b.device:
        raise ValueError("VICReg embedding views must be on the same device")
    if not torch.isfinite(embedding_a).all() or not torch.isfinite(embedding_b).all():
        raise ValueError("VICReg embeddings must be finite")

    invariance = functional.mse_loss(embedding_a, embedding_b)
    centered_a = embedding_a - embedding_a.mean(dim=0)
    centered_b = embedding_b - embedding_b.mean(dim=0)
    std_a = torch.sqrt(embedding_a.var(dim=0, unbiased=False) + float(config.eps))
    std_b = torch.sqrt(embedding_b.var(dim=0, unbiased=False) + float(config.eps))
    variance_a = torch.relu(float(config.target_std) - std_a).mean()
    variance_b = torch.relu(float(config.target_std) - std_b).mean()
    variance = 0.5 * (variance_a + variance_b)

    denominator = float(int(embedding_a.shape[0]) - 1)
    covariance_a_matrix = centered_a.T @ centered_a / denominator
    covariance_b_matrix = centered_b.T @ centered_b / denominator
    feature_count = float(int(embedding_a.shape[1]))
    offdiag_a = _off_diagonal(covariance_a_matrix)
    offdiag_b = _off_diagonal(covariance_b_matrix)
    covariance_a = offdiag_a.square().sum() / feature_count
    covariance_b = offdiag_b.square().sum() / feature_count
    covariance = 0.5 * (covariance_a + covariance_b)

    weighted_invariance = float(config.invariance_weight) * invariance
    weighted_variance = float(config.variance_weight) * variance
    weighted_covariance = float(config.covariance_weight) * covariance
    total = weighted_invariance + weighted_variance + weighted_covariance
    diagnostics = {
        "loss": total,
        "invariance": invariance,
        "variance": variance,
        "covariance": covariance,
        "weighted_invariance": weighted_invariance,
        "weighted_variance": weighted_variance,
        "weighted_covariance": weighted_covariance,
        "std_mean_a": std_a.mean(),
        "std_mean_b": std_b.mean(),
        "variance_a": variance_a,
        "variance_b": variance_b,
        "covariance_a": covariance_a,
        "covariance_b": covariance_b,
    }
    return total, diagnostics


def build_ssl_cache_identity(
    *,
    training_table_sha256: str,
    native_registry_sha256: str,
    split_registry_sha256: str,
    selected_rows: pd.DataFrame,
    profile: str,
    held_out_fold: int,
    seed: int,
    model_config: Mapping[str, Any],
    augmentation_config: EventPreservingAugmentationConfig = EventPreservingAugmentationConfig(),
    vicreg_config: VICRegConfig = VICRegConfig(),
    code_revision: str,
) -> SSLCacheIdentity:
    """Build the exact cache identity for one fold-local SSL encoder."""

    selected_rows_sha256, selected_tics_sha256 = _selection_hashes(selected_rows)
    selected_tics = _integer_series(selected_rows["tic"], name="tic", positive=True)
    return SSLCacheIdentity(
        training_table_sha256=str(training_table_sha256),
        native_registry_sha256=str(native_registry_sha256),
        split_registry_sha256=str(split_registry_sha256),
        selected_rows_sha256=selected_rows_sha256,
        selected_tics_sha256=selected_tics_sha256,
        model_config_sha256=_canonical_sha256(model_config),
        augmentation_config_sha256=_canonical_sha256(augmentation_config),
        vicreg_config_sha256=_canonical_sha256(vicreg_config),
        profile=str(profile),
        held_out_fold=int(held_out_fold),
        seed=int(seed),
        code_revision=str(code_revision),
        n_selected_rows=int(len(selected_rows)),
        n_selected_tics=int(selected_tics.nunique()),
    )


def validate_ssl_cache_identity(
    expected: SSLCacheIdentity,
    observed: SSLCacheIdentity | Mapping[str, Any],
) -> None:
    """Fail closed unless a cached encoder has exactly the expected identity."""

    if not isinstance(expected, SSLCacheIdentity):
        raise TypeError("expected must be an SSLCacheIdentity")
    if isinstance(observed, SSLCacheIdentity):
        observed_identity = observed
        observed_digest = observed.digest()
    elif isinstance(observed, Mapping):
        identity_keys = set(expected.as_dict())
        required_keys = identity_keys | {"identity_sha256"}
        observed_keys = set(observed)
        missing = sorted(required_keys - observed_keys)
        extra = sorted(observed_keys - required_keys)
        if missing or extra:
            raise ValueError(
                f"SSL cache identity keys differ: missing={missing}, extra={extra}"
            )
        payload = {key: observed[key] for key in identity_keys}
        observed_identity = SSLCacheIdentity(**payload)
        observed_digest = str(observed["identity_sha256"])
        if _SHA256_PATTERN.fullmatch(observed_digest) is None:
            raise ValueError("identity_sha256 must be a lowercase 64-character SHA-256")
        if observed_digest != observed_identity.digest():
            raise ValueError("SSL cache identity digest does not match its payload")
    else:
        raise TypeError("observed must be an SSLCacheIdentity or strict manifest mapping")

    if observed_identity != expected or observed_digest != expected.digest():
        mismatches = [
            key
            for key, value in expected.as_dict().items()
            if observed_identity.as_dict().get(key) != value
        ]
        raise ValueError(f"SSL cache identity mismatch: {mismatches}")


__all__ = [
    "HARMONIC_SSL_CONTRACT_VERSION",
    "EventPreservingAugmentationConfig",
    "SSLCacheIdentity",
    "VICRegConfig",
    "augment_ssl_batch",
    "build_ssl_cache_identity",
    "select_ssl_fold_rows",
    "validate_ssl_cache_identity",
    "vicreg_loss",
]
