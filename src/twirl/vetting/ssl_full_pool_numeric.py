"""Fail-closed numerical audit for Teacher v4-SSL model-facing tensors.

The native-input release validator proves file, schema, provenance, and
coverage contracts.  This module supplies the separate model-facing numerical
contract: it collates one exact dataset sample with the production collator and
audits only values that the model mask marks as active.  It never clips a
value, changes a mask, excludes a row, or otherwise repairs an input.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Mapping

import numpy as np
import pandas as pd

from twirl.vetting.harmonic_dataset import collate_native_samples
from twirl.vetting.harmonic_inputs import (
    HARMONIC_FACTORS,
    HARMONIC_VIEW_CHANNELS,
    MODEL_INPUT_CONTRACT_VERSION,
    MODEL_INPUT_NUMERIC_ABS_MAX,
    PERIODOGRAM_CHANNELS,
)
from twirl.vetting.ssl_full_pool_eligibility import (
    PRODUCTION_ELIGIBLE_IDENTITY_SHA256,
    PRODUCTION_ELIGIBLE_OBSERVATIONS,
    PRODUCTION_EXCLUDED_IDENTITY_SHA256,
    PRODUCTION_EXCLUDED_OBSERVATIONS,
    PRODUCTION_FULL_IDENTITY_SHA256,
    PRODUCTION_FULL_OBSERVATIONS,
    observation_identity_sha256,
)


TEACHER_SSL_NUMERIC_ENVELOPE_V1 = (
    "twirl_teacher_ssl_fullpool_numeric_envelope_v1"
)
MODEL_INPUT_NUMERIC_AUDIT_SCHEMA = (
    "twirl_teacher_ssl_fullpool_model_input_numeric_audit_v1"
)
MODEL_INPUT_NUMERIC_RELEASE_SCHEMA = (
    "twirl_teacher_ssl_fullpool_model_input_numeric_release_v1"
)
MODEL_INPUT_NUMERIC_ENVELOPE_CONTRACT = TEACHER_SSL_NUMERIC_ENVELOPE_V1
MODEL_INPUT_NUMERIC_AUTHORITY_NAMES: tuple[str, ...] = (
    "ssl_registry",
    "ssl_registry_summary",
    "native_registry",
    "native_registry_summary",
    "native_release_summary",
)
_MODEL_INPUT_NUMERIC_AUDIT_COLUMNS: tuple[str, ...] = (
    "ssl_observation_id",
    "sector",
    "tic",
    "ssl_pool_include",
    "numeric_status",
    "model_input_numeric_passed",
    "n_failures",
    "failure_codes",
    "failures_json",
    "harmonic_max_abs",
    "local_max_abs",
    "periodogram_max_abs",
)
_MODEL_INPUT_NUMERIC_AUDIT_DTYPES: tuple[str, ...] = (
    "object",
    "int64",
    "int64",
    "bool",
    "object",
    "boolean",
    "int64",
    "object",
    "object",
    "float64",
    "float64",
    "float64",
)
_MODEL_INPUT_NUMERIC_AUDIT_TEXT_COLUMNS: tuple[str, ...] = (
    "ssl_observation_id",
    "numeric_status",
    "failure_codes",
    "failures_json",
)
_MODEL_INPUT_NUMERIC_MAX_COLUMNS: tuple[str, ...] = (
    "harmonic_max_abs",
    "local_max_abs",
    "periodogram_max_abs",
)
_MODEL_INPUT_NUMERIC_REPORT_ROW_FIELDS = frozenset(
    {
        "ssl_observation_id",
        "sector",
        "tic",
        "passed",
        "n_failures",
        "failure_codes",
        "failures",
        "harmonic_max_abs",
        "local_max_abs",
        "periodogram_max_abs",
    }
)

_FLOAT32_MAX = np.finfo(np.float32).max
_FLOAT32_SQRT_MAX = np.float32(np.sqrt(_FLOAT32_MAX))
FLOAT32_SQUARE_SAFE_ABS_MAX = float(
    np.nextafter(_FLOAT32_SQRT_MAX, np.float32(0.0))
)


@dataclass(frozen=True)
class FullPoolNumericEnvelope:
    """Explicit domains for the v1 Teacher v4-SSL numerical audit."""

    contract_version: str = TEACHER_SSL_NUMERIC_ENVELOPE_V1
    unbounded_model_abs_max: float = MODEL_INPUT_NUMERIC_ABS_MAX
    float32_square_safe_abs_max: float = FLOAT32_SQUARE_SAFE_ABS_MAX
    phase_min: float = -0.5
    phase_max: float = 0.5
    binary_values: tuple[float, float] = (0.0, 1.0)
    harmonic_view_count: int = len(HARMONIC_FACTORS)
    local_view_count: int = 2 * len(HARMONIC_FACTORS)

    def __post_init__(self) -> None:
        if self.contract_version != TEACHER_SSL_NUMERIC_ENVELOPE_V1:
            raise ValueError(
                "unsupported Teacher SSL numerical-envelope contract"
            )
        if (
            not np.isfinite(self.unbounded_model_abs_max)
            or self.unbounded_model_abs_max <= 0
        ):
            raise ValueError("unbounded_model_abs_max must be finite and positive")
        if self.unbounded_model_abs_max != MODEL_INPUT_NUMERIC_ABS_MAX:
            raise ValueError(
                "unbounded_model_abs_max differs from the locked v1 envelope"
            )
        if (
            not np.isfinite(self.float32_square_safe_abs_max)
            or self.float32_square_safe_abs_max <= 0
        ):
            raise ValueError(
                "float32_square_safe_abs_max must be finite and positive"
            )
        if (
            not np.isfinite(self.phase_min)
            or not np.isfinite(self.phase_max)
            or self.phase_min >= self.phase_max
        ):
            raise ValueError("phase bounds must be finite and increasing")
        if tuple(float(value) for value in self.binary_values) != (0.0, 1.0):
            raise ValueError("binary_values must be exactly (0.0, 1.0)")
        if self.harmonic_view_count != len(HARMONIC_FACTORS):
            raise ValueError("harmonic_view_count differs from the input contract")
        if self.local_view_count != 2 * len(HARMONIC_FACTORS):
            raise ValueError("local_view_count differs from the input contract")

    def as_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable statement of the numerical envelope."""

        return asdict(self)


FULL_POOL_NUMERIC_ENVELOPE_V1 = FullPoolNumericEnvelope()


def _canonical_sha256(value: Mapping[str, Any]) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


MODEL_INPUT_NUMERIC_ENVELOPE_SHA256 = _canonical_sha256(
    FULL_POOL_NUMERIC_ENVELOPE_V1.as_dict()
)


def _as_numpy(value: Any) -> np.ndarray:
    """Return a detached CPU NumPy view without importing torch eagerly."""

    detached = value.detach() if hasattr(value, "detach") else value
    cpu_value = detached.cpu() if hasattr(detached, "cpu") else detached
    return np.asarray(cpu_value)


def _json_number(value: Any) -> int | float | str:
    item = value.item() if isinstance(value, np.generic) else value
    if isinstance(item, (int, np.integer)):
        return int(item)
    number = float(item)
    if np.isnan(number):
        return "nan"
    if np.isposinf(number):
        return "inf"
    if np.isneginf(number):
        return "-inf"
    return number


def _first_index(mask: np.ndarray) -> list[int] | None:
    indices = np.argwhere(mask)
    return [int(value) for value in indices[0]] if len(indices) else None


def _failure(
    failures: list[dict[str, Any]],
    *,
    code: str,
    tensor: str,
    count: int,
    first_index: list[int] | None = None,
    channel_index: int | None = None,
    channel_name: str | None = None,
    **details: Any,
) -> None:
    item: dict[str, Any] = {
        "code": str(code),
        "tensor": str(tensor),
        "count": int(count),
    }
    if first_index is not None:
        item["first_index"] = first_index
    if channel_index is not None:
        item["channel_index"] = int(channel_index)
    if channel_name is not None:
        item["channel_name"] = str(channel_name)
    item.update(
        {
            key: _json_number(value)
            if isinstance(value, (int, float, np.number))
            else value
            for key, value in details.items()
        }
    )
    failures.append(item)


def _channel_statistics(
    values: np.ndarray,
    mask: np.ndarray,
    *,
    channel_axis: int,
    channel_names: tuple[str, ...],
) -> dict[str, dict[str, int | float | None]]:
    values_by_channel = np.moveaxis(values, channel_axis, 0)
    masks_by_channel = np.moveaxis(mask, channel_axis, 0)
    result: dict[str, dict[str, int | float | None]] = {}
    for channel_index, channel_name in enumerate(channel_names):
        channel_values = values_by_channel[channel_index]
        channel_mask = masks_by_channel[channel_index]
        selected = channel_values[channel_mask & np.isfinite(channel_values)]
        result[channel_name] = {
            "n_model_values": int(channel_mask.sum()),
            "n_masked_values": int(channel_mask.size - channel_mask.sum()),
            "min": float(selected.min()) if selected.size else None,
            "max": float(selected.max()) if selected.size else None,
            "max_abs": float(np.abs(selected).max()) if selected.size else None,
        }
    return result


def _audit_masked_sequence(
    *,
    tensor_name: str,
    values: np.ndarray,
    mask: np.ndarray,
    channel_names: tuple[str, ...],
    expected_views: int | None,
    envelope: FullPoolNumericEnvelope,
    failures: list[dict[str, Any]],
) -> dict[str, dict[str, int | float | None]]:
    if values.ndim not in (2, 3):
        _failure(
            failures,
            code="rank_mismatch",
            tensor=tensor_name,
            count=1,
            observed_rank=values.ndim,
        )
        return {}
    if values.shape != mask.shape:
        _failure(
            failures,
            code="mask_shape_mismatch",
            tensor=tensor_name,
            count=1,
            values_shape=list(values.shape),
            mask_shape=list(mask.shape),
        )
        return {}
    if mask.dtype != np.bool_:
        _failure(
            failures,
            code="mask_dtype_mismatch",
            tensor=tensor_name,
            count=1,
            observed_dtype=str(mask.dtype),
        )
        return {}
    if values.dtype != np.float32:
        _failure(
            failures,
            code="value_dtype_mismatch",
            tensor=tensor_name,
            count=1,
            observed_dtype=str(values.dtype),
        )
    if expected_views is not None and (
        values.ndim != 3 or values.shape[0] != expected_views
    ):
        _failure(
            failures,
            code="view_count_mismatch",
            tensor=tensor_name,
            count=1,
            observed_views=int(values.shape[0]) if values.ndim else 0,
            expected_views=expected_views,
        )
    channel_axis = values.ndim - 2
    if values.shape[channel_axis] != len(channel_names):
        _failure(
            failures,
            code="channel_count_mismatch",
            tensor=tensor_name,
            count=1,
            observed_channels=int(values.shape[channel_axis]),
            expected_channels=len(channel_names),
        )
        return {}

    work = values.astype(np.float32, copy=False)
    nonfinite_payload = ~np.isfinite(work)
    if np.any(nonfinite_payload):
        _failure(
            failures,
            code="tensor_payload_nonfinite",
            tensor=tensor_name,
            count=int(nonfinite_payload.sum()),
            first_index=_first_index(nonfinite_payload),
        )

    model_nonfinite = mask & nonfinite_payload
    if np.any(model_nonfinite):
        _failure(
            failures,
            code="model_value_nonfinite",
            tensor=tensor_name,
            count=int(model_nonfinite.sum()),
            first_index=_first_index(model_nonfinite),
        )

    envelope_exceeded = mask & np.isfinite(work) & (
        np.abs(work) > envelope.unbounded_model_abs_max
    )
    if np.any(envelope_exceeded):
        _failure(
            failures,
            code="numeric_envelope_exceeded",
            tensor=tensor_name,
            count=int(envelope_exceeded.sum()),
            first_index=_first_index(envelope_exceeded),
            maximum_allowed_abs=envelope.unbounded_model_abs_max,
        )

    # LayerNorm and variance-like operations square float32 activations.  Audit
    # the input-domain operation directly; a masked finite extreme is harmless
    # because the model multiplies it by a false channel mask before the stem.
    square_overflow = mask & np.isfinite(work) & (
        np.abs(work) > envelope.float32_square_safe_abs_max
    )
    with np.errstate(over="ignore", invalid="ignore"):
        squared = np.square(np.where(mask, work, np.float32(0.0)), dtype=np.float32)
    square_overflow |= mask & ~np.isfinite(squared)
    if np.any(square_overflow):
        _failure(
            failures,
            code="float32_square_overflow",
            tensor=tensor_name,
            count=int(square_overflow.sum()),
            first_index=_first_index(square_overflow),
            maximum_allowed_abs=envelope.float32_square_safe_abs_max,
        )

    values_by_channel = np.moveaxis(work, channel_axis, 0)
    masks_by_channel = np.moveaxis(mask, channel_axis, 0)
    if tuple(channel_names) == tuple(HARMONIC_VIEW_CHANNELS):
        phase_mask = masks_by_channel[5]
        quality_mask = masks_by_channel[6]
        coordinate_mask_mismatch = phase_mask != quality_mask
        if np.any(coordinate_mask_mismatch):
            _failure(
                failures,
                code="phase_quality_mask_mismatch",
                tensor=tensor_name,
                count=int(coordinate_mask_mismatch.sum()),
                first_index=_first_index(coordinate_mask_mismatch),
            )
        photometry_without_coordinates = (
            np.any(masks_by_channel[:5], axis=0)
            & ~(phase_mask & quality_mask)
        )
        if np.any(photometry_without_coordinates):
            _failure(
                failures,
                code="photometry_without_coordinates",
                tensor=tensor_name,
                count=int(photometry_without_coordinates.sum()),
                first_index=_first_index(photometry_without_coordinates),
            )
        quality_values = values_by_channel[6]
        quality_active = (
            quality_mask
            & np.isfinite(quality_values)
            & (quality_values == 1.0)
        )
        photometry_unmasked = (
            masks_by_channel[:5]
            & np.expand_dims(quality_active, axis=0)
        )
        if np.any(photometry_unmasked):
            original_axes = np.moveaxis(
                photometry_unmasked, 0, channel_axis
            )
            _failure(
                failures,
                code="quality_bad_photometry_unmasked",
                tensor=tensor_name,
                count=int(photometry_unmasked.sum()),
                first_index=_first_index(original_axes),
            )
    for channel_index, channel_name in enumerate(channel_names):
        channel_values = values_by_channel[channel_index]
        channel_mask = masks_by_channel[channel_index] & np.isfinite(channel_values)
        if channel_name in {
            "raw_error_small_scaled",
            "raw_error_primary_scaled",
        }:
            invalid = channel_mask & (channel_values < 0)
            if np.any(invalid):
                _failure(
                    failures,
                    code="negative_error",
                    tensor=tensor_name,
                    count=int(invalid.sum()),
                    first_index=_first_index(invalid),
                    channel_index=channel_index,
                    channel_name=channel_name,
                    minimum_allowed=0.0,
                )
        elif channel_name == "orbital_phase":
            invalid = channel_mask & (
                (channel_values < envelope.phase_min)
                | (channel_values > envelope.phase_max)
            )
            if np.any(invalid):
                _failure(
                    failures,
                    code="coordinate_out_of_range",
                    tensor=tensor_name,
                    count=int(invalid.sum()),
                    first_index=_first_index(invalid),
                    channel_index=channel_index,
                    channel_name=channel_name,
                    minimum_allowed=envelope.phase_min,
                    maximum_allowed=envelope.phase_max,
                )
        elif channel_name == "quality_nonzero":
            invalid = channel_mask & ~np.isin(
                channel_values, envelope.binary_values
            )
            if np.any(invalid):
                _failure(
                    failures,
                    code="binary_domain_violation",
                    tensor=tensor_name,
                    count=int(invalid.sum()),
                    first_index=_first_index(invalid),
                    channel_index=channel_index,
                    channel_name=channel_name,
                    allowed_values=list(envelope.binary_values),
                )

    return _channel_statistics(
        work,
        mask,
        channel_axis=channel_axis,
        channel_names=channel_names,
    )


def _audit_unmasked_values(
    *,
    tensor_name: str,
    values: np.ndarray,
    envelope: FullPoolNumericEnvelope,
    failures: list[dict[str, Any]],
) -> dict[str, int | float | None]:
    if values.dtype != np.float32:
        _failure(
            failures,
            code="value_dtype_mismatch",
            tensor=tensor_name,
            count=1,
            observed_dtype=str(values.dtype),
        )
    work = values.astype(np.float32, copy=False)
    nonfinite = ~np.isfinite(work)
    if np.any(nonfinite):
        _failure(
            failures,
            code="model_value_nonfinite",
            tensor=tensor_name,
            count=int(nonfinite.sum()),
            first_index=_first_index(nonfinite),
        )
    envelope_exceeded = np.isfinite(work) & (
        np.abs(work) > envelope.unbounded_model_abs_max
    )
    if np.any(envelope_exceeded):
        _failure(
            failures,
            code="numeric_envelope_exceeded",
            tensor=tensor_name,
            count=int(envelope_exceeded.sum()),
            first_index=_first_index(envelope_exceeded),
            maximum_allowed_abs=envelope.unbounded_model_abs_max,
        )
    square_overflow = np.isfinite(work) & (
        np.abs(work) > envelope.float32_square_safe_abs_max
    )
    with np.errstate(over="ignore", invalid="ignore"):
        squared = np.square(work, dtype=np.float32)
    square_overflow |= ~np.isfinite(squared)
    if np.any(square_overflow):
        _failure(
            failures,
            code="float32_square_overflow",
            tensor=tensor_name,
            count=int(square_overflow.sum()),
            first_index=_first_index(square_overflow),
            maximum_allowed_abs=envelope.float32_square_safe_abs_max,
        )
    finite = work[np.isfinite(work)]
    return {
        "n_model_values": int(work.size),
        "min": float(finite.min()) if finite.size else None,
        "max": float(finite.max()) if finite.size else None,
        "max_abs": float(np.abs(finite).max()) if finite.size else None,
    }


def audit_collated_sample(
    batch: Mapping[str, Any],
    *,
    sample_index: int = 0,
    envelope: FullPoolNumericEnvelope = FULL_POOL_NUMERIC_ENVELOPE_V1,
) -> dict[str, Any]:
    """Audit one row of an already-collated exact model-facing batch.

    The returned mapping is JSON serializable.  Numerical violations are
    accumulated in ``failures`` rather than repaired or silently omitted.
    """

    required = {
        "review_id",
        "tic",
        "period_d",
        "duration_min",
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
        raise KeyError(f"collated numeric audit is missing keys: {missing}")
    batch_size = int(_as_numpy(batch["period_d"]).shape[0])
    index = int(sample_index)
    if index < 0 or index >= batch_size:
        raise IndexError(
            f"sample_index {index} is outside collated batch size {batch_size}"
        )

    failures: list[dict[str, Any]] = []
    maxima: dict[str, Any] = {}
    for tensor_name, mask_name, channel_names, expected_views in (
        (
            "harmonic_values",
            "harmonic_mask",
            HARMONIC_VIEW_CHANNELS,
            envelope.harmonic_view_count,
        ),
        (
            "local_values",
            "local_mask",
            HARMONIC_VIEW_CHANNELS,
            envelope.local_view_count,
        ),
        (
            "periodogram_values",
            "periodogram_mask",
            PERIODOGRAM_CHANNELS,
            None,
        ),
    ):
        maxima[tensor_name] = _audit_masked_sequence(
            tensor_name=tensor_name,
            values=_as_numpy(batch[tensor_name])[index],
            mask=_as_numpy(batch[mask_name])[index],
            channel_names=tuple(channel_names),
            expected_views=expected_views,
            envelope=envelope,
            failures=failures,
        )

    maxima["metadata"] = _audit_unmasked_values(
        tensor_name="metadata",
        values=_as_numpy(batch["metadata"])[index],
        envelope=envelope,
        failures=failures,
    )

    scalar_values: dict[str, float | None] = {}
    for tensor_name in ("period_d", "duration_min"):
        value_array = _as_numpy(batch[tensor_name])
        value = np.asarray(value_array[index])
        maxima[tensor_name] = _audit_unmasked_values(
            tensor_name=tensor_name,
            values=value,
            envelope=envelope,
            failures=failures,
        )
        scalar = float(value)
        scalar_values[tensor_name] = scalar if np.isfinite(scalar) else None
        if not np.isfinite(scalar) or scalar <= 0:
            _failure(
                failures,
                code="nonpositive_scalar",
                tensor=tensor_name,
                count=1,
                observed_value=scalar,
            )

    review_ids = batch["review_id"]
    review_id = str(review_ids[index])
    tic = int(_as_numpy(batch["tic"])[index])
    return {
        "schema_version": MODEL_INPUT_NUMERIC_AUDIT_SCHEMA,
        "contract_version": envelope.contract_version,
        "envelope": envelope.as_dict(),
        "review_id": review_id,
        "tic": tic,
        **scalar_values,
        "passed": not failures,
        "n_failures": len(failures),
        "failures": failures,
        "maxima": maxima,
        "action": "audit_only_no_clip_no_exclusion",
    }


def audit_model_facing_sample(
    sample: Mapping[str, Any],
    *,
    envelope: FullPoolNumericEnvelope = FULL_POOL_NUMERIC_ENVELOPE_V1,
) -> dict[str, Any]:
    """Collate and audit one exact ``HarmonicNativeDataset`` sample."""

    batch = collate_native_samples([dict(sample)])
    return audit_collated_sample(batch, sample_index=0, envelope=envelope)


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _normalized_file_binding(
    record: Any,
    *,
    context: str,
    verify_file: bool,
) -> dict[str, Any]:
    if not isinstance(record, Mapping):
        raise ValueError(f"{context} must be a file-binding mapping")
    required = {"path", "sha256", "size_bytes"}
    if set(record) != required:
        raise ValueError(
            f"{context} must contain exactly {sorted(required)}"
        )
    path = Path(str(record["path"])).expanduser()
    if not path.is_absolute():
        raise ValueError(f"{context}.path must be absolute")
    digest = str(record["sha256"]).strip()
    if re.fullmatch(r"[0-9a-f]{64}", digest) is None:
        raise ValueError(f"{context}.sha256 is not a lowercase SHA-256")
    size = record["size_bytes"]
    if type(size) is not int or size <= 0:
        raise ValueError(f"{context}.size_bytes must be a positive integer")
    normalized = {
        "path": str(path),
        "sha256": digest,
        "size_bytes": int(size),
    }
    if verify_file:
        resolved = path.resolve(strict=True)
        if not resolved.is_file():
            raise ValueError(f"{context}.path is not a regular file")
        if int(resolved.stat().st_size) != int(size):
            raise ValueError(f"{context} file size differs from its binding")
        if _file_sha256(resolved) != digest:
            raise ValueError(f"{context} file SHA-256 differs from its binding")
    return normalized


def _verify_sha256_sidecar(
    path: Path,
    *,
    digest: str,
    context: str,
) -> None:
    sidecar = Path(str(path) + ".sha256")
    expected = f"{digest}  {path.name}\n"
    try:
        observed = sidecar.read_text(encoding="ascii")
    except (OSError, UnicodeDecodeError) as exc:
        raise ValueError(f"{context} SHA-256 sidecar is missing or invalid") from exc
    if observed != expected:
        raise ValueError(f"{context} SHA-256 sidecar differs")


def _strict_json_object(path: Path, *, context: str) -> dict[str, Any]:
    def no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        value: dict[str, Any] = {}
        for key, item in pairs:
            if key in value:
                raise ValueError(f"{context} contains duplicate key {key!r}")
            value[key] = item
        return value

    def reject_constant(value: str) -> None:
        raise ValueError(f"{context} contains non-finite constant {value!r}")

    try:
        payload = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=no_duplicates,
            parse_constant=reject_constant,
        )
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"{context} is not strict UTF-8 JSON") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"{context} root must be a mapping")
    return payload


def _validate_numeric_audit(path: Path) -> pd.DataFrame:
    try:
        frame = pd.read_parquet(path)
    except Exception as exc:
        raise ValueError("numeric audit is not a readable Parquet table") from exc
    if tuple(frame.columns) != _MODEL_INPUT_NUMERIC_AUDIT_COLUMNS:
        raise ValueError("numeric audit table has the wrong exact schema")
    observed_dtypes = tuple(
        str(frame[name].dtype) for name in _MODEL_INPUT_NUMERIC_AUDIT_COLUMNS
    )
    if observed_dtypes != _MODEL_INPUT_NUMERIC_AUDIT_DTYPES:
        raise ValueError("numeric audit table has the wrong exact dtypes")
    for name in _MODEL_INPUT_NUMERIC_AUDIT_TEXT_COLUMNS:
        if frame[name].isna().any() or not frame[name].map(
            lambda value: type(value) is str
        ).all():
            raise ValueError(f"numeric audit column {name} must contain strings")
    if (
        frame["ssl_observation_id"].map(
            lambda value: bool(value) and value.strip() == value
        ).eq(False).any()
        or (frame["sector"] <= 0).any()
        or (frame["tic"] <= 0).any()
        or frame["ssl_observation_id"].duplicated().any()
        or frame[["sector", "tic"]].duplicated().any()
        or len(frame) != PRODUCTION_FULL_OBSERVATIONS
        or observation_identity_sha256(frame)
        != PRODUCTION_FULL_IDENTITY_SHA256
    ):
        raise ValueError("numeric audit table differs from the full pool")
    maxima = frame.loc[:, list(_MODEL_INPUT_NUMERIC_MAX_COLUMNS)].to_numpy(
        dtype=np.float64,
        copy=False,
    )
    if np.isinf(maxima).any() or (maxima[np.isfinite(maxima)] < 0).any():
        raise ValueError("numeric audit maxima must be nonnegative finite or null")

    include = frame["ssl_pool_include"]
    included = frame.loc[include]
    excluded = frame.loc[~include]
    if (
        len(included) != PRODUCTION_ELIGIBLE_OBSERVATIONS
        or observation_identity_sha256(included)
        != PRODUCTION_ELIGIBLE_IDENTITY_SHA256
        or not included["numeric_status"].eq("passed").all()
        or included["model_input_numeric_passed"].isna().any()
        or not included["model_input_numeric_passed"].eq(True).all()
        or not included["n_failures"].eq(0).all()
        or not included["failure_codes"].eq("[]").all()
        or not included["failures_json"].eq("[]").all()
    ):
        raise ValueError("numeric audit eligible partition is invalid")
    if (
        len(excluded) != PRODUCTION_EXCLUDED_OBSERVATIONS
        or observation_identity_sha256(excluded)
        != PRODUCTION_EXCLUDED_IDENTITY_SHA256
        or not excluded["numeric_status"].eq("not_model_eligible").all()
        or not excluded["model_input_numeric_passed"].isna().all()
        or not excluded["n_failures"].eq(0).all()
        or not excluded["failure_codes"].eq("[]").all()
        or not excluded["failures_json"].eq("[]").all()
        or not excluded[
            list(_MODEL_INPUT_NUMERIC_MAX_COLUMNS)
        ].isna().all().all()
    ):
        raise ValueError("numeric audit excluded partition is invalid")
    return frame


def _validate_numeric_report_rows(
    rows: Any,
    *,
    sector: int,
    scanned_observations: int,
    observation_identity: str,
    context: str,
) -> tuple[list[tuple[str, int, int]], int, int]:
    if type(rows) is not list or len(rows) != scanned_observations:
        raise ValueError(f"{context} rows differ from the declared count")
    identities: list[tuple[str, int, int]] = []
    passed_observations = 0
    failed_observations = 0
    for row_index, row in enumerate(rows):
        row_context = f"{context}.rows[{row_index}]"
        if type(row) is not dict or set(row) != (
            _MODEL_INPUT_NUMERIC_REPORT_ROW_FIELDS
        ):
            raise ValueError(f"{row_context} has an invalid exact shape")
        ssl_observation_id = row["ssl_observation_id"]
        row_sector = row["sector"]
        tic = row["tic"]
        passed = row["passed"]
        n_failures = row["n_failures"]
        failure_codes = row["failure_codes"]
        failures = row["failures"]
        if (
            type(ssl_observation_id) is not str
            or not ssl_observation_id
            or ssl_observation_id.strip() != ssl_observation_id
            or type(row_sector) is not int
            or row_sector != sector
            or type(tic) is not int
            or tic <= 0
            or type(passed) is not bool
            or type(n_failures) is not int
            or n_failures < 0
        ):
            raise ValueError(f"{row_context} has invalid strict field types")
        if (
            type(failure_codes) is not list
            or any(type(value) is not str or not value for value in failure_codes)
            or failure_codes != sorted(set(failure_codes))
            or type(failures) is not list
            or any(type(value) is not dict for value in failures)
            or n_failures != len(failures)
        ):
            raise ValueError(f"{row_context} has invalid failure evidence")
        for name in _MODEL_INPUT_NUMERIC_MAX_COLUMNS:
            value = row[name]
            if value is not None and (
                type(value) is not float
                or not np.isfinite(value)
                or value < 0
            ):
                raise ValueError(
                    f"{row_context}.{name} must be a nonnegative float or null"
                )
        if passed:
            passed_observations += 1
            if n_failures != 0 or failure_codes or failures:
                raise ValueError(
                    f"{row_context} passed with nonempty failure evidence"
                )
        else:
            failed_observations += 1
            if n_failures <= 0 or not failure_codes:
                raise ValueError(
                    f"{row_context} failed without complete failure evidence"
                )
        identities.append((ssl_observation_id, row_sector, tic))

    if (
        len({identity[0] for identity in identities}) != len(identities)
        or len({identity[1:] for identity in identities}) != len(identities)
        or observation_identity_sha256(
            [(identity[1], identity[2]) for identity in identities]
        )
        != observation_identity
    ):
        raise ValueError(f"{context} row identity coverage is invalid")
    return identities, passed_observations, failed_observations


def _validate_numeric_evidence(
    payload: Mapping[str, Any],
    *,
    release_path: Path,
    code_revision: str,
    normalized_authorities: Mapping[str, Mapping[str, Any]],
) -> None:
    for name, expected in (
        ("quality_bad_photometry_policy_verified", True),
        ("float32_conversion_verified", True),
        ("real_only", True),
        ("labels_consumed", False),
        ("injections_consumed", False),
    ):
        if payload.get(name) is not expected:
            raise ValueError(f"numeric-gate release {name} differs from contract")
    if payload.get("action") != "audit_only_no_clip_no_exclusion":
        raise ValueError("numeric-gate release action differs from contract")

    outputs = payload.get("outputs")
    if not isinstance(outputs, Mapping) or set(outputs) != {"numeric_audit"}:
        raise ValueError("numeric-gate release output inventory is invalid")
    audit_binding = _normalized_file_binding(
        outputs["numeric_audit"],
        context="outputs.numeric_audit",
        verify_file=True,
    )
    _verify_sha256_sidecar(
        Path(audit_binding["path"]),
        digest=str(audit_binding["sha256"]),
        context="numeric audit",
    )
    audit = _validate_numeric_audit(Path(audit_binding["path"]))

    shard_reports = payload.get("shard_reports")
    if not isinstance(shard_reports, list) or len(shard_reports) != 112:
        raise ValueError("numeric-gate release requires exactly 112 shard reports")
    expected_inventory = {
        (sector, shard_index)
        for sector in range(56, 63)
        for shard_index in range(16)
    }
    observed_inventory: set[tuple[int, int]] = set()
    observed_report_paths: set[str] = set()
    observed_native_paths: set[str] = set()
    reported_eligible_rows: list[tuple[str, int, int]] = []
    scanned_observations = 0
    for index, record in enumerate(shard_reports):
        context = f"shard_reports[{index}]"
        if not isinstance(record, Mapping) or set(record) != {
            "sector",
            "shard_index",
            "path",
            "size_bytes",
            "sha256",
            "native_h5",
            "observation_identity_sha256",
            "scanned_observations",
        }:
            raise ValueError(f"{context} has an invalid record shape")
        sector = record.get("sector")
        shard_index = record.get("shard_index")
        scanned = record.get("scanned_observations")
        if (
            type(sector) is not int
            or type(shard_index) is not int
            or type(scanned) is not int
            or scanned <= 0
        ):
            raise ValueError(f"{context} has invalid integer fields")
        key = (int(sector), int(shard_index))
        if key not in expected_inventory or key in observed_inventory:
            raise ValueError(f"{context} has an invalid sector/shard key")
        observed_inventory.add(key)
        scanned_observations += int(scanned)
        identity = record.get("observation_identity_sha256")
        if not isinstance(identity, str) or re.fullmatch(
            r"[0-9a-f]{64}",
            identity,
        ) is None:
            raise ValueError(f"{context} has an invalid observation identity")

        report_binding = _normalized_file_binding(
            {
                "path": record["path"],
                "size_bytes": record["size_bytes"],
                "sha256": record["sha256"],
            },
            context=context,
            verify_file=True,
        )
        if report_binding["path"] in observed_report_paths:
            raise ValueError("numeric-gate release reuses a shard report path")
        observed_report_paths.add(str(report_binding["path"]))
        report_path = Path(report_binding["path"])
        _verify_sha256_sidecar(
            report_path,
            digest=str(report_binding["sha256"]),
            context=context,
        )

        native_binding = _normalized_file_binding(
            record["native_h5"],
            context=f"{context}.native_h5",
            verify_file=False,
        )
        if native_binding["path"] in observed_native_paths:
            raise ValueError("numeric-gate release reuses a native HDF5 path")
        observed_native_paths.add(str(native_binding["path"]))

        report = _strict_json_object(
            report_path,
            context=f"{context} report",
        )
        counts = report.get("counts")
        if (
            report.get("schema_version") != MODEL_INPUT_NUMERIC_AUDIT_SCHEMA
            or report.get("code_revision") != code_revision
            or report.get("model_input_contract_version")
            != MODEL_INPUT_CONTRACT_VERSION
            or report.get("envelope_contract")
            != MODEL_INPUT_NUMERIC_ENVELOPE_CONTRACT
            or report.get("envelope_canonical_sha256")
            != MODEL_INPUT_NUMERIC_ENVELOPE_SHA256
            or not isinstance(report.get("envelope"), Mapping)
            or _canonical_sha256(report["envelope"])
            != MODEL_INPUT_NUMERIC_ENVELOPE_SHA256
            or report.get("passed") is not True
            or type(report.get("sector")) is not int
            or report.get("sector") != sector
            or type(report.get("shard_index")) is not int
            or report.get("shard_index") != shard_index
            or type(report.get("n_shards")) is not int
            or report.get("n_shards") != 16
            or report.get("observation_identity_sha256") != identity
            or report.get("action") != "audit_only_no_clip_no_exclusion"
            or not isinstance(counts, Mapping)
            or set(counts)
            != {
                "scanned_observations",
                "passed_observations",
                "failed_observations",
            }
        ):
            raise ValueError(f"{context} report header differs from release")
        report_counts = {
            name: counts.get(name)
            for name in (
                "scanned_observations",
                "passed_observations",
                "failed_observations",
            )
        }
        if any(type(value) is not int for value in report_counts.values()):
            raise ValueError(f"{context} report counts must be strict integers")
        row_identities, passed_observations, failed_observations = (
            _validate_numeric_report_rows(
                report.get("rows"),
                sector=sector,
                scanned_observations=scanned,
                observation_identity=identity,
                context=f"{context} report",
            )
        )
        if report_counts != {
            "scanned_observations": len(row_identities),
            "passed_observations": passed_observations,
            "failed_observations": failed_observations,
        } or report_counts != {
            "scanned_observations": scanned,
            "passed_observations": scanned,
            "failed_observations": 0,
        }:
            raise ValueError(f"{context} report row counts differ from release")
        reported_eligible_rows.extend(row_identities)
        report_native = _normalized_file_binding(
            report.get("native_h5"),
            context=f"{context} report native_h5",
            verify_file=False,
        )
        if report_native != native_binding:
            raise ValueError(f"{context} report native binding differs")
        report_authorities = report.get("authority_bindings")
        if not isinstance(report_authorities, Mapping) or set(
            report_authorities
        ) != set(MODEL_INPUT_NUMERIC_AUTHORITY_NAMES):
            raise ValueError(f"{context} report authorities are invalid")
        normalized_report_authorities = {
            name: _normalized_file_binding(
                report_authorities[name],
                context=f"{context}.authority_bindings.{name}",
                verify_file=False,
            )
            for name in MODEL_INPUT_NUMERIC_AUTHORITY_NAMES
        }
        if normalized_report_authorities != dict(normalized_authorities):
            raise ValueError(f"{context} report authorities differ")

    if observed_inventory != expected_inventory:
        raise ValueError("numeric-gate shard report inventory is not exact")
    if scanned_observations != PRODUCTION_ELIGIBLE_OBSERVATIONS:
        raise ValueError(
            "numeric-gate shard reports do not cover all eligible observations"
        )
    if (
        len(reported_eligible_rows) != PRODUCTION_ELIGIBLE_OBSERVATIONS
        or len({row[0] for row in reported_eligible_rows})
        != len(reported_eligible_rows)
        or len({row[1:] for row in reported_eligible_rows})
        != len(reported_eligible_rows)
        or observation_identity_sha256(
            [(row[1], row[2]) for row in reported_eligible_rows]
        )
        != PRODUCTION_ELIGIBLE_IDENTITY_SHA256
    ):
        raise ValueError(
            "numeric-gate shard report union differs from eligible production"
        )
    audit_eligible_rows = {
        (str(row.ssl_observation_id), int(row.sector), int(row.tic))
        for row in audit.loc[audit["ssl_pool_include"]].itertuples(index=False)
    }
    if audit_eligible_rows != set(reported_eligible_rows):
        raise ValueError(
            "numeric audit eligible union differs from shard report rows"
        )
    _verify_sha256_sidecar(
        release_path,
        digest=_file_sha256(release_path),
        context="numeric-gate release",
    )


def _exact_int(mapping: Mapping[str, Any], name: str) -> int:
    value = mapping.get(name)
    if type(value) is not int:
        raise ValueError(f"numeric release count {name} must be an integer")
    return int(value)


def validate_numeric_gate_release(
    path: Path,
    *,
    expected_code_revision: str | None = None,
    expected_authority_bindings: Mapping[str, Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    """Load and strictly validate a full-pool model-input numerical release.

    The release is a training authority, not merely a summary.  Validation
    therefore re-hashes every bound registry/native authority before returning
    the payload.  Optional expected bindings let a caller require exact
    agreement with authorities it has already loaded independently.
    """

    release_path = Path(path).expanduser().resolve(strict=True)
    if not release_path.is_file():
        raise ValueError("numeric-gate release path is not a regular file")
    payload = _strict_json_object(
        release_path,
        context="numeric-gate release",
    )
    if payload.get("schema_version") != MODEL_INPUT_NUMERIC_RELEASE_SCHEMA:
        raise ValueError("numeric-gate release schema_version mismatch")
    if payload.get("passed") is not True:
        raise ValueError("numeric-gate release did not pass")
    if (
        payload.get("model_input_contract_version")
        != MODEL_INPUT_CONTRACT_VERSION
    ):
        raise ValueError("numeric-gate model-input contract mismatch")
    if (
        payload.get("envelope_contract")
        != MODEL_INPUT_NUMERIC_ENVELOPE_CONTRACT
    ):
        raise ValueError("numeric-gate envelope contract mismatch")
    envelope = payload.get("envelope")
    if not isinstance(envelope, Mapping):
        raise ValueError("numeric-gate envelope content mismatch")
    if (
        payload.get("envelope_canonical_sha256")
        != MODEL_INPUT_NUMERIC_ENVELOPE_SHA256
        or _canonical_sha256(envelope) != MODEL_INPUT_NUMERIC_ENVELOPE_SHA256
    ):
        raise ValueError("numeric-gate envelope canonical SHA-256 mismatch")

    counts = payload.get("counts")
    if not isinstance(counts, Mapping):
        raise ValueError("numeric-gate release counts must be a mapping")
    expected_counts = {
        "full_observations": PRODUCTION_FULL_OBSERVATIONS,
        "eligible_observations": PRODUCTION_ELIGIBLE_OBSERVATIONS,
        "excluded_observations": PRODUCTION_EXCLUDED_OBSERVATIONS,
        "scanned_observations": PRODUCTION_ELIGIBLE_OBSERVATIONS,
        "failed_observations": 0,
        "native_shards": 112,
    }
    if set(counts) != set(expected_counts):
        raise ValueError("numeric-gate release count keys differ from contract")
    observed_counts = {
        name: _exact_int(counts, name) for name in expected_counts
    }
    if observed_counts != expected_counts:
        raise ValueError("numeric-gate release counts differ from production")

    identities = payload.get("identity_hashes")
    expected_identities = {
        "full": PRODUCTION_FULL_IDENTITY_SHA256,
        "eligible": PRODUCTION_ELIGIBLE_IDENTITY_SHA256,
        "excluded": PRODUCTION_EXCLUDED_IDENTITY_SHA256,
    }
    if not isinstance(identities, Mapping) or dict(identities) != (
        expected_identities
    ):
        raise ValueError(
            "numeric-gate release identities differ from production"
        )

    code_revision = str(payload.get("code_revision", "")).strip()
    if re.fullmatch(r"[0-9a-f]{40}", code_revision) is None:
        raise ValueError("numeric-gate code_revision is not a Git SHA")
    if (
        expected_code_revision is not None
        and code_revision != str(expected_code_revision).strip()
    ):
        raise ValueError("numeric-gate code_revision differs from expected")

    bindings = payload.get("authority_bindings")
    if not isinstance(bindings, Mapping) or set(bindings) != set(
        MODEL_INPUT_NUMERIC_AUTHORITY_NAMES
    ):
        raise ValueError(
            "numeric-gate release authority binding inventory is invalid"
        )
    normalized_bindings = {
        name: _normalized_file_binding(
            bindings[name],
            context=f"authority_bindings.{name}",
            verify_file=True,
        )
        for name in MODEL_INPUT_NUMERIC_AUTHORITY_NAMES
    }
    if expected_authority_bindings is not None:
        if set(expected_authority_bindings) != set(
            MODEL_INPUT_NUMERIC_AUTHORITY_NAMES
        ):
            raise ValueError("expected numeric authority inventory is invalid")
        normalized_expected = {
            name: _normalized_file_binding(
                expected_authority_bindings[name],
                context=f"expected_authority_bindings.{name}",
                verify_file=False,
            )
            for name in MODEL_INPUT_NUMERIC_AUTHORITY_NAMES
        }
        if normalized_bindings != normalized_expected:
            raise ValueError(
                "numeric-gate release authorities differ from expected"
            )
    _validate_numeric_evidence(
        payload,
        release_path=release_path,
        code_revision=code_revision,
        normalized_authorities=normalized_bindings,
    )
    return dict(payload)


__all__ = [
    "FLOAT32_SQUARE_SAFE_ABS_MAX",
    "FULL_POOL_NUMERIC_ENVELOPE_V1",
    "FullPoolNumericEnvelope",
    "MODEL_INPUT_CONTRACT_VERSION",
    "MODEL_INPUT_NUMERIC_ABS_MAX",
    "MODEL_INPUT_NUMERIC_AUDIT_SCHEMA",
    "MODEL_INPUT_NUMERIC_AUTHORITY_NAMES",
    "MODEL_INPUT_NUMERIC_ENVELOPE_CONTRACT",
    "MODEL_INPUT_NUMERIC_ENVELOPE_SHA256",
    "MODEL_INPUT_NUMERIC_RELEASE_SCHEMA",
    "TEACHER_SSL_NUMERIC_ENVELOPE_V1",
    "audit_collated_sample",
    "audit_model_facing_sample",
    "validate_numeric_gate_release",
]
