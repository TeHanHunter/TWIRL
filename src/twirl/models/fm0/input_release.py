"""Model-safe six-view input release for the frozen TWIRL-FM0.1 POC.

This module applies the frozen quality and detrending rules to an already
verified raw observation array, emits deterministic NPZ shards, and exposes
window tensors through an exact allowlist.  The strict-NPZ fixture and the
checksum-bound A2v1 HDF5 adapter deliberately remain separate input paths.
"""
from __future__ import annotations

from dataclasses import dataclass
from io import BytesIO
import csv
import hashlib
import json
from pathlib import Path
import re
import zipfile
from typing import Any, Mapping, Sequence

import numpy as np

from .registry import (
    FM0ContractError,
    FrozenContract,
    assert_no_search_columns,
    load_frozen_contract,
    publish_immutable,
    read_rows,
    sha256_file,
)


INPUT_RELEASE_SCHEMA_VERSION = "twirl_fm0_1_input_release_v1"
BUILD_SUMMARY_SCHEMA_VERSION = "twirl_fm0_1_input_release_build_summary_v1"
FIXTURE_ADAPTER_NAME = "strict_npz_fixture_v1"
FLUX_VIEW_NAMES = (
    "raw_relative_1x1",
    "raw_relative_3x3",
    "adp_1x1",
    "adp_3x3",
    "adp015_1x1",
    "adp015_3x3",
)
ERROR_VIEW_NAMES = (
    "raw_flux_error_1x1_over_robust_scale",
    "raw_flux_error_3x3_over_robust_scale",
)
RAW_ARRAY_KEYS = frozenset(
    {
        "time",
        "cadence",
        "orbit",
        "internal_quality",
        "spoc_quality",
        "qlp_quality",
        "authority_excluded",
        "raw_flux_1x1",
        "raw_flux_error_1x1",
        "raw_flux_3x3",
        "raw_flux_error_3x3",
    }
)
SHARD_ARRAY_KEYS = frozenset(
    {
        "flux",
        "flux_valid",
        "flux_error",
        "error_valid",
        "local_time_cadences",
        "delta_time_cadences",
        "time_valid",
        "segment_boundary",
        "segment_id",
        "view_present",
    }
)
FORWARD_TENSOR_KEYS = frozenset(
    {
        "flux",
        "flux_valid",
        "flux_error",
        "error_valid",
        "local_time_cadences",
        "delta_time_cadences",
        "time_valid",
        "segment_boundary",
        "view_present",
    }
)
MANIFEST_COLUMNS = (
    "input_release_schema_version",
    "observation_key",
    "product_instance_id",
    "source_sha256",
    "leakage_component_id",
    "source_partition",
    "product_state",
    "relative_path",
    "sha256",
    "input_source_sha256",
    "n_cadences",
    "n_segments",
    "view_present_json",
    "host_visit_offset_cadences",
    "host_visit_gap_cadences",
    "host_visit_overlaps_previous",
    "input_adapter",
    "scientific_training_eligible",
)
RAW_MANIFEST_COLUMNS = frozenset(
    {
        "observation_key",
        "product_instance_id",
        "source_sha256",
        "raw_npz_path",
        "raw_npz_sha256",
    }
)
WINDOW_LENGTH = 2048
EVALUATION_STRIDE = 1024
CADENCE_SECONDS = 200.0
GAP_DAYS = 0.2


@dataclass(frozen=True)
class ObservationRelease:
    flux: np.ndarray
    flux_valid: np.ndarray
    flux_error: np.ndarray
    error_valid: np.ndarray
    local_time_cadences: np.ndarray
    delta_time_cadences: np.ndarray
    time_valid: np.ndarray
    segment_boundary: np.ndarray
    segment_id: np.ndarray
    view_present: np.ndarray
    audit: Mapping[str, Any] | None = None

    @property
    def n_cadences(self) -> int:
        return int(self.flux.shape[0])

    @property
    def n_segments(self) -> int:
        return int(np.unique(self.segment_id).size) if self.segment_id.size else 0


@dataclass(frozen=True)
class WindowSpec:
    segment_id: int
    start_offset: int
    n_observed: int
    n_padded: int


def _as_1d(raw: Mapping[str, Any], key: str, dtype: Any) -> np.ndarray:
    value = np.asarray(raw[key], dtype=dtype)
    if value.ndim != 1:
        raise FM0ContractError(f"raw array {key!r} must be one-dimensional")
    return value


def _absolute_visit_bounds(time: np.ndarray) -> tuple[float, float]:
    """Return finite visit bounds after checking cadence-ordered absolute time.

    Individual non-finite timestamps remain permitted by the frozen
    time-validity-mask contract.  The represented visit must nevertheless have
    a finite start/end, and its finite timestamps must be strictly increasing
    in authoritative cadence order.
    """

    finite_time = np.asarray(time, dtype=np.float64)[np.isfinite(time)]
    if finite_time.size == 0:
        raise FM0ContractError("raw observation has no finite absolute time")
    if finite_time.size > 1 and np.any(np.diff(finite_time) <= 0):
        raise FM0ContractError(
            "finite absolute times must be strictly increasing in cadence order"
        )
    start = float(finite_time[0])
    end = float(finite_time[-1])
    if not np.isfinite(start) or not np.isfinite(end) or end < start:
        raise FM0ContractError("raw observation has invalid absolute visit bounds")
    return start, end


def _strict_raw_arrays(raw: Mapping[str, Any]) -> dict[str, np.ndarray]:
    keys = frozenset(raw)
    if keys != RAW_ARRAY_KEYS:
        missing = sorted(RAW_ARRAY_KEYS - keys)
        extra = sorted(keys - RAW_ARRAY_KEYS)
        raise FM0ContractError(f"raw observation keys mismatch; missing={missing}, extra={extra}")
    assert_no_search_columns(keys, context="raw observation arrays")
    typed = {
        "time": _as_1d(raw, "time", np.float64),
        "cadence": _as_1d(raw, "cadence", np.int64),
        "orbit": _as_1d(raw, "orbit", np.int64),
        "internal_quality": _as_1d(raw, "internal_quality", np.uint64),
        "spoc_quality": _as_1d(raw, "spoc_quality", np.uint64),
        "qlp_quality": _as_1d(raw, "qlp_quality", np.uint64),
        "authority_excluded": _as_1d(raw, "authority_excluded", np.bool_),
        "raw_flux_1x1": _as_1d(raw, "raw_flux_1x1", np.float64),
        "raw_flux_error_1x1": _as_1d(raw, "raw_flux_error_1x1", np.float64),
        "raw_flux_3x3": _as_1d(raw, "raw_flux_3x3", np.float64),
        "raw_flux_error_3x3": _as_1d(raw, "raw_flux_error_3x3", np.float64),
    }
    lengths = {array.size for array in typed.values()}
    if len(lengths) != 1 or not lengths or next(iter(lengths)) == 0:
        raise FM0ContractError("all raw arrays must have the same positive length")
    if np.unique(typed["cadence"]).size != typed["cadence"].size:
        raise FM0ContractError("duplicate cadence identifiers are not admissible")
    order = np.lexsort(
        (np.nan_to_num(typed["time"], nan=np.inf), typed["cadence"])
    )
    ordered = {key: value[order] for key, value in typed.items()}
    _absolute_visit_bounds(ordered["time"])
    return ordered


def _segments(time: np.ndarray, orbit: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = time.size
    boundary = np.zeros(n, dtype=np.bool_)
    boundary[0] = True
    last_finite_time: float | None = None
    previous_orbit: int | None = None
    for index in range(n):
        current_orbit = int(orbit[index])
        if index and previous_orbit is not None and current_orbit != previous_orbit:
            boundary[index] = True
        if np.isfinite(time[index]):
            value = float(time[index])
            if index and last_finite_time is not None and value - last_finite_time > GAP_DAYS:
                boundary[index] = True
            last_finite_time = value
        previous_orbit = current_orbit
    return np.cumsum(boundary, dtype=np.int64).astype(np.int32) - 1, boundary


def _derive_aperture(
    fit_time: np.ndarray,
    flux: np.ndarray,
    error: np.ndarray,
    effective_good: np.ndarray,
) -> dict[str, Any]:
    # Release readers and GPU training need only NumPy.  Keep SciPy-backed
    # detrending imports local to construction so an immutable finished shard
    # remains loadable on compute nodes without importing the build stack.
    from twirl.lightcurves.detrend_presets import adp015q_config, adp03q_config
    from twirl.lightcurves.flux_detrend import flux_space_detrend_result

    aperture_good = effective_good & np.isfinite(flux)
    error_good = aperture_good & np.isfinite(error) & (error > 0)
    if not np.any(aperture_good):
        raise FM0ContractError("aperture has no effective-good finite flux")
    median = float(np.median(flux[aperture_good]))
    fit_quality = (~error_good).astype(np.uint8)
    adp = flux_space_detrend_result(
        fit_time, flux, quality=fit_quality, flux_err=error, cfg=adp03q_config()
    )
    adp015 = flux_space_detrend_result(
        fit_time, flux, quality=fit_quality, flux_err=error, cfg=adp015q_config()
    )
    scale = float(adp.scale)
    if not np.isfinite(median) or not np.isfinite(scale) or scale <= 0:
        raise FM0ContractError(f"invalid aperture median/ADP scale: m={median}, s={scale}")
    raw_relative = (flux - median) / scale
    error_relative = error / scale
    values = (raw_relative, adp.det_flux - 1.0, adp015.det_flux - 1.0)
    validity = tuple(aperture_good & np.isfinite(value) for value in values)
    return {
        "values": values,
        "validity": validity,
        "error": error_relative,
        "error_good": error_good & np.isfinite(error_relative),
        "median": median,
        "scale": scale,
        "scale_source": adp.scale_source,
        "fit_count": adp.fit_count,
    }


def build_observation_release(
    raw: Mapping[str, Any],
    *,
    input_adapter: str = FIXTURE_ADAPTER_NAME,
    scientific_training_eligible: bool = False,
) -> ObservationRelease:
    """Apply the frozen quality, six-view, error, and timing contract."""

    arrays = _strict_raw_arrays(raw)
    time = arrays["time"]
    visit_start, visit_end = _absolute_visit_bounds(time)
    external_quality = arrays["spoc_quality"] | (arrays["qlp_quality"] << np.uint64(30))
    time_valid = np.isfinite(time)
    effective_good = (
        time_valid
        & (arrays["internal_quality"] == 0)
        & (external_quality == 0)
        & ~arrays["authority_excluded"]
    )
    segment_id, segment_boundary = _segments(time, arrays["orbit"])
    # ADP/ADP015 are time-domain detrenders.  Use the measured timestamps for
    # their knot spacing and gap segmentation, rather than an inferred uniform
    # cadence grid.  Non-finite timestamps are excluded from their fit and
    # remain false in the corresponding model-valid masks.
    fit_time = arrays["time"]
    small = _derive_aperture(
        fit_time,
        arrays["raw_flux_1x1"],
        arrays["raw_flux_error_1x1"],
        effective_good,
    )
    large = _derive_aperture(
        fit_time,
        arrays["raw_flux_3x3"],
        arrays["raw_flux_error_3x3"],
        effective_good,
    )
    # Frozen interleaving is raw(1,3), ADP(1,3), ADP015(1,3).
    flux_columns = (
        small["values"][0],
        large["values"][0],
        small["values"][1],
        large["values"][1],
        small["values"][2],
        large["values"][2],
    )
    valid_columns = (
        small["validity"][0],
        large["validity"][0],
        small["validity"][1],
        large["validity"][1],
        small["validity"][2],
        large["validity"][2],
    )
    flux_valid = np.column_stack(valid_columns).astype(np.bool_)
    flux = np.column_stack(flux_columns).astype(np.float32)
    flux[~flux_valid] = 0.0
    error_valid = np.column_stack((small["error_good"], large["error_good"])).astype(np.bool_)
    flux_error = np.column_stack((small["error"], large["error"])).astype(np.float32)
    flux_error[~error_valid] = 0.0

    local_time = np.zeros(time.size, dtype=np.float32)
    if np.any(time_valid):
        reference = float(time[np.flatnonzero(time_valid)[0]])
        local_time[time_valid] = ((time[time_valid] - reference) * 86400.0 / CADENCE_SECONDS).astype(
            np.float32
        )
    delta = np.zeros(time.size, dtype=np.float32)
    adjacent = time_valid[1:] & time_valid[:-1] & ~segment_boundary[1:]
    delta[1:][adjacent] = (
        (time[1:][adjacent] - time[:-1][adjacent]) * 86400.0 / CADENCE_SECONDS
    ).astype(np.float32)
    view_present = np.any(flux_valid, axis=0).astype(np.bool_)
    release = ObservationRelease(
        flux=flux,
        flux_valid=flux_valid,
        flux_error=flux_error,
        error_valid=error_valid,
        local_time_cadences=local_time,
        delta_time_cadences=delta,
        time_valid=time_valid.astype(np.bool_),
        segment_boundary=segment_boundary,
        segment_id=segment_id,
        view_present=view_present,
        audit={
            "external_quality_formula": "spoc_quality | (qlp_quality << 30)",
            "n_effective_good": int(np.sum(effective_good)),
            "median_1x1": small["median"],
            "median_3x3": large["median"],
            "scale_1x1": small["scale"],
            "scale_3x3": large["scale"],
            "scale_source_1x1": small["scale_source"],
            "scale_source_3x3": large["scale_source"],
            "absolute_visit_start": visit_start,
            "absolute_visit_end": visit_end,
            "input_adapter": str(input_adapter),
            "scientific_training_eligible": bool(scientific_training_eligible),
        },
    )
    validate_observation_release(release)
    return release


def validate_observation_release(release: ObservationRelease) -> None:
    """Fail closed on shape, dtype, mask, and neutral-fill violations."""

    n = release.n_cadences
    shapes = {
        "flux": (n, 6),
        "flux_valid": (n, 6),
        "flux_error": (n, 2),
        "error_valid": (n, 2),
        "local_time_cadences": (n,),
        "delta_time_cadences": (n,),
        "time_valid": (n,),
        "segment_boundary": (n,),
        "segment_id": (n,),
        "view_present": (6,),
    }
    for name, expected in shapes.items():
        actual = np.asarray(getattr(release, name)).shape
        if actual != expected:
            raise FM0ContractError(f"{name} shape {actual} != {expected}")
    for name in ("flux_valid", "error_valid", "time_valid", "segment_boundary", "view_present"):
        if np.asarray(getattr(release, name)).dtype != np.bool_:
            raise FM0ContractError(f"{name} must have boolean dtype")
    if release.flux.dtype != np.float32 or release.flux_error.dtype != np.float32:
        raise FM0ContractError("flux and flux_error must have float32 dtype")
    if release.local_time_cadences.dtype != np.float32 or release.delta_time_cadences.dtype != np.float32:
        raise FM0ContractError("stored timing arrays must have float32 dtype")
    if not np.all(np.isfinite(release.flux)) or not np.all(np.isfinite(release.flux_error)):
        raise FM0ContractError("release flux/error arrays must use finite neutral fill")
    if np.any(release.flux[~release.flux_valid] != 0) or np.any(
        release.flux_error[~release.error_valid] != 0
    ):
        raise FM0ContractError("invalid flux/error positions must be neutral zero")
    if n and (not release.segment_boundary[0] or release.segment_id[0] != 0):
        raise FM0ContractError("first cadence must begin segment zero")
    if n and not np.array_equal(
        release.segment_id, np.cumsum(release.segment_boundary, dtype=np.int64) - 1
    ):
        raise FM0ContractError("segment_id and segment_boundary disagree")
    if not np.array_equal(release.view_present, np.any(release.flux_valid, axis=0)):
        raise FM0ContractError("view_present disagrees with flux-valid masks")


def _release_arrays(release: ObservationRelease) -> dict[str, np.ndarray]:
    return {name: np.asarray(getattr(release, name)) for name in sorted(SHARD_ARRAY_KEYS)}


def deterministic_npz_bytes(release: ObservationRelease) -> bytes:
    """Serialize a shard with sorted entries and fixed ZIP metadata."""

    validate_observation_release(release)
    output = BytesIO()
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_STORED) as archive:
        for name, array in _release_arrays(release).items():
            array_bytes = BytesIO()
            np.lib.format.write_array(array_bytes, np.asarray(array), allow_pickle=False)
            info = zipfile.ZipInfo(f"{name}.npy", date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_STORED
            info.create_system = 3
            info.external_attr = 0o600 << 16
            archive.writestr(info, array_bytes.getvalue())
    return output.getvalue()


def load_input_release(path: str | Path) -> ObservationRelease:
    """Load and validate one deterministic release shard."""

    with np.load(Path(path), allow_pickle=False) as archive:
        keys = frozenset(archive.files)
        if keys != SHARD_ARRAY_KEYS:
            raise FM0ContractError(
                f"release shard keys mismatch; missing={sorted(SHARD_ARRAY_KEYS-keys)}, "
                f"extra={sorted(keys-SHARD_ARRAY_KEYS)}"
            )
        release = ObservationRelease(**{name: archive[name] for name in SHARD_ARRAY_KEYS})
    validate_observation_release(release)
    return release


def evaluation_windows(release: ObservationRelease) -> tuple[WindowSpec, ...]:
    """Return deterministic stride-1024 windows, never crossing a segment."""

    specs: list[WindowSpec] = []
    for segment in np.unique(release.segment_id):
        size = int(np.sum(release.segment_id == segment))
        for start in range(0, size, EVALUATION_STRIDE):
            observed = min(WINDOW_LENGTH, size - start)
            specs.append(WindowSpec(int(segment), start, observed, WINDOW_LENGTH - observed))
    return tuple(specs)


def deterministic_training_window(
    release: ObservationRelease,
    *,
    observation_key: str,
    epoch: int,
    draw_index: int,
) -> WindowSpec:
    """Choose a reproducible segment and random start without data-dependent scores."""

    segments = np.unique(release.segment_id)
    if not segments.size:
        raise FM0ContractError("cannot sample an empty observation")
    seed = hashlib.sha256(
        f"twirl_fm0_1_window_v1:{observation_key}:{epoch}:{draw_index}".encode("utf-8")
    ).digest()
    value = int.from_bytes(seed, "big")
    segment = int(segments[value % segments.size])
    size = int(np.sum(release.segment_id == segment))
    maximum = max(0, size - WINDOW_LENGTH)
    start = (value // max(1, segments.size)) % (maximum + 1)
    observed = min(WINDOW_LENGTH, size - start)
    return WindowSpec(segment, int(start), observed, WINDOW_LENGTH - observed)


def extract_window(
    release: ObservationRelease, *, segment_id: int, start_offset: int
) -> dict[str, np.ndarray]:
    """Extract one padded, model-visible window through the exact allowlist."""

    indices = np.flatnonzero(release.segment_id == int(segment_id))
    if start_offset < 0 or start_offset >= indices.size:
        raise FM0ContractError("window start is outside the requested segment")
    selected = indices[start_offset : start_offset + WINDOW_LENGTH]
    count = selected.size

    def padded(source: np.ndarray, shape: tuple[int, ...], dtype: Any) -> np.ndarray:
        result = np.zeros(shape, dtype=dtype)
        result[:count] = source[selected]
        return result

    time_valid = padded(release.time_valid, (WINDOW_LENGTH,), np.bool_)
    stored_local_time = padded(release.local_time_cadences, (WINDOW_LENGTH,), np.float32)
    local_time = np.zeros(WINDOW_LENGTH, dtype=np.float32)
    valid_indices = np.flatnonzero(time_valid)
    if valid_indices.size:
        reference = stored_local_time[valid_indices[0]]
        local_time[time_valid] = stored_local_time[time_valid] - reference
    tensors = {
        "flux": padded(release.flux, (WINDOW_LENGTH, 6), np.float32),
        "flux_valid": padded(release.flux_valid, (WINDOW_LENGTH, 6), np.bool_),
        "flux_error": padded(release.flux_error, (WINDOW_LENGTH, 2), np.float32),
        "error_valid": padded(release.error_valid, (WINDOW_LENGTH, 2), np.bool_),
        "local_time_cadences": local_time,
        "delta_time_cadences": padded(
            release.delta_time_cadences, (WINDOW_LENGTH,), np.float32
        ),
        "time_valid": time_valid,
        "segment_boundary": padded(release.segment_boundary, (WINDOW_LENGTH,), np.bool_),
        "view_present": release.view_present.copy(),
    }
    validate_model_tensors(tensors)
    return tensors


def validate_model_tensors(tensors: Mapping[str, Any]) -> None:
    """Enforce that no identity, detector, BLS, or audit field reaches forward."""

    keys = frozenset(tensors)
    if keys != FORWARD_TENSOR_KEYS:
        raise FM0ContractError(
            f"forward tensor allowlist mismatch; missing={sorted(FORWARD_TENSOR_KEYS-keys)}, "
            f"extra={sorted(keys-FORWARD_TENSOR_KEYS)}"
        )
    assert_no_search_columns(keys, context="forward tensors")
    expected = {
        "flux": (WINDOW_LENGTH, 6),
        "flux_valid": (WINDOW_LENGTH, 6),
        "flux_error": (WINDOW_LENGTH, 2),
        "error_valid": (WINDOW_LENGTH, 2),
        "local_time_cadences": (WINDOW_LENGTH,),
        "delta_time_cadences": (WINDOW_LENGTH,),
        "time_valid": (WINDOW_LENGTH,),
        "segment_boundary": (WINDOW_LENGTH,),
        "view_present": (6,),
    }
    for name, shape in expected.items():
        if np.asarray(tensors[name]).shape != shape:
            raise FM0ContractError(f"forward tensor {name} has wrong shape")


def _csv_bytes(rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> bytes:
    from io import StringIO

    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(fields), lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row.get(field, "") for field in fields})
    return stream.getvalue().encode("utf-8")


def _raw_npz(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as archive:
        return {name: archive[name] for name in archive.files}


def write_input_release(
    *,
    registry_dir: str | Path,
    raw_manifest_path: str | Path,
    out_dir: str | Path,
    contract: FrozenContract | None = None,
) -> dict[str, Any]:
    """Build an immutable, non-scientific fixture release from admitted rows.

    This build summary is not an ORCD admission receipt and does not certify
    that the full S56--S67 campaign is present.  The strict-NPZ adapter exists
    for contract and numerical tests only; its outputs are explicitly
    ineligible for scientific training until a separately reviewed immutable
    A2v1 HDF5 adapter exists.
    """

    contract = contract or load_frozen_contract()
    registry_path = Path(registry_dir)
    observations_path = registry_path / "observations.csv"
    if not observations_path.is_file():
        raise FM0ContractError("registry release has no observations.csv")
    observations = {row["observation_key"]: row for row in read_rows(observations_path)}
    raw_manifest = read_rows(raw_manifest_path)
    if not raw_manifest:
        raise FM0ContractError("raw manifest is empty")
    assert_no_search_columns(raw_manifest[0].keys(), context="raw manifest")
    output = Path(out_dir)
    manifest_rows: list[dict[str, Any]] = []
    shard_payloads: dict[str, bytes] = {}
    seen: set[str] = set()
    raw_hash_components: dict[str, str] = {}
    visit_timing: list[tuple[dict[str, Any], str, float, float]] = []
    for raw_row in raw_manifest:
        assert_no_search_columns(raw_row.keys(), context="raw manifest")
        if frozenset(raw_row) != RAW_MANIFEST_COLUMNS:
            missing = sorted(RAW_MANIFEST_COLUMNS - frozenset(raw_row))
            extra = sorted(frozenset(raw_row) - RAW_MANIFEST_COLUMNS)
            raise FM0ContractError(
                f"raw manifest columns mismatch; missing={missing}, extra={extra}"
            )
        key = str(raw_row.get("observation_key", "")).strip()
        if not key or key in seen:
            raise FM0ContractError(f"missing or duplicate observation_key: {key!r}")
        seen.add(key)
        observation = observations.get(key)
        if observation is None:
            raise FM0ContractError(f"raw manifest observation is absent from registry: {key}")
        if str(observation.get("quarantined", "")).lower() in {"true", "1"}:
            raise FM0ContractError(f"quarantined observation cannot enter release: {key}")
        manifest_product_id = str(raw_row.get("product_instance_id", "")).strip()
        if manifest_product_id != str(observation["product_instance_id"]):
            raise FM0ContractError(
                f"raw manifest product_instance_id does not match registry for {key}"
            )
        manifest_source_hash = str(raw_row.get("source_sha256", "")).strip().lower()
        if not re.fullmatch(r"[0-9a-f]{64}", manifest_source_hash):
            raise FM0ContractError(
                f"raw manifest source_sha256 is not a lowercase SHA-256 for {key}"
            )
        if manifest_source_hash != str(observation["source_sha256"]):
            raise FM0ContractError(
                f"raw manifest source_sha256 does not match registry for {key}"
            )
        raw_path = Path(str(raw_row.get("raw_npz_path", ""))).expanduser()
        if not raw_path.is_absolute():
            raw_path = Path(raw_manifest_path).resolve().parent / raw_path
        expected_hash = str(raw_row.get("raw_npz_sha256", "")).strip().lower()
        if not re.fullmatch(r"[0-9a-f]{64}", expected_hash):
            raise FM0ContractError(f"invalid raw_npz_sha256 for {key}")
        if not raw_path.is_file() or sha256_file(raw_path) != expected_hash:
            raise FM0ContractError(f"raw NPZ missing or hash mismatch for {key}")
        component = str(observation["leakage_component_id"])
        previous_component = raw_hash_components.get(expected_hash)
        if previous_component is not None and previous_component != component:
            raise FM0ContractError(
                "identical raw NPZ hashes cannot represent different leakage "
                f"components: {previous_component} and {component}"
            )
        raw_hash_components[expected_hash] = component
        release = build_observation_release(_raw_npz(raw_path))
        if release.audit is None:
            raise FM0ContractError(f"fixture adapter omitted timing audit for {key}")
        visit_start = float(release.audit["absolute_visit_start"])
        visit_end = float(release.audit["absolute_visit_end"])
        if not np.isfinite(visit_start) or not np.isfinite(visit_end):
            raise FM0ContractError(f"non-finite visit timing for {key}")
        relative_path = f"shards/{key}.npz"
        payload = deterministic_npz_bytes(release)
        digest = hashlib.sha256(payload).hexdigest()
        shard_payloads[relative_path] = payload
        manifest_rows.append(
            {
                "input_release_schema_version": INPUT_RELEASE_SCHEMA_VERSION,
                "observation_key": key,
                "product_instance_id": observation["product_instance_id"],
                "source_sha256": observation["source_sha256"],
                "leakage_component_id": observation["leakage_component_id"],
                "source_partition": observation["source_partition"],
                "product_state": observation["product_state"],
                "relative_path": relative_path,
                "sha256": digest,
                "input_source_sha256": expected_hash,
                "n_cadences": release.n_cadences,
                "n_segments": release.n_segments,
                "view_present_json": json.dumps(
                    release.view_present.astype(int).tolist(), separators=(",", ":")
                ),
                "host_visit_offset_cadences": "",
                "host_visit_gap_cadences": "",
                "host_visit_overlaps_previous": "",
                "input_adapter": FIXTURE_ADAPTER_NAME,
                "scientific_training_eligible": False,
            }
        )
        visit_timing.append(
            (
                manifest_rows[-1],
                str(observation["physical_source_id"]),
                visit_start,
                visit_end,
            )
        )

    visits_by_source: dict[str, list[tuple[dict[str, Any], float, float]]] = {}
    for row, physical_source_id, visit_start, visit_end in visit_timing:
        visits_by_source.setdefault(physical_source_id, []).append(
            (row, visit_start, visit_end)
        )
    cadence_units_per_day = 86400.0 / CADENCE_SECONDS
    for physical_source_id, visits in visits_by_source.items():
        visits.sort(key=lambda item: (item[1], str(item[0]["observation_key"])))
        first_start = visits[0][1]
        previous_start: float | None = None
        previous_end: float | None = None
        for row, visit_start, visit_end in visits:
            if previous_start is not None and visit_start <= previous_start:
                raise FM0ContractError(
                    f"visit starts are not strictly monotonic for {physical_source_id}"
                )
            offset = (visit_start - first_start) * cadence_units_per_day
            gap = (
                0.0
                if previous_end is None
                else (visit_start - previous_end) * cadence_units_per_day
            )
            if not np.isfinite(offset) or not np.isfinite(gap) or offset < 0:
                raise FM0ContractError(
                    f"invalid derived visit offset/gap for {physical_source_id}"
                )
            row["host_visit_offset_cadences"] = float(offset)
            row["host_visit_gap_cadences"] = float(gap)
            row["host_visit_overlaps_previous"] = bool(
                previous_end is not None and gap < 0
            )
            previous_start = visit_start
            previous_end = visit_end
    manifest_rows.sort(key=lambda row: row["observation_key"])
    manifest_payload = _csv_bytes(manifest_rows, MANIFEST_COLUMNS)
    summary = {
        "summary_schema_version": BUILD_SUMMARY_SCHEMA_VERSION,
        "input_release_schema_version": INPUT_RELEASE_SCHEMA_VERSION,
        "campaign_id": contract.config["campaign_id"],
        "design_sha256": contract.design_sha256,
        "config_sha256": contract.config_sha256,
        "freeze_receipt_sha256": contract.freeze_receipt_sha256,
        "registry_observations_sha256": sha256_file(observations_path),
        "raw_manifest_sha256": sha256_file(raw_manifest_path),
        "manifest_sha256": hashlib.sha256(manifest_payload).hexdigest(),
        "n_observations": len(manifest_rows),
        "n_cadences": sum(int(row["n_cadences"]) for row in manifest_rows),
        "flux_view_names": list(FLUX_VIEW_NAMES),
        "error_view_names": list(ERROR_VIEW_NAMES),
        "input_adapter": FIXTURE_ADAPTER_NAME,
        "scientific_training_eligible": False,
        "host_visit_timing_derivation": (
            "absolute_raw_time_grouped_by_physical_source_id"
        ),
        "host_visit_gap_definition": "signed_current_start_minus_previous_end",
        "partial_release": len(manifest_rows) != len(observations),
        "certifies_full_campaign": False,
    }
    summary_payload = (json.dumps(summary, indent=2, sort_keys=True) + "\n").encode("utf-8")
    for relative_path, payload in shard_payloads.items():
        publish_immutable(output / relative_path, payload)
    publish_immutable(output / "manifest.csv", manifest_payload)
    publish_immutable(output / "summary.json", summary_payload)
    return summary


def write_a2v1_hdf5_input_release(
    *,
    registry_dir: str | Path,
    hdf5_manifest_path: str | Path,
    out_dir: str | Path,
    contract: FrozenContract | None = None,
) -> dict[str, Any]:
    """Build a scientific-eligible release from bound A2v1 HDF5 inputs.

    This function is the only FM0.1 path that may mark an input release
    eligible for scientific training.  It still certifies *only* the named
    source products, not survey completeness, a trained checkpoint, or a
    foundation-model claim.
    """

    from .a2v1_adapter import (
        A2V1_HDF5_ADAPTER_NAME,
        A2V1_HDF5_MANIFEST_COLUMNS,
        A2V1AdapterCache,
        load_a2v1_hdf5_observation,
    )

    contract = contract or load_frozen_contract()
    registry_path = Path(registry_dir)
    observations_path = registry_path / "observations.csv"
    if not observations_path.is_file():
        raise FM0ContractError("registry release has no observations.csv")
    observations = {row["observation_key"]: row for row in read_rows(observations_path)}
    hdf5_manifest = read_rows(hdf5_manifest_path)
    if not hdf5_manifest:
        raise FM0ContractError("A2v1 HDF5 manifest is empty")
    assert_no_search_columns(hdf5_manifest[0].keys(), context="A2v1 HDF5 manifest")
    output = Path(out_dir)
    manifest_rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    source_components: dict[str, str] = {}
    visit_timing: list[tuple[dict[str, Any], str, float, float]] = []
    hdf5_manifest_path = Path(hdf5_manifest_path)
    adapter_cache = A2V1AdapterCache()
    for source_row in hdf5_manifest:
        assert_no_search_columns(source_row.keys(), context="A2v1 HDF5 manifest")
        if frozenset(source_row) != A2V1_HDF5_MANIFEST_COLUMNS:
            missing = sorted(A2V1_HDF5_MANIFEST_COLUMNS - frozenset(source_row))
            extra = sorted(frozenset(source_row) - A2V1_HDF5_MANIFEST_COLUMNS)
            raise FM0ContractError(
                f"A2v1 HDF5 manifest columns mismatch; missing={missing}, extra={extra}"
            )
        key = str(source_row["observation_key"]).strip()
        if not key or key in seen:
            raise FM0ContractError(f"missing or duplicate observation_key: {key!r}")
        seen.add(key)
        observation = observations.get(key)
        if observation is None:
            raise FM0ContractError(f"A2v1 source is absent from registry: {key}")
        if str(observation.get("quarantined", "")).lower() in {"true", "1"}:
            raise FM0ContractError(f"quarantined observation cannot enter release: {key}")
        if str(source_row["product_instance_id"]).strip() != str(
            observation["product_instance_id"]
        ):
            raise FM0ContractError(
                f"A2v1 source product_instance_id does not match registry for {key}"
            )
        if str(source_row["source_sha256"]).strip().lower() != str(
            observation["source_sha256"]
        ):
            raise FM0ContractError(
                f"A2v1 source source_sha256 does not match registry for {key}"
            )
        source = load_a2v1_hdf5_observation(
            source_row,
            manifest_dir=hdf5_manifest_path.resolve().parent,
            sector=int(observation["sector"]),
            tic_id=str(observation["tic_id"]),
            cache=adapter_cache,
        )
        if source.source_sha256 != str(observation["source_sha256"]):
            raise FM0ContractError(f"A2v1 source identity does not match registry for {key}")
        component = str(observation["leakage_component_id"])
        prior_component = source_components.get(source.source_sha256)
        if prior_component is not None and prior_component != component:
            raise FM0ContractError(
                "identical A2v1 source bindings cannot cross leakage components: "
                f"{prior_component} and {component}"
            )
        source_components[source.source_sha256] = component
        release = build_observation_release(
            source.raw_arrays,
            input_adapter=A2V1_HDF5_ADAPTER_NAME,
            scientific_training_eligible=True,
        )
        if release.audit is None:
            raise FM0ContractError(f"A2v1 adapter omitted visit timing for {key}")
        visit_start = float(release.audit["absolute_visit_start"])
        visit_end = float(release.audit["absolute_visit_end"])
        relative_path = f"shards/{key}.npz"
        payload = deterministic_npz_bytes(release)
        # A full multi-sector release can contain hundreds of thousands of
        # observations.  Publish each immutable shard immediately instead of
        # retaining every encoded payload in RAM until the final manifest is
        # assembled.  A failed build therefore leaves only checksum-safe,
        # idempotently reusable shards; without manifest.csv it is not a
        # consumable release.
        publish_immutable(output / relative_path, payload)
        manifest_rows.append(
            {
                "input_release_schema_version": INPUT_RELEASE_SCHEMA_VERSION,
                "observation_key": key,
                "product_instance_id": observation["product_instance_id"],
                "source_sha256": observation["source_sha256"],
                "leakage_component_id": observation["leakage_component_id"],
                "source_partition": observation["source_partition"],
                "product_state": observation["product_state"],
                "relative_path": relative_path,
                "sha256": hashlib.sha256(payload).hexdigest(),
                "input_source_sha256": source.source_sha256,
                "n_cadences": release.n_cadences,
                "n_segments": release.n_segments,
                "view_present_json": json.dumps(
                    release.view_present.astype(int).tolist(), separators=(",", ":")
                ),
                "host_visit_offset_cadences": "",
                "host_visit_gap_cadences": "",
                "host_visit_overlaps_previous": "",
                "input_adapter": A2V1_HDF5_ADAPTER_NAME,
                "scientific_training_eligible": True,
            }
        )
        visit_timing.append(
            (
                manifest_rows[-1],
                str(observation["physical_source_id"]),
                visit_start,
                visit_end,
            )
        )
        if len(manifest_rows) % 1000 == 0:
            print(
                "FM0_INPUT_PROGRESS "
                f"observations={len(manifest_rows)} "
                f"cadences={sum(int(row['n_cadences']) for row in manifest_rows)}",
                flush=True,
            )

    try:
        adapter_cache.assert_unchanged()
    except Exception as exc:
        raise FM0ContractError(
            "an external-quality authority changed while building the release"
        ) from exc

    visits_by_source: dict[str, list[tuple[dict[str, Any], float, float]]] = {}
    for row, physical_source_id, visit_start, visit_end in visit_timing:
        visits_by_source.setdefault(physical_source_id, []).append(
            (row, visit_start, visit_end)
        )
    cadence_units_per_day = 86400.0 / CADENCE_SECONDS
    for physical_source_id, visits in visits_by_source.items():
        visits.sort(key=lambda item: (item[1], str(item[0]["observation_key"])))
        first_start = visits[0][1]
        previous_start: float | None = None
        previous_end: float | None = None
        for row, visit_start, visit_end in visits:
            if previous_start is not None and visit_start <= previous_start:
                raise FM0ContractError(
                    f"visit starts are not strictly monotonic for {physical_source_id}"
                )
            offset = (visit_start - first_start) * cadence_units_per_day
            gap = 0.0 if previous_end is None else (visit_start - previous_end) * cadence_units_per_day
            if not np.isfinite(offset) or not np.isfinite(gap) or offset < 0:
                raise FM0ContractError(
                    f"invalid derived visit offset/gap for {physical_source_id}"
                )
            row["host_visit_offset_cadences"] = float(offset)
            row["host_visit_gap_cadences"] = float(gap)
            row["host_visit_overlaps_previous"] = bool(
                previous_end is not None and gap < 0
            )
            previous_start = visit_start
            previous_end = visit_end

    manifest_rows.sort(key=lambda row: row["observation_key"])
    manifest_payload = _csv_bytes(manifest_rows, MANIFEST_COLUMNS)
    summary = {
        "summary_schema_version": BUILD_SUMMARY_SCHEMA_VERSION,
        "input_release_schema_version": INPUT_RELEASE_SCHEMA_VERSION,
        "campaign_id": contract.config["campaign_id"],
        "design_sha256": contract.design_sha256,
        "config_sha256": contract.config_sha256,
        "freeze_receipt_sha256": contract.freeze_receipt_sha256,
        "registry_observations_sha256": sha256_file(observations_path),
        "hdf5_manifest_sha256": sha256_file(hdf5_manifest_path),
        "manifest_sha256": hashlib.sha256(manifest_payload).hexdigest(),
        "n_observations": len(manifest_rows),
        "n_cadences": sum(int(row["n_cadences"]) for row in manifest_rows),
        "flux_view_names": list(FLUX_VIEW_NAMES),
        "error_view_names": list(ERROR_VIEW_NAMES),
        "input_adapter": A2V1_HDF5_ADAPTER_NAME,
        "scientific_training_eligible": True,
        "host_visit_timing_derivation": "absolute_raw_time_grouped_by_physical_source_id",
        "host_visit_gap_definition": "signed_current_start_minus_previous_end",
        "partial_release": len(manifest_rows) != len(observations),
        "certifies_full_campaign": False,
    }
    summary_payload = (json.dumps(summary, indent=2, sort_keys=True) + "\n").encode("utf-8")
    publish_immutable(output / "manifest.csv", manifest_payload)
    publish_immutable(output / "summary.json", summary_payload)
    return summary
