"""Payload-screened crop schedule for the matched TWIRL-FM0.3 canary.

This is the deliberately small second stage after the immutable identity-only
matched-canary plan.  It opens exactly the 1,440 shards named by that plan,
checks each source checksum, and freezes one deterministic, unpadded
128-cadence crop per observation.  It does not inject an event, load a
checkpoint, fit a probe, compute a metric, or train a model.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import shutil
import tempfile
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from io import StringIO
from pathlib import Path, PurePosixPath
from typing import Any

import numpy as np

from .centered_event_context_diagnostic import centered_trapezoid
from .composite_release import SOURCE_BINDING_FIELDS
from .input_release import (
    SHARD_ARRAY_KEYS,
    ObservationRelease,
    load_input_release_bytes,
)
from .later_sector_release import LATER_SIX_VIEW_MANIFEST_FIELDS
from .matched_canary_plan import (
    ARCHITECTURES,
    CONTEXT_CADENCES,
    EVENT_CENTER_INDEX_ZERO_BASED,
    EVENT_DEPTHS_TEXT,
    EVENT_DURATIONS,
    MATCHED_CANARY_PLAN_SCHEMA_VERSION,
    REQUIRED_VIEW_INDICES,
    SCHEDULE_FIELDS,
    SPLIT_ORDER,
    TARGET_COMPONENTS_PER_COHORT,
    TEMPORAL_COHORTS,
    _load_composite_identity_authority,
    validate_matched_canary_evaluation_plan,
)
from .registry import FM0ContractError, sha256_file
from .temporal_zero_shot import load_temporal_panel

PAYLOAD_PLAN_SCHEMA_VERSION = "twirl_fm0_3_matched_canary_payload_plan_v1"
PAYLOAD_PLAN_RECEIPT_SCHEMA_VERSION = (
    "twirl_fm0_3_matched_canary_payload_plan_receipt_v1"
)
PAYLOAD_PLAN_READY_STATE = "FM0_3_MATCHED_CANARY_PAYLOAD_PLAN_READY"
PAYLOAD_PLAN_SALT = "twirl_fm0_3_matched_canary_payload_plan_v1"
CROP_PAYLOAD_HASH_SCHEMA = "twirl_fm0_3_crop_payload_sha256_v1"

EXPECTED_ROWS = (
    len(SPLIT_ORDER) * len(TEMPORAL_COHORTS) * TARGET_COMPONENTS_PER_COHORT
)
MINIMUM_VALID_FRACTION = 0.8
MINIMUM_VALID_CADENCES = math.ceil(MINIMUM_VALID_FRACTION * CONTEXT_CADENCES)
CENTER_SUPPORT_START = 60
CENTER_SUPPORT_STOP_EXCLUSIVE = 69
PROBE_EPOCHS = 400
PROBE_LEARNING_RATE = 0.02
PROBE_L2 = 0.001
PROBE_SEED = 560306
BOOTSTRAP_REPLICATES = 1000
BOOTSTRAP_SEED = 560305

RAW_FEATURE_CHANNELS = (
    "adp_1x1_flux",
    "adp_3x3_flux",
    "adp_1x1_flux_valid",
    "adp_3x3_flux_valid",
    "raw_flux_error_1x1_over_robust_scale",
    "raw_flux_error_3x3_over_robust_scale",
    "error_valid_1x1",
    "error_valid_3x3",
    "local_time_cadences",
    "delta_time_cadences",
    "time_valid",
    "segment_boundary",
)
QUALITY_FEATURE_CHANNELS = (
    "adp_1x1_flux_valid",
    "adp_3x3_flux_valid",
    "error_valid_1x1",
    "error_valid_3x3",
    "local_time_cadences",
    "delta_time_cadences",
    "time_valid",
    "segment_boundary",
)

PAYLOAD_SCHEDULE_FIELDS = (
    "payload_schema_version",
    *SCHEDULE_FIELDS,
    "identity_row_sha256",
    "source_release_root",
    "source_relative_path",
    "source_shard_sha256",
    "source_n_cadences",
    "source_n_segments",
    "crop_segment_id",
    "segment_start_index_zero_based",
    "segment_stop_index_exclusive",
    "crop_start_index_zero_based",
    "crop_stop_index_exclusive",
    "crop_start_offset_within_segment",
    "crop_stop_offset_within_segment_exclusive",
    "n_time_valid",
    "n_adp_1x1_valid",
    "n_adp_3x3_valid",
    "n_joint_valid",
    "joint_valid_fraction",
    "center_support_start_index_zero_based",
    "center_support_stop_index_exclusive",
    "n_center_support_joint_valid",
    "n_eligible_crops",
    "selected_eligible_crop_rank_zero_based",
    "crop_payload_sha256",
)

_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_GIT_SHA = re.compile(r"^[0-9a-f]{40}$")


@dataclass(frozen=True, slots=True)
class SourceBinding:
    """One identity-plan observation bound to one immutable NPZ shard."""

    release_root: Path
    relative_path: str
    shard_sha256: str
    n_cadences: int
    n_segments: int
    view_present: tuple[bool, ...]


@dataclass(frozen=True, slots=True)
class CropChoice:
    """One exact eligible crop within one source segment."""

    segment_id: int
    segment_start: int
    segment_stop: int
    crop_start: int
    crop_stop: int
    crop_start_offset: int
    crop_stop_offset: int
    time_valid: int
    adp_1x1_valid: int
    adp_3x3_valid: int
    joint_valid: int
    eligible_count: int
    selected_rank: int


@dataclass(frozen=True, slots=True)
class MatchedCanaryPayloadPlanResult:
    root: Path
    schedule_path: Path
    receipt_path: Path
    ready_path: Path
    schedule_sha256: str
    receipt_sha256: str
    receipt: Mapping[str, Any]


def _digest(value: Any, *, label: str) -> str:
    digest = str(value).strip().lower()
    if _SHA256.fullmatch(digest) is None:
        raise FM0ContractError(f"{label} must be a lowercase SHA-256")
    return digest


def _git_sha(value: Any) -> str:
    sha = str(value).strip().lower()
    if _GIT_SHA.fullmatch(sha) is None:
        raise FM0ContractError("producer_git_sha must be a full lowercase Git SHA")
    return sha


def _canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _stable_index(count: int, *parts: object) -> int:
    if count <= 0:
        raise FM0ContractError("cannot select from an empty crop set")
    digest = hashlib.sha256(
        "\x1f".join((PAYLOAD_PLAN_SALT, *(str(part) for part in parts))).encode()
    ).digest()
    return int.from_bytes(digest[:8], "big", signed=False) % count


def _read_exact_csv(
    path: Path, fields: Sequence[str], *, label: str
) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise FM0ContractError(f"{label} columns drifted")
        rows = [dict(row) for row in reader]
    if not rows:
        raise FM0ContractError(f"{label} is empty")
    return rows


def _identity_rows(path: Path) -> list[dict[str, str]]:
    return _read_exact_csv(path, SCHEDULE_FIELDS, label="identity schedule")


def _safe_relative_shard(value: Any, *, observation_key: str) -> str:
    text = str(value).strip()
    relative = PurePosixPath(text)
    if (
        not text
        or relative.is_absolute()
        or ".." in relative.parts
        or relative.parts != ("shards", f"{observation_key}.npz")
    ):
        raise FM0ContractError("source shard path escaped its immutable release")
    return relative.as_posix()


def _source_binding_from_row(
    row: Mapping[str, str], *, release_root: str | Path, require_read_only: bool
) -> SourceBinding:
    root_raw = Path(release_root).expanduser()
    if root_raw.is_symlink():
        raise FM0ContractError("source release root must be materialized")
    root = root_raw.resolve(strict=True)
    if not root.is_dir():
        raise FM0ContractError("source release root is not a directory")
    if require_read_only and root.stat().st_mode & 0o222:
        raise FM0ContractError("source release root must be read-only")
    try:
        n_cadences = int(row["n_cadences"])
        n_segments = int(row["n_segments"])
    except (KeyError, TypeError, ValueError) as exc:
        raise FM0ContractError("source shard dimensions are invalid") from exc
    if n_cadences <= 0 or n_segments <= 0:
        raise FM0ContractError("source shard dimensions are invalid")
    try:
        present_raw = json.loads(row["view_present_json"])
    except (KeyError, TypeError, json.JSONDecodeError) as exc:
        raise FM0ContractError("source view-presence metadata is invalid") from exc
    if (
        not isinstance(present_raw, list)
        or len(present_raw) != 6
        or any(type(value) is not int or value not in (0, 1) for value in present_raw)
    ):
        raise FM0ContractError("source view-presence metadata is invalid")
    return SourceBinding(
        release_root=root,
        relative_path=_safe_relative_shard(
            row["relative_path"], observation_key=row["observation_key"]
        ),
        shard_sha256=_digest(row["sha256"], label="source shard hash"),
        n_cadences=n_cadences,
        n_segments=n_segments,
        view_present=tuple(bool(value) for value in present_raw),
    )


def _require_authority(value: Any, *, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise FM0ContractError(f"{label} authority is invalid")
    return value


def _resolve_sources(
    identity_rows: Sequence[Mapping[str, str]],
    *,
    identity_receipt: Mapping[str, Any],
    require_read_only: bool,
) -> dict[str, SourceBinding]:
    authorities = _require_authority(
        identity_receipt.get("input_authorities"), label="identity-plan input"
    )
    panel = _require_authority(
        authorities.get("temporal_panel"), label="temporal-panel"
    )
    composite = _require_authority(
        authorities.get("composite_release"), label="composite"
    )

    try:
        panel_rows, panel_loaded = load_temporal_panel(
            panel["root"], receipt_sha256=panel["receipt_sha256"]
        )
    except (KeyError, OSError, TypeError, ValueError) as exc:
        raise FM0ContractError("bound temporal-panel authority failed validation") from exc
    if (
        panel_loaded["panel_sha256"] != panel.get("panel_sha256")
        or panel_loaded["sector_bindings_sha256"]
        != panel.get("sector_bindings_sha256")
    ):
        raise FM0ContractError("temporal-panel authority hashes drifted")
    panel_by_observation = {row["observation_key"]: row for row in panel_rows}

    try:
        composite_loaded = _load_composite_identity_authority(
            composite["root"],
            receipt_sha256=composite["receipt_sha256"],
            source_bindings_sha256=composite["source_bindings_sha256"],
            role_index_sha256=composite["role_index_sha256"],
            require_read_only=require_read_only,
        )
    except (KeyError, OSError, TypeError, ValueError) as exc:
        raise FM0ContractError("bound composite authority failed validation") from exc

    composite_root = Path(composite["root"]).expanduser().resolve(strict=True)
    source_bindings_path = composite_root / "source_bindings.csv"
    bindings = _read_exact_csv(
        source_bindings_path,
        SOURCE_BINDING_FIELDS,
        label="composite source bindings",
    )
    s77_bindings = [row for row in bindings if row["source_id"] == "later_s0077"]
    if len(s77_bindings) != 1:
        raise FM0ContractError("composite authority has no unique S77 source")
    s77_binding = s77_bindings[0]
    s77_root = Path(s77_binding["release_root"]).expanduser().resolve(strict=True)
    manifest_path = s77_root / "manifest.csv"
    if (
        manifest_path.is_symlink()
        or not manifest_path.is_file()
        or sha256_file(manifest_path) != s77_binding["manifest_sha256"]
    ):
        raise FM0ContractError("composite S77 manifest binding drifted")
    if require_read_only and manifest_path.stat().st_mode & 0o222:
        raise FM0ContractError("composite S77 manifest must be read-only")
    s77_manifest_rows = _read_exact_csv(
        manifest_path,
        LATER_SIX_VIEW_MANIFEST_FIELDS,
        label="composite S77 manifest",
    )
    s77_by_observation = {
        row["observation_key"]: row for row in s77_manifest_rows
    }
    composite_holdout = {
        item.observation_key: item for item in composite_loaded.holdout
    }

    resolved: dict[str, SourceBinding] = {}
    for identity in identity_rows:
        observation = identity["observation_key"]
        authority = identity["source_authority"]
        if authority == "temporal_panel":
            try:
                source = panel_by_observation[observation]
            except KeyError as exc:
                raise FM0ContractError(
                    "identity-plan observation is absent from temporal panel"
                ) from exc
            if (
                source["leakage_component_id"]
                != identity["leakage_component_id"]
                or source["sector"] != identity["sector"]
                or source["source_partition"] != identity["source_partition"]
                or identity["source_id"] != f"later_s{int(source['sector']):04d}"
            ):
                raise FM0ContractError("temporal-panel observation identity drifted")
            binding = _source_binding_from_row(
                source,
                release_root=source["sector_release_root"],
                require_read_only=require_read_only,
            )
        elif authority == "composite_release":
            try:
                source = s77_by_observation[observation]
                holdout = composite_holdout[observation]
            except KeyError as exc:
                raise FM0ContractError(
                    "identity-plan observation is absent from composite S77 holdout"
                ) from exc
            if (
                holdout.component != identity["leakage_component_id"]
                or source["leakage_component_id"]
                != identity["leakage_component_id"]
                or source["sector"] != identity["sector"]
                or identity["source_id"] != "later_s0077"
                or source["source_partition"] != identity["source_partition"]
            ):
                raise FM0ContractError("composite S77 observation identity drifted")
            binding = _source_binding_from_row(
                source, release_root=s77_root, require_read_only=require_read_only
            )
            if (
                binding.relative_path != holdout.relative_path
                or binding.shard_sha256 != holdout.shard_sha256
                or binding.n_cadences != holdout.n_cadences
            ):
                raise FM0ContractError("composite S77 shard metadata drifted")
        else:
            raise FM0ContractError("identity plan names an unknown source authority")
        if observation in resolved:
            raise FM0ContractError("identity plan repeats one source observation")
        resolved[observation] = binding
    if len(resolved) != EXPECTED_ROWS:
        raise FM0ContractError("payload-plan source count is not exactly 1,440")
    return resolved


def _joint_valid(release: ObservationRelease) -> np.ndarray:
    return (
        np.asarray(release.time_valid, dtype=bool)
        & np.asarray(release.flux_valid[:, REQUIRED_VIEW_INDICES[0]], dtype=bool)
        & np.asarray(release.flux_valid[:, REQUIRED_VIEW_INDICES[1]], dtype=bool)
    )


def eligible_crop_starts(release: ObservationRelease) -> tuple[int, ...]:
    """Return every exact 128-cadence start that passes the frozen screen."""

    joint = _joint_valid(release)
    cumulative = np.concatenate(
        (np.zeros(1, dtype=np.int64), np.cumsum(joint, dtype=np.int64))
    )
    starts: list[int] = []
    for segment in np.unique(release.segment_id):
        indices = np.flatnonzero(release.segment_id == segment)
        if indices.size < CONTEXT_CADENCES:
            continue
        segment_start = int(indices[0])
        segment_stop = int(indices[-1]) + 1
        if not np.array_equal(indices, np.arange(segment_start, segment_stop)):
            raise FM0ContractError("source segment is not contiguous")
        candidates = np.arange(
            segment_start,
            segment_stop - CONTEXT_CADENCES + 1,
            dtype=np.int64,
        )
        valid_counts = (
            cumulative[candidates + CONTEXT_CADENCES] - cumulative[candidates]
        )
        support_counts = (
            cumulative[candidates + CENTER_SUPPORT_STOP_EXCLUSIVE]
            - cumulative[candidates + CENTER_SUPPORT_START]
        )
        eligible = candidates[
            (valid_counts >= MINIMUM_VALID_CADENCES)
            & (
                support_counts
                == CENTER_SUPPORT_STOP_EXCLUSIVE - CENTER_SUPPORT_START
            )
        ]
        starts.extend(int(value) for value in eligible)
    return tuple(starts)


def choose_crop(
    release: ObservationRelease,
    *,
    identity_plan_receipt_sha256: str,
    identity_schedule_sha256: str,
    event_pair_id: str,
    observation_key: str,
) -> CropChoice:
    """Choose one deterministic eligible crop without replacing the identity."""

    _digest(
        identity_plan_receipt_sha256,
        label="identity-plan receipt hash",
    )
    schedule_hash = _digest(
        identity_schedule_sha256,
        label="identity-plan schedule hash",
    )
    starts = eligible_crop_starts(release)
    if not starts:
        raise FM0ContractError(
            f"preselected observation {observation_key} has no eligible 128-cadence crop"
        )
    rank = _stable_index(
        len(starts),
        schedule_hash,
        event_pair_id,
        observation_key,
    )
    start = starts[rank]
    stop = start + CONTEXT_CADENCES
    segment_id = int(release.segment_id[start])
    segment_indices = np.flatnonzero(release.segment_id == segment_id)
    joint = _joint_valid(release)[start:stop]
    return CropChoice(
        segment_id=segment_id,
        segment_start=int(segment_indices[0]),
        segment_stop=int(segment_indices[-1]) + 1,
        crop_start=start,
        crop_stop=stop,
        crop_start_offset=start - int(segment_indices[0]),
        crop_stop_offset=stop - int(segment_indices[0]),
        time_valid=int(np.count_nonzero(release.time_valid[start:stop])),
        adp_1x1_valid=int(
            np.count_nonzero(release.flux_valid[start:stop, REQUIRED_VIEW_INDICES[0]])
        ),
        adp_3x3_valid=int(
            np.count_nonzero(release.flux_valid[start:stop, REQUIRED_VIEW_INDICES[1]])
        ),
        joint_valid=int(np.count_nonzero(joint)),
        eligible_count=len(starts),
        selected_rank=rank,
    )


def crop_arrays(
    release: ObservationRelease, *, start: int
) -> dict[str, np.ndarray]:
    """Materialize the exact source tensors represented by one schedule crop."""

    stop = int(start) + CONTEXT_CADENCES
    if start < 0 or stop > release.n_cadences:
        raise FM0ContractError("crop bounds escaped the source shard")
    arrays: dict[str, np.ndarray] = {}
    for name in sorted(SHARD_ARRAY_KEYS):
        value = np.asarray(getattr(release, name))
        arrays[name] = value.copy() if name == "view_present" else value[start:stop].copy()
    return arrays


def crop_payload_sha256(arrays: Mapping[str, np.ndarray]) -> str:
    """Hash exact cropped source arrays with dtype and shape framing."""

    if frozenset(arrays) != SHARD_ARRAY_KEYS:
        raise FM0ContractError("crop payload differs from the shard-array allowlist")
    digest = hashlib.sha256(CROP_PAYLOAD_HASH_SCHEMA.encode("ascii") + b"\0")
    for name in sorted(arrays):
        array = np.ascontiguousarray(np.asarray(arrays[name]))
        digest.update(name.encode("utf-8") + b"\0")
        digest.update(array.dtype.str.encode("ascii") + b"\0")
        digest.update(
            json.dumps(array.shape, separators=(",", ":")).encode("ascii") + b"\0"
        )
        digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def inject_centered_event(
    clean_crop: Mapping[str, np.ndarray],
    *,
    duration_cadences: int,
    fractional_depth: float,
) -> dict[str, np.ndarray]:
    """Apply the frozen event formula to the two ADP flux views only.

    This helper records the later evaluator's exact mechanic.  The payload
    freezer itself never calls it.
    """

    injected = {name: np.asarray(value).copy() for name, value in clean_crop.items()}
    if frozenset(injected) != SHARD_ARRAY_KEYS:
        raise FM0ContractError("event crop differs from the shard-array allowlist")
    flux = injected["flux"]
    if flux.shape != (CONTEXT_CADENCES, 6):
        raise FM0ContractError("event crop flux has the wrong shape")
    duration = int(duration_cadences)
    if duration not in EVENT_DURATIONS:
        raise FM0ContractError("event duration differs from the frozen grid")
    depth_text = f"{float(fractional_depth):g}"
    if depth_text not in EVENT_DEPTHS_TEXT:
        raise FM0ContractError("event depth differs from the frozen grid")
    profile = centered_trapezoid(duration, float(fractional_depth))
    start = EVENT_CENTER_INDEX_ZERO_BASED - duration // 2
    stop = start + duration
    support_valid = (
        np.asarray(injected["time_valid"], dtype=bool)[start:stop]
        & np.asarray(injected["flux_valid"], dtype=bool)[start:stop, 2]
        & np.asarray(injected["flux_valid"], dtype=bool)[start:stop, 3]
    )
    if not np.all(support_valid):
        raise FM0ContractError("event support is not valid in time and both ADP views")
    for view in REQUIRED_VIEW_INDICES:
        original = flux[start:stop, view].copy()
        flux[start:stop, view] = (1.0 + original) * (1.0 - profile) - 1.0
    for name in SHARD_ARRAY_KEYS - {"flux"}:
        if not np.array_equal(injected[name], np.asarray(clean_crop[name])):
            raise FM0ContractError("event injection changed a non-flux tensor")
    unchanged_views = sorted(set(range(6)) - set(REQUIRED_VIEW_INDICES))
    if not np.array_equal(
        injected["flux"][:, unchanged_views],
        np.asarray(clean_crop["flux"])[:, unchanged_views],
    ):
        raise FM0ContractError("event injection changed a non-ADP flux view")
    return injected


def _screen_source(
    identity: Mapping[str, str],
    binding: SourceBinding,
    *,
    identity_plan_receipt_sha256: str,
    identity_schedule_sha256: str,
    require_read_only: bool,
) -> dict[str, Any]:
    relative = PurePosixPath(binding.relative_path)
    raw_path = binding.release_root.joinpath(*relative.parts)
    if raw_path.is_symlink():
        raise FM0ContractError("source shard must be materialized")
    path = raw_path.resolve(strict=True)
    try:
        path.relative_to(binding.release_root)
    except ValueError as exc:
        raise FM0ContractError("source shard escaped its release root") from exc
    if not path.is_file() or path.stat().st_size <= 0:
        raise FM0ContractError("source shard is missing or empty")
    if require_read_only and path.stat().st_mode & 0o222:
        raise FM0ContractError("source shard must be read-only")
    payload = path.read_bytes()
    if hashlib.sha256(payload).hexdigest() != binding.shard_sha256:
        raise FM0ContractError("source shard SHA-256 drifted")
    release = load_input_release_bytes(payload)
    if (
        release.n_cadences != binding.n_cadences
        or release.n_segments != binding.n_segments
        or tuple(bool(value) for value in release.view_present)
        != binding.view_present
    ):
        raise FM0ContractError("source shard dimensions differ from its manifest")
    choice = choose_crop(
        release,
        identity_plan_receipt_sha256=identity_plan_receipt_sha256,
        identity_schedule_sha256=identity_schedule_sha256,
        event_pair_id=identity["event_pair_id"],
        observation_key=identity["observation_key"],
    )
    crop_hash = crop_payload_sha256(crop_arrays(release, start=choice.crop_start))
    identity_row = {field: identity[field] for field in SCHEDULE_FIELDS}
    row: dict[str, Any] = {
        "payload_schema_version": PAYLOAD_PLAN_SCHEMA_VERSION,
        **identity_row,
    }
    row.update(
        {
            "identity_row_sha256": _canonical_sha256(identity_row),
            "source_release_root": str(binding.release_root),
            "source_relative_path": binding.relative_path,
            "source_shard_sha256": binding.shard_sha256,
            "source_n_cadences": str(binding.n_cadences),
            "source_n_segments": str(binding.n_segments),
            "crop_segment_id": str(choice.segment_id),
            "segment_start_index_zero_based": str(choice.segment_start),
            "segment_stop_index_exclusive": str(choice.segment_stop),
            "crop_start_index_zero_based": str(choice.crop_start),
            "crop_stop_index_exclusive": str(choice.crop_stop),
            "crop_start_offset_within_segment": str(choice.crop_start_offset),
            "crop_stop_offset_within_segment_exclusive": str(
                choice.crop_stop_offset
            ),
            "n_time_valid": str(choice.time_valid),
            "n_adp_1x1_valid": str(choice.adp_1x1_valid),
            "n_adp_3x3_valid": str(choice.adp_3x3_valid),
            "n_joint_valid": str(choice.joint_valid),
            "joint_valid_fraction": f"{choice.joint_valid / CONTEXT_CADENCES:.8f}",
            "center_support_start_index_zero_based": str(CENTER_SUPPORT_START),
            "center_support_stop_index_exclusive": str(
                CENTER_SUPPORT_STOP_EXCLUSIVE
            ),
            "n_center_support_joint_valid": str(
                CENTER_SUPPORT_STOP_EXCLUSIVE - CENTER_SUPPORT_START
            ),
            "n_eligible_crops": str(choice.eligible_count),
            "selected_eligible_crop_rank_zero_based": str(choice.selected_rank),
            "crop_payload_sha256": crop_hash,
        }
    )
    if tuple(row) != PAYLOAD_SCHEDULE_FIELDS:
        raise FM0ContractError("payload schedule row columns drifted")
    return row


def _csv_payload(rows: Sequence[Mapping[str, Any]]) -> bytes:
    stream = StringIO(newline="")
    writer = csv.DictWriter(
        stream, fieldnames=list(PAYLOAD_SCHEDULE_FIELDS), lineterminator="\n"
    )
    writer.writeheader()
    for row in rows:
        if tuple(row) != PAYLOAD_SCHEDULE_FIELDS:
            raise FM0ContractError("payload schedule columns drifted")
        writer.writerow(row)
    return stream.getvalue().encode("utf-8")


def _binding_hash(
    rows: Sequence[Mapping[str, str]], *, crop: bool
) -> str:
    records: list[Mapping[str, str]] = []
    for row in rows:
        if crop:
            records.append(
                {
                    "event_pair_id": row["event_pair_id"],
                    "observation_key": row["observation_key"],
                    "crop_payload_sha256": row["crop_payload_sha256"],
                }
            )
        else:
            records.append(
                {
                    "observation_key": row["observation_key"],
                    "source_release_root": row["source_release_root"],
                    "source_relative_path": row["source_relative_path"],
                    "source_shard_sha256": row["source_shard_sha256"],
                }
            )
    return _canonical_sha256(records)


def science_mechanics() -> dict[str, Any]:
    """Return the complete fixed evaluator contract bound by the receipt."""

    return {
        "scientific_question": (
            "does either short-context cadence-preserving backbone encode a "
            "single centered transit-like event for reusable downstream "
            "classification and triage"
        ),
        "bls_or_periodic_search": False,
        "matched_architectures": list(ARCHITECTURES),
        "cadence": {
            "context_cadences": CONTEXT_CADENCES,
            "nominal_seconds": 200,
            "one_encoder_token_per_native_cadence": True,
            "patch_stride": 1,
            "averaging_merging_downsampling_or_pooling": False,
        },
        "event": {
            "profile": "symmetric_trapezoid_sampled_at_cadence_centers",
            "ingress_egress_fraction_of_duration": 1.0 / 3.0,
            "center_index_zero_based": EVENT_CENTER_INDEX_ZERO_BASED,
            "duration_cadences": list(EVENT_DURATIONS),
            "fractional_depths": [float(value) for value in EVENT_DEPTHS_TEXT],
            "support_indices_for_longest_event_inclusive": [60, 68],
            "apply_to_flux_views_only": ["adp_1x1", "adp_3x3"],
            "formula": "injected=(1+flux)*(1-profile)-1",
            "all_other_tensors_and_flux_views_bitwise_identical": True,
            "one_adjacent_clean_injected_pair_per_schedule_row": True,
            "period_defined": False,
        },
        "inference": {
            "temporal_mask_all_128_cadences": "zero_false",
            "reconstruction_mask_all_two_adp_views_by_128_cadences": "zero_false",
            "mask_token_inserted": False,
            "readout": "linear_logit_from_h_cadence_token_64_only",
            "window_or_token_pooling_before_readout": False,
        },
        "features": {
            "raw_adp_validity_error_exact_center": list(RAW_FEATURE_CHANNELS),
            "quality_only_exact_center": list(QUALITY_FEATURE_CHANNELS),
            "model_arms": [
                "TWIRL-FM0.3.1_step0_h_cadence_token64",
                "TWIRL-FM0.3.1_step2000_h_cadence_token64",
                "TWIRL-FM0.3.2_step0_h_cadence_token64",
                "TWIRL-FM0.3.2_step2000_h_cadence_token64",
            ],
            "event_duration_depth_support_period_bls_or_search_score": "forbidden",
        },
        "probe": {
            "family": "linear_logit_binary_classifier",
            "one_separate_fitted_weight_vector_per_feature_arm": True,
            "identical_fixed_algorithm_for_every_arm": True,
            "input": "exact_center_token_or_exact_center_raw_control_only",
            "standardization": "dimensionwise_mean_and_std_fit_on_probe_train_only",
            "fit_split": "probe_train",
            "validation_role": "diagnostic_only_no_tuning_or_selection",
            "fresh_s77_test_used_for_fit_tuning_threshold_or_selection": False,
            "optimizer": "adam",
            "epochs": PROBE_EPOCHS,
            "learning_rate": PROBE_LEARNING_RATE,
            "l2_weight": PROBE_L2,
            "initialization_seed": PROBE_SEED,
            "hyperparameter_tuning": False,
        },
        "metrics": {
            "primary": "sample_roc_auc",
            "confidence_interval_level": 0.95,
            "bootstrap": {
                "method": "component_pair_clustered_bootstrap",
                "resampling_unit": (
                    "whole_leakage_component_with_its_clean_injected_pair"
                ),
                "replicates": BOOTSTRAP_REPLICATES,
                "seed": BOOTSTRAP_SEED,
                "same_resamples_across_arms_and_paired_deltas": True,
            },
        },
    }


def screen_contract() -> dict[str, Any]:
    """Return the exact direct-128 payload-eligibility contract."""

    return {
        "context_cadences": CONTEXT_CADENCES,
        "unpadded": True,
        "within_exactly_one_segment": True,
        "minimum_joint_time_and_both_adp_valid_fraction": MINIMUM_VALID_FRACTION,
        "minimum_joint_valid_cadences": MINIMUM_VALID_CADENCES,
        "complete_joint_valid_support_indices_inclusive": [60, 68],
        "deterministic_selection_salt": PAYLOAD_PLAN_SALT,
        "ineligible_preselected_identity_policy": "fail_closed_no_replacement",
    }


def _make_tree_read_only(root: Path) -> None:
    for path in root.rglob("*"):
        path.chmod(0o555 if path.is_dir() else 0o444)
    root.chmod(0o555)


def _validate_row(row: Mapping[str, str]) -> None:
    identity = {field: row[field] for field in SCHEDULE_FIELDS}
    integers = (
        "source_n_cadences",
        "source_n_segments",
        "crop_segment_id",
        "segment_start_index_zero_based",
        "segment_stop_index_exclusive",
        "crop_start_index_zero_based",
        "crop_stop_index_exclusive",
        "crop_start_offset_within_segment",
        "crop_stop_offset_within_segment_exclusive",
        "n_time_valid",
        "n_adp_1x1_valid",
        "n_adp_3x3_valid",
        "n_joint_valid",
        "center_support_start_index_zero_based",
        "center_support_stop_index_exclusive",
        "n_center_support_joint_valid",
        "n_eligible_crops",
        "selected_eligible_crop_rank_zero_based",
    )
    try:
        values = {field: int(row[field]) for field in integers}
    except ValueError as exc:
        raise FM0ContractError("payload schedule integer is invalid") from exc
    start = values["crop_start_index_zero_based"]
    stop = values["crop_stop_index_exclusive"]
    segment_start = values["segment_start_index_zero_based"]
    segment_stop = values["segment_stop_index_exclusive"]
    offset_start = values["crop_start_offset_within_segment"]
    offset_stop = values["crop_stop_offset_within_segment_exclusive"]
    if (
        row["payload_schema_version"] != PAYLOAD_PLAN_SCHEMA_VERSION
        or row["schema_version"] != MATCHED_CANARY_PLAN_SCHEMA_VERSION
        or row["identity_row_sha256"] != _canonical_sha256(identity)
        or _SHA256.fullmatch(row["source_shard_sha256"]) is None
        or _SHA256.fullmatch(row["crop_payload_sha256"]) is None
        or stop - start != CONTEXT_CADENCES
        or not (segment_start <= start < stop <= segment_stop)
        or offset_start != start - segment_start
        or offset_stop != stop - segment_start
        or offset_stop - offset_start != CONTEXT_CADENCES
        or values["n_joint_valid"] < MINIMUM_VALID_CADENCES
        or values["n_joint_valid"] > CONTEXT_CADENCES
        or any(
            not 0 <= values[field] <= CONTEXT_CADENCES
            for field in (
                "n_time_valid",
                "n_adp_1x1_valid",
                "n_adp_3x3_valid",
            )
        )
        or row["joint_valid_fraction"]
        != f"{values['n_joint_valid'] / CONTEXT_CADENCES:.8f}"
        or values["center_support_start_index_zero_based"]
        != CENTER_SUPPORT_START
        or values["center_support_stop_index_exclusive"]
        != CENTER_SUPPORT_STOP_EXCLUSIVE
        or values["n_center_support_joint_valid"]
        != CENTER_SUPPORT_STOP_EXCLUSIVE - CENTER_SUPPORT_START
        or values["n_eligible_crops"] <= 0
        or not 0
        <= values["selected_eligible_crop_rank_zero_based"]
        < values["n_eligible_crops"]
        or values["source_n_cadences"] < stop
        or values["crop_segment_id"] < 0
        or values["source_n_segments"] <= values["crop_segment_id"]
        or _safe_relative_shard(
            row["source_relative_path"], observation_key=row["observation_key"]
        )
        != row["source_relative_path"]
    ):
        raise FM0ContractError("payload schedule row contract drifted")


def validate_matched_canary_payload_plan(
    root: str | Path,
    *,
    expected_receipt_sha256: str,
    require_read_only: bool = True,
) -> MatchedCanaryPayloadPlanResult:
    """Validate the frozen schedule without reopening any source payload."""

    output_raw = Path(root).expanduser()
    if output_raw.is_symlink():
        raise FM0ContractError("payload-plan root must be materialized")
    output = output_raw.resolve(strict=True)
    if (
        not output.is_dir()
        or {path.name for path in output.iterdir()}
        != {"schedule.csv", "receipt.json", "READY"}
    ):
        raise FM0ContractError("payload-plan root closure drifted")
    receipt_path = output / "receipt.json"
    schedule_path = output / "schedule.csv"
    ready_path = output / "READY"
    for path in (receipt_path, schedule_path, ready_path):
        if path.is_symlink() or not path.is_file():
            raise FM0ContractError("payload-plan artifact must be materialized")
        if require_read_only and path.stat().st_mode & 0o222:
            raise FM0ContractError("payload-plan artifact must be read-only")
    if require_read_only and output.stat().st_mode & 0o222:
        raise FM0ContractError("payload-plan root must be read-only")
    receipt_hash = sha256_file(receipt_path)
    if receipt_hash != _digest(
        expected_receipt_sha256, label="payload-plan receipt hash"
    ):
        raise FM0ContractError("payload-plan receipt SHA-256 drifted")
    if ready_path.read_text(encoding="utf-8").strip() != receipt_hash:
        raise FM0ContractError("payload-plan READY/receipt binding drifted")
    try:
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FM0ContractError("payload-plan receipt is invalid") from exc
    if (
        not isinstance(receipt, Mapping)
        or receipt.get("schema_version") != PAYLOAD_PLAN_RECEIPT_SCHEMA_VERSION
        or receipt.get("ready_state") != PAYLOAD_PLAN_READY_STATE
        or receipt.get("passed") is not True
        or _GIT_SHA.fullmatch(str(receipt.get("producer_git_sha", ""))) is None
        or receipt.get("screen") != screen_contract()
        or receipt.get("science_mechanics") != science_mechanics()
    ):
        raise FM0ContractError("payload-plan receipt contract drifted")
    schedule_hash = sha256_file(schedule_path)
    schedule = receipt.get("schedule")
    if (
        not isinstance(schedule, Mapping)
        or schedule.get("path") != "schedule.csv"
        or schedule.get("sha256") != schedule_hash
        or schedule.get("n_rows") != EXPECTED_ROWS
    ):
        raise FM0ContractError("payload-plan schedule binding drifted")
    rows = _read_exact_csv(
        schedule_path, PAYLOAD_SCHEDULE_FIELDS, label="payload schedule"
    )
    if len(rows) != EXPECTED_ROWS:
        raise FM0ContractError("payload-plan schedule is not exactly 1,440 rows")
    counts: dict[str, dict[str, int]] = {
        split: defaultdict(int) for split in SPLIT_ORDER
    }
    observations: set[str] = set()
    pairs: set[str] = set()
    components: dict[str, set[str]] = defaultdict(set)
    sources: set[tuple[str, str]] = set()
    for row in rows:
        _validate_row(row)
        if row["split"] not in SPLIT_ORDER or row["cohort"] not in TEMPORAL_COHORTS:
            raise FM0ContractError("payload schedule split/cohort drifted")
        if (
            row["observation_key"] in observations
            or row["event_pair_id"] in pairs
            or (row["source_release_root"], row["source_relative_path"]) in sources
        ):
            raise FM0ContractError("payload schedule identities or sources duplicate")
        observations.add(row["observation_key"])
        pairs.add(row["event_pair_id"])
        sources.add((row["source_release_root"], row["source_relative_path"]))
        components[row["split"]].add(row["leakage_component_id"])
        counts[row["split"]][row["cohort"]] += 1
    if any(
        counts[split][cohort] != TARGET_COMPONENTS_PER_COHORT
        for split in SPLIT_ORDER
        for cohort in TEMPORAL_COHORTS
    ):
        raise FM0ContractError("payload schedule split/cohort counts drifted")
    for left_index, left in enumerate(SPLIT_ORDER):
        for right in SPLIT_ORDER[left_index + 1 :]:
            if components[left] & components[right]:
                raise FM0ContractError("payload schedule leaks components across splits")
    bindings = receipt.get("payload_bindings")
    if (
        not isinstance(bindings, Mapping)
        or set(bindings)
        != {
            "n_source_shards_opened",
            "n_source_sha256_verified",
            "n_crops_frozen",
            "source_shard_bindings_sha256",
            "crop_payload_hash_schema",
            "crop_payload_bindings_sha256",
        }
        or bindings.get("source_shard_bindings_sha256")
        != _binding_hash(rows, crop=False)
        or bindings.get("crop_payload_bindings_sha256")
        != _binding_hash(rows, crop=True)
        or bindings.get("crop_payload_hash_schema") != CROP_PAYLOAD_HASH_SCHEMA
        or bindings.get("n_source_shards_opened") != EXPECTED_ROWS
        or bindings.get("n_source_sha256_verified") != EXPECTED_ROWS
        or bindings.get("n_crops_frozen") != EXPECTED_ROWS
    ):
        raise FM0ContractError("payload source/crop aggregate binding drifted")
    limits = receipt.get("limits")
    if limits != {
        "identity_only": False,
        "identity_only_predecessor_validated": True,
        "payload_screening_executed": True,
        "bound_source_payloads_opened": True,
        "events_injected": False,
        "checkpoints_loaded": False,
        "probe_fitted": False,
        "metrics_computed": False,
        "evaluation_executed": False,
        "model_training_executed": False,
        "sealed_test_opened": False,
        "bls_or_search_executed": False,
    }:
        raise FM0ContractError("payload-plan execution boundary drifted")
    identity = receipt.get("identity_plan")
    if (
        not isinstance(identity, Mapping)
        or set(identity)
        != {
            "root",
            "receipt_sha256",
            "schedule_sha256",
            "producer_git_sha",
            "input_authorities",
        }
        or _SHA256.fullmatch(str(identity.get("receipt_sha256", ""))) is None
        or _SHA256.fullmatch(str(identity.get("schedule_sha256", ""))) is None
    ):
        raise FM0ContractError("payload-plan predecessor binding drifted")
    try:
        predecessor = validate_matched_canary_evaluation_plan(
            identity["root"],
            expected_receipt_sha256=identity["receipt_sha256"],
            require_read_only=require_read_only,
        )
    except (OSError, TypeError, ValueError) as exc:
        raise FM0ContractError(
            "payload-plan predecessor failed transitive validation"
        ) from exc
    if (
        predecessor.schedule_sha256 != identity["schedule_sha256"]
        or predecessor.receipt.get("input_authorities")
        != identity["input_authorities"]
        or predecessor.receipt.get("producer_git_sha")
        != identity["producer_git_sha"]
        or identity["producer_git_sha"] != receipt["producer_git_sha"]
    ):
        raise FM0ContractError("payload-plan predecessor authority drifted")
    predecessor_rows = _identity_rows(predecessor.schedule_path)
    if len(predecessor_rows) != len(rows):
        raise FM0ContractError("payload-plan predecessor row count drifted")
    for payload_row, identity_row in zip(rows, predecessor_rows, strict=True):
        if any(payload_row[field] != identity_row[field] for field in SCHEDULE_FIELDS):
            raise FM0ContractError(
                "payload schedule differs from its ordered identity predecessor"
            )
    return MatchedCanaryPayloadPlanResult(
        root=output,
        schedule_path=schedule_path,
        receipt_path=receipt_path,
        ready_path=ready_path,
        schedule_sha256=schedule_hash,
        receipt_sha256=receipt_hash,
        receipt=receipt,
    )


def freeze_matched_canary_payload_plan(
    output_dir: str | Path,
    *,
    identity_plan_root: str | Path,
    identity_plan_receipt_sha256: str,
    producer_git_sha: str,
    require_read_only: bool = True,
) -> MatchedCanaryPayloadPlanResult:
    """Open and screen the exact preselected shards, then freeze their crops."""

    identity_receipt_hash = _digest(
        identity_plan_receipt_sha256, label="identity-plan receipt hash"
    )
    identity_plan = validate_matched_canary_evaluation_plan(
        identity_plan_root,
        expected_receipt_sha256=identity_receipt_hash,
        require_read_only=require_read_only,
    )
    producer = _git_sha(producer_git_sha)
    if identity_plan.receipt.get("producer_git_sha") != producer:
        raise FM0ContractError(
            "identity plan was not frozen by the requested exact Git revision"
        )
    output = Path(output_dir).expanduser().absolute()
    if output.exists():
        ready = output / "READY"
        if not ready.is_file():
            raise FM0ContractError("existing payload plan has no READY binding")
        result = validate_matched_canary_payload_plan(
            output,
            expected_receipt_sha256=ready.read_text(encoding="utf-8").strip(),
            require_read_only=require_read_only,
        )
        if (
            result.receipt["identity_plan"]["receipt_sha256"]
            != identity_plan.receipt_sha256
            or result.receipt["identity_plan"]["schedule_sha256"]
            != identity_plan.schedule_sha256
            or result.receipt["producer_git_sha"] != producer
        ):
            raise FM0ContractError("existing payload plan binds another identity plan")
        return result
    if output.is_symlink() or output.name in {"", ".", ".."}:
        raise FM0ContractError("payload-plan output must be a named directory")

    identity_rows = _identity_rows(identity_plan.schedule_path)
    if len(identity_rows) != EXPECTED_ROWS:
        raise FM0ContractError("identity plan is not exactly 1,440 rows")
    sources = _resolve_sources(
        identity_rows,
        identity_receipt=identity_plan.receipt,
        require_read_only=require_read_only,
    )
    screened: list[dict[str, Any]] = []
    for index, identity in enumerate(identity_rows, start=1):
        screened.append(
            _screen_source(
                identity,
                sources[identity["observation_key"]],
                identity_plan_receipt_sha256=identity_plan.receipt_sha256,
                identity_schedule_sha256=identity_plan.schedule_sha256,
                require_read_only=require_read_only,
            )
        )
        if index == 1 or index % 25 == 0 or index == EXPECTED_ROWS:
            print(
                f"FM0_3_PAYLOAD_PLAN_PROGRESS screened={index}/{EXPECTED_ROWS}",
                flush=True,
            )
    schedule_payload = _csv_payload(screened)
    schedule_hash = hashlib.sha256(schedule_payload).hexdigest()
    receipt = {
        "schema_version": PAYLOAD_PLAN_RECEIPT_SCHEMA_VERSION,
        "ready_state": PAYLOAD_PLAN_READY_STATE,
        "passed": True,
        "producer_git_sha": producer,
        "identity_plan": {
            "root": str(identity_plan.root),
            "receipt_sha256": identity_plan.receipt_sha256,
            "schedule_sha256": identity_plan.schedule_sha256,
            "producer_git_sha": identity_plan.receipt["producer_git_sha"],
            "input_authorities": identity_plan.receipt["input_authorities"],
        },
        "screen": screen_contract(),
        "payload_bindings": {
            "n_source_shards_opened": EXPECTED_ROWS,
            "n_source_sha256_verified": EXPECTED_ROWS,
            "n_crops_frozen": EXPECTED_ROWS,
            "source_shard_bindings_sha256": _binding_hash(screened, crop=False),
            "crop_payload_hash_schema": CROP_PAYLOAD_HASH_SCHEMA,
            "crop_payload_bindings_sha256": _binding_hash(screened, crop=True),
        },
        "science_mechanics": science_mechanics(),
        "schedule": {
            "path": "schedule.csv",
            "sha256": schedule_hash,
            "n_rows": EXPECTED_ROWS,
        },
        "limits": {
            "identity_only": False,
            "identity_only_predecessor_validated": True,
            "payload_screening_executed": True,
            "bound_source_payloads_opened": True,
            "events_injected": False,
            "checkpoints_loaded": False,
            "probe_fitted": False,
            "metrics_computed": False,
            "evaluation_executed": False,
            "model_training_executed": False,
            "sealed_test_opened": False,
            "bls_or_search_executed": False,
        },
        "claim_limit": (
            "Payload-screened pre-checkpoint crop schedule only; source payloads "
            "were opened and hashed, but no event injection, checkpoint inference, "
            "probe fitting, metric evaluation, model training, BLS/search, sealed "
            "access, architecture result, or promotion occurred."
        ),
    }
    receipt_payload = (
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    receipt_hash = hashlib.sha256(receipt_payload).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = Path(
        tempfile.mkdtemp(prefix=f".{output.name}.partial.", dir=output.parent)
    )
    try:
        (partial / "schedule.csv").write_bytes(schedule_payload)
        (partial / "receipt.json").write_bytes(receipt_payload)
        (partial / "READY").write_text(receipt_hash + "\n", encoding="utf-8")
        _make_tree_read_only(partial)
        os.replace(partial, output)
    except BaseException:
        if partial.exists():
            partial.chmod(0o755)
            for path in partial.rglob("*"):
                path.chmod(0o755 if path.is_dir() else 0o644)
            shutil.rmtree(partial)
        raise
    return validate_matched_canary_payload_plan(
        output,
        expected_receipt_sha256=receipt_hash,
        require_read_only=require_read_only,
    )
