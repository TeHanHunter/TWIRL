"""Mission-provider-neutral cadence quality for later TWIRL-FM inputs.

The frozen A2v1 v1 quality reference names its mission column
``spoc_quality``.  That is scientifically correct for the original sectors,
but it cannot describe the QLP provider transition at Sector 67.  This module
defines a separate v2 contract with provider-neutral ``mission_quality`` and
``qlp_quality`` columns.  It does not modify or reinterpret the v1 contract.

Each table row is authorized by a QLP camera quaternion cadence.  Detector QLP
qflags must cover that cadence exactly, while the sector's mission authority
is SPOC before S67 and TICA from S67 onward.  Mission-only rows are retained
outside the table as immutable authority exclusions: they may be recognized
and masked in a light curve, but they are never promoted into the quaternion
cadence authority.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from twirl.models.fm0.mission_quality_admission import (
    MISSION_QUALITY_EXCLUSION_POLICY,
    MISSION_QUALITY_POLICY,
    MISSION_QUALITY_READY_STATE,
    MISSION_QUALITY_RECEIPT_SCHEMA_VERSION,
    _read_flag_file,
    _read_quat,
    mission_quality_type,
    verify_mission_quality_sources,
)
from twirl.models.fm0.mission_quality_transfer import (
    MISSION_QUALITY_TRANSFER_SCHEMA_VERSION,
)
from twirl.models.fm0.registry import FM0ContractError, sha256_file

MISSION_QUALITY_REFERENCE_SCHEMA_VERSION = "twirl_fm0_mission_quality_reference_v2"
MISSION_QUALITY_REFERENCE_BUILDER_VERSION = (
    "twirl_fm0_mission_quality_reference_builder_v2"
)
MISSION_QUALITY_REFERENCE_TABLE = "cadence_quality.csv"
MISSION_QUALITY_REFERENCE_MANIFEST = "manifest.json"
MISSION_QUALITY_REFERENCE_READY = "READY"
MISSION_QUALITY_REFERENCE_COLUMNS = (
    "sector",
    "orbit",
    "camera",
    "ccd",
    "cadence",
    "mission_quality",
    "qlp_quality",
)
MISSION_QUALITY_CADENCE_AUTHORITY = "qlp_camera_quaternion"
QLP_QUALITY_EXTERNAL_BIT = 30
MISSION_QUALITY_COMPOSITION = {
    "external_quality": "mission_quality | (qlp_quality << 30)",
    "qlp_quality_raw_values": [0, 1],
    "qlp_quality_external_bit": QLP_QUALITY_EXTERNAL_BIT,
    "authority_excluded": (
        "mission-only detector cadence absent from the QLP camera quaternion; "
        "mask independently"
    ),
}

_SHA40 = re.compile(r"^[0-9a-f]{40}$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_EXPECTED_DETECTORS = tuple(
    (camera, ccd) for camera in range(1, 5) for ccd in range(1, 5)
)
_REFERENCE_MANIFEST_FIELDS = frozenset(
    {
        "schema_version",
        "builder_version",
        "generated_utc",
        "producer_git_sha",
        "sector",
        "expected_orbits",
        "mission_quality_provider",
        "quality_policy",
        "cadence_authority",
        "quality_composition",
        "table_file",
        "table_sha256",
        "table_columns",
        "n_rows",
        "n_rows_by_cell",
        "n_nonzero_mission_quality",
        "n_nonzero_qlp_quality",
        "detectors",
        "authority_exclusions",
        "authority_exclusions_sha256",
        "quality_transfer_manifest_path",
        "quality_transfer_manifest_sha256",
        "quality_transfer_schema_version",
        "quality_transfer_producer_git_sha",
        "mission_quality_receipt_path",
        "mission_quality_receipt_sha256",
        "mission_quality_receipt_schema_version",
        "source_bindings",
        "source_declaration_sha256",
        "hdf5_quality_join_required",
        "six_view_shards_verified",
        "panel_admission_authorized",
        "claim_limit",
    }
)


def _canonical_json_sha256(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _make_output_tree_read_only(root: Path) -> None:
    """Remove every write bit before an output tree is published."""

    if root.is_symlink() or not root.is_dir():
        raise FM0ContractError(
            f"quality-reference output is not a materialized directory: {root}"
        )
    for path in sorted(root.rglob("*"), key=lambda value: len(value.parts), reverse=True):
        if path.is_symlink():
            raise FM0ContractError(f"quality-reference output contains a symlink: {path}")
        path.chmod(0o550 if path.is_dir() else 0o440)
    root.chmod(0o550)


def _assert_output_tree_read_only(root: Path) -> None:
    if root.is_symlink() or not root.is_dir():
        raise FM0ContractError(
            f"quality-reference output is not a materialized directory: {root}"
        )
    for path in (root, *root.rglob("*")):
        if path.is_symlink() or path.stat().st_mode & 0o222:
            raise FM0ContractError(f"quality-reference output is not read-only: {path}")


def _remove_owned_output_tree(root: Path) -> None:
    """Remove exactly one builder-owned tree without following symlinks."""

    if not root.exists() and not root.is_symlink():
        return
    if root.is_symlink() or not root.is_dir():
        raise FM0ContractError(f"refusing to clean unsafe output path: {root}")
    # Unlinking a read-only file only requires a writable parent. Restore
    # directory modes, not every file, so scheduler cleanup stays bounded.
    for directory, names, _files in os.walk(root, topdown=True, followlinks=False):
        directory_path = Path(directory)
        directory_path.chmod(0o700)
        for name in names:
            child = directory_path / name
            if not child.is_symlink():
                child.chmod(0o700)
    shutil.rmtree(root)


def _output_leaf_path(value: str | Path, *, label: str) -> Path:
    """Resolve the parent while preserving the output leaf for symlink checks."""

    requested = Path(value).expanduser()
    if requested.name in {"", ".", ".."}:
        raise FM0ContractError(f"{label} must name a concrete output directory")
    return requested.parent.resolve() / requested.name


def _load_json_object(path: Path, *, label: str) -> Mapping[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FM0ContractError(f"invalid {label}: {path}") from exc
    if not isinstance(payload, Mapping):
        raise FM0ContractError(f"{label} must be a JSON object: {path}")
    return payload


def _positive_int(value: Any, *, label: str, allow_zero: bool = False) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise FM0ContractError(f"{label} must be an integer")
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise FM0ContractError(f"{label} must be an integer") from exc
    if str(value).strip() != str(result) and not isinstance(value, (int, np.integer)):
        raise FM0ContractError(f"{label} must be an exact integer")
    if result < (0 if allow_zero else 1):
        raise FM0ContractError(f"{label} is outside its allowed range")
    return result


def _digest(value: Any, *, label: str, length: int = 64) -> str:
    result = str(value).strip().lower()
    pattern = _SHA40 if length == 40 else _SHA256
    if pattern.fullmatch(result) is None:
        raise FM0ContractError(
            f"{label} must be a lowercase {length}-character hexadecimal digest"
        )
    return result


def _expected_orbits(
    sector: int, supplied: Sequence[int] | None = None
) -> tuple[int, int]:
    expected = (2 * int(sector) + 7, 2 * int(sector) + 8)
    observed = (
        tuple(
            _positive_int(value, label=f"S{sector} supplied orbit")
            for value in supplied
        )
        if supplied is not None
        else None
    )
    if observed is not None and observed != expected:
        raise FM0ContractError(
            f"S{sector} expected orbits are {expected}, observed {observed}"
        )
    return expected


def _staged_path(root: Path, value: Any, *, label: str) -> Path:
    relative = Path(str(value))
    if relative.is_absolute() or ".." in relative.parts or str(relative) in {"", "."}:
        raise FM0ContractError(f"{label} must be a safe relative path")
    resolved = (root / relative).resolve()
    if not resolved.is_relative_to(root):
        raise FM0ContractError(f"{label} escapes the transfer root")
    return resolved


def _validate_transfer_manifest(
    root: Path,
) -> tuple[Mapping[str, Any], Path, str]:
    manifest_path = root / "manifest.json"
    ready_path = root / "READY"
    manifest = _load_json_object(
        manifest_path, label="mission-quality transfer manifest"
    )
    producer = _digest(
        manifest.get("producer_git_sha"),
        label="quality transfer producer_git_sha",
        length=40,
    )
    if manifest.get("schema_version") != MISSION_QUALITY_TRANSFER_SCHEMA_VERSION:
        raise FM0ContractError("mission-quality transfer schema mismatch")
    if (
        not ready_path.is_file()
        or ready_path.read_text(encoding="utf-8").strip() != producer
    ):
        raise FM0ContractError("mission-quality transfer READY marker mismatch")

    sectors_raw = manifest.get("sectors")
    if not isinstance(sectors_raw, list):
        raise FM0ContractError("mission-quality transfer sectors must be a list")
    sectors = [_positive_int(value, label="transfer sector") for value in sectors_raw]
    if sectors != sorted(set(sectors)) or any(value < 56 for value in sectors):
        raise FM0ContractError("mission-quality transfer sectors are invalid")

    receipts = manifest.get("receipts")
    sources = manifest.get("sources")
    if not isinstance(receipts, list) or not isinstance(sources, list):
        raise FM0ContractError("mission-quality transfer inventories are invalid")
    if (
        _positive_int(manifest.get("n_sectors"), label="transfer n_sectors")
        != len(sectors)
        or _positive_int(manifest.get("n_receipts"), label="transfer n_receipts")
        != len(receipts)
        or _positive_int(
            manifest.get("n_source_files"), label="transfer n_source_files"
        )
        != len(sources)
    ):
        raise FM0ContractError("mission-quality transfer inventory count mismatch")

    observed_paths: set[Path] = set()
    receipt_sectors: list[int] = []
    observed_bytes = 0
    for kind, rows in (("receipt", receipts), ("source", sources)):
        for index, raw_row in enumerate(rows):
            if not isinstance(raw_row, Mapping):
                raise FM0ContractError(
                    f"quality transfer {kind} row {index} is invalid"
                )
            row_sector = _positive_int(
                raw_row.get("sector"), label=f"quality transfer {kind} sector"
            )
            if row_sector not in sectors:
                raise FM0ContractError(
                    f"quality transfer {kind} references undeclared S{row_sector}"
                )
            if kind == "receipt":
                receipt_sectors.append(row_sector)
            else:
                if not str(raw_row.get("role", "")).strip():
                    raise FM0ContractError("quality transfer source role is empty")
                if not Path(str(raw_row.get("original_path", ""))).is_absolute():
                    raise FM0ContractError(
                        "quality transfer source original_path must be absolute"
                    )
            path = _staged_path(
                root,
                raw_row.get("staged_path"),
                label=f"quality transfer {kind} staged_path",
            )
            if path in observed_paths:
                raise FM0ContractError(f"duplicate quality transfer path: {path}")
            observed_paths.add(path)
            expected_hash = _digest(
                raw_row.get("sha256"), label=f"quality transfer {kind} sha256"
            )
            expected_bytes = _positive_int(
                raw_row.get("bytes"),
                label=f"quality transfer {kind} bytes",
            )
            if (
                not path.is_file()
                or path.stat().st_size != expected_bytes
                or sha256_file(path) != expected_hash
            ):
                raise FM0ContractError(f"quality transfer {kind} drifted: {path}")
            observed_bytes += expected_bytes
    if receipt_sectors != sectors:
        raise FM0ContractError(
            "mission-quality transfer must contain one sector-sorted receipt per sector"
        )
    if (
        _positive_int(manifest.get("n_bytes"), label="transfer n_bytes")
        != observed_bytes
    ):
        raise FM0ContractError("mission-quality transfer byte count mismatch")
    return manifest, manifest_path, sha256_file(manifest_path)


def _binding_identity(binding: Mapping[str, Any]) -> tuple[Any, ...]:
    role = str(binding.get("role", ""))
    expected_fields: set[str]
    if role == "qlp_quaternion":
        expected_fields = {"role", "orbit", "camera", "path", "sha256", "n_rows"}
        identity = (
            role,
            _positive_int(binding.get("orbit"), label="quaternion orbit"),
            _positive_int(binding.get("camera"), label="quaternion camera"),
        )
    elif role == "qlp_detector_qflag":
        expected_fields = {
            "role",
            "orbit",
            "camera",
            "ccd",
            "path",
            "sha256",
            "n_rows",
        }
        identity = (
            role,
            _positive_int(binding.get("orbit"), label="qflag orbit"),
            _positive_int(binding.get("camera"), label="qflag camera"),
            _positive_int(binding.get("ccd"), label="qflag ccd"),
        )
    elif role in {"spoc_mission_quality", "tica_mission_quality"}:
        expected_fields = {"role", "camera", "ccd", "path", "sha256", "n_rows"}
        identity = (
            role,
            _positive_int(binding.get("camera"), label="mission camera"),
            _positive_int(binding.get("ccd"), label="mission ccd"),
        )
    elif role in {"tica_materialization_summary", "tica_materialization_ready"}:
        expected_fields = {"role", "path", "sha256"}
        identity = (role,)
    else:
        raise FM0ContractError(f"unsupported mission-quality source role: {role!r}")
    if set(binding) != expected_fields:
        raise FM0ContractError(f"{role} source binding fields drifted")
    path = Path(str(binding.get("path", ""))).expanduser()
    if not path.is_absolute():
        raise FM0ContractError(f"{role} source binding path must be absolute")
    _digest(binding.get("sha256"), label=f"{role} source binding sha256")
    if "n_rows" in binding:
        _positive_int(binding["n_rows"], label=f"{role} source binding n_rows")
    return identity


def _expected_staged_relative(*, sector: int, binding: Mapping[str, Any]) -> Path:
    role = str(binding["role"])
    name = Path(str(binding.get("path", ""))).name
    if not name:
        raise FM0ContractError("mission-quality source binding has no filename")
    if role in {"qlp_quaternion", "qlp_detector_qflag"}:
        orbit = _positive_int(binding.get("orbit"), label="quality source orbit")
        return Path(f"sources/s{sector:04d}/qlp/orbit-{orbit}/ffi/run/{name}")
    return Path(f"sources/s{sector:04d}/mission/{name}")


def _normalized_source_binding(
    *, binding: Mapping[str, Any], staged_path: Path, regenerated: Mapping[str, Any]
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "role": str(binding["role"]),
        "path": str(staged_path),
        "sha256": _digest(binding.get("sha256"), label="source binding sha256"),
    }
    identity = _binding_identity(binding)
    role = identity[0]
    if role == "qlp_quaternion":
        result.update(orbit=identity[1], camera=identity[2])
    elif role == "qlp_detector_qflag":
        result.update(orbit=identity[1], camera=identity[2], ccd=identity[3])
    elif role in {"spoc_mission_quality", "tica_mission_quality"}:
        result.update(camera=identity[1], ccd=identity[2])
    if "n_rows" in regenerated:
        result["n_rows"] = _positive_int(
            regenerated["n_rows"], label="source binding n_rows"
        )
    return result


@dataclass(frozen=True)
class _VerifiedMissionSources:
    sector: int
    provider: str
    orbits: tuple[int, int]
    transfer_root: Path
    transfer_manifest: Mapping[str, Any]
    transfer_manifest_path: Path
    transfer_manifest_sha256: str
    quality_receipt_path: Path
    quality_receipt_sha256: str
    source_bindings: tuple[Mapping[str, Any], ...]
    quaternion: Mapping[tuple[int, int], tuple[int, ...]]
    qflag: Mapping[tuple[int, int, int], Mapping[int, int]]
    mission: Mapping[tuple[int, int], Mapping[int, int]]
    authority_exclusions: Mapping[str, Any]
    authority_exclusions_sha256: str


def _verify_sector_transfer(
    *, quality_transfer_root: str | Path, sector: int
) -> _VerifiedMissionSources:
    sector = int(sector)
    provider = mission_quality_type(sector)
    orbits = _expected_orbits(sector)
    root = Path(quality_transfer_root).expanduser().resolve()
    if not root.is_dir():
        raise FM0ContractError(f"mission-quality transfer root is missing: {root}")
    transfer, transfer_path, transfer_hash = _validate_transfer_manifest(root)
    if sector not in [int(value) for value in transfer["sectors"]]:
        raise FM0ContractError(f"S{sector} is absent from the mission-quality transfer")

    receipt_rows = [
        row
        for row in transfer["receipts"]
        if isinstance(row, Mapping) and int(row.get("sector", -1)) == sector
    ]
    if len(receipt_rows) != 1:
        raise FM0ContractError(f"S{sector} transfer lacks exactly one quality receipt")
    receipt_path = _staged_path(
        root, receipt_rows[0].get("staged_path"), label="quality receipt staged_path"
    )
    receipt_hash = _digest(
        receipt_rows[0].get("sha256"), label="quality receipt sha256"
    )
    receipt = _load_json_object(
        receipt_path, label=f"S{sector} mission-quality receipt"
    )
    if (
        receipt.get("schema_version") != MISSION_QUALITY_RECEIPT_SCHEMA_VERSION
        or receipt.get("quality_state") != MISSION_QUALITY_READY_STATE
        or receipt.get("quality_policy") != MISSION_QUALITY_POLICY
        or receipt.get("mission_quality_type") != provider
        or receipt.get("passed") is not True
        or int(receipt.get("sector", -1)) != sector
        or tuple(int(value) for value in receipt.get("expected_orbits", ())) != orbits
        or receipt.get("hdf5_quality_join_verified") is not False
        or receipt.get("six_view_shards_verified") is not False
        or receipt.get("panel_admission_authorized") is not False
    ):
        raise FM0ContractError(f"S{sector} mission-quality receipt violates v2")

    source_rows = [
        row
        for row in transfer["sources"]
        if isinstance(row, Mapping) and int(row.get("sector", -1)) == sector
    ]
    bindings = receipt.get("source_bindings")
    if not isinstance(bindings, list) or not all(
        isinstance(value, Mapping) for value in bindings
    ):
        raise FM0ContractError(f"S{sector} quality receipt source bindings are invalid")
    source_by_original: dict[str, Mapping[str, Any]] = {}
    for row in source_rows:
        original = str(row.get("original_path", ""))
        if original in source_by_original:
            raise FM0ContractError(f"S{sector} transfer has duplicate original sources")
        source_by_original[original] = row
    if len(bindings) != len(source_rows):
        raise FM0ContractError(f"S{sector} transfer/source binding count mismatch")

    binding_by_identity: dict[tuple[Any, ...], Mapping[str, Any]] = {}
    staged_by_identity: dict[tuple[Any, ...], Path] = {}
    for binding in bindings:
        identity = _binding_identity(binding)
        if identity in binding_by_identity:
            raise FM0ContractError(f"S{sector} source binding identity is duplicated")
        binding_by_identity[identity] = binding
        original = str(Path(str(binding.get("path", ""))).expanduser().resolve())
        row = source_by_original.get(original)
        if row is None:
            raise FM0ContractError(
                f"S{sector} source binding is absent from transfer: {original}"
            )
        if (
            str(row.get("role", "")) != identity[0]
            or _digest(row.get("sha256"), label="transfer source sha256")
            != _digest(binding.get("sha256"), label="receipt source sha256")
            or Path(str(row.get("staged_path", "")))
            != _expected_staged_relative(sector=sector, binding=binding)
        ):
            raise FM0ContractError(f"S{sector} transfer source role/path/hash drifted")
        staged_by_identity[identity] = _staged_path(
            root, row["staged_path"], label="quality source staged_path"
        )

    mission_role = f"{provider}_mission_quality"
    expected_identities: set[tuple[Any, ...]] = {
        ("qlp_quaternion", orbit, camera) for orbit in orbits for camera in range(1, 5)
    }
    expected_identities.update(
        ("qlp_detector_qflag", orbit, camera, ccd)
        for orbit in orbits
        for camera, ccd in _EXPECTED_DETECTORS
    )
    expected_identities.update(
        (mission_role, camera, ccd) for camera, ccd in _EXPECTED_DETECTORS
    )
    if provider == "tica":
        expected_identities.update(
            (("tica_materialization_summary",), ("tica_materialization_ready",))
        )
    if set(binding_by_identity) != expected_identities:
        missing = sorted(expected_identities - set(binding_by_identity))
        extra = sorted(set(binding_by_identity) - expected_identities)
        raise FM0ContractError(
            f"S{sector} source-role inventory mismatch; missing={missing[:5]}, "
            f"extra={extra[:5]}"
        )

    qlp_root = root / f"sources/s{sector:04d}/qlp"
    mission_root = root / f"sources/s{sector:04d}/mission"
    regenerated = verify_mission_quality_sources(
        sector=sector,
        expected_orbits=orbits,
        qlp_root=qlp_root,
        mission_flag_root=mission_root,
    )
    comparison_fields = (
        "quality_policy",
        "mission_quality_type",
        "n_detectors",
        "n_quaternion_files",
        "n_qflag_files",
        "n_mission_quality_files",
        "n_quaternion_rows",
        "n_qflag_rows",
        "n_mission_quality_rows",
        "n_mission_quality_rows_retained",
        "n_mission_quality_rows_excluded_by_quat",
        "n_nonzero_qflag_rows",
        "n_nonzero_mission_quality_rows",
        "n_nonzero_mission_quality_rows_retained",
        "mission_quality_authority_exclusions",
        "mission_quality_authority_exclusions_sha256",
    )
    for field in comparison_fields:
        if receipt.get(field) != regenerated.get(field):
            raise FM0ContractError(
                f"S{sector} relocated mission-quality evidence disagrees on {field}"
            )

    regenerated_by_identity = {
        _binding_identity(value): value for value in regenerated["source_bindings"]
    }
    if set(regenerated_by_identity) != expected_identities:
        raise FM0ContractError(f"S{sector} regenerated source inventory drifted")
    normalized_bindings: list[Mapping[str, Any]] = []
    for identity in sorted(
        expected_identities, key=lambda value: tuple(map(str, value))
    ):
        original = binding_by_identity[identity]
        regenerated_binding = regenerated_by_identity[identity]
        if _digest(original.get("sha256"), label="receipt source sha256") != _digest(
            regenerated_binding.get("sha256"), label="regenerated source sha256"
        ):
            raise FM0ContractError(f"S{sector} regenerated source hash drifted")
        if original.get("n_rows") != regenerated_binding.get("n_rows"):
            raise FM0ContractError(f"S{sector} regenerated source row count drifted")
        normalized_bindings.append(
            _normalized_source_binding(
                binding=original,
                staged_path=staged_by_identity[identity],
                regenerated=regenerated_binding,
            )
        )

    quaternion: dict[tuple[int, int], tuple[int, ...]] = {}
    qflag: dict[tuple[int, int, int], Mapping[int, int]] = {}
    mission: dict[tuple[int, int], Mapping[int, int]] = {}
    for orbit in orbits:
        for camera in range(1, 5):
            run = qlp_root / f"orbit-{orbit}/ffi/run"
            quaternion[(orbit, camera)] = tuple(
                sorted(_read_quat(run / f"cam{camera}_quat.txt"))
            )
            for ccd in range(1, 5):
                qflag[(orbit, camera, ccd)] = _read_flag_file(
                    run / f"cam{camera}ccd{ccd}_qflag.txt",
                    label="transferred QLP qflag authority",
                )
    prefix = "spocffiflag" if provider == "spoc" else "ticaffiflag"
    for camera, ccd in _EXPECTED_DETECTORS:
        mission[(camera, ccd)] = _read_flag_file(
            mission_root / f"{prefix}_s{sector}_cam{camera}_ccd{ccd}.txt",
            label=f"transferred {provider.upper()} mission-quality authority",
        )

    exclusions = receipt.get("mission_quality_authority_exclusions")
    exclusion_hash = _digest(
        receipt.get("mission_quality_authority_exclusions_sha256"),
        label="mission-quality authority exclusions sha256",
    )
    if (
        not isinstance(exclusions, Mapping)
        or exclusions.get("policy") != MISSION_QUALITY_EXCLUSION_POLICY
        or _canonical_json_sha256(exclusions) != exclusion_hash
    ):
        raise FM0ContractError(f"S{sector} mission-quality exclusions are invalid")
    return _VerifiedMissionSources(
        sector=sector,
        provider=provider,
        orbits=orbits,
        transfer_root=root,
        transfer_manifest=transfer,
        transfer_manifest_path=transfer_path,
        transfer_manifest_sha256=transfer_hash,
        quality_receipt_path=receipt_path,
        quality_receipt_sha256=receipt_hash,
        source_bindings=tuple(normalized_bindings),
        quaternion=quaternion,
        qflag=qflag,
        mission=mission,
        authority_exclusions=dict(exclusions),
        authority_exclusions_sha256=exclusion_hash,
    )


def _quality_frame(sources: _VerifiedMissionSources) -> pd.DataFrame:
    rows: list[tuple[int, int, int, int, int, int, int]] = []
    for orbit in sources.orbits:
        for camera in range(1, 5):
            cadences = sources.quaternion[(orbit, camera)]
            for ccd in range(1, 5):
                qflag = sources.qflag[(orbit, camera, ccd)]
                mission = sources.mission[(camera, ccd)]
                for cadence in cadences:
                    rows.append(
                        (
                            sources.sector,
                            orbit,
                            camera,
                            ccd,
                            cadence,
                            int(mission[cadence]),
                            int(qflag[cadence]),
                        )
                    )
    frame = pd.DataFrame(rows, columns=MISSION_QUALITY_REFERENCE_COLUMNS)
    for column in MISSION_QUALITY_REFERENCE_COLUMNS:
        frame[column] = frame[column].astype(np.int64)
    if frame.empty:
        raise FM0ContractError(f"S{sources.sector} normalized quality table is empty")
    if frame.duplicated(["sector", "camera", "ccd", "cadence"]).any():
        raise FM0ContractError(
            f"S{sources.sector} quaternion cadences overlap between orbits"
        )
    if not set(frame["qlp_quality"]).issubset({0, 1}):
        raise FM0ContractError("normalized QLP quality must be binary")
    if (frame["mission_quality"] < 0).any() or (
        frame["mission_quality"] >= (1 << QLP_QUALITY_EXTERNAL_BIT)
    ).any():
        raise FM0ContractError(
            "mission quality collides with the reserved QLP external bit"
        )
    return frame


def _validate_exclusions(
    *, payload: Any, sector: int, frame: pd.DataFrame
) -> dict[tuple[int, int], tuple[np.ndarray, np.ndarray]]:
    if not isinstance(payload, Mapping) or set(payload) != {
        "policy",
        "n_rows",
        "by_detector",
    }:
        raise FM0ContractError("mission-quality authority exclusions are malformed")
    if payload["policy"] != MISSION_QUALITY_EXCLUSION_POLICY:
        raise FM0ContractError("mission-quality authority exclusion policy mismatch")
    by_detector = payload["by_detector"]
    if not isinstance(by_detector, Mapping) or set(by_detector) != {
        f"cam{camera}_ccd{ccd}" for camera, ccd in _EXPECTED_DETECTORS
    }:
        raise FM0ContractError("mission-quality exclusion detector inventory mismatch")
    table_keys = {
        (int(row.camera), int(row.ccd), int(row.cadence))
        for row in frame.itertuples(index=False)
    }
    result: dict[tuple[int, int], tuple[np.ndarray, np.ndarray]] = {}
    total = 0
    for camera, ccd in _EXPECTED_DETECTORS:
        name = f"cam{camera}_ccd{ccd}"
        entry = by_detector[name]
        if not isinstance(entry, Mapping) or set(entry) != {"n_rows", "rows"}:
            raise FM0ContractError(f"S{sector} {name} exclusion entry is malformed")
        rows = entry["rows"]
        if not isinstance(rows, list) or _positive_int(
            entry["n_rows"], label=f"{name} exclusion n_rows", allow_zero=True
        ) != len(rows):
            raise FM0ContractError(f"S{sector} {name} exclusion count mismatch")
        values: list[tuple[int, int]] = []
        for row in rows:
            if not isinstance(row, Mapping) or set(row) != {
                "cadence",
                "mission_quality",
            }:
                raise FM0ContractError(f"S{sector} {name} exclusion row is malformed")
            cadence = _positive_int(row["cadence"], label="excluded cadence")
            quality = _positive_int(
                row["mission_quality"],
                label="excluded mission quality",
                allow_zero=True,
            )
            if quality >= (1 << QLP_QUALITY_EXTERNAL_BIT):
                raise FM0ContractError("excluded mission quality exceeds bit budget")
            if (camera, ccd, cadence) in table_keys:
                raise FM0ContractError(
                    f"S{sector} {name} exclusion is present in quaternion table"
                )
            values.append((cadence, quality))
        ordered = sorted(values)
        if len({cadence for cadence, _ in ordered}) != len(ordered):
            raise FM0ContractError(f"S{sector} {name} exclusions are duplicated")
        cadences = np.asarray([value[0] for value in ordered], dtype=np.int64)
        qualities = np.asarray([value[1] for value in ordered], dtype=np.uint64)
        cadences.setflags(write=False)
        qualities.setflags(write=False)
        result[(camera, ccd)] = (cadences, qualities)
        total += len(ordered)
    if (
        _positive_int(
            payload["n_rows"], label="authority exclusion n_rows", allow_zero=True
        )
        != total
    ):
        raise FM0ContractError("mission-quality authority exclusion total mismatch")
    return result


def build_mission_quality_reference(
    *,
    quality_transfer_root: str | Path,
    expected_quality_transfer_manifest_sha256: str,
    sector: int,
    output_dir: str | Path,
    producer_git_sha: str,
) -> dict[str, Any]:
    """Publish one immutable provider-neutral cadence-quality reference."""

    producer = _digest(producer_git_sha, label="reference producer_git_sha", length=40)
    sector = _positive_int(sector, label="mission-quality reference sector")
    sources = _verify_sector_transfer(
        quality_transfer_root=quality_transfer_root, sector=sector
    )
    expected_transfer_hash = _digest(
        expected_quality_transfer_manifest_sha256,
        label="expected quality-transfer manifest sha256",
    )
    if sources.transfer_manifest_sha256 != expected_transfer_hash:
        raise FM0ContractError("mission-quality transfer manifest hash drifted")
    frame = _quality_frame(sources)
    final = _output_leaf_path(output_dir, label="quality-reference output")
    partial = final.with_name(final.name + ".partial")
    if (
        final.exists()
        or final.is_symlink()
        or partial.exists()
        or partial.is_symlink()
    ):
        raise FM0ContractError(f"refusing to overwrite quality reference: {final}")
    owns_partial = False
    owns_final = False
    try:
        partial.mkdir(parents=True)
        owns_partial = True
        table_path = partial / MISSION_QUALITY_REFERENCE_TABLE
        frame.to_csv(table_path, index=False, lineterminator="\n")
        table_hash = sha256_file(table_path)
        n_rows_by_cell = {
            f"o{int(orbit)}_cam{int(camera)}_ccd{int(ccd)}": len(group)
            for (orbit, camera, ccd), group in frame.groupby(
                ["orbit", "camera", "ccd"], sort=True
            )
        }
        manifest = {
            "schema_version": MISSION_QUALITY_REFERENCE_SCHEMA_VERSION,
            "builder_version": MISSION_QUALITY_REFERENCE_BUILDER_VERSION,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "producer_git_sha": producer,
            "sector": sources.sector,
            "expected_orbits": list(sources.orbits),
            "mission_quality_provider": sources.provider,
            "quality_policy": MISSION_QUALITY_POLICY,
            "cadence_authority": MISSION_QUALITY_CADENCE_AUTHORITY,
            "quality_composition": MISSION_QUALITY_COMPOSITION,
            "table_file": MISSION_QUALITY_REFERENCE_TABLE,
            "table_sha256": table_hash,
            "table_columns": list(MISSION_QUALITY_REFERENCE_COLUMNS),
            "n_rows": len(frame),
            "n_rows_by_cell": n_rows_by_cell,
            "n_nonzero_mission_quality": int(
                np.count_nonzero(frame["mission_quality"].to_numpy())
            ),
            "n_nonzero_qlp_quality": int(
                np.count_nonzero(frame["qlp_quality"].to_numpy())
            ),
            "detectors": [
                f"cam{camera}_ccd{ccd}" for camera, ccd in _EXPECTED_DETECTORS
            ],
            "authority_exclusions": sources.authority_exclusions,
            "authority_exclusions_sha256": sources.authority_exclusions_sha256,
            "quality_transfer_manifest_path": str(sources.transfer_manifest_path),
            "quality_transfer_manifest_sha256": sources.transfer_manifest_sha256,
            "quality_transfer_schema_version": MISSION_QUALITY_TRANSFER_SCHEMA_VERSION,
            "quality_transfer_producer_git_sha": sources.transfer_manifest[
                "producer_git_sha"
            ],
            "mission_quality_receipt_path": str(sources.quality_receipt_path),
            "mission_quality_receipt_sha256": sources.quality_receipt_sha256,
            "mission_quality_receipt_schema_version": (
                MISSION_QUALITY_RECEIPT_SCHEMA_VERSION
            ),
            "source_bindings": list(sources.source_bindings),
            "source_declaration_sha256": _canonical_json_sha256(
                list(sources.source_bindings)
            ),
            "hdf5_quality_join_required": True,
            "six_view_shards_verified": False,
            "panel_admission_authorized": False,
            "claim_limit": (
                "normalized mission/QLP cadence quality only; HDF5 provenance, "
                "six-view shards, identity, temporal-panel admission, and A2v1 "
                "acceptance remain separate gates"
            ),
        }
        manifest_path = partial / MISSION_QUALITY_REFERENCE_MANIFEST
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        (partial / MISSION_QUALITY_REFERENCE_READY).write_text(
            producer + "\n", encoding="utf-8"
        )
        load_mission_quality_reference(
            reference_dir=partial,
            sector=sources.sector,
            expected_orbits=sources.orbits,
        )
        _make_output_tree_read_only(partial)
        _assert_output_tree_read_only(partial)
        os.replace(partial, final)
        owns_partial = False
        owns_final = True
        _assert_output_tree_read_only(final)
        return manifest
    except BaseException:
        # Both paths were absent before this attempt. Clean only a tree whose
        # ownership was established here, including the narrow replace window.
        if owns_partial and not partial.exists() and final.exists():
            owns_partial = False
            owns_final = True
        cleanup_error: Exception | None = None
        for path, owned in ((partial, owns_partial), (final, owns_final)):
            if not owned:
                continue
            try:
                _remove_owned_output_tree(path)
            except (FM0ContractError, OSError) as exc:
                cleanup_error = exc
        if cleanup_error is not None:
            raise FM0ContractError(
                "quality-reference build failed and cleanup also failed: "
                f"{cleanup_error}"
            ) from cleanup_error
        raise


def _integer_column(frame: pd.DataFrame, column: str) -> np.ndarray:
    try:
        numeric = pd.to_numeric(frame[column], errors="raise").to_numpy(
            dtype=np.float64
        )
    except (TypeError, ValueError) as exc:
        raise FM0ContractError(f"quality reference {column} is not numeric") from exc
    if (
        numeric.ndim != 1
        or not np.all(np.isfinite(numeric))
        or not np.all(numeric == np.floor(numeric))
    ):
        raise FM0ContractError(f"quality reference {column} must contain integers")
    return numeric.astype(np.int64)


def _input_integer_array(values: Any, *, label: str) -> np.ndarray:
    array = np.asarray(values)
    if array.ndim != 1 or array.size == 0 or array.dtype == np.bool_:
        raise FM0ContractError(f"{label} must be a nonempty one-dimensional array")
    try:
        numeric = array.astype(np.float64)
    except (TypeError, ValueError) as exc:
        raise FM0ContractError(f"{label} must contain integers") from exc
    if not np.all(np.isfinite(numeric)) or not np.all(numeric == np.floor(numeric)):
        raise FM0ContractError(f"{label} must contain finite integers")
    return numeric.astype(np.int64)


@dataclass(frozen=True)
class MissionQualityLookupResult:
    """Provider-neutral external quality arrays for one HDF5 cadence vector."""

    mission_quality: np.ndarray
    qlp_quality: np.ndarray
    authority_excluded: np.ndarray

    @property
    def external_quality(self) -> np.ndarray:
        return np.bitwise_or(
            self.mission_quality.astype(np.uint64, copy=False),
            np.left_shift(
                self.qlp_quality.astype(np.uint64, copy=False),
                np.uint64(QLP_QUALITY_EXTERNAL_BIT),
            ),
        )

    @property
    def external_bad(self) -> np.ndarray:
        return (self.external_quality != 0) | self.authority_excluded


@dataclass(frozen=True)
class MissionQualityReference:
    """Validated immutable later-sector detector/cadence quality lookup."""

    sector: int
    provider: str
    expected_orbits: tuple[int, int]
    reference_dir: Path
    table_path: Path
    manifest_path: Path
    ready_path: Path
    table_sha256: str
    manifest_sha256: str
    ready_sha256: str
    transfer_manifest_path: Path
    transfer_manifest_sha256: str
    quality_receipt_path: Path
    quality_receipt_sha256: str
    source_file_sha256: Mapping[Path, str]
    source_declaration_sha256: str
    by_cell: Mapping[tuple[int, int, int], tuple[np.ndarray, np.ndarray, np.ndarray]]
    authority_exclusions: Mapping[tuple[int, int], tuple[np.ndarray, np.ndarray]]
    authority_exclusions_sha256: str

    @property
    def provenance(self) -> dict[str, Any]:
        return {
            "schema_version": MISSION_QUALITY_REFERENCE_SCHEMA_VERSION,
            "sector": self.sector,
            "mission_quality_provider": self.provider,
            "quality_policy": MISSION_QUALITY_POLICY,
            "cadence_authority": MISSION_QUALITY_CADENCE_AUTHORITY,
            "quality_composition": MISSION_QUALITY_COMPOSITION,
            "table_path": str(self.table_path),
            "table_sha256": self.table_sha256,
            "manifest_path": str(self.manifest_path),
            "manifest_sha256": self.manifest_sha256,
            "source_declaration_sha256": self.source_declaration_sha256,
            "authority_exclusions_sha256": self.authority_exclusions_sha256,
            "n_authority_exclusions": int(
                sum(len(values[0]) for values in self.authority_exclusions.values())
            ),
        }

    def assert_unchanged(self) -> None:
        """Fail if the reference or any bound upstream authority changed."""

        checks = {
            self.table_path: self.table_sha256,
            self.manifest_path: self.manifest_sha256,
            self.ready_path: self.ready_sha256,
            self.transfer_manifest_path: self.transfer_manifest_sha256,
            self.quality_receipt_path: self.quality_receipt_sha256,
            **self.source_file_sha256,
        }
        for path, expected in checks.items():
            if not path.is_file() or sha256_file(path) != expected:
                raise RuntimeError(f"mission-quality reference input changed: {path}")

    def lookup(
        self,
        *,
        orbit: int,
        camera: int,
        ccd: int,
        cadence: Any,
        context: str = "light curve",
    ) -> MissionQualityLookupResult:
        """Resolve provider-neutral quality, accepting only declared exclusions."""

        orbit = _positive_int(orbit, label=f"{context} orbit")
        camera = _positive_int(camera, label=f"{context} camera")
        ccd = _positive_int(ccd, label=f"{context} ccd")
        key = (orbit, camera, ccd)
        if orbit not in self.expected_orbits:
            raise FM0ContractError(f"{context}: orbit {orbit} is not in S{self.sector}")
        reference = self.by_cell.get(key)
        if reference is None:
            raise FM0ContractError(f"{context}: no quality mapping for cell {key}")
        cadences = _input_integer_array(cadence, label=f"{context} cadence")
        if np.any(cadences <= 0) or np.unique(cadences).size != cadences.size:
            raise FM0ContractError(f"{context}: cadences must be positive and unique")
        reference_cadence, reference_mission, reference_qlp = reference
        positions = np.searchsorted(reference_cadence, cadences)
        in_bounds = positions < len(reference_cadence)
        matched = np.zeros(len(cadences), dtype=bool)
        matched[in_bounds] = (
            reference_cadence[positions[in_bounds]] == cadences[in_bounds]
        )
        excluded = np.zeros(len(cadences), dtype=bool)
        mission = np.zeros(len(cadences), dtype=np.uint64)
        qlp = np.zeros(len(cadences), dtype=np.uint64)
        mission[matched] = reference_mission[positions[matched]]
        qlp[matched] = reference_qlp[positions[matched]]
        if not matched.all():
            exclusion_cadence, exclusion_mission = self.authority_exclusions[
                (camera, ccd)
            ]
            indices = np.flatnonzero(~matched)
            exclusion_positions = np.searchsorted(exclusion_cadence, cadences[indices])
            exclusion_in_bounds = exclusion_positions < len(exclusion_cadence)
            exclusion_matched = np.zeros(len(indices), dtype=bool)
            exclusion_matched[exclusion_in_bounds] = (
                exclusion_cadence[exclusion_positions[exclusion_in_bounds]]
                == cadences[indices[exclusion_in_bounds]]
            )
            accepted_indices = indices[exclusion_matched]
            excluded[accepted_indices] = True
            mission[accepted_indices] = exclusion_mission[
                exclusion_positions[exclusion_matched]
            ]
        uncovered = ~matched & ~excluded
        if uncovered.any():
            examples = cadences[uncovered][:20].astype(int).tolist()
            raise FM0ContractError(
                f"{context}: unexplained cadence coverage gap for {key}: {examples}"
            )
        return MissionQualityLookupResult(
            mission_quality=mission,
            qlp_quality=qlp,
            authority_excluded=excluded,
        )


def load_mission_quality_reference(
    *,
    reference_dir: str | Path,
    sector: int,
    expected_orbits: Sequence[int] | None = None,
) -> MissionQualityReference:
    """Load a v2 bundle and re-prove its transferred source authorities."""

    sector = _positive_int(sector, label="mission-quality reference sector")
    orbits = _expected_orbits(sector, expected_orbits)
    root = Path(reference_dir).expanduser().resolve()
    table_path = root / MISSION_QUALITY_REFERENCE_TABLE
    manifest_path = root / MISSION_QUALITY_REFERENCE_MANIFEST
    ready_path = root / MISSION_QUALITY_REFERENCE_READY
    if (
        not table_path.is_file()
        or not manifest_path.is_file()
        or not ready_path.is_file()
    ):
        raise FM0ContractError(
            f"mission-quality reference bundle is incomplete: {root}"
        )
    table_hash = sha256_file(table_path)
    manifest_hash = sha256_file(manifest_path)
    ready_hash = sha256_file(ready_path)
    manifest = _load_json_object(
        manifest_path, label="mission-quality reference manifest"
    )
    missing = sorted(_REFERENCE_MANIFEST_FIELDS - set(manifest))
    unknown = sorted(set(manifest) - _REFERENCE_MANIFEST_FIELDS)
    if missing or unknown:
        raise FM0ContractError(
            f"mission-quality reference manifest fields drifted; "
            f"missing={missing}, unknown={unknown}"
        )
    producer = _digest(
        manifest["producer_git_sha"], label="reference producer_git_sha", length=40
    )
    try:
        generated = datetime.fromisoformat(str(manifest["generated_utc"]))
    except (TypeError, ValueError) as exc:
        raise FM0ContractError("reference generated_utc is invalid") from exc
    if generated.tzinfo is None:
        raise FM0ContractError("reference generated_utc must be timezone-aware")
    if (
        manifest["schema_version"] != MISSION_QUALITY_REFERENCE_SCHEMA_VERSION
        or manifest["builder_version"] != MISSION_QUALITY_REFERENCE_BUILDER_VERSION
        or int(manifest["sector"]) != sector
        or tuple(int(value) for value in manifest["expected_orbits"]) != orbits
        or manifest["mission_quality_provider"] != mission_quality_type(sector)
        or manifest["quality_policy"] != MISSION_QUALITY_POLICY
        or manifest["cadence_authority"] != MISSION_QUALITY_CADENCE_AUTHORITY
        or manifest["quality_composition"] != MISSION_QUALITY_COMPOSITION
        or manifest["table_file"] != MISSION_QUALITY_REFERENCE_TABLE
        or manifest["table_columns"] != list(MISSION_QUALITY_REFERENCE_COLUMNS)
        or manifest["table_sha256"] != table_hash
        or ready_path.read_text(encoding="utf-8").strip() != producer
        or manifest["hdf5_quality_join_required"] is not True
        or manifest["six_view_shards_verified"] is not False
        or manifest["panel_admission_authorized"] is not False
    ):
        raise FM0ContractError("mission-quality reference manifest violates v2")

    transfer_declared = Path(
        str(manifest["quality_transfer_manifest_path"])
    ).expanduser()
    if not transfer_declared.is_absolute():
        raise FM0ContractError("reference transfer manifest path must be absolute")
    transfer_path = transfer_declared.resolve()
    if not transfer_path.is_file() or transfer_path.name != "manifest.json":
        raise FM0ContractError("reference transfer manifest path is invalid")
    if manifest[
        "quality_transfer_schema_version"
    ] != MISSION_QUALITY_TRANSFER_SCHEMA_VERSION or sha256_file(
        transfer_path
    ) != _digest(
        manifest["quality_transfer_manifest_sha256"],
        label="reference transfer manifest sha256",
    ):
        raise FM0ContractError("reference transfer manifest binding drifted")
    sources = _verify_sector_transfer(
        quality_transfer_root=transfer_path.parent, sector=sector
    )
    receipt_declared = Path(str(manifest["mission_quality_receipt_path"])).expanduser()
    if not receipt_declared.is_absolute():
        raise FM0ContractError(
            "reference mission-quality receipt path must be absolute"
        )
    receipt_path = receipt_declared.resolve()
    if (
        receipt_path != sources.quality_receipt_path
        or manifest["mission_quality_receipt_schema_version"]
        != MISSION_QUALITY_RECEIPT_SCHEMA_VERSION
        or _digest(
            manifest["mission_quality_receipt_sha256"],
            label="reference mission-quality receipt sha256",
        )
        != sources.quality_receipt_sha256
        or manifest["quality_transfer_producer_git_sha"]
        != sources.transfer_manifest["producer_git_sha"]
    ):
        raise FM0ContractError("reference mission-quality upstream binding drifted")
    source_bindings = list(sources.source_bindings)
    if manifest["source_bindings"] != source_bindings or manifest[
        "source_declaration_sha256"
    ] != _canonical_json_sha256(source_bindings):
        raise FM0ContractError("reference source declaration drifted")

    frame = pd.read_csv(table_path, low_memory=False)
    if tuple(frame.columns) != MISSION_QUALITY_REFERENCE_COLUMNS:
        raise FM0ContractError("mission-quality reference table columns drifted")
    for column in MISSION_QUALITY_REFERENCE_COLUMNS:
        frame[column] = _integer_column(frame, column)
    expected_frame = _quality_frame(sources)
    if frame.shape != expected_frame.shape or not np.array_equal(
        frame.to_numpy(dtype=np.int64), expected_frame.to_numpy(dtype=np.int64)
    ):
        raise FM0ContractError(
            "mission-quality reference table disagrees with transferred authorities"
        )
    if _positive_int(manifest["n_rows"], label="reference n_rows") != len(frame):
        raise FM0ContractError("mission-quality reference row count mismatch")
    rows_by_cell = {
        f"o{int(orbit)}_cam{int(camera)}_ccd{int(ccd)}": len(group)
        for (orbit, camera, ccd), group in frame.groupby(
            ["orbit", "camera", "ccd"], sort=True
        )
    }
    if manifest["n_rows_by_cell"] != rows_by_cell:
        raise FM0ContractError("mission-quality reference cell counts mismatch")
    if manifest["detectors"] != [
        f"cam{camera}_ccd{ccd}" for camera, ccd in _EXPECTED_DETECTORS
    ]:
        raise FM0ContractError("mission-quality reference detectors drifted")
    for field, column in (
        ("n_nonzero_mission_quality", "mission_quality"),
        ("n_nonzero_qlp_quality", "qlp_quality"),
    ):
        if _positive_int(
            manifest[field], label=f"reference {field}", allow_zero=True
        ) != int(np.count_nonzero(frame[column].to_numpy())):
            raise FM0ContractError(f"mission-quality reference {field} mismatch")

    if (
        manifest["authority_exclusions"] != sources.authority_exclusions
        or manifest["authority_exclusions_sha256"]
        != sources.authority_exclusions_sha256
        or _canonical_json_sha256(manifest["authority_exclusions"])
        != sources.authority_exclusions_sha256
    ):
        raise FM0ContractError("reference authority exclusion binding drifted")
    exclusions = _validate_exclusions(
        payload=manifest["authority_exclusions"], sector=sector, frame=frame
    )

    by_cell: dict[tuple[int, int, int], tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for (orbit, camera, ccd), group in frame.groupby(
        ["orbit", "camera", "ccd"], sort=True
    ):
        ordered = group.sort_values("cadence", kind="stable")
        arrays = (
            ordered["cadence"].to_numpy(dtype=np.int64),
            ordered["mission_quality"].to_numpy(dtype=np.uint64),
            ordered["qlp_quality"].to_numpy(dtype=np.uint64),
        )
        for array in arrays:
            array.setflags(write=False)
        by_cell[(int(orbit), int(camera), int(ccd))] = arrays
    if set(by_cell) != {
        (orbit, camera, ccd) for orbit in orbits for camera, ccd in _EXPECTED_DETECTORS
    }:
        raise FM0ContractError("mission-quality reference cell inventory is incomplete")

    if (
        sha256_file(table_path) != table_hash
        or sha256_file(manifest_path) != manifest_hash
    ):
        raise RuntimeError("mission-quality reference changed while loading")
    source_files = {
        Path(str(binding["path"])): str(binding["sha256"])
        for binding in source_bindings
    }
    reference = MissionQualityReference(
        sector=sector,
        provider=str(manifest["mission_quality_provider"]),
        expected_orbits=orbits,
        reference_dir=root,
        table_path=table_path,
        manifest_path=manifest_path,
        ready_path=ready_path,
        table_sha256=table_hash,
        manifest_sha256=manifest_hash,
        ready_sha256=ready_hash,
        transfer_manifest_path=transfer_path,
        transfer_manifest_sha256=sources.transfer_manifest_sha256,
        quality_receipt_path=receipt_path,
        quality_receipt_sha256=sources.quality_receipt_sha256,
        source_file_sha256=source_files,
        source_declaration_sha256=str(manifest["source_declaration_sha256"]),
        by_cell=by_cell,
        authority_exclusions=exclusions,
        authority_exclusions_sha256=sources.authority_exclusions_sha256,
    )
    reference.assert_unchanged()
    return reference


__all__ = [
    "MISSION_QUALITY_CADENCE_AUTHORITY",
    "MISSION_QUALITY_COMPOSITION",
    "MISSION_QUALITY_REFERENCE_BUILDER_VERSION",
    "MISSION_QUALITY_REFERENCE_COLUMNS",
    "MISSION_QUALITY_REFERENCE_SCHEMA_VERSION",
    "MissionQualityLookupResult",
    "MissionQualityReference",
    "build_mission_quality_reference",
    "load_mission_quality_reference",
]
