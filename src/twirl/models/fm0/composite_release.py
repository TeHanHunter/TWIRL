"""Identity-only FM0.3 composite release over existing immutable shards.

The freezer binds the S56--S64 release plus the twelve S66--S77 sector
releases. It writes only a manifest-of-manifests, a compact observation-role
index, and a receipt. It never reads or rewrites shard payloads.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import tempfile
from collections import OrderedDict, defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from io import StringIO
from pathlib import Path, PurePosixPath
from typing import Any, Literal

import numpy as np

from .dataset import (
    FM0ReleaseDataset,
    _manifest_view_presence,
    variant_view_indices,
)
from .input_release import INPUT_RELEASE_SCHEMA_VERSION, MANIFEST_COLUMNS
from .later_sector_admission_v2 import (
    POOL_RECEIPT_SCHEMA_VERSION,
    PREPARATION_SECTORS,
)
from .later_sector_release import (
    LATER_SIX_VIEW_MANIFEST_FIELDS,
    validate_later_sector_release,
)
from .registry import FM0ContractError, sha256_file

COMPOSITE_RELEASE_SCHEMA_VERSION = "twirl_fm0_3_composite_release_v1"
COMPOSITE_RELEASE_STATE = "FM0_3_COMPOSITE_IDENTITY_READY"
SOURCE_BINDING_SCHEMA_VERSION = "twirl_fm0_3_source_manifest_binding_v1"
ROLE_INDEX_SCHEMA_VERSION = "twirl_fm0_3_observation_role_index_v1"

LEGACY_SECTORS = tuple(range(56, 65))
LATER_SECTORS = tuple(PREPARATION_SECTORS)
INCLUDED_SECTORS = (*LEGACY_SECTORS, *LATER_SECTORS)
TEMPORAL_HOLDOUT_SECTOR = 77

TRAIN_ROLE = "fm03_train"
HOLDOUT_ROLE = "temporal_holdout"
EXCLUDED_OVERLAP_ROLE = "excluded_temporal_overlap"
ROLE_ORDER = (TRAIN_ROLE, HOLDOUT_ROLE, EXCLUDED_OVERLAP_ROLE)

SOURCE_BINDING_FIELDS = (
    "schema_version",
    "source_id",
    "source_kind",
    "sector_min",
    "sector_max",
    "release_root",
    "manifest_sha256",
    "receipt_path",
    "receipt_sha256",
    "n_rows",
    "n_cadences",
    "rows_by_partition_json",
)
ROLE_INDEX_FIELDS = (
    "schema_version",
    "source_id",
    "observation_key",
    "leakage_component_id",
    "role",
)

_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_GIT_SHA = re.compile(r"^[0-9a-f]{40}$")
_COMPONENT = re.compile(r"^leakage_[0-9a-f]{64}$")
_OBSERVATION = re.compile(r"^observation_[0-9a-f]{64}$")
_PARTITIONS = ("poc_development", "poc_sealed_test", "poc_train")
_FM03_VARIANTS = frozenset(("TWIRL-FM0.3.1", "TWIRL-FM0.3.2"))


@dataclass(frozen=True, slots=True)
class SourceManifestBinding:
    """One exact source release binding supplied by the remote freeze step."""

    source_kind: Literal["legacy", "later"]
    sector_min: int
    sector_max: int
    release_root: Path
    manifest_sha256: str
    receipt_path: Path
    receipt_sha256: str

    @property
    def source_id(self) -> str:
        return (
            "legacy_s56_s64"
            if self.source_kind == "legacy"
            else f"later_s{self.sector_min:04d}"
        )


@dataclass(frozen=True, slots=True)
class CompositeObservation:
    source_id: str
    sector: int | None
    observation_key: str
    component: str
    relative_path: str
    shard_sha256: str
    n_cadences: int
    view_present_json: str


@dataclass(frozen=True, slots=True)
class CompositeReleaseResult:
    root: Path
    source_bindings_sha256: str
    role_index_sha256: str
    receipt_sha256: str
    receipt: Mapping[str, Any]


@dataclass(frozen=True, slots=True)
class CompositeReleaseValidation:
    result: CompositeReleaseResult
    observations_by_role: Mapping[str, tuple[CompositeObservation, ...]]


@dataclass(frozen=True, slots=True)
class _SourceScan:
    binding: SourceManifestBinding
    binding_row: Mapping[str, str]
    train_rows: tuple[CompositeObservation, ...]
    n_rows: int
    n_cadences: int
    rows_by_partition: Mapping[str, int]


def _digest(value: Any, label: str) -> str:
    text = str(value).strip().lower()
    if _SHA256.fullmatch(text) is None:
        raise FM0ContractError(f"{label} must be a lowercase SHA-256")
    return text


def _exact_positive_int(value: Any, label: str) -> int:
    text = str(value).strip()
    try:
        result = int(text)
    except ValueError as exc:
        raise FM0ContractError(f"{label} must be a positive integer") from exc
    if text != str(result) or result <= 0:
        raise FM0ContractError(f"{label} must be a positive integer")
    return result


def _bound_file(
    path: str | Path,
    digest: str,
    label: str,
    *,
    require_read_only: bool,
) -> Path:
    raw = Path(path).expanduser()
    expected = _digest(digest, f"{label} hash")
    if raw.is_symlink():
        raise FM0ContractError(f"{label} must not be a symlink")
    resolved = raw.resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise FM0ContractError(f"{label} is missing or empty: {resolved}")
    if require_read_only and resolved.stat().st_mode & 0o222:
        raise FM0ContractError(f"{label} is writable: {resolved}")
    if sha256_file(resolved) != expected:
        raise FM0ContractError(f"{label} hash drifted")
    return resolved


def _json(path: Path, label: str) -> Mapping[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FM0ContractError(f"invalid {label}: {path}") from exc
    if not isinstance(payload, Mapping):
        raise FM0ContractError(f"{label} must be an object")
    return payload


def _csv_payload(rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> bytes:
    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(fields), lineterminator="\n")
    writer.writeheader()
    for row in rows:
        if tuple(row) != tuple(fields):
            raise FM0ContractError("composite CSV columns drifted")
        writer.writerow(row)
    return stream.getvalue().encode()


def _read_csv(path: Path, fields: Sequence[str], label: str) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise FM0ContractError(f"{label} columns drifted")
        return [dict(row) for row in reader]


def _view_json(value: Any) -> str:
    try:
        parsed = json.loads(str(value))
    except json.JSONDecodeError as exc:
        raise FM0ContractError("invalid view_present_json") from exc
    if (
        not isinstance(parsed, list)
        or len(parsed) != 6
        or any(type(item) is not int or item not in (0, 1) for item in parsed)
    ):
        raise FM0ContractError("invalid view_present_json")
    return json.dumps(parsed, separators=(",", ":"))


def _relative_path(value: Any, observation_key: str) -> str:
    path = PurePosixPath(str(value))
    if (
        path.is_absolute()
        or ".." in path.parts
        or path.parts != ("shards", f"{observation_key}.npz")
    ):
        raise FM0ContractError("source shard path is noncanonical")
    return path.as_posix()


def _observation(
    row: Mapping[str, str], *, source_id: str, sector: int | None
) -> CompositeObservation:
    observation_key = str(row["observation_key"])
    component = str(row["leakage_component_id"])
    if (
        _OBSERVATION.fullmatch(observation_key) is None
        or _COMPONENT.fullmatch(component) is None
    ):
        raise FM0ContractError("invalid observation/component identity")
    return CompositeObservation(
        source_id=source_id,
        sector=sector,
        observation_key=observation_key,
        component=component,
        relative_path=_relative_path(row["relative_path"], observation_key),
        shard_sha256=_digest(row["sha256"], "source shard hash"),
        n_cadences=_exact_positive_int(row["n_cadences"], "cadence count"),
        view_present_json=_view_json(row["view_present_json"]),
    )


def _binding_shape(binding: SourceManifestBinding) -> None:
    if binding.source_kind == "legacy":
        valid = (binding.sector_min, binding.sector_max) == (56, 64)
    else:
        valid = (
            binding.source_kind == "later"
            and binding.sector_min == binding.sector_max
            and binding.sector_min in LATER_SECTORS
        )
    if not valid:
        raise FM0ContractError("source binding sector/kind drifted")
    _digest(binding.manifest_sha256, "source manifest hash")
    _digest(binding.receipt_sha256, "source receipt hash")


def _scan_legacy(
    binding: SourceManifestBinding, *, require_read_only: bool
) -> tuple[list[dict[str, str]], Path, Path]:
    root = Path(binding.release_root).expanduser().resolve()
    manifest = _bound_file(
        root / "manifest.csv",
        binding.manifest_sha256,
        "legacy manifest",
        require_read_only=require_read_only,
    )
    receipt_path = _bound_file(
        binding.receipt_path,
        binding.receipt_sha256,
        "legacy receipt",
        require_read_only=require_read_only,
    )
    rows = _read_csv(manifest, MANIFEST_COLUMNS, "legacy manifest")
    receipt = _json(receipt_path, "legacy receipt")
    release = receipt.get("release")
    if (
        receipt.get("schema_version") != "twirl_fm0_1_input_release_receipt_v1"
        or receipt.get("passed") is not True
        or receipt.get("scientific_training_eligible") is not True
        or not isinstance(release, Mapping)
        or release.get("manifest_sha256") != binding.manifest_sha256
    ):
        raise FM0ContractError("legacy receipt contract drifted")
    return rows, manifest, receipt_path


def _scan_later(
    binding: SourceManifestBinding, *, require_read_only: bool
) -> tuple[list[dict[str, str]], Path, Path]:
    root = Path(binding.release_root).expanduser().resolve()
    receipt_path, receipt, bundle = validate_later_sector_release(
        root,
        expected_receipt_sha256=binding.receipt_sha256,
        require_read_only=require_read_only,
        verify_shard_payloads=False,
    )
    manifest = root / "manifest.csv"
    if (
        int(bundle["sector"]) != binding.sector_min
        or receipt["manifest_sha256"] != binding.manifest_sha256
        or sha256_file(manifest) != binding.manifest_sha256
    ):
        raise FM0ContractError("later source manifest binding drifted")
    return [dict(row) for row in bundle["manifest_rows"]], manifest, receipt_path


def _scan_source(
    binding: SourceManifestBinding, *, require_read_only: bool
) -> _SourceScan:
    _binding_shape(binding)
    rows, manifest, receipt_path = (
        _scan_legacy(binding, require_read_only=require_read_only)
        if binding.source_kind == "legacy"
        else _scan_later(binding, require_read_only=require_read_only)
    )
    if not rows:
        raise FM0ContractError("source manifest is empty")
    counts = {partition: 0 for partition in _PARTITIONS}
    train_rows: list[CompositeObservation] = []
    n_cadences = 0
    seen: set[str] = set()
    expected_fields = (
        MANIFEST_COLUMNS
        if binding.source_kind == "legacy"
        else LATER_SIX_VIEW_MANIFEST_FIELDS
    )
    for row in rows:
        if tuple(row) != tuple(expected_fields):
            raise FM0ContractError("source manifest row columns drifted")
        partition = str(row["source_partition"])
        if partition not in counts:
            raise FM0ContractError("source partition drifted")
        if binding.source_kind == "legacy":
            if (
                row["input_release_schema_version"] != INPUT_RELEASE_SCHEMA_VERSION
                or row["product_state"] != "A2V1_ACCEPTED"
                or row["input_adapter"] != "a2v1_hdf5_quality_aware_v1"
            ):
                raise FM0ContractError("legacy source row contract drifted")
            sector = None
        else:
            if (
                int(row["sector"]) != binding.sector_min
                or row["scientific_training_eligible"] != "false"
                or row["panel_admission_authorized"] != "false"
            ):
                raise FM0ContractError("later source row contract drifted")
            sector = binding.sector_min
        item = _observation(row, source_id=binding.source_id, sector=sector)
        if item.observation_key in seen:
            raise FM0ContractError("duplicate observation within source manifest")
        seen.add(item.observation_key)
        counts[partition] += 1
        n_cadences += item.n_cadences
        if partition == "poc_train":
            if (
                binding.source_kind == "legacy"
                and row["scientific_training_eligible"] != "True"
            ):
                raise FM0ContractError("legacy poc_train row is ineligible")
            train_rows.append(item)
    if binding.source_kind == "later" and counts["poc_sealed_test"]:
        raise FM0ContractError("later six-view manifest contains sealed rows")
    canonical = SourceManifestBinding(
        source_kind=binding.source_kind,
        sector_min=binding.sector_min,
        sector_max=binding.sector_max,
        release_root=manifest.parent,
        manifest_sha256=binding.manifest_sha256,
        receipt_path=receipt_path,
        receipt_sha256=binding.receipt_sha256,
    )
    binding_row = {
        "schema_version": SOURCE_BINDING_SCHEMA_VERSION,
        "source_id": canonical.source_id,
        "source_kind": canonical.source_kind,
        "sector_min": str(canonical.sector_min),
        "sector_max": str(canonical.sector_max),
        "release_root": str(canonical.release_root),
        "manifest_sha256": canonical.manifest_sha256,
        "receipt_path": str(canonical.receipt_path),
        "receipt_sha256": canonical.receipt_sha256,
        "n_rows": str(len(rows)),
        "n_cadences": str(n_cadences),
        "rows_by_partition_json": json.dumps(
            counts, sort_keys=True, separators=(",", ":")
        ),
    }
    return _SourceScan(
        binding=canonical,
        binding_row=binding_row,
        train_rows=tuple(train_rows),
        n_rows=len(rows),
        n_cadences=n_cadences,
        rows_by_partition=counts,
    )


def _ordered_bindings(
    legacy: SourceManifestBinding, later: Sequence[SourceManifestBinding]
) -> tuple[SourceManifestBinding, ...]:
    if legacy.source_kind != "legacy":
        raise FM0ContractError("first source must be legacy S56--S64")
    _binding_shape(legacy)
    later = tuple(later)
    if tuple(item.sector_min for item in later) != LATER_SECTORS:
        raise FM0ContractError("later sources must be chronological S66--S77")
    for item in later:
        _binding_shape(item)
    return (legacy, *later)


def _assign_roles(
    scans: Sequence[_SourceScan],
) -> tuple[list[dict[str, str]], dict[str, tuple[CompositeObservation, ...]]]:
    observations = [item for scan in scans for item in scan.train_rows]
    keys = [item.observation_key for item in observations]
    if len(keys) != len(set(keys)):
        raise FM0ContractError("observation keys collide across source manifests")
    holdout_components = {
        item.component
        for item in observations
        if item.sector == TEMPORAL_HOLDOUT_SECTOR
    }
    if not holdout_components:
        raise FM0ContractError("S77 has no poc_train holdout components")
    grouped: dict[str, list[CompositeObservation]] = {role: [] for role in ROLE_ORDER}
    for item in observations:
        if item.sector == TEMPORAL_HOLDOUT_SECTOR:
            role = HOLDOUT_ROLE
        elif item.component in holdout_components:
            role = EXCLUDED_OVERLAP_ROLE
        else:
            role = TRAIN_ROLE
        grouped[role].append(item)
    if {item.component for item in grouped[TRAIN_ROLE]} & holdout_components:
        raise FM0ContractError("training and temporal holdout components overlap")
    role_rows = [
        {
            "schema_version": ROLE_INDEX_SCHEMA_VERSION,
            "source_id": item.source_id,
            "observation_key": item.observation_key,
            "leakage_component_id": item.component,
            "role": role,
        }
        for role in ROLE_ORDER
        for item in sorted(
            grouped[role], key=lambda row: (row.source_id, row.observation_key)
        )
    ]
    return role_rows, {role: tuple(grouped[role]) for role in ROLE_ORDER}


def _admission(
    path: str | Path,
    digest: str,
    scans: Sequence[_SourceScan],
    *,
    require_read_only: bool,
) -> tuple[Path, str]:
    receipt_path = _bound_file(
        path, digest, "later admission receipt", require_read_only=require_read_only
    )
    receipt = _json(receipt_path, "later admission receipt")
    later = [scan for scan in scans if scan.binding.source_kind == "later"]
    n_rows = sum(scan.n_rows for scan in later)
    if (
        receipt.get("schema_version") != POOL_RECEIPT_SCHEMA_VERSION
        or receipt.get("preparation_pool_sectors") != list(LATER_SECTORS)
        or receipt.get("excluded_sectors") != [65]
        or receipt.get("preparation_pool_admitted") is not True
        or receipt.get("n_nonsealed_preparation_rows") != n_rows
        or receipt.get("n_six_view_shards") != n_rows
        or receipt.get("scientific_training_eligible") is not False
        or receipt.get("model_training_authorized") is not False
        or receipt.get("sealed_test_access_authorized") is not False
    ):
        raise FM0ContractError("later admission receipt contract drifted")
    sector_receipts = receipt.get("sector_receipts")
    if not isinstance(sector_receipts, list) or len(sector_receipts) != 12:
        raise FM0ContractError("later admission sector closure drifted")
    by_sector = {scan.binding.sector_min: scan for scan in later}
    for expected_sector, item in zip(LATER_SECTORS, sector_receipts, strict=True):
        if not isinstance(item, Mapping) or item.get("sector") != expected_sector:
            raise FM0ContractError("later admission sector order drifted")
        evidence = item.get("evidence")
        six_view = (
            evidence.get("six_view_release") if isinstance(evidence, Mapping) else None
        )
        scan = by_sector[expected_sector]
        if (
            not isinstance(six_view, Mapping)
            or six_view.get("sha256") != scan.binding.receipt_sha256
            or six_view.get("manifest_sha256") != scan.binding.manifest_sha256
        ):
            raise FM0ContractError("later admission source binding drifted")
    return receipt_path, sha256_file(receipt_path)


def _products(
    *,
    bindings: Sequence[SourceManifestBinding],
    admission_path: str | Path,
    admission_sha256: str,
    producer_git_sha: str,
    require_source_read_only: bool,
) -> tuple[
    bytes, bytes, dict[str, Any], Mapping[str, tuple[CompositeObservation, ...]]
]:
    if _GIT_SHA.fullmatch(producer_git_sha) is None:
        raise FM0ContractError("producer_git_sha must be a full Git SHA")
    scans = [
        _scan_source(binding, require_read_only=require_source_read_only)
        for binding in bindings
    ]
    admission_path, admission_hash = _admission(
        admission_path,
        admission_sha256,
        scans,
        require_read_only=require_source_read_only,
    )
    source_payload = _csv_payload(
        [scan.binding_row for scan in scans], SOURCE_BINDING_FIELDS
    )
    role_rows, by_role = _assign_roles(scans)
    role_payload = _csv_payload(role_rows, ROLE_INDEX_FIELDS)
    role_counts = {role: len(by_role[role]) for role in ROLE_ORDER}
    role_cadences = {
        role: sum(item.n_cadences for item in by_role[role]) for role in ROLE_ORDER
    }
    receipt = {
        "schema_version": COMPOSITE_RELEASE_SCHEMA_VERSION,
        "release_state": COMPOSITE_RELEASE_STATE,
        "passed": True,
        "producer_git_sha": producer_git_sha,
        "sources": {
            "included_sectors": list(INCLUDED_SECTORS),
            "excluded_sectors": [65],
            "n_manifests": 13,
            "n_rows": sum(scan.n_rows for scan in scans),
            "n_cadences": sum(scan.n_cadences for scan in scans),
            "rows_by_partition": {
                partition: sum(scan.rows_by_partition[partition] for scan in scans)
                for partition in _PARTITIONS
            },
            "source_bindings_sha256": hashlib.sha256(source_payload).hexdigest(),
            "later_admission_receipt_path": str(admission_path),
            "later_admission_receipt_sha256": admission_hash,
        },
        "selection": {
            "temporal_holdout_sector": TEMPORAL_HOLDOUT_SECTOR,
            "role_index_sha256": hashlib.sha256(role_payload).hexdigest(),
            "role_counts": role_counts,
            "role_cadences": role_cadences,
            "pre_quarantine_training_rows": (
                role_counts[TRAIN_ROLE] + role_counts[EXCLUDED_OVERLAP_ROLE]
            ),
            "pre_quarantine_training_cadences": (
                role_cadences[TRAIN_ROLE] + role_cadences[EXCLUDED_OVERLAP_ROLE]
            ),
            "n_training_components": len(
                {item.component for item in by_role[TRAIN_ROLE]}
            ),
            "n_holdout_components": len(
                {item.component for item in by_role[HOLDOUT_ROLE]}
            ),
            "n_excluded_overlap_components": len(
                {item.component for item in by_role[EXCLUDED_OVERLAP_ROLE]}
            ),
        },
        "limits": {
            "identity_only": True,
            "source_shards_opened": False,
            "shards_rewritten": False,
            "development_rows_selected": 0,
            "sealed_rows_selected": 0,
            "scientific_training_eligible": True,
            "model_training_authorized": False,
            "real_cli_training_enabled": False,
            "sealed_test_access_authorized": False,
        },
        "claim_limit": (
            "identity-frozen FM0.3 selection only; not a training run, model "
            "result, Stage-1 acceptance, or sealed test"
        ),
    }
    return source_payload, role_payload, receipt, by_role


def freeze_composite_release(
    output_dir: str | Path,
    *,
    legacy_binding: SourceManifestBinding,
    later_bindings: Sequence[SourceManifestBinding],
    admission_receipt_path: str | Path,
    admission_receipt_sha256: str,
    producer_git_sha: str,
    require_source_read_only: bool = True,
) -> CompositeReleaseResult:
    """Freeze the composite identity selection without opening any shard."""

    bindings = _ordered_bindings(legacy_binding, later_bindings)
    source_payload, role_payload, receipt, _by_role = _products(
        bindings=bindings,
        admission_path=admission_receipt_path,
        admission_sha256=admission_receipt_sha256,
        producer_git_sha=producer_git_sha,
        require_source_read_only=require_source_read_only,
    )
    receipt_payload = (
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode()
    receipt_hash = hashlib.sha256(receipt_payload).hexdigest()
    output = Path(output_dir).expanduser().resolve()
    if output.exists():
        return validate_composite_release(
            output,
            expected_receipt_sha256=receipt_hash,
            expected_source_bindings_sha256=receipt["sources"][
                "source_bindings_sha256"
            ],
            expected_role_index_sha256=receipt["selection"]["role_index_sha256"],
            require_source_read_only=require_source_read_only,
        ).result
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = Path(
        tempfile.mkdtemp(prefix=f".{output.name}.partial.", dir=output.parent)
    )
    try:
        for name, payload in {
            "source_bindings.csv": source_payload,
            "role_index.csv": role_payload,
            "receipt.json": receipt_payload,
            "READY": f"{receipt_hash}\n".encode(),
        }.items():
            path = partial / name
            path.write_bytes(payload)
            path.chmod(0o444)
        partial.chmod(0o555)
        os.replace(partial, output)
    except BaseException:
        if partial.exists():
            partial.chmod(0o755)
            shutil.rmtree(partial)
        raise
    return CompositeReleaseResult(
        root=output,
        source_bindings_sha256=receipt["sources"]["source_bindings_sha256"],
        role_index_sha256=receipt["selection"]["role_index_sha256"],
        receipt_sha256=receipt_hash,
        receipt=receipt,
    )


def _binding_from_row(row: Mapping[str, str]) -> SourceManifestBinding:
    if (
        tuple(row) != SOURCE_BINDING_FIELDS
        or row["schema_version"] != SOURCE_BINDING_SCHEMA_VERSION
    ):
        raise FM0ContractError("source-binding row drifted")
    binding = SourceManifestBinding(
        source_kind=row["source_kind"],
        sector_min=int(row["sector_min"]),
        sector_max=int(row["sector_max"]),
        release_root=Path(row["release_root"]),
        manifest_sha256=row["manifest_sha256"],
        receipt_path=Path(row["receipt_path"]),
        receipt_sha256=row["receipt_sha256"],
    )
    _binding_shape(binding)
    if row["source_id"] != binding.source_id:
        raise FM0ContractError("source-binding identity drifted")
    return binding


def validate_composite_release(
    root: str | Path,
    *,
    expected_receipt_sha256: str,
    expected_source_bindings_sha256: str,
    expected_role_index_sha256: str,
    require_read_only: bool = True,
    require_source_read_only: bool = True,
) -> CompositeReleaseValidation:
    """Re-derive the role index from the 13 bound source manifests."""

    root = Path(root).expanduser().resolve()
    if (
        root.is_symlink()
        or not root.is_dir()
        or {path.name for path in root.iterdir()}
        != {"source_bindings.csv", "role_index.csv", "receipt.json", "READY"}
    ):
        raise FM0ContractError("composite root closure drifted")
    receipt_path = _bound_file(
        root / "receipt.json",
        expected_receipt_sha256,
        "composite receipt",
        require_read_only=require_read_only,
    )
    source_path = _bound_file(
        root / "source_bindings.csv",
        expected_source_bindings_sha256,
        "source bindings",
        require_read_only=require_read_only,
    )
    role_path = _bound_file(
        root / "role_index.csv",
        expected_role_index_sha256,
        "role index",
        require_read_only=require_read_only,
    )
    receipt = _json(receipt_path, "composite receipt")
    if (
        set(receipt)
        != {
            "schema_version",
            "release_state",
            "passed",
            "producer_git_sha",
            "sources",
            "selection",
            "limits",
            "claim_limit",
        }
        or receipt.get("schema_version") != COMPOSITE_RELEASE_SCHEMA_VERSION
        or receipt.get("release_state") != COMPOSITE_RELEASE_STATE
        or receipt.get("passed") is not True
        or receipt.get("limits")
        != {
            "identity_only": True,
            "source_shards_opened": False,
            "shards_rewritten": False,
            "development_rows_selected": 0,
            "sealed_rows_selected": 0,
            "scientific_training_eligible": True,
            "model_training_authorized": False,
            "real_cli_training_enabled": False,
            "sealed_test_access_authorized": False,
        }
    ):
        raise FM0ContractError("composite receipt contract drifted")
    if (
        receipt["sources"]["source_bindings_sha256"] != expected_source_bindings_sha256
        or receipt["selection"]["role_index_sha256"] != expected_role_index_sha256
        or (root / "READY").read_text(encoding="utf-8").strip()
        != expected_receipt_sha256
        or (require_read_only and root.stat().st_mode & 0o222)
    ):
        raise FM0ContractError("composite artifact binding drifted")
    source_rows = _read_csv(source_path, SOURCE_BINDING_FIELDS, "source bindings")
    if len(source_rows) != 13:
        raise FM0ContractError("composite release must bind 13 source manifests")
    bindings = tuple(_binding_from_row(row) for row in source_rows)
    _ordered_bindings(bindings[0], bindings[1:])
    source_payload, role_payload, expected_receipt, by_role = _products(
        bindings=bindings,
        admission_path=receipt["sources"]["later_admission_receipt_path"],
        admission_sha256=receipt["sources"]["later_admission_receipt_sha256"],
        producer_git_sha=receipt["producer_git_sha"],
        require_source_read_only=require_source_read_only,
    )
    if source_payload != source_path.read_bytes():
        raise FM0ContractError("source-binding content drifted")
    if role_payload != role_path.read_bytes():
        raise FM0ContractError("role-index selection drifted")
    if expected_receipt != receipt:
        raise FM0ContractError("composite aggregate receipt drifted")
    return CompositeReleaseValidation(
        result=CompositeReleaseResult(
            root=root,
            source_bindings_sha256=expected_source_bindings_sha256,
            role_index_sha256=expected_role_index_sha256,
            receipt_sha256=expected_receipt_sha256,
            receipt=receipt,
        ),
        observations_by_role=by_role,
    )


@dataclass(frozen=True, slots=True)
class FM0CompositeDatasetConfig:
    composite_root: str
    receipt_sha256: str
    source_bindings_sha256: str
    role_index_sha256: str
    variant: str
    role: str = TRAIN_ROLE
    seed: int = 560067
    windows_per_epoch: int = 1_280_000
    window_length: int = 128
    shard_cache_size: int = 8
    mask_target_fraction: float = 0.15
    mask_span_range: tuple[int, int] = (1, 4)
    require_read_only: bool = True

    def __post_init__(self) -> None:
        if self.variant not in _FM03_VARIANTS:
            raise ValueError("composite dataset supports only FM0.3.1/0.3.2")
        if self.role not in {TRAIN_ROLE, HOLDOUT_ROLE}:
            raise ValueError("composite dataset role must be train or temporal_holdout")
        if (
            self.window_length != 128
            or self.mask_target_fraction != 0.15
            or self.mask_span_range != (1, 4)
        ):
            raise ValueError("FM0.3 composite geometry must be L128/.15/spans1--4")
        if self.windows_per_epoch <= 0 or self.shard_cache_size <= 0:
            raise ValueError("dataset sizes must be positive")
        for value in (
            self.receipt_sha256,
            self.source_bindings_sha256,
            self.role_index_sha256,
        ):
            if _SHA256.fullmatch(value) is None:
                raise ValueError("composite dataset hashes must be lowercase SHA-256")


class FM0CompositeReleaseDataset(FM0ReleaseDataset):
    """Resolve selected observations at their original immutable roots."""

    def __init__(self, config: FM0CompositeDatasetConfig) -> None:
        self.config = config
        self.epoch = 0
        self.release_root = Path(config.composite_root).resolve(strict=True)
        validated = validate_composite_release(
            self.release_root,
            expected_receipt_sha256=config.receipt_sha256,
            expected_source_bindings_sha256=config.source_bindings_sha256,
            expected_role_index_sha256=config.role_index_sha256,
            require_read_only=config.require_read_only,
            require_source_read_only=config.require_read_only,
        )
        source_rows = _read_csv(
            self.release_root / "source_bindings.csv",
            SOURCE_BINDING_FIELDS,
            "source bindings",
        )
        roots = {
            row["source_id"]: Path(row["release_root"]).resolve(strict=True)
            for row in source_rows
        }
        required = np.asarray(variant_view_indices(config.variant), dtype=int)
        grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
        excluded_missing_views = 0
        for item in validated.observations_by_role[config.role]:
            row = {
                "observation_key": item.observation_key,
                "leakage_component_id": item.component,
                "view_present_json": item.view_present_json,
                "sha256": item.shard_sha256,
            }
            if not np.all(_manifest_view_presence(row)[required]):
                excluded_missing_views += 1
                continue
            source_root = roots[item.source_id]
            shard = (source_root / item.relative_path).resolve()
            if (
                source_root not in shard.parents
                or shard.is_symlink()
                or not shard.is_file()
                or (config.require_read_only and shard.stat().st_mode & 0o222)
            ):
                raise ValueError("composite source shard is unavailable or mutable")
            row["_resolved_shard"] = str(shard)
            grouped[item.component].append(row)
        if not grouped:
            raise ValueError(f"no usable composite observations for {config.role}")
        self._source_ids = tuple(sorted(grouped))
        self._visits = {
            source: tuple(sorted(rows, key=lambda row: row["observation_key"]))
            for source, rows in grouped.items()
        }
        self._n_excluded_missing_required_views = excluded_missing_views
        self._cache: OrderedDict[str, Any] = OrderedDict()

    @property
    def contract(self) -> dict[str, Any]:
        return {
            "kind": "fm0_3_composite_release",
            "composite_root": str(self.release_root),
            "receipt_sha256": self.config.receipt_sha256,
            "source_bindings_sha256": self.config.source_bindings_sha256,
            "role_index_sha256": self.config.role_index_sha256,
            "variant": self.config.variant,
            "role": self.config.role,
            "seed": self.config.seed,
            "windows_per_epoch": self.config.windows_per_epoch,
            "window_length": self.config.window_length,
            "mask_target_fraction": self.config.mask_target_fraction,
            "mask_span_range": list(self.config.mask_span_range),
            "n_sources": len(self._source_ids),
            "n_observations": sum(len(rows) for rows in self._visits.values()),
            "n_excluded_missing_required_views": self._n_excluded_missing_required_views,
        }


__all__ = [
    "COMPOSITE_RELEASE_SCHEMA_VERSION",
    "EXCLUDED_OVERLAP_ROLE",
    "HOLDOUT_ROLE",
    "ROLE_INDEX_FIELDS",
    "SOURCE_BINDING_FIELDS",
    "TRAIN_ROLE",
    "FM0CompositeDatasetConfig",
    "FM0CompositeReleaseDataset",
    "SourceManifestBinding",
    "freeze_composite_release",
    "validate_composite_release",
]
