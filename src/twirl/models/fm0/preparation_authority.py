"""Exact authorities for the deferred S66--S77 FM full-visit preparation.

This module is deliberately limited to CPU-side data preparation.  It binds
the pre-model S65 exclusion, validates Phase-A mission/source/HDF5 evidence,
and supplies exact receipt bindings to the separately frozen admission-v2
contract.  It cannot freeze a temporal panel or authorize model training.
"""

from __future__ import annotations

import hashlib
import json
import re
import stat
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .later_sector_admission_v2 import (
    PREPARATION_SECTORS,
    ReceiptBinding,
    SectorPreparationBindings,
    _verify_hdf5_quality,
    _verify_mission_quality,
    _verify_source_inventory,
    admit_preparation_pool,
)
from .registry import FM0ContractError, publish_immutable, sha256_file

AUTHORITY_MAP_SCHEMA_VERSION = (
    "twirl_fm0_s66_s77_full_visit_preparation_authority_v1"
)
AUTHORITY_MAP_CAMPAIGN_ID = "twirl_fm0_s66_s77_full_visit_preparation_v1"
PHASE_A_RECORD_SCHEMA_VERSION = "twirl_fm0_s66_s77_phase_a_authority_record_v1"
PHASE_A_READY_STATE = "FM_PHASE_A_AUTHORITIES_READY"

_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_GIT_SHA = re.compile(r"^[0-9a-f]{40}$")
_PROJECT_ROOT = Path("/orcd/data/mki_aryeh/001/twirl")


@dataclass(frozen=True)
class PreparationAuthorityMap:
    """Validated immutable preparation authority map."""

    path: Path
    sha256: str
    payload: Mapping[str, Any]
    hdf5_receipt_sha256: Mapping[int, str]


def _mapping(value: Any, *, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise FM0ContractError(f"{label} must be a mapping")
    return value


def _sequence(value: Any, *, label: str) -> Sequence[Any]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise FM0ContractError(f"{label} must be a sequence")
    return value


def _digest(value: Any, *, label: str) -> str:
    digest = str(value).strip().lower()
    if _SHA256.fullmatch(digest) is None:
        raise FM0ContractError(f"{label} must be a lowercase SHA-256 digest")
    return digest


def _absolute_project_path(value: Any, *, label: str) -> Path:
    text = str(value).strip()
    if not text:
        raise FM0ContractError(f"{label} path is missing")
    path = Path(text)
    if (
        not path.is_absolute()
        or ".." in path.parts
        or not path.is_relative_to(_PROJECT_ROOT)
    ):
        raise FM0ContractError(f"{label} must stay under {_PROJECT_ROOT}")
    return path


def _repository_path(value: Any, *, label: str) -> Path:
    path = Path(str(value).strip())
    if not str(path) or path.is_absolute() or ".." in path.parts:
        raise FM0ContractError(f"{label} must be a safe repository-relative path")
    return path


def _load_json(path: Path, *, label: str) -> Mapping[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FM0ContractError(f"invalid {label}: {path}") from exc
    return _mapping(payload, label=label)


def _validate_exact_sector_map(value: Any, *, label: str) -> Mapping[str, Any]:
    mapping = _mapping(value, label=label)
    if tuple(sorted(int(key) for key in mapping)) != PREPARATION_SECTORS:
        raise FM0ContractError(f"{label} must cover exactly S66--S77")
    return mapping


def load_preparation_authority_map(
    path: str | Path,
    *,
    expected_sha256: str,
    s66_hdf5_receipt_sha256: str | None = None,
    require_s66_hdf5: bool = True,
) -> PreparationAuthorityMap:
    """Load and strictly validate the exact S66--S77 authority map."""

    source = Path(path).expanduser().resolve()
    expected = _digest(expected_sha256, label="authority-map SHA-256")
    if source.is_symlink() or not source.is_file() or sha256_file(source) != expected:
        raise FM0ContractError("preparation authority-map path or hash drifted")
    payload = _load_json(source, label="preparation authority map")
    if (
        payload.get("schema_version") != AUTHORITY_MAP_SCHEMA_VERSION
        or payload.get("campaign_id") != AUTHORITY_MAP_CAMPAIGN_ID
    ):
        raise FM0ContractError("preparation authority-map identity drifted")

    scope = _mapping(payload.get("scope"), label="authority-map scope")
    sectors = tuple(int(value) for value in scope.get("ordered_sectors", ()))
    if sectors != PREPARATION_SECTORS:
        raise FM0ContractError("authority map must select chronological S66--S77")
    if (
        tuple(int(value) for value in scope.get("excluded_sectors", ())) != (65,)
        or tuple(int(value) for value in scope.get("blocked_sectors", ())) != (78,)
        or tuple(int(value) for value in scope.get("untouched_sectors", ()))
        != (79, 80)
        or scope.get("full_visit_six_view_shards") is not True
        or scope.get("label_blind_sector_wide_hdf5_quality_qa_allowed") is not True
        or scope.get(
            "sealed_aperture_photometry_or_derived_shard_access_forbidden"
        )
        is not True
        or scope.get("temporal_panel_freeze_authorized") is not False
        or scope.get("model_training_authorized") is not False
        or scope.get("h200_use_authorized") is not False
    ):
        raise FM0ContractError("authority-map safety boundary drifted")

    orcd = _mapping(payload.get("orcd"), label="authority-map ORCD settings")
    if (
        orcd.get("host") != "tehan@orcd-login.mit.edu"
        or orcd.get("control_socket")
        != "/Users/tehan/.ssh/cm/tehan@orcd-login.mit.edu:22"
        or orcd.get("partition") != "pg_mki_aryeh"
        or orcd.get("excluded_node") != "node4900"
    ):
        raise FM0ContractError("authority-map ORCD endpoint/resource drifted")
    for field in (
        "project_root",
        "source_repository",
        "detached_repository_prefix",
        "python",
    ):
        _absolute_project_path(orcd.get(field), label=f"ORCD {field}")

    authorities = _mapping(payload.get("authorities"), label="authorities")
    for name in (
        "mission_quality_transfer",
        "corpus_selection",
        "target_observations",
        "gaia_tic_aliases",
    ):
        authority = _mapping(authorities.get(name), label=name)
        for field in ("root", "manifest_path") if name == "mission_quality_transfer" else ("path",):
            _absolute_project_path(authority.get(field), label=f"{name} {field}")
        hash_fields = (
            ("manifest_sha256",)
            if name == "mission_quality_transfer"
            else ("sha256",)
        )
        for field in hash_fields:
            _digest(authority.get(field), label=f"{name} {field}")
    corpus = _mapping(authorities.get("corpus_selection"), label="corpus selection")
    _absolute_project_path(corpus.get("summary_path"), label="corpus summary")
    _digest(corpus.get("summary_sha256"), label="corpus-summary SHA-256")
    for name in ("admission_policy", "exclusion_ledger"):
        authority = _mapping(authorities.get(name), label=name)
        _repository_path(authority.get("repository_path"), label=name)
        _digest(authority.get("sha256"), label=f"{name} SHA-256")

    source_receipts = _validate_exact_sector_map(
        authorities.get("source_receipts"), label="source receipts"
    )
    hdf5_receipts = _validate_exact_sector_map(
        authorities.get("hdf5_quality_receipts"), label="HDF5-quality receipts"
    )
    resolved_hdf5: dict[int, str] = {}
    for sector in PREPARATION_SECTORS:
        source_receipt = _mapping(
            source_receipts[str(sector)], label=f"S{sector} source receipt"
        )
        hdf5_receipt = _mapping(
            hdf5_receipts[str(sector)], label=f"S{sector} HDF5 receipt"
        )
        _absolute_project_path(
            source_receipt.get("path"), label=f"S{sector} source receipt"
        )
        _absolute_project_path(
            hdf5_receipt.get("path"), label=f"S{sector} HDF5 receipt"
        )
        _digest(source_receipt.get("sha256"), label=f"S{sector} source receipt")
        raw_hdf5_sha = hdf5_receipt.get("sha256")
        if sector == 66 and raw_hdf5_sha is None:
            if (
                hdf5_receipt.get("required_runtime_env")
                != "TWIRL_FM0_S66_HDF5_QUALITY_RECEIPT_SHA256"
            ):
                raise FM0ContractError("S66 runtime HDF5 hash requirement drifted")
            if s66_hdf5_receipt_sha256 is None:
                if require_s66_hdf5:
                    raise FM0ContractError(
                        "S66 HDF5-quality receipt SHA-256 is required at runtime"
                    )
                continue
            raw_hdf5_sha = s66_hdf5_receipt_sha256
        resolved_hdf5[sector] = _digest(
            raw_hdf5_sha, label=f"S{sector} HDF5-quality receipt SHA-256"
        )

    outputs = _mapping(payload.get("outputs"), label="authority-map outputs")
    output_paths = {
        field: _absolute_project_path(outputs.get(field), label=f"output {field}")
        for field in (
            "run_root",
            "quality_reference_root",
            "source_inventory_root",
            "phase_a_authority_record",
            "six_view_root",
            "admission_receipt",
            "launch_root",
        )
    }
    run_root = output_paths["run_root"]
    metadata_fields = (
        "quality_reference_root",
        "source_inventory_root",
        "phase_a_authority_record",
        "admission_receipt",
        "launch_root",
    )
    if any(
        not output_paths[field].is_relative_to(run_root) for field in metadata_fields
    ):
        raise FM0ContractError("preparation metadata escapes the scoped run root")
    expected_export_root = (
        _PROJECT_ROOT
        / "exports/fm0/later_sector_releases/"
        "twirl_fm0_s66_s77_full_visit_v1"
    )
    if output_paths["six_view_root"] != expected_export_root:
        raise FM0ContractError("six-view release root drifted from stable FM exports")
    return PreparationAuthorityMap(source, expected, payload, resolved_hdf5)


def _assert_tree_read_only(root: Path) -> None:
    if root.is_symlink() or not root.is_dir():
        raise FM0ContractError(f"immutable authority directory is missing: {root}")
    for path in (root, *root.rglob("*")):
        if path.is_symlink() or path.stat().st_mode & (
            stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH
        ):
            raise FM0ContractError(f"authority output is not materialized/read-only: {path}")


def freeze_phase_a_authorities(
    authority: PreparationAuthorityMap,
    *,
    producer_git_sha: str,
    output_path: str | Path,
) -> tuple[dict[str, Any], str]:
    """Fully validate and freeze all Phase-A plus HDF5 receipt bindings."""

    if _GIT_SHA.fullmatch(str(producer_git_sha)) is None:
        raise FM0ContractError("producer Git SHA must be full lowercase 40-hex")
    payload = authority.payload
    authorities = _mapping(payload["authorities"], label="authorities")
    outputs = _mapping(payload["outputs"], label="outputs")
    quality_root = Path(str(outputs["quality_reference_root"]))
    source_root = Path(str(outputs["source_inventory_root"]))
    hdf5_map = _mapping(authorities["hdf5_quality_receipts"], label="HDF5 map")
    rows: list[dict[str, Any]] = []
    for sector in PREPARATION_SECTORS:
        quality_dir = quality_root / f"s{sector:04d}"
        source_dir = source_root / f"s{sector:04d}"
        _assert_tree_read_only(quality_dir)
        _assert_tree_read_only(source_dir)
        mission_path = quality_dir / "manifest.json"
        source_path = source_dir / "summary.json"
        hdf5_path = Path(str(_mapping(hdf5_map[str(sector)], label="HDF5 row")["path"]))
        mission_binding = ReceiptBinding(mission_path, sha256_file(mission_path))
        source_binding = ReceiptBinding(source_path, sha256_file(source_path))
        hdf5_binding = ReceiptBinding(
            hdf5_path, authority.hdf5_receipt_sha256[sector]
        )
        _mission_path, reference, mission = _verify_mission_quality(
            sector=sector, binding=mission_binding
        )
        (
            _source_path,
            inventory_summary,
            _selected_rows,
            source,
        ) = _verify_source_inventory(sector=sector, binding=source_binding)
        _hdf5_path, _hdf5_receipt, hdf5 = _verify_hdf5_quality(
            sector=sector,
            binding=hdf5_binding,
            inventory_summary=inventory_summary,
            reference=reference,
        )
        rows.append(
            {
                "sector": sector,
                "mission_quality_reference": mission,
                "source_inventory": source,
                "hdf5_quality": hdf5,
            }
        )
    record = {
        "schema_version": PHASE_A_RECORD_SCHEMA_VERSION,
        "authority_state": PHASE_A_READY_STATE,
        "campaign_id": AUTHORITY_MAP_CAMPAIGN_ID,
        "authority_map_path": str(authority.path),
        "authority_map_sha256": authority.sha256,
        "producer_git_sha": producer_git_sha,
        "ordered_sectors": list(PREPARATION_SECTORS),
        "sectors": rows,
        "phase_a_fully_validated": True,
        "hdf5_quality_fully_validated": True,
        "label_blind_sector_wide_hdf5_quality_qa_performed": True,
        "sealed_aperture_photometry_or_derived_shard_accessed": False,
        "six_view_shards_verified": False,
        "temporal_panel_frozen": False,
        "model_training_authorized": False,
        "h200_used": False,
    }
    encoded = (
        json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    destination = Path(output_path).expanduser().resolve()
    expected_output = Path(str(outputs["phase_a_authority_record"])).resolve()
    if destination != expected_output:
        raise FM0ContractError("Phase-A authority-record output path drifted")
    publish_immutable(destination, encoded)
    destination.chmod(0o440)
    return record, hashlib.sha256(encoded).hexdigest()


def load_phase_a_authority_record(
    authority: PreparationAuthorityMap,
    *,
    expected_producer_git_sha: str,
    expected_sha256: str,
) -> tuple[Mapping[str, Any], Path, str]:
    """Load the immutable Phase-A record and recheck exact sector bindings."""

    output = _mapping(authority.payload["outputs"], label="outputs")
    path = Path(str(output["phase_a_authority_record"])).resolve()
    expected = _digest(expected_sha256, label="Phase-A authority record SHA-256")
    if (
        path.is_symlink()
        or not path.is_file()
        or path.stat().st_mode & 0o222
        or sha256_file(path) != expected
    ):
        raise FM0ContractError(
            "Phase-A authority record is missing, writable, or hash-drifted"
        )
    record = _load_json(path, label="Phase-A authority record")
    if (
        record.get("schema_version") != PHASE_A_RECORD_SCHEMA_VERSION
        or record.get("authority_state") != PHASE_A_READY_STATE
        or record.get("campaign_id") != AUTHORITY_MAP_CAMPAIGN_ID
        or record.get("authority_map_sha256") != authority.sha256
        or record.get("producer_git_sha") != expected_producer_git_sha
        or tuple(int(value) for value in record.get("ordered_sectors", ()))
        != PREPARATION_SECTORS
        or record.get("phase_a_fully_validated") is not True
        or record.get("hdf5_quality_fully_validated") is not True
        or record.get("label_blind_sector_wide_hdf5_quality_qa_performed")
        is not True
        or record.get(
            "sealed_aperture_photometry_or_derived_shard_accessed"
        )
        is not False
        or record.get("six_view_shards_verified") is not False
        or record.get("model_training_authorized") is not False
        or record.get("h200_used") is not False
    ):
        raise FM0ContractError("Phase-A authority record contract drifted")
    rows = _sequence(record.get("sectors"), label="Phase-A sector rows")
    if tuple(int(_mapping(row, label="Phase-A row")["sector"]) for row in rows) != (
        PREPARATION_SECTORS
    ):
        raise FM0ContractError("Phase-A authority record sector order drifted")
    authorities = _mapping(authority.payload["authorities"], label="authorities")
    outputs = _mapping(authority.payload["outputs"], label="outputs")
    hdf5_map = _mapping(
        authorities["hdf5_quality_receipts"], label="HDF5-quality map"
    )
    for raw_row in rows:
        row = _mapping(raw_row, label="Phase-A row")
        sector = int(row["sector"])
        mission = _mapping(
            row.get("mission_quality_reference"), label="Phase-A mission row"
        )
        source = _mapping(
            row.get("source_inventory"), label="Phase-A source row"
        )
        hdf5 = _mapping(row.get("hdf5_quality"), label="Phase-A HDF5 row")
        expected_bindings = (
            (
                mission,
                Path(str(outputs["quality_reference_root"]))
                / f"s{sector:04d}"
                / "manifest.json",
                str(mission.get("sha256", "")),
            ),
            (
                source,
                Path(str(outputs["source_inventory_root"]))
                / f"s{sector:04d}"
                / "summary.json",
                str(source.get("sha256", "")),
            ),
            (
                hdf5,
                Path(
                    str(
                        _mapping(
                            hdf5_map[str(sector)], label="HDF5-quality map row"
                        )["path"]
                    )
                ),
                authority.hdf5_receipt_sha256[sector],
            ),
        )
        for evidence, expected_path, expected_sha in expected_bindings:
            bound_path = Path(str(evidence.get("path", ""))).expanduser().resolve()
            digest = _digest(
                evidence.get("sha256"), label=f"S{sector} Phase-A evidence SHA-256"
            )
            if (
                bound_path != expected_path.resolve()
                or digest != expected_sha
                or bound_path.is_symlink()
                or not bound_path.is_file()
                or sha256_file(bound_path) != digest
            ):
                raise FM0ContractError(f"S{sector} Phase-A evidence binding drifted")
    return record, path, expected


def admit_from_preparation_authorities(
    authority: PreparationAuthorityMap,
    *,
    phase_a_producer_git_sha: str,
    phase_a_record_sha256: str,
    producer_git_sha: str,
    output_path: str | Path,
) -> tuple[dict[str, Any], str]:
    """Fully validate six-view releases and run the frozen admission-v2 gate."""

    record, _record_path, _record_sha = load_phase_a_authority_record(
        authority,
        expected_producer_git_sha=phase_a_producer_git_sha,
        expected_sha256=phase_a_record_sha256,
    )
    authorities = _mapping(authority.payload["authorities"], label="authorities")
    outputs = _mapping(authority.payload["outputs"], label="outputs")
    phase_rows = {
        int(_mapping(row, label="Phase-A row")["sector"]): _mapping(
            row, label="Phase-A row"
        )
        for row in _sequence(record["sectors"], label="Phase-A rows")
    }
    bindings: list[SectorPreparationBindings] = []
    six_root = Path(str(outputs["six_view_root"]))
    for sector in PREPARATION_SECTORS:
        row = phase_rows[sector]
        mission = _mapping(row["mission_quality_reference"], label="mission row")
        hdf5 = _mapping(row["hdf5_quality"], label="HDF5 row")
        source = _mapping(row["source_inventory"], label="source row")
        six_path = six_root / f"s{sector:04d}" / "receipt.json"
        if six_path.is_symlink() or not six_path.is_file():
            raise FM0ContractError(
                f"S{sector} six-view receipt is not a materialized exact path"
            )
        receipt = _load_json(six_path, label=f"S{sector} six-view receipt")
        if int(receipt.get("sector", -1)) != sector:
            raise FM0ContractError(f"S{sector} six-view receipt identity drifted")
        six_sha = sha256_file(six_path)
        bindings.append(
            SectorPreparationBindings(
                sector=sector,
                mission_quality_reference=ReceiptBinding(
                    Path(str(mission["path"])), str(mission["sha256"])
                ),
                hdf5_quality=ReceiptBinding(
                    Path(str(hdf5["path"])), str(hdf5["sha256"])
                ),
                source_inventory=ReceiptBinding(
                    Path(str(source["path"])), str(source["sha256"])
                ),
                six_view_release=ReceiptBinding(six_path, six_sha),
            )
        )
    policy = _mapping(authorities["admission_policy"], label="admission policy")
    ledger = _mapping(authorities["exclusion_ledger"], label="exclusion ledger")
    # The actual detached repository is supplied through the map path, which is
    # itself inside the exact checked-out repository used by the Slurm wrapper.
    repo = authority.path.parents[2]
    destination = Path(output_path).expanduser().resolve()
    if destination != Path(str(outputs["admission_receipt"])).resolve():
        raise FM0ContractError("admission-v2 output path drifted")
    return admit_preparation_pool(
        config_path=repo / _repository_path(
            policy["repository_path"], label="admission policy"
        ),
        expected_config_sha256=str(policy["sha256"]),
        exclusion_ledger_path=repo / _repository_path(
            ledger["repository_path"], label="exclusion ledger"
        ),
        expected_exclusion_ledger_sha256=str(ledger["sha256"]),
        ordered_sector_bindings=bindings,
        producer_git_sha=producer_git_sha,
        output_path=destination,
    )


__all__ = [
    "AUTHORITY_MAP_CAMPAIGN_ID",
    "AUTHORITY_MAP_SCHEMA_VERSION",
    "PHASE_A_READY_STATE",
    "PHASE_A_RECORD_SCHEMA_VERSION",
    "PreparationAuthorityMap",
    "admit_from_preparation_authorities",
    "freeze_phase_a_authorities",
    "load_phase_a_authority_record",
    "load_preparation_authority_map",
]
