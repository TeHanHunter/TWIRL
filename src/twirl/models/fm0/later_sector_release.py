"""Immutable full-visit six-view preparation for later TWIRL-FM sectors.

The release produced here is a deferred ORCD development checkpoint.  It is
not accepted A2v1 and is not, by itself, authorized for model training or a
temporal panel.  The builder consumes the exact checksum-bound corpus
selection through the later source-inventory bundle, the provider-neutral
mission-quality reference, and the all-HDF5 cadence/openability receipt.

Sealed-test source rows remain in inventory counts, but are filtered before
any HDF5 path is resolved, hashed, or opened.  Shards retain every cadence in
the two-orbit source visit; model context length is intentionally a later
loader-level decision.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import shutil
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path, PurePosixPath
from typing import Any

from twirl.lightcurves.mission_quality_reference import (
    MISSION_QUALITY_REFERENCE_SCHEMA_VERSION,
    MissionQualityReference,
    load_mission_quality_reference,
)

from .corpus import (
    CORPUS_SELECTION_FIELDS,
    CORPUS_SELECTION_SCHEMA_VERSION,
)
from .hdf5_quality_admission import (
    HDF5_QUALITY_READY_STATE,
    HDF5_QUALITY_RECEIPT_SCHEMA_VERSION,
)
from .input_release import (
    ERROR_VIEW_NAMES,
    FLUX_VIEW_NAMES,
    build_mission_quality_observation_release,
    deterministic_npz_bytes,
    load_input_release,
)
from .later_hdf5_adapter import (
    LATER_ALLOWED_SOURCE_PARTITIONS,
    LATER_FORBIDDEN_SOURCE_PARTITION,
    LATER_HDF5_ADAPTER_NAME,
    build_later_adapter_cache,
    load_later_hdf5_observation,
)
from .later_source_inventory import (
    IDENTITY_SOURCE_CORPUS_SELECTION,
    SOURCE_FIELDS,
    SOURCE_ROW_SCHEMA_VERSION,
)
from .later_source_inventory import (
    SUMMARY_SCHEMA_VERSION as SOURCE_INVENTORY_SUMMARY_SCHEMA_VERSION,
)
from .registry import (
    FM0ContractError,
    assert_no_search_columns,
    deterministic_source_partition,
    publish_immutable,
    read_rows,
    sha256_file,
)

LATER_SIX_VIEW_RECEIPT_SCHEMA_VERSION = (
    "twirl_fm0_later_six_view_sector_release_v1"
)
LATER_SIX_VIEW_READY_STATE = "FM_SIX_VIEW_DEFERRED_READY"
LATER_SIX_VIEW_MANIFEST_SCHEMA_VERSION = "twirl_fm0_later_six_view_manifest_v1"
LATER_SIX_VIEW_SOURCE_MANIFEST_SCHEMA_VERSION = (
    "twirl_fm0_later_six_view_source_manifest_v1"
)
LATER_SIX_VIEW_TIMING_SCHEMA_VERSION = "twirl_fm0_later_six_view_timing_v1"
LATER_SIX_VIEW_PRODUCT_STATE = "ORCD_COMPLETE_DEFERRED"
LATER_SIX_VIEW_ALLOWED_PARTITIONS = tuple(sorted(LATER_ALLOWED_SOURCE_PARTITIONS))

LATER_SIX_VIEW_RECEIPT_FIELDS = frozenset(
    {
        "schema_version",
        "release_state",
        "passed",
        "created_utc",
        "producer_git_sha",
        "sector",
        "expected_orbits",
        "product_state",
        "n_observations",
        "n_shards",
        "n_cadences",
        "n_segments",
        "source_rows_by_partition",
        "n_sealed_source_rows_excluded",
        "flux_view_names",
        "error_view_names",
        "source_inventory_summary_path",
        "source_inventory_summary_sha256",
        "source_inventory_sources_path",
        "source_inventory_sources_sha256",
        "corpus_selection_path",
        "corpus_selection_sha256",
        "mission_quality_reference_schema_version",
        "mission_quality_provider",
        "mission_quality_composition",
        "mission_quality_reference_manifest_path",
        "mission_quality_reference_manifest_sha256",
        "mission_quality_reference_table_path",
        "mission_quality_reference_table_sha256",
        "hdf5_quality_receipt_path",
        "hdf5_quality_receipt_sha256",
        "manifest_sha256",
        "source_manifest_sha256",
        "visit_timing_sha256",
        "input_adapter",
        "full_visit_shards",
        "model_context_length_bound",
        "n_internal_bad_cadences",
        "n_external_bad_cadences",
        "n_authority_excluded_cadences",
        "sealed_hdf5_content_opened",
        "sealed_shards_written",
        "model_outcome_blind",
        "six_view_shards_verified",
        "scientific_training_eligible",
        "panel_admission_authorized",
        "claim_limit",
    }
)

LATER_SIX_VIEW_MANIFEST_FIELDS = (
    "manifest_schema_version",
    "observation_key",
    "product_instance_id",
    "physical_source_id",
    "gaia_dr3_source_id",
    "tic_id",
    "sector",
    "leakage_component_id",
    "source_partition",
    "product_state",
    "relative_path",
    "sha256",
    "input_source_sha256",
    "source_row_sha256",
    "n_cadences",
    "n_segments",
    "view_present_json",
    "input_adapter",
    "mission_quality_provider",
    "mission_quality_reference_manifest_sha256",
    "hdf5_quality_receipt_sha256",
    "full_visit_shard",
    "model_context_length_bound",
    "scientific_training_eligible",
    "panel_admission_authorized",
)
LATER_SIX_VIEW_TIMING_FIELDS = (
    "timing_schema_version",
    "observation_key",
    "physical_source_id",
    "absolute_visit_start",
    "absolute_visit_end",
)

_SHA40 = re.compile(r"^[0-9a-f]{40}$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_LEAKAGE_COMPONENT = re.compile(r"^leakage_[0-9a-f]{64}$")
_ALL_PARTITIONS = frozenset(
    {*LATER_ALLOWED_SOURCE_PARTITIONS, LATER_FORBIDDEN_SOURCE_PARTITION}
)
_EVIDENCE_FIELDS = frozenset(
    {
        "orbit",
        "camera",
        "ccd",
        "cell",
        "hdf5_sha256",
        "cell_manifest_sha256",
        "output_manifest_sha256",
        "hdf5_path",
        "output_manifest_path",
    }
)


def _digest(value: Any, *, label: str, length: int = 64) -> str:
    result = str(value).strip().lower()
    pattern = _SHA40 if length == 40 else _SHA256
    if pattern.fullmatch(result) is None:
        raise FM0ContractError(
            f"{label} must be a lowercase {length}-character hexadecimal digest"
        )
    return result


def _read_json(path: Path, *, label: str) -> Mapping[str, Any]:
    if path.is_symlink() or not path.is_file() or path.stat().st_size <= 0:
        raise FM0ContractError(f"{label} must be a nonempty materialized file: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FM0ContractError(f"invalid {label}: {path}") from exc
    if not isinstance(value, Mapping):
        raise FM0ContractError(f"{label} must be a JSON object: {path}")
    return value


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()


def _stable_id(prefix: str, payload: Any) -> str:
    return f"{prefix}_{_canonical_sha256(payload)}"


def _csv_bytes(rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> bytes:
    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(fields), lineterminator="\n")
    writer.writeheader()
    for row in rows:
        if tuple(row) != tuple(fields):
            raise FM0ContractError("later six-view CSV columns drifted")
        writer.writerow({field: row[field] for field in fields})
    return stream.getvalue().encode("utf-8")


def _json_list(value: Any, *, label: str) -> list[Any]:
    try:
        result = json.loads(str(value))
    except (TypeError, json.JSONDecodeError) as exc:
        raise FM0ContractError(f"{label} must be a JSON list") from exc
    if not isinstance(result, list) or not result:
        raise FM0ContractError(f"{label} must be a nonempty JSON list")
    return result


def _bool_text(value: Any, *, label: str) -> bool:
    normalized = str(value).strip().lower()
    if normalized in {"true", "1"}:
        return True
    if normalized in {"false", "0"}:
        return False
    raise FM0ContractError(f"{label} must be a serialized boolean")


def _source_row_portable_evidence(row: Mapping[str, Any]) -> list[dict[str, Any]]:
    values = _json_list(row["hdf5_evidence_json"], label="hdf5_evidence_json")
    portable: list[dict[str, Any]] = []
    for value in values:
        if not isinstance(value, Mapping) or set(value) != _EVIDENCE_FIELDS:
            raise FM0ContractError("later source HDF5 evidence fields drifted")
        portable.append(
            {
                "orbit": int(value["orbit"]),
                "camera": int(value["camera"]),
                "ccd": int(value["ccd"]),
                "cell": str(value["cell"]),
                "hdf5_sha256": _digest(
                    value["hdf5_sha256"], label="source evidence HDF5 SHA-256"
                ),
                "cell_manifest_sha256": _digest(
                    value["cell_manifest_sha256"],
                    label="source evidence cell-manifest SHA-256",
                ),
                "output_manifest_sha256": _digest(
                    value["output_manifest_sha256"],
                    label="source evidence output-manifest SHA-256",
                ),
            }
        )
    return portable


def _validate_source_row(
    row: Mapping[str, Any],
    *,
    sector: int,
    target_authority_sha256: str,
    identity_binding_sha256: str,
    corpus_summary_sha256: str,
    gaia_tic_alias_authority_sha256: str,
    source_receipt_path: str,
    source_receipt_sha256: str,
) -> None:
    if tuple(row) != tuple(SOURCE_FIELDS):
        raise FM0ContractError("later source-inventory columns drifted")
    assert_no_search_columns(row.keys(), context="later source inventory")
    if (
        row["source_row_schema_version"] != SOURCE_ROW_SCHEMA_VERSION
        or int(row["sector"]) != sector
        or row["identity_source"] != IDENTITY_SOURCE_CORPUS_SELECTION
        or row["target_authority_sha256"] != target_authority_sha256
        or row["corpus_selection_sha256"] != identity_binding_sha256
        or row["corpus_summary_sha256"] != corpus_summary_sha256
        or row["gaia_tic_alias_authority_sha256"]
        != gaia_tic_alias_authority_sha256
        or row["product_state"] != LATER_SIX_VIEW_PRODUCT_STATE
        or _bool_text(row["identity_ambiguous"], label="identity_ambiguous")
        is not False
        or int(row["n_hdf5_products"]) != 2
        or row["source_receipt_path"] != source_receipt_path
        or row["source_receipt_sha256"] != source_receipt_sha256
        or _bool_text(row["model_outcome_blind"], label="model_outcome_blind")
        is not True
        or _bool_text(
            row["panel_admission_authorized"],
            label="panel_admission_authorized",
        )
        is not False
    ):
        raise FM0ContractError("later source-inventory row contract drifted")
    gaia = str(row["gaia_dr3_source_id"])
    tic = str(row["tic_id"])
    if not gaia.isdigit() or int(gaia) <= 0 or not tic.isdigit() or int(tic) <= 0:
        raise FM0ContractError("later source identifiers are invalid")
    component = str(row["leakage_component_id"])
    partition = str(row["source_partition"])
    if _LEAKAGE_COMPONENT.fullmatch(component) is None or partition not in _ALL_PARTITIONS:
        raise FM0ContractError("later source leakage component or partition is invalid")
    expected_partition, _bucket = deterministic_source_partition(component)
    if partition != expected_partition:
        raise FM0ContractError("later source partition differs from frozen split")
    orbits = tuple(int(value) for value in _json_list(row["orbits_json"], label="orbits_json"))
    if orbits != (2 * sector + 7, 2 * sector + 8):
        raise FM0ContractError("later source orbit identity drifted")
    paths = _json_list(row["hdf5_paths_json"], label="hdf5_paths_json")
    hashes = _json_list(row["hdf5_sha256_json"], label="hdf5_sha256_json")
    portable = _source_row_portable_evidence(row)
    if len(paths) != 2 or len(hashes) != 2 or len(portable) != 2:
        raise FM0ContractError("later source does not bind two HDF5 products")
    normalized_hashes = [
        _digest(value, label="source HDF5 SHA-256") for value in hashes
    ]
    if normalized_hashes != [value["hdf5_sha256"] for value in portable]:
        raise FM0ContractError("later source HDF5/evidence hashes disagree")
    retained_source_sha = _canonical_sha256(
        {
            "schema_version": "twirl_fm0_later_retained_source_inventory_v1",
            "sector": sector,
            "tic_id": str(int(tic)),
            "target_authority_sha256": target_authority_sha256,
            "evidence": portable,
        }
    )
    if row["retained_source_sha256"] != retained_source_sha:
        raise FM0ContractError("later retained_source_sha256 drifted")
    source_row_sha = _canonical_sha256(
        {
            "retained_source_sha256": retained_source_sha,
            "gaia_dr3_source_id": str(int(gaia)),
            "tic_id": str(int(tic)),
            "leakage_component_id": component,
            "source_partition": partition,
            "identity_source": IDENTITY_SOURCE_CORPUS_SELECTION,
            "identity_binding_sha256": identity_binding_sha256,
            "corpus_summary_sha256": corpus_summary_sha256,
            "gaia_tic_alias_authority_sha256": gaia_tic_alias_authority_sha256,
        }
    )
    if row["source_row_sha256"] != source_row_sha:
        raise FM0ContractError("later source_row_sha256 drifted")


def _load_inventory_bundle(
    *,
    source_inventory_dir: str | Path,
    sector: int,
    expected_summary_sha256: str,
) -> tuple[Mapping[str, Any], Path, str, list[dict[str, Any]], Path, str]:
    root = Path(source_inventory_dir).expanduser().resolve()
    summary_path = root / "summary.json"
    sources_path = root / "sources.csv"
    ready_path = root / "READY"
    expected_summary_hash = _digest(
        expected_summary_sha256,
        label="expected source-inventory summary SHA-256",
    )
    if summary_path.is_symlink() or not summary_path.is_file():
        raise FM0ContractError("later source-inventory summary is not materialized")
    summary_hash = sha256_file(summary_path)
    if summary_hash != expected_summary_hash:
        raise FM0ContractError("later source-inventory summary caller binding drifted")
    summary = _read_json(summary_path, label="later source-inventory summary")
    if (
        not sources_path.is_file()
        or sources_path.is_symlink()
        or not ready_path.is_file()
        or ready_path.is_symlink()
        or ready_path.read_text(encoding="utf-8").strip() != summary_hash
    ):
        raise FM0ContractError("later source-inventory bundle is incomplete")
    if (
        summary.get("summary_schema_version")
        != SOURCE_INVENTORY_SUMMARY_SCHEMA_VERSION
        or int(summary.get("sector", -1)) != sector
        or summary.get("selection_identity_source")
        != IDENTITY_SOURCE_CORPUS_SELECTION
        or summary.get("hdf5_content_opened") is not False
        or summary.get("sealed_hdf5_content_opened") is not False
        or summary.get("model_outcome_blind") is not True
        or summary.get("six_view_shards_verified") is not False
        or summary.get("panel_admission_authorized") is not False
    ):
        raise FM0ContractError("later source-inventory summary contract drifted")
    sources_hash = sha256_file(sources_path)
    if sources_hash != _digest(
        summary.get("sources_sha256"), label="source-inventory sources SHA-256"
    ):
        raise FM0ContractError("later source-inventory sources hash drifted")
    rows = read_rows(sources_path)
    if not rows or any(tuple(row) != tuple(SOURCE_FIELDS) for row in rows):
        raise FM0ContractError("later source-inventory sources columns drifted")
    if int(summary.get("n_source_rows", -1)) != len(rows):
        raise FM0ContractError("later source-inventory source count drifted")
    partitions = {
        partition: sum(row["source_partition"] == partition for row in rows)
        for partition in sorted(_ALL_PARTITIONS)
    }
    if summary.get("source_rows_by_partition") != partitions:
        raise FM0ContractError("later source-inventory partition counts drifted")

    target_authority_sha = _digest(
        summary.get("target_authority_sha256"),
        label="source-inventory target-authority SHA-256",
    )
    identity_binding_sha = _digest(
        summary.get("identity_binding_sha256"),
        label="source-inventory identity-binding SHA-256",
    )
    identity_binding_path = Path(
        str(summary.get("identity_binding_path", ""))
    ).expanduser()
    if not identity_binding_path.is_absolute():
        raise FM0ContractError("corpus-selection identity binding must be absolute")
    identity_binding_path = identity_binding_path.resolve()
    if (
        identity_binding_path.is_symlink()
        or not identity_binding_path.is_file()
        or sha256_file(identity_binding_path) != identity_binding_sha
    ):
        raise FM0ContractError("corpus-selection identity binding drifted")

    corpus_summary_sha = _digest(
        summary.get("corpus_summary_sha256"),
        label="source-inventory corpus-summary SHA-256",
    )
    corpus_summary_path = Path(
        str(summary.get("corpus_summary_path", ""))
    ).expanduser()
    if not corpus_summary_path.is_absolute():
        raise FM0ContractError("corpus-summary binding must be absolute")
    corpus_summary_path = corpus_summary_path.resolve()
    corpus_summary = _read_json(corpus_summary_path, label="corpus-selection summary")
    if sha256_file(corpus_summary_path) != corpus_summary_sha:
        raise FM0ContractError("corpus-selection summary binding drifted")
    try:
        corpus_sectors = tuple(int(value) for value in corpus_summary["sectors"])
        corpus_sector_count = int(corpus_summary["observations_by_sector"][str(sector)])
        corpus_selection_sha = str(corpus_summary["outputs_sha256"]["selection.csv"])
        corpus_authorities = corpus_summary["input_authorities"]
        corpus_observation_sha = str(corpus_authorities["observation_fits"]["sha256"])
        corpus_alias_sha = str(corpus_authorities["gaia_tic_alias_table"]["sha256"])
    except (KeyError, TypeError, ValueError) as exc:
        raise FM0ContractError("corpus-selection summary fields drifted") from exc
    if (
        corpus_summary.get("schema_version") != CORPUS_SELECTION_SCHEMA_VERSION
        or corpus_sectors != tuple(range(66, 78))
        or corpus_sector_count != len(rows)
        or corpus_selection_sha != identity_binding_sha
        or corpus_observation_sha != target_authority_sha
    ):
        raise FM0ContractError("corpus-selection summary authority drifted")

    alias_authority_sha = _digest(
        summary.get("gaia_tic_alias_authority_sha256"),
        label="source-inventory Gaia--TIC alias authority SHA-256",
    )
    alias_authority_path = Path(
        str(summary.get("gaia_tic_alias_authority_path", ""))
    ).expanduser()
    if not alias_authority_path.is_absolute():
        raise FM0ContractError("Gaia--TIC alias authority path must be absolute")
    alias_authority_path = alias_authority_path.resolve()
    if (
        summary.get("gaia_tic_alias_authority_file_verified") is not True
        or corpus_alias_sha != alias_authority_sha
        or alias_authority_path.is_symlink()
        or not alias_authority_path.is_file()
        or sha256_file(alias_authority_path) != alias_authority_sha
    ):
        raise FM0ContractError("Gaia--TIC alias authority binding drifted")
    selection_rows = read_rows(identity_binding_path)
    if not selection_rows or any(
        tuple(row) != tuple(CORPUS_SELECTION_FIELDS) for row in selection_rows
    ):
        raise FM0ContractError("corpus-selection columns drifted")
    selected = {
        (str(int(row["gaia_dr3_source_id"])), str(int(row["tic_id"]))): row
        for row in selection_rows
        if int(row["sector"]) == sector
        and row["schema_version"] == CORPUS_SELECTION_SCHEMA_VERSION
    }
    source_keys = {
        (str(int(row["gaia_dr3_source_id"])), str(int(row["tic_id"])))
        for row in rows
    }
    if not selected or set(selected) != source_keys or len(source_keys) != len(rows):
        raise FM0ContractError("source inventory differs from exact corpus selection")

    source_receipt_path = str(summary.get("source_receipt_path", ""))
    receipt_path = Path(source_receipt_path).expanduser()
    if not receipt_path.is_absolute():
        raise FM0ContractError("source receipt path must be absolute")
    receipt_path = receipt_path.resolve()
    source_receipt_sha = _digest(
        summary.get("source_receipt_sha256"), label="source receipt SHA-256"
    )
    if (
        receipt_path.is_symlink()
        or not receipt_path.is_file()
        or sha256_file(receipt_path) != source_receipt_sha
    ):
        raise FM0ContractError("source receipt binding drifted")

    for row in rows:
        key = (str(int(row["gaia_dr3_source_id"])), str(int(row["tic_id"])))
        selection = selected[key]
        camera = int(selection["camera"])
        ccd = int(selection["ccd"])
        expected_selection_paths = [
            f"orbit-{orbit}/ffi/cam{camera}/ccd{ccd}/LC/{key[1]}.h5"
            for orbit in (2 * sector + 7, 2 * sector + 8)
        ]
        portable_evidence = _source_row_portable_evidence(row)
        if (
            selection["leakage_component_id"] != row["leakage_component_id"]
            or selection["source_partition"] != row["source_partition"]
            or tuple(int(value) for value in json.loads(selection["orbits_json"]))
            != tuple(int(value) for value in json.loads(row["orbits_json"]))
            or json.loads(selection["hdf5_paths_json"])
            != expected_selection_paths
            or any(
                int(value["camera"]) != camera or int(value["ccd"]) != ccd
                for value in portable_evidence
            )
        ):
            raise FM0ContractError("source row differs from corpus-selection identity")
        _validate_source_row(
            row,
            sector=sector,
            target_authority_sha256=target_authority_sha,
            identity_binding_sha256=identity_binding_sha,
            corpus_summary_sha256=corpus_summary_sha,
            gaia_tic_alias_authority_sha256=alias_authority_sha,
            source_receipt_path=str(receipt_path),
            source_receipt_sha256=source_receipt_sha,
        )
    return summary, summary_path, summary_hash, rows, sources_path, sources_hash


def _load_hdf5_quality_receipt(
    *,
    path: str | Path,
    sector: int,
    inventory_summary: Mapping[str, Any],
    reference: MissionQualityReference,
    expected_receipt_sha256: str,
) -> tuple[Mapping[str, Any], Path, str]:
    receipt_path = Path(path).expanduser().resolve()
    expected_receipt_hash = _digest(
        expected_receipt_sha256,
        label="expected HDF5-quality receipt SHA-256",
    )
    if receipt_path.is_symlink() or not receipt_path.is_file():
        raise FM0ContractError("later HDF5-quality receipt is not materialized")
    receipt_hash = sha256_file(receipt_path)
    if receipt_hash != expected_receipt_hash:
        raise FM0ContractError("later HDF5-quality receipt caller binding drifted")
    receipt = _read_json(receipt_path, label="later HDF5-quality receipt")
    expected_orbits = (2 * sector + 7, 2 * sector + 8)
    if (
        receipt.get("schema_version") != HDF5_QUALITY_RECEIPT_SCHEMA_VERSION
        or receipt.get("quality_state") != HDF5_QUALITY_READY_STATE
        or receipt.get("passed") is not True
        or int(receipt.get("sector", -1)) != sector
        or tuple(int(value) for value in receipt.get("expected_orbits", ()))
        != expected_orbits
        or receipt.get("mission_quality_type") != reference.provider
        or receipt.get("hdf5_openability_verified") is not True
        or receipt.get("internal_cadence_quality_verified") is not True
        or receipt.get("external_cadence_quality_verified") is not True
        or receipt.get("six_view_shards_verified") is not False
        or receipt.get("panel_admission_authorized") is not False
        or int(receipt.get("n_unreadable_hdf5", -1)) != 0
        or int(receipt.get("n_hdf5_opened", -1))
        != int(receipt.get("n_hdf5_products", -2))
        or int(receipt.get("n_hdf5_products", -1))
        != int(inventory_summary.get("n_hdf5_products_declared", -2))
        or receipt.get("source_receipt_sha256")
        != inventory_summary.get("source_receipt_sha256")
        or str(Path(str(receipt.get("source_receipt_path", ""))).resolve())
        != str(Path(str(inventory_summary.get("source_receipt_path", ""))).resolve())
        or receipt.get("quality_transfer_manifest_sha256")
        != reference.transfer_manifest_sha256
        or receipt.get("mission_quality_authority_exclusions_sha256")
        != reference.authority_exclusions_sha256
    ):
        raise FM0ContractError("later HDF5-quality receipt contract drifted")
    return receipt, receipt_path, receipt_hash


def _validate_shard(path: Path, row: Mapping[str, Any]) -> None:
    if not path.is_file() or sha256_file(path) != row["sha256"]:
        raise FM0ContractError(f"later six-view shard hash drifted: {path}")
    release = load_input_release(path)
    if (
        release.n_cadences != int(row["n_cadences"])
        or release.n_segments != int(row["n_segments"])
        or release.view_present.astype(int).tolist()
        != json.loads(str(row["view_present_json"]))
    ):
        raise FM0ContractError(f"later six-view shard metadata drifted: {path}")


def _read_csv_exact(
    path: Path, *, fields: Sequence[str], label: str
) -> list[dict[str, Any]]:
    if path.is_symlink() or not path.is_file() or path.stat().st_size <= 0:
        raise FM0ContractError(f"{label} must be a nonempty materialized CSV")
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if tuple(reader.fieldnames or ()) != tuple(fields):
                raise FM0ContractError(f"{label} columns drifted")
            rows = [dict(row) for row in reader]
    except (OSError, UnicodeError, csv.Error) as exc:
        raise FM0ContractError(f"cannot read {label}: {path}") from exc
    if not rows or any(tuple(row) != tuple(fields) for row in rows):
        raise FM0ContractError(f"{label} rows drifted or are empty")
    return rows


def _make_tree_read_only(root: Path) -> None:
    if root.is_symlink() or not root.is_dir():
        raise FM0ContractError(f"release root is not a materialized directory: {root}")
    for path in sorted(root.rglob("*"), key=lambda value: len(value.parts), reverse=True):
        if path.is_symlink():
            raise FM0ContractError(f"release contains a symlink: {path}")
        path.chmod(0o550 if path.is_dir() else 0o440)
    root.chmod(0o550)


def _remove_build_tree(root: Path) -> None:
    """Remove only a builder-owned partial/final tree after a failed attempt."""

    if not root.exists():
        return
    if root.is_symlink() or not root.is_dir():
        raise FM0ContractError(f"refusing to clean unsafe build path: {root}")
    # Unlinking a read-only shard only requires write access to its parent.
    # Restore directory modes only; chmod-ing every one of ~200k shards would
    # duplicate the most expensive metadata traversal during a time-limit exit.
    for directory, names, _files in os.walk(root, topdown=True, followlinks=False):
        directory_path = Path(directory)
        directory_path.chmod(0o700)
        for name in names:
            child = directory_path / name
            if not child.is_symlink():
                child.chmod(0o700)
    shutil.rmtree(root)


def _output_leaf_path(value: str | Path) -> Path:
    """Resolve the parent while preserving the output leaf for symlink checks."""

    requested = Path(value).expanduser()
    if requested.name in {"", ".", ".."}:
        raise FM0ContractError(
            "later six-view output must name a concrete output directory"
        )
    return requested.parent.resolve() / requested.name


def validate_later_sector_release(
    root: str | Path,
    *,
    expected_receipt_sha256: str | None = None,
    require_read_only: bool = True,
) -> tuple[Path, Mapping[str, Any], dict[str, Any]]:
    """Revalidate a complete immutable later-sector six-view bundle."""

    release_root = Path(root).expanduser().resolve()
    if release_root.is_symlink() or not release_root.is_dir():
        raise FM0ContractError(f"six-view release root is missing: {release_root}")
    expected_root_entries = {
        "manifest.csv",
        "source_manifest.csv",
        "visit_timing.csv",
        "receipt.json",
        "READY",
        "shards",
    }
    observed_root_entries = {path.name for path in release_root.iterdir()}
    if observed_root_entries != expected_root_entries:
        raise FM0ContractError(
            "six-view release root closure drifted; "
            f"missing={sorted(expected_root_entries - observed_root_entries)}, "
            f"extra={sorted(observed_root_entries - expected_root_entries)}"
        )
    receipt_path = release_root / "receipt.json"
    if receipt_path.is_symlink() or not receipt_path.is_file():
        raise FM0ContractError("six-view receipt is not materialized")
    receipt_hash = sha256_file(receipt_path)
    if expected_receipt_sha256 is not None and receipt_hash != _digest(
        expected_receipt_sha256, label="expected six-view receipt SHA-256"
    ):
        raise FM0ContractError("six-view receipt caller binding drifted")
    receipt = _read_json(receipt_path, label="later six-view receipt")
    if set(receipt) != LATER_SIX_VIEW_RECEIPT_FIELDS:
        raise FM0ContractError("later six-view receipt fields drifted")
    try:
        sector = int(receipt["sector"])
        expected_orbits = tuple(int(value) for value in receipt["expected_orbits"])
        n_observations = int(receipt["n_observations"])
        n_shards = int(receipt["n_shards"])
        n_cadences = int(receipt["n_cadences"])
        n_segments = int(receipt["n_segments"])
        datetime.fromisoformat(str(receipt["created_utc"]))
    except (TypeError, ValueError) as exc:
        raise FM0ContractError("later six-view receipt numeric/time fields drifted") from exc
    provider = "spoc" if sector == 66 else "tica"
    expected_provider_or_none = provider if sector in range(66, 78) else None
    if (
        expected_provider_or_none is None
        or expected_orbits != (2 * sector + 7, 2 * sector + 8)
        or receipt["schema_version"] != LATER_SIX_VIEW_RECEIPT_SCHEMA_VERSION
        or receipt["release_state"] != LATER_SIX_VIEW_READY_STATE
        or receipt["passed"] is not True
        or receipt["product_state"] != LATER_SIX_VIEW_PRODUCT_STATE
        or _SHA40.fullmatch(str(receipt["producer_git_sha"])) is None
        or n_observations <= 0
        or n_shards != n_observations
        or n_cadences <= 0
        or n_segments <= 0
        or receipt["flux_view_names"] != list(FLUX_VIEW_NAMES)
        or receipt["error_view_names"] != list(ERROR_VIEW_NAMES)
        or receipt["mission_quality_reference_schema_version"]
        != MISSION_QUALITY_REFERENCE_SCHEMA_VERSION
        or receipt["mission_quality_provider"] != provider
        or receipt["mission_quality_composition"]
        != "mission_quality | (qlp_quality << 30)"
        or receipt["input_adapter"] != LATER_HDF5_ADAPTER_NAME
        or receipt["full_visit_shards"] is not True
        or receipt["model_context_length_bound"] is not False
        or receipt["sealed_hdf5_content_opened"] is not False
        or receipt["sealed_shards_written"] is not False
        or receipt["model_outcome_blind"] is not True
        or receipt["six_view_shards_verified"] is not True
        or receipt["scientific_training_eligible"] is not False
        or receipt["panel_admission_authorized"] is not False
    ):
        raise FM0ContractError("later six-view receipt contract drifted")
    for path_field, sha_field in (
        ("source_inventory_summary_path", "source_inventory_summary_sha256"),
        ("source_inventory_sources_path", "source_inventory_sources_sha256"),
        ("corpus_selection_path", "corpus_selection_sha256"),
        (
            "mission_quality_reference_manifest_path",
            "mission_quality_reference_manifest_sha256",
        ),
        (
            "mission_quality_reference_table_path",
            "mission_quality_reference_table_sha256",
        ),
        ("hdf5_quality_receipt_path", "hdf5_quality_receipt_sha256"),
    ):
        if not Path(str(receipt[path_field])).is_absolute():
            raise FM0ContractError(f"six-view receipt path is not absolute: {path_field}")
        _digest(receipt[sha_field], label=f"six-view {sha_field}")

    ready_path = release_root / "READY"
    if (
        ready_path.is_symlink()
        or not ready_path.is_file()
        or ready_path.read_text(encoding="utf-8").strip() != receipt_hash
    ):
        raise FM0ContractError("later six-view READY marker drifted")
    manifest_path = release_root / "manifest.csv"
    source_manifest_path = release_root / "source_manifest.csv"
    timing_path = release_root / "visit_timing.csv"
    for path, receipt_field in (
        (manifest_path, "manifest_sha256"),
        (source_manifest_path, "source_manifest_sha256"),
        (timing_path, "visit_timing_sha256"),
    ):
        if sha256_file(path) != _digest(receipt[receipt_field], label=receipt_field):
            raise FM0ContractError(f"later six-view {path.name} hash drifted")
    manifest_rows = _read_csv_exact(
        manifest_path,
        fields=LATER_SIX_VIEW_MANIFEST_FIELDS,
        label="later six-view manifest",
    )
    source_rows = _read_csv_exact(
        source_manifest_path,
        fields=SOURCE_FIELDS,
        label="later six-view source manifest",
    )
    timing_rows = _read_csv_exact(
        timing_path,
        fields=LATER_SIX_VIEW_TIMING_FIELDS,
        label="later six-view visit timing",
    )
    if not (
        len(manifest_rows) == len(source_rows) == len(timing_rows) == n_observations
    ):
        raise FM0ContractError("later six-view CSV row counts drifted")

    observation_keys: set[str] = set()
    physical_source_by_observation: dict[str, str] = {}
    source_row_hashes: set[str] = set()
    expected_shards: set[str] = set()
    manifest_total_cadences = 0
    manifest_total_segments = 0
    for row in manifest_rows:
        observation_key = str(row["observation_key"])
        gaia = str(row["gaia_dr3_source_id"])
        relative = PurePosixPath(str(row["relative_path"]))
        if (
            row["manifest_schema_version"] != LATER_SIX_VIEW_MANIFEST_SCHEMA_VERSION
            or not observation_key.startswith("observation_")
            or observation_key in observation_keys
            or row["physical_source_id"] != f"gaia_dr3:{gaia}"
            or int(row["sector"]) != sector
            or row["source_partition"] not in LATER_ALLOWED_SOURCE_PARTITIONS
            or row["product_state"] != LATER_SIX_VIEW_PRODUCT_STATE
            or relative.is_absolute()
            or ".." in relative.parts
            or relative.parts != ("shards", f"{observation_key}.npz")
            or row["input_adapter"] != LATER_HDF5_ADAPTER_NAME
            or row["mission_quality_provider"] != provider
            or row["mission_quality_reference_manifest_sha256"]
            != receipt["mission_quality_reference_manifest_sha256"]
            or row["hdf5_quality_receipt_sha256"]
            != receipt["hdf5_quality_receipt_sha256"]
            or row["full_visit_shard"] != "true"
            or row["model_context_length_bound"] != "false"
            or row["scientific_training_eligible"] != "false"
            or row["panel_admission_authorized"] != "false"
        ):
            raise FM0ContractError("later six-view manifest row contract drifted")
        _digest(row["sha256"], label="six-view shard SHA-256")
        source_row_hash = _digest(
            row["source_row_sha256"], label="six-view source-row SHA-256"
        )
        if source_row_hash in source_row_hashes:
            raise FM0ContractError("later six-view source-row binding is duplicated")
        source_row_hashes.add(source_row_hash)
        observation_keys.add(observation_key)
        physical_source_by_observation[observation_key] = row["physical_source_id"]
        expected_shards.add(relative.name)
        try:
            row_cadences = int(row["n_cadences"])
            row_segments = int(row["n_segments"])
        except ValueError as exc:
            raise FM0ContractError("six-view manifest counts are invalid") from exc
        if row_cadences <= 0 or row_segments <= 0:
            raise FM0ContractError("six-view manifest counts must be positive")
        manifest_total_cadences += row_cadences
        manifest_total_segments += row_segments
        _validate_shard(release_root / relative.as_posix(), row)
    if (
        manifest_total_cadences != n_cadences
        or manifest_total_segments != n_segments
    ):
        raise FM0ContractError("later six-view aggregate shard counts drifted")

    source_by_hash: dict[str, Mapping[str, Any]] = {}
    for row in source_rows:
        row_hash = _digest(row["source_row_sha256"], label="source-row SHA-256")
        if (
            row_hash in source_by_hash
            or int(row["sector"]) != sector
            or row["source_partition"] not in LATER_ALLOWED_SOURCE_PARTITIONS
            or row["product_state"] != LATER_SIX_VIEW_PRODUCT_STATE
            or row["model_outcome_blind"] != "true"
            or row["panel_admission_authorized"] != "false"
        ):
            raise FM0ContractError("later six-view source-manifest row drifted")
        source_by_hash[row_hash] = row
    if set(source_by_hash) != source_row_hashes:
        raise FM0ContractError("manifest/source-manifest closure drifted")
    manifest_by_source = {row["source_row_sha256"]: row for row in manifest_rows}
    for row_hash, source_row in source_by_hash.items():
        manifest_row = manifest_by_source[row_hash]
        for field in (
            "gaia_dr3_source_id",
            "tic_id",
            "leakage_component_id",
            "source_partition",
        ):
            if source_row[field] != manifest_row[field]:
                raise FM0ContractError("manifest/source identity binding drifted")

    timing_keys: set[str] = set()
    for row in timing_rows:
        try:
            start = float(row["absolute_visit_start"])
            end = float(row["absolute_visit_end"])
        except ValueError as exc:
            raise FM0ContractError("later visit timing is not numeric") from exc
        key = row["observation_key"]
        if (
            row["timing_schema_version"] != LATER_SIX_VIEW_TIMING_SCHEMA_VERSION
            or key in timing_keys
            or key not in observation_keys
            or row["physical_source_id"]
            != physical_source_by_observation.get(key)
            or not math.isfinite(start)
            or not math.isfinite(end)
            or not start < end
        ):
            raise FM0ContractError("later visit-timing row contract drifted")
        timing_keys.add(key)
    if timing_keys != observation_keys:
        raise FM0ContractError("manifest/timing observation closure drifted")

    shards_root = release_root / "shards"
    if shards_root.is_symlink() or not shards_root.is_dir():
        raise FM0ContractError("later six-view shards directory is invalid")
    actual_shards = {
        path.name
        for path in shards_root.iterdir()
        if path.is_file() and not path.is_symlink()
    }
    if actual_shards != expected_shards or any(
        not path.is_file() or path.is_symlink() for path in shards_root.iterdir()
    ):
        raise FM0ContractError("later six-view shard closure drifted")
    rows_by_partition = {
        partition: sum(row["source_partition"] == partition for row in manifest_rows)
        for partition in LATER_SIX_VIEW_ALLOWED_PARTITIONS
    }
    if receipt["source_rows_by_partition"] != rows_by_partition:
        raise FM0ContractError("later six-view receipt partition counts drifted")
    for field in (
        "n_internal_bad_cadences",
        "n_external_bad_cadences",
        "n_authority_excluded_cadences",
        "n_sealed_source_rows_excluded",
    ):
        try:
            if int(receipt[field]) < 0:
                raise ValueError
        except (TypeError, ValueError) as exc:
            raise FM0ContractError(f"later six-view receipt {field} is invalid") from exc

    if require_read_only:
        for path in (release_root, *release_root.rglob("*")):
            if path.is_symlink() or path.stat().st_mode & 0o222:
                raise FM0ContractError(f"six-view release is not read-only: {path}")
    return receipt_path, receipt, {
        "receipt_sha256": receipt_hash,
        "sector": sector,
        "mission_quality_provider": provider,
        "n_observations": n_observations,
        "n_shards": n_shards,
        "n_cadences": n_cadences,
        "manifest_rows": tuple(manifest_rows),
        "source_rows": tuple(source_rows),
        "timing_rows": tuple(timing_rows),
    }


def _build_later_sector_release(
    *,
    sector: int,
    source_inventory_dir: str | Path,
    expected_source_inventory_summary_sha256: str,
    mission_quality_reference_dir: str | Path,
    expected_mission_quality_reference_manifest_sha256: str,
    hdf5_quality_receipt_path: str | Path,
    expected_hdf5_quality_receipt_sha256: str,
    output_dir: str | Path,
    producer_git_sha: str,
) -> dict[str, Any]:
    """Publish one immutable later-sector, full-visit six-view bundle."""

    sector = int(sector)
    producer = _digest(producer_git_sha, label="producer_git_sha", length=40)
    if sector not in range(66, 78):
        raise FM0ContractError("later six-view release is restricted to S66--S77")
    expected_orbits = (2 * sector + 7, 2 * sector + 8)
    (
        inventory_summary,
        inventory_summary_path,
        inventory_summary_hash,
        all_source_rows,
        inventory_sources_path,
        inventory_sources_hash,
    ) = _load_inventory_bundle(
        source_inventory_dir=source_inventory_dir,
        sector=sector,
        expected_summary_sha256=expected_source_inventory_summary_sha256,
    )

    # Partition filtering precedes reference loading and, critically, every
    # retained-HDF5 loader invocation.  Sealed paths are never passed onward.
    selected_rows = [
        row
        for row in all_source_rows
        if row["source_partition"] in LATER_ALLOWED_SOURCE_PARTITIONS
    ]
    sealed_rows = [
        row
        for row in all_source_rows
        if row["source_partition"] == LATER_FORBIDDEN_SOURCE_PARTITION
    ]
    if len(selected_rows) + len(sealed_rows) != len(all_source_rows) or not selected_rows:
        raise FM0ContractError("later source partitions are incomplete or unusable")

    expected_reference_manifest_hash = _digest(
        expected_mission_quality_reference_manifest_sha256,
        label="expected mission-quality reference manifest SHA-256",
    )
    reference_root = Path(mission_quality_reference_dir).expanduser().resolve()
    reference_manifest_path = reference_root / "manifest.json"
    if (
        reference_manifest_path.is_symlink()
        or not reference_manifest_path.is_file()
        or sha256_file(reference_manifest_path) != expected_reference_manifest_hash
    ):
        raise FM0ContractError(
            "mission-quality reference manifest caller binding drifted"
        )
    reference = load_mission_quality_reference(
        reference_dir=mission_quality_reference_dir,
        sector=sector,
        expected_orbits=expected_orbits,
    )
    expected_provider = "spoc" if sector == 66 else "tica"
    if (
        reference.provider != expected_provider
        or reference.manifest_sha256 != expected_reference_manifest_hash
    ):
        raise FM0ContractError("mission-quality reference provider/binding is invalid")
    _hdf5_receipt, hdf5_receipt_path, hdf5_receipt_hash = (
        _load_hdf5_quality_receipt(
            path=hdf5_quality_receipt_path,
            sector=sector,
            inventory_summary=inventory_summary,
            reference=reference,
            expected_receipt_sha256=expected_hdf5_quality_receipt_sha256,
        )
    )

    cache = build_later_adapter_cache(
        source_receipt_path=inventory_summary["source_receipt_path"],
        expected_source_receipt_sha256=inventory_summary["source_receipt_sha256"],
        sector=sector,
        expected_orbits=expected_orbits,
        expected_target_authority_sha256=inventory_summary[
            "target_authority_sha256"
        ],
        reference=reference,
    )
    cache.track_file(inventory_summary_path, inventory_summary_hash)
    cache.track_file(inventory_sources_path, inventory_sources_hash)
    cache.track_file(
        inventory_summary["identity_binding_path"],
        inventory_summary["identity_binding_sha256"],
    )
    cache.track_file(
        inventory_summary["corpus_summary_path"],
        inventory_summary["corpus_summary_sha256"],
    )
    cache.track_file(
        inventory_summary["gaia_tic_alias_authority_path"],
        inventory_summary["gaia_tic_alias_authority_sha256"],
    )
    cache.track_marker(inventory_summary_path.parent / "READY", inventory_summary_hash)
    cache.track_file(hdf5_receipt_path, hdf5_receipt_hash)

    final = _output_leaf_path(output_dir)
    partial = final.with_name(final.name + ".partial")
    if (
        final.exists()
        or final.is_symlink()
        or partial.exists()
        or partial.is_symlink()
    ):
        raise FM0ContractError(f"refusing to overwrite later six-view release: {final}")
    partial.mkdir(parents=True)
    (partial / "shards").mkdir()

    manifest_rows: list[dict[str, Any]] = []
    timing_rows: list[dict[str, Any]] = []
    selected_source_rows: list[dict[str, Any]] = []
    quality_totals = {
        "n_internal_bad_cadences": 0,
        "n_external_bad_cadences": 0,
        "n_authority_excluded_cadences": 0,
    }
    seen_observations: set[str] = set()
    selected_rows.sort(
        key=lambda row: (int(row["gaia_dr3_source_id"]), int(row["tic_id"]))
    )
    for source_row in selected_rows:
        source = load_later_hdf5_observation(
            source_row,
            reference=reference,
            hdf5_quality_receipt_sha256=hdf5_receipt_hash,
            cache=cache,
        )
        release = build_mission_quality_observation_release(
            source.raw_arrays,
            mission_quality_provider=reference.provider,
            input_adapter=LATER_HDF5_ADAPTER_NAME,
            scientific_training_eligible=False,
        )
        if release.audit is None:
            raise FM0ContractError("later six-view builder omitted visit audit")
        gaia = str(source_row["gaia_dr3_source_id"])
        tic = str(source_row["tic_id"])
        physical_source_id = f"gaia_dr3:{gaia}"
        observation_key = _stable_id(
            "observation",
            {"physical_source_id": physical_source_id, "sector": sector},
        )
        if observation_key in seen_observations:
            raise FM0ContractError("later release has duplicate physical-source visit")
        seen_observations.add(observation_key)
        product_instance_id = _stable_id(
            "product",
            {
                "observation_key": observation_key,
                "tic_id": tic,
                "a2v1_product_version": "A2v1_ORCD_DEFERRED",
                "source_sha256": source.source_sha256,
            },
        )
        relative_path = f"shards/{observation_key}.npz"
        shard_path = partial / relative_path
        payload = deterministic_npz_bytes(release)
        publish_immutable(shard_path, payload)
        shard_hash = hashlib.sha256(payload).hexdigest()
        manifest_row = {
            "manifest_schema_version": LATER_SIX_VIEW_MANIFEST_SCHEMA_VERSION,
            "observation_key": observation_key,
            "product_instance_id": product_instance_id,
            "physical_source_id": physical_source_id,
            "gaia_dr3_source_id": gaia,
            "tic_id": tic,
            "sector": sector,
            "leakage_component_id": source_row["leakage_component_id"],
            "source_partition": source_row["source_partition"],
            "product_state": LATER_SIX_VIEW_PRODUCT_STATE,
            "relative_path": relative_path,
            "sha256": shard_hash,
            "input_source_sha256": source.source_sha256,
            "source_row_sha256": source_row["source_row_sha256"],
            "n_cadences": release.n_cadences,
            "n_segments": release.n_segments,
            "view_present_json": json.dumps(
                release.view_present.astype(int).tolist(), separators=(",", ":")
            ),
            "input_adapter": LATER_HDF5_ADAPTER_NAME,
            "mission_quality_provider": reference.provider,
            "mission_quality_reference_manifest_sha256": reference.manifest_sha256,
            "hdf5_quality_receipt_sha256": hdf5_receipt_hash,
            "full_visit_shard": "true",
            "model_context_length_bound": "false",
            "scientific_training_eligible": "false",
            "panel_admission_authorized": "false",
        }
        manifest_rows.append(manifest_row)
        timing_rows.append(
            {
                "timing_schema_version": LATER_SIX_VIEW_TIMING_SCHEMA_VERSION,
                "observation_key": observation_key,
                "physical_source_id": physical_source_id,
                "absolute_visit_start": repr(
                    float(release.audit["absolute_visit_start"])
                ),
                "absolute_visit_end": repr(float(release.audit["absolute_visit_end"])),
            }
        )
        selected_source_rows.append(dict(source_row))
        counts = source.audit["quality_counts"]
        quality_totals["n_internal_bad_cadences"] += int(counts["n_internal_bad"])
        quality_totals["n_external_bad_cadences"] += int(counts["n_external_bad"])
        quality_totals["n_authority_excluded_cadences"] += int(
            counts["n_authority_excluded"]
        )
        if len(manifest_rows) % 1000 == 0:
            print(
                "FM_LATER_SIX_VIEW_PROGRESS "
                f"sector={sector} observations={len(manifest_rows)} "
                f"cadences={sum(int(row['n_cadences']) for row in manifest_rows)}",
                flush=True,
            )

    cache.assert_unchanged()
    manifest_rows.sort(key=lambda row: row["observation_key"])
    timing_rows.sort(key=lambda row: row["observation_key"])
    selected_source_rows.sort(
        key=lambda row: (int(row["gaia_dr3_source_id"]), int(row["tic_id"]))
    )
    manifest_payload = _csv_bytes(manifest_rows, LATER_SIX_VIEW_MANIFEST_FIELDS)
    timing_payload = _csv_bytes(timing_rows, LATER_SIX_VIEW_TIMING_FIELDS)
    source_payload = _csv_bytes(selected_source_rows, SOURCE_FIELDS)
    publish_immutable(partial / "manifest.csv", manifest_payload)
    publish_immutable(partial / "visit_timing.csv", timing_payload)
    publish_immutable(partial / "source_manifest.csv", source_payload)

    n_cadences = sum(int(row["n_cadences"]) for row in manifest_rows)
    rows_by_partition = {
        partition: sum(row["source_partition"] == partition for row in manifest_rows)
        for partition in LATER_SIX_VIEW_ALLOWED_PARTITIONS
    }
    receipt = {
        "schema_version": LATER_SIX_VIEW_RECEIPT_SCHEMA_VERSION,
        "release_state": LATER_SIX_VIEW_READY_STATE,
        "passed": True,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "producer_git_sha": producer,
        "sector": sector,
        "expected_orbits": list(expected_orbits),
        "product_state": LATER_SIX_VIEW_PRODUCT_STATE,
        "n_observations": len(manifest_rows),
        "n_shards": len(manifest_rows),
        "n_cadences": n_cadences,
        "n_segments": sum(int(row["n_segments"]) for row in manifest_rows),
        "source_rows_by_partition": rows_by_partition,
        "n_sealed_source_rows_excluded": len(sealed_rows),
        "flux_view_names": list(FLUX_VIEW_NAMES),
        "error_view_names": list(ERROR_VIEW_NAMES),
        "source_inventory_summary_path": str(inventory_summary_path),
        "source_inventory_summary_sha256": inventory_summary_hash,
        "source_inventory_sources_path": str(inventory_sources_path),
        "source_inventory_sources_sha256": inventory_sources_hash,
        "corpus_selection_path": str(inventory_summary["identity_binding_path"]),
        "corpus_selection_sha256": inventory_summary["identity_binding_sha256"],
        "mission_quality_reference_schema_version": (
            MISSION_QUALITY_REFERENCE_SCHEMA_VERSION
        ),
        "mission_quality_provider": reference.provider,
        "mission_quality_composition": "mission_quality | (qlp_quality << 30)",
        "mission_quality_reference_manifest_path": str(reference.manifest_path),
        "mission_quality_reference_manifest_sha256": reference.manifest_sha256,
        "mission_quality_reference_table_path": str(reference.table_path),
        "mission_quality_reference_table_sha256": reference.table_sha256,
        "hdf5_quality_receipt_path": str(hdf5_receipt_path),
        "hdf5_quality_receipt_sha256": hdf5_receipt_hash,
        "manifest_sha256": hashlib.sha256(manifest_payload).hexdigest(),
        "source_manifest_sha256": hashlib.sha256(source_payload).hexdigest(),
        "visit_timing_sha256": hashlib.sha256(timing_payload).hexdigest(),
        "input_adapter": LATER_HDF5_ADAPTER_NAME,
        "full_visit_shards": True,
        "model_context_length_bound": False,
        **quality_totals,
        "sealed_hdf5_content_opened": False,
        "sealed_shards_written": False,
        "model_outcome_blind": True,
        "six_view_shards_verified": True,
        "scientific_training_eligible": False,
        "panel_admission_authorized": False,
        "claim_limit": (
            "deferred ORCD full-visit six-view preparation only; not accepted "
            "A2v1, a temporal panel, training authorization, or a model result"
        ),
    }
    receipt_payload = (
        json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    publish_immutable(partial / "receipt.json", receipt_payload)
    publish_immutable(
        partial / "READY", hashlib.sha256(receipt_payload).hexdigest().encode("ascii") + b"\n"
    )
    _make_tree_read_only(partial)
    os.replace(partial, final)
    return receipt


def build_later_sector_release(
    *,
    sector: int,
    source_inventory_dir: str | Path,
    expected_source_inventory_summary_sha256: str,
    mission_quality_reference_dir: str | Path,
    expected_mission_quality_reference_manifest_sha256: str,
    hdf5_quality_receipt_path: str | Path,
    expected_hdf5_quality_receipt_sha256: str,
    output_dir: str | Path,
    producer_git_sha: str,
) -> dict[str, Any]:
    """Build, fully validate, and publish one read-only six-view release."""

    sector = int(sector)
    if sector not in range(66, 78):
        raise FM0ContractError("later six-view release is restricted to S66--S77")
    final = _output_leaf_path(output_dir)
    partial = final.with_name(final.name + ".partial")
    if (
        final.exists()
        or final.is_symlink()
        or partial.exists()
        or partial.is_symlink()
    ):
        raise FM0ContractError(f"refusing to overwrite later six-view release: {final}")
    try:
        receipt = _build_later_sector_release(
            sector=sector,
            source_inventory_dir=source_inventory_dir,
            expected_source_inventory_summary_sha256=(
                expected_source_inventory_summary_sha256
            ),
            mission_quality_reference_dir=mission_quality_reference_dir,
            expected_mission_quality_reference_manifest_sha256=(
                expected_mission_quality_reference_manifest_sha256
            ),
            hdf5_quality_receipt_path=hdf5_quality_receipt_path,
            expected_hdf5_quality_receipt_sha256=(
                expected_hdf5_quality_receipt_sha256
            ),
            output_dir=final,
            producer_git_sha=producer_git_sha,
        )
        receipt_hash = sha256_file(final / "receipt.json")
        _receipt_path, validated_receipt, _summary = validate_later_sector_release(
            final,
            expected_receipt_sha256=receipt_hash,
            require_read_only=True,
        )
        if dict(validated_receipt) != receipt:
            raise FM0ContractError("published six-view receipt differs from build result")
        return receipt
    except Exception:
        # Both paths were proven absent above, so any tree now present belongs
        # to this failed attempt and may be safely removed for a clean retry.
        cleanup_error: Exception | None = None
        for path in (partial, final):
            try:
                _remove_build_tree(path)
            except (FM0ContractError, OSError) as exc:
                cleanup_error = exc
        if cleanup_error is not None:
            raise FM0ContractError(
                f"six-view build failed and cleanup also failed: {cleanup_error}"
            ) from cleanup_error
        raise


__all__ = [
    "LATER_SIX_VIEW_ALLOWED_PARTITIONS",
    "LATER_SIX_VIEW_MANIFEST_FIELDS",
    "LATER_SIX_VIEW_MANIFEST_SCHEMA_VERSION",
    "LATER_SIX_VIEW_PRODUCT_STATE",
    "LATER_SIX_VIEW_READY_STATE",
    "LATER_SIX_VIEW_RECEIPT_FIELDS",
    "LATER_SIX_VIEW_RECEIPT_SCHEMA_VERSION",
    "LATER_SIX_VIEW_TIMING_FIELDS",
    "LATER_SIX_VIEW_TIMING_SCHEMA_VERSION",
    "build_later_sector_release",
    "validate_later_sector_release",
]
