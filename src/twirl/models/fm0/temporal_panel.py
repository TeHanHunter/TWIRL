"""Freeze the label-blind S66--S77 FM development-evaluation panel.

The freezer consumes only checksum-bound JSON/CSV authorities.  It never
resolves, stats, hashes, or opens a six-view NPZ shard.  The emitted table is
therefore an identity/provenance panel, not a training release or a model
result.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path, PurePosixPath
from typing import Any

import yaml

from .input_release import (
    FLUX_VIEW_NAMES,
    INPUT_RELEASE_SCHEMA_VERSION,
    MANIFEST_COLUMNS,
)
from .later_sector_admission_v2 import (
    DEFERRED_PRODUCT_STATE,
    EXCLUDED_SECTORS,
    POOL_RECEIPT_SCHEMA_VERSION,
)
from .later_sector_release import (
    LATER_SIX_VIEW_MANIFEST_FIELDS,
    LATER_SIX_VIEW_MANIFEST_SCHEMA_VERSION,
    LATER_SIX_VIEW_READY_STATE,
    LATER_SIX_VIEW_RECEIPT_SCHEMA_VERSION,
)
from .registry import (
    FM0ContractError,
    assert_no_search_columns,
    publish_immutable,
    sha256_file,
)

TEMPORAL_PANEL_CONFIG_SCHEMA_VERSION = "twirl_fm0_temporal_panel_config_v1"
TEMPORAL_PANEL_CAMPAIGN_ID = "twirl_fm0_2_s66_s77_temporal_panel_v1"
TEMPORAL_PANEL_FROZEN_POLICY_SHA256 = (
    "1927d1fa66a87e2ccb42a5661d56b53dc8e250ab8ac0ccaf6c67b01e8632ea9f"
)
TEMPORAL_PANEL_SCHEMA_VERSION = "twirl_fm0_s66_s77_temporal_panel_v1"
TEMPORAL_PANEL_RECEIPT_SCHEMA_VERSION = "twirl_fm0_s66_s77_temporal_panel_receipt_v1"
TEMPORAL_PANEL_SECTOR_BINDING_SCHEMA_VERSION = (
    "twirl_fm0_s66_s77_temporal_panel_sector_binding_v1"
)
TEMPORAL_PANEL_READY_STATE = "FM_TEMPORAL_PANEL_READY"
TEMPORAL_PANEL_SECTORS = tuple(range(66, 78))
BASELINE_SECTORS = tuple(range(56, 65))
DEVELOPMENT_PARTITION = "poc_development"
TRAIN_PARTITION = "poc_train"
SEALED_PARTITION = "poc_sealed_test"

TEMPORAL_PANEL_FIELDS = (
    "panel_schema_version",
    "observation_key",
    "product_instance_id",
    "physical_source_id",
    "gaia_dr3_source_id",
    "tic_id",
    "sector",
    "leakage_component_id",
    "source_partition",
    "temporal_cohort",
    "sector_release_root",
    "sector_receipt_sha256",
    "sector_manifest_sha256",
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
    "development_evaluation_eligible",
)

TEMPORAL_PANEL_SECTOR_BINDING_FIELDS = (
    "binding_schema_version",
    "sector",
    "sector_release_root",
    "sector_receipt_path",
    "sector_receipt_sha256",
    "sector_manifest_path",
    "sector_manifest_sha256",
    "n_manifest_rows",
    "n_development_rows",
    "n_train_rows_excluded",
    "n_sealed_rows_emitted",
    "scientific_training_eligible",
    "development_evaluation_eligible",
)

_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_GIT_SHA = re.compile(r"^[0-9a-f]{40}$")
_POSITIVE_ID = re.compile(r"^[1-9][0-9]*$")
_LEAKAGE_COMPONENT = re.compile(r"^leakage_[0-9a-f]{64}$")
_ALLOWED_RELEASE_PARTITIONS = frozenset({TRAIN_PARTITION, DEVELOPMENT_PARTITION})
_ALLOWED_BASELINE_PARTITIONS = frozenset(
    {TRAIN_PARTITION, DEVELOPMENT_PARTITION, SEALED_PARTITION}
)


@dataclass(frozen=True)
class TemporalPanelFreeze:
    """Paths and receipt for one immutable temporal-panel freeze."""

    output_dir: Path
    panel_path: Path
    sector_bindings_path: Path
    receipt_path: Path
    ready_path: Path
    receipt_sha256: str
    receipt: Mapping[str, Any]


def _mapping(value: Any, *, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise FM0ContractError(f"{label} must be a mapping")
    return value


def _sequence(value: Any, *, label: str) -> Sequence[Any]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise FM0ContractError(f"{label} must be a sequence")
    return value


def _digest(value: Any, *, label: str) -> str:
    result = str(value).strip().lower()
    if _SHA256.fullmatch(result) is None:
        raise FM0ContractError(f"{label} must be a lowercase SHA-256 digest")
    return result


def _integer(value: Any, *, label: str, allow_zero: bool = False) -> int:
    if isinstance(value, bool):
        raise FM0ContractError(f"{label} must be an integer")
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise FM0ContractError(f"{label} must be an integer") from exc
    if isinstance(value, str) and value.strip() != str(result):
        raise FM0ContractError(f"{label} must be an exact integer")
    if result < (0 if allow_zero else 1):
        raise FM0ContractError(f"{label} is outside its allowed range")
    return result


def _bool_text(value: Any, *, label: str) -> bool:
    normalized = str(value).strip().lower()
    if normalized in {"true", "1"}:
        return True
    if normalized in {"false", "0"}:
        return False
    raise FM0ContractError(f"{label} must be a serialized boolean")


def _require_read_only(path: Path, *, label: str) -> None:
    if path.stat().st_mode & 0o222:
        raise FM0ContractError(f"{label} is not read-only: {path}")


def _verify_file(
    path: str | Path,
    expected_sha256: str,
    *,
    label: str,
    require_read_only: bool = False,
) -> tuple[Path, str]:
    expected = _digest(expected_sha256, label=f"{label} SHA-256")
    raw = Path(path).expanduser()
    if raw.is_symlink():
        raise FM0ContractError(f"{label} must be materialized, not a symlink")
    source = raw.resolve()
    if not source.is_file() or source.stat().st_size <= 0:
        raise FM0ContractError(f"missing materialized {label}: {source}")
    if require_read_only:
        _require_read_only(source, label=label)
    observed = sha256_file(source)
    if observed != expected:
        raise FM0ContractError(
            f"{label} hash mismatch: expected {expected}, observed {observed}"
        )
    return source, observed


def _read_json(path: Path, *, label: str) -> Mapping[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FM0ContractError(f"invalid {label}: {path}") from exc
    return _mapping(value, label=label)


def _read_csv(path: Path, fields: Sequence[str], *, label: str) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if tuple(reader.fieldnames or ()) != tuple(fields):
                raise FM0ContractError(f"{label} columns drifted")
            rows = [dict(row) for row in reader]
    except (OSError, UnicodeError, csv.Error) as exc:
        raise FM0ContractError(f"invalid {label}: {path}") from exc
    if not rows:
        raise FM0ContractError(f"{label} is empty")
    return rows


def _csv_bytes(rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> bytes:
    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(fields), lineterminator="\n")
    writer.writeheader()
    for row in rows:
        if tuple(row) != tuple(fields):
            raise FM0ContractError("temporal-panel CSV columns drifted")
        writer.writerow({field: row[field] for field in fields})
    return stream.getvalue().encode("utf-8")


def _validate_policy(config: Mapping[str, Any]) -> str:
    if config.get("schema_version") != TEMPORAL_PANEL_CONFIG_SCHEMA_VERSION:
        raise FM0ContractError("temporal-panel config schema mismatch")
    scope = _mapping(config.get("scope"), label="temporal-panel scope")
    admission = _mapping(config.get("admission"), label="temporal-panel admission")
    release = _mapping(
        config.get("sector_release"), label="temporal-panel sector release"
    )
    selection = _mapping(config.get("selection"), label="temporal-panel selection")

    if (
        scope.get("label_blind") is not True
        or scope.get("identity_only_freeze") is not True
        or scope.get("shard_payload_access_authorized") is not False
        or scope.get("model_training_authorized") is not False
        or scope.get("sealed_test_access_authorized") is not False
    ):
        raise FM0ContractError("temporal-panel scope drifted")
    for field in (
        "event_retention_authorized",
        "formal_model_gate_authorized",
        "production_model_claim",
        "prospective_test_claim",
    ):
        if scope.get(field) is not False:
            raise FM0ContractError(f"temporal-panel scope improperly enables {field}")

    if (
        admission.get("required_schema_version") != POOL_RECEIPT_SCHEMA_VERSION
        or admission.get("required_product_state") != DEFERRED_PRODUCT_STATE
        or admission.get("required_preparation_pool_admitted") is not True
        or tuple(int(value) for value in admission.get("required_sectors", ()))
        != TEMPORAL_PANEL_SECTORS
        or tuple(int(value) for value in admission.get("excluded_sectors", ()))
        != EXCLUDED_SECTORS
        or admission.get("later_sector_policy") != "fail"
    ):
        raise FM0ContractError("temporal-panel admission policy drifted")

    if (
        release.get("required_receipt_schema_version")
        != LATER_SIX_VIEW_RECEIPT_SCHEMA_VERSION
        or release.get("required_release_state") != LATER_SIX_VIEW_READY_STATE
        or release.get("required_manifest_schema_version")
        != LATER_SIX_VIEW_MANIFEST_SCHEMA_VERSION
        or release.get("exact_receipt_and_manifest_sha256_required") is not True
        or release.get("read_only_materialized_authorities_required") is not True
        or release.get("full_visit_shards_required") is not True
        or release.get("model_context_length_bound") is not False
        or tuple(release.get("source_partitions_allowed_in_release", ()))
        != (TRAIN_PARTITION, DEVELOPMENT_PARTITION)
        or release.get("sealed_rows_allowed_in_release") is not False
    ):
        raise FM0ContractError("temporal-panel sector-release policy drifted")

    baseline_hash = _digest(
        selection.get("baseline_manifest_sha256"),
        label="baseline manifest policy SHA-256",
    )
    if (
        selection.get("emitted_source_partition") != DEVELOPMENT_PARTITION
        or tuple(int(value) for value in selection.get("baseline_manifest_sectors", ()))
        != BASELINE_SECTORS
        or selection.get("repeated_definition")
        != "leakage_component_present_in_frozen_baseline_manifest"
        or selection.get("new_definition")
        != "leakage_component_absent_from_frozen_baseline_manifest"
        or selection.get("retain_all_eligible_visits") is not True
        or selection.get("duplicate_observation_policy") != "fail"
        or selection.get("component_partition_drift_policy") != "fail"
        or selection.get("scientific_training_eligible") is not False
        or selection.get("development_evaluation_eligible") is not True
    ):
        raise FM0ContractError("temporal-panel selection policy drifted")

    forbidden = _sequence(
        config.get("forbidden_input_tokens"), label="forbidden input tokens"
    )
    if not forbidden or any(not str(value).strip() for value in forbidden):
        raise FM0ContractError("temporal-panel forbidden-input policy is empty")
    if not str(config.get("claim_limit", "")).strip():
        raise FM0ContractError("temporal-panel claim limit is empty")
    return baseline_hash


def _load_policy(
    path: str | Path, *, expected_sha256: str
) -> tuple[Mapping[str, Any], Path, str, str]:
    source, observed = _verify_file(
        path, expected_sha256, label="temporal-panel config"
    )
    try:
        value = yaml.safe_load(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        raise FM0ContractError(f"invalid temporal-panel config: {source}") from exc
    config = _mapping(value, label="temporal-panel config")
    if (
        config.get("campaign_id") == TEMPORAL_PANEL_CAMPAIGN_ID
        and observed != TEMPORAL_PANEL_FROZEN_POLICY_SHA256
    ):
        raise FM0ContractError("production temporal-panel config hash drifted")
    baseline_hash = _validate_policy(config)
    return config, source, observed, baseline_hash


def _load_baseline_components(
    path: str | Path, *, expected_sha256: str
) -> tuple[Path, str, dict[str, str], int]:
    source, observed = _verify_file(
        path, expected_sha256, label="frozen S56--S64 manifest"
    )
    rows = _read_csv(source, MANIFEST_COLUMNS, label="frozen S56--S64 manifest")
    assert_no_search_columns(MANIFEST_COLUMNS, context="baseline manifest")
    component_partition: dict[str, str] = {}
    observation_keys: set[str] = set()
    product_ids: set[str] = set()
    for row in rows:
        component = str(row["leakage_component_id"]).strip()
        partition = str(row["source_partition"]).strip()
        observation_key = str(row["observation_key"]).strip()
        product_id = str(row["product_instance_id"]).strip()
        if (
            row["input_release_schema_version"] != INPUT_RELEASE_SCHEMA_VERSION
            or _LEAKAGE_COMPONENT.fullmatch(component) is None
            or partition not in _ALLOWED_BASELINE_PARTITIONS
            or not observation_key
            or not product_id
        ):
            raise FM0ContractError("baseline manifest identity contract drifted")
        prior = component_partition.setdefault(component, partition)
        if prior != partition:
            raise FM0ContractError(
                "baseline leakage component crosses source partitions"
            )
        if observation_key in observation_keys or product_id in product_ids:
            raise FM0ContractError("baseline manifest contains duplicate identity rows")
        observation_keys.add(observation_key)
        product_ids.add(product_id)
    return source, observed, component_partition, len(rows)


def _safe_shard_reference(value: Any, *, label: str) -> str:
    text = str(value).strip()
    path = PurePosixPath(text)
    if (
        not text
        or path.is_absolute()
        or len(path.parts) != 2
        or path.parts[0] != "shards"
        or path.parts[1] in {"", ".", ".."}
        or path.suffix != ".npz"
        or ".." in path.parts
    ):
        raise FM0ContractError(f"{label} is not a safe sector-relative NPZ path")
    return text


def _validate_manifest_row(row: Mapping[str, str], *, sector: int) -> None:
    component = str(row["leakage_component_id"]).strip()
    partition = str(row["source_partition"]).strip()
    gaia = str(row["gaia_dr3_source_id"]).strip()
    tic = str(row["tic_id"]).strip()
    if (
        row["manifest_schema_version"] != LATER_SIX_VIEW_MANIFEST_SCHEMA_VERSION
        or _integer(row["sector"], label="manifest sector") != sector
        or row["product_state"] != DEFERRED_PRODUCT_STATE
        or partition not in _ALLOWED_RELEASE_PARTITIONS
        or _LEAKAGE_COMPONENT.fullmatch(component) is None
        or _POSITIVE_ID.fullmatch(gaia) is None
        or _POSITIVE_ID.fullmatch(tic) is None
        or row["physical_source_id"] != f"gaia_dr3:{gaia}"
        or not row["observation_key"].strip()
        or not row["product_instance_id"].strip()
        or not row["input_adapter"].strip()
        or not row["mission_quality_provider"].strip()
    ):
        raise FM0ContractError(f"S{sector} six-view manifest identity drifted")
    _safe_shard_reference(row["relative_path"], label=f"S{sector} shard reference")
    for field in (
        "sha256",
        "input_source_sha256",
        "source_row_sha256",
        "mission_quality_reference_manifest_sha256",
        "hdf5_quality_receipt_sha256",
    ):
        _digest(row[field], label=f"S{sector} manifest {field}")
    _integer(row["n_cadences"], label=f"S{sector} manifest n_cadences")
    _integer(row["n_segments"], label=f"S{sector} manifest n_segments")
    try:
        present = json.loads(row["view_present_json"])
    except json.JSONDecodeError as exc:
        raise FM0ContractError(f"S{sector} view_present_json is invalid") from exc
    if (
        not isinstance(present, list)
        or len(present) != len(FLUX_VIEW_NAMES)
        or any(value not in (0, 1, False, True) for value in present)
        or _bool_text(row["full_visit_shard"], label="full_visit_shard") is not True
        or _bool_text(
            row["model_context_length_bound"], label="model_context_length_bound"
        )
        is not False
        or _bool_text(
            row["scientific_training_eligible"],
            label="scientific_training_eligible",
        )
        is not False
        or _bool_text(
            row["panel_admission_authorized"],
            label="panel_admission_authorized",
        )
        is not False
    ):
        raise FM0ContractError(f"S{sector} six-view manifest boundary drifted")


def _load_sector_manifest(
    *,
    admission_sector: Mapping[str, Any],
    sector: int,
    seen_observations: set[str],
    seen_products: set[str],
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    if (
        _integer(admission_sector.get("sector"), label="admission sector") != sector
        or tuple(
            int(value)
            for value in _sequence(
                admission_sector.get("expected_orbits"), label="expected orbits"
            )
        )
        != (2 * sector + 7, 2 * sector + 8)
        or admission_sector.get("product_state") != DEFERRED_PRODUCT_STATE
        or admission_sector.get("preparation_admitted") is not True
        or admission_sector.get("a2v1_accepted") is not False
        or admission_sector.get("scientific_training_eligible") is not False
    ):
        raise FM0ContractError(f"S{sector} admission-sector contract drifted")
    evidence = _mapping(admission_sector.get("evidence"), label=f"S{sector} evidence")
    summary = _mapping(
        evidence.get("six_view_release"), label=f"S{sector} six-view summary"
    )
    if (
        summary.get("schema_version") != LATER_SIX_VIEW_RECEIPT_SCHEMA_VERSION
        or summary.get("release_state") != LATER_SIX_VIEW_READY_STATE
    ):
        raise FM0ContractError(f"S{sector} six-view admission summary drifted")
    receipt_sha = _digest(summary.get("sha256"), label=f"S{sector} receipt SHA-256")
    receipt_path_text = str(summary.get("path", "")).strip()
    if not receipt_path_text or not Path(receipt_path_text).expanduser().is_absolute():
        raise FM0ContractError(f"S{sector} receipt path must be absolute")
    receipt_path, _ = _verify_file(
        receipt_path_text,
        receipt_sha,
        label=f"S{sector} six-view receipt",
        require_read_only=True,
    )
    root = receipt_path.parent
    if root.is_symlink() or not root.is_dir():
        raise FM0ContractError(f"S{sector} release root is invalid")
    _require_read_only(root, label=f"S{sector} release root")
    receipt = _read_json(receipt_path, label=f"S{sector} six-view receipt")
    if (
        receipt.get("schema_version") != LATER_SIX_VIEW_RECEIPT_SCHEMA_VERSION
        or receipt.get("release_state") != LATER_SIX_VIEW_READY_STATE
        or receipt.get("passed") is not True
        or _integer(receipt.get("sector"), label=f"S{sector} receipt sector") != sector
        or tuple(int(value) for value in receipt.get("expected_orbits", ()))
        != (2 * sector + 7, 2 * sector + 8)
        or receipt.get("product_state") != DEFERRED_PRODUCT_STATE
        or receipt.get("full_visit_shards") is not True
        or receipt.get("model_context_length_bound") is not False
        or receipt.get("sealed_hdf5_content_opened") is not False
        or receipt.get("sealed_shards_written") is not False
        or receipt.get("model_outcome_blind") is not True
        or receipt.get("six_view_shards_verified") is not True
        or receipt.get("scientific_training_eligible") is not False
        or receipt.get("panel_admission_authorized") is not False
    ):
        raise FM0ContractError(f"S{sector} six-view receipt contract drifted")

    ready = root / "READY"
    if ready.is_symlink() or not ready.is_file():
        raise FM0ContractError(f"S{sector} READY is not materialized")
    _require_read_only(ready, label=f"S{sector} READY")
    if ready.read_text(encoding="utf-8").strip() != receipt_sha:
        raise FM0ContractError(f"S{sector} READY/receipt binding drifted")

    manifest_sha = _digest(
        receipt.get("manifest_sha256"), label=f"S{sector} manifest SHA-256"
    )
    if manifest_sha != _digest(
        summary.get("manifest_sha256"), label=f"S{sector} admission manifest SHA-256"
    ):
        raise FM0ContractError(f"S{sector} admission/receipt manifest binding drifted")
    manifest_path, _ = _verify_file(
        root / "manifest.csv",
        manifest_sha,
        label=f"S{sector} six-view manifest",
        require_read_only=True,
    )
    rows = _read_csv(
        manifest_path,
        LATER_SIX_VIEW_MANIFEST_FIELDS,
        label=f"S{sector} six-view manifest",
    )
    assert_no_search_columns(
        LATER_SIX_VIEW_MANIFEST_FIELDS, context=f"S{sector} six-view manifest"
    )
    counts = {partition: 0 for partition in _ALLOWED_RELEASE_PARTITIONS}
    n_cadences = 0
    n_segments = 0
    for row in rows:
        _validate_manifest_row(row, sector=sector)
        observation = row["observation_key"]
        product = row["product_instance_id"]
        if observation in seen_observations or product in seen_products:
            raise FM0ContractError("S66--S77 manifest identities are not unique")
        seen_observations.add(observation)
        seen_products.add(product)
        counts[row["source_partition"]] += 1
        n_cadences += int(row["n_cadences"])
        n_segments += int(row["n_segments"])

    n_rows = len(rows)
    n_observations = _integer(receipt.get("n_observations"), label="n_observations")
    n_shards = _integer(receipt.get("n_shards"), label="n_shards")
    if (
        n_observations != n_rows
        or n_shards != n_rows
        or _integer(receipt.get("n_cadences"), label="n_cadences") != n_cadences
        or _integer(receipt.get("n_segments"), label="n_segments") != n_segments
        or receipt.get("source_rows_by_partition") != counts
        or _integer(summary.get("n_observations"), label="summary n_observations")
        != n_rows
        or _integer(summary.get("n_shards"), label="summary n_shards") != n_rows
        or _integer(summary.get("n_cadences"), label="summary n_cadences") != n_cadences
    ):
        raise FM0ContractError(f"S{sector} six-view count binding drifted")
    for field in ("source_manifest_sha256", "visit_timing_sha256"):
        if _digest(receipt.get(field), label=f"S{sector} receipt {field}") != _digest(
            summary.get(field), label=f"S{sector} admission {field}"
        ):
            raise FM0ContractError(f"S{sector} admission/receipt {field} drifted")

    binding = {
        "binding_schema_version": TEMPORAL_PANEL_SECTOR_BINDING_SCHEMA_VERSION,
        "sector": sector,
        "sector_release_root": str(root),
        "sector_receipt_path": str(receipt_path),
        "sector_receipt_sha256": receipt_sha,
        "sector_manifest_path": str(manifest_path),
        "sector_manifest_sha256": manifest_sha,
        "n_manifest_rows": n_rows,
        "n_development_rows": counts[DEVELOPMENT_PARTITION],
        "n_train_rows_excluded": counts[TRAIN_PARTITION],
        "n_sealed_rows_emitted": 0,
        "scientific_training_eligible": "false",
        "development_evaluation_eligible": "true",
    }
    return rows, binding


def _validate_admission(
    receipt: Mapping[str, Any],
) -> Sequence[Mapping[str, Any]]:
    if (
        receipt.get("schema_version") != POOL_RECEIPT_SCHEMA_VERSION
        or tuple(int(value) for value in receipt.get("preparation_pool_sectors", ()))
        != TEMPORAL_PANEL_SECTORS
        or tuple(int(value) for value in receipt.get("excluded_sectors", ()))
        != EXCLUDED_SECTORS
        or receipt.get("product_state") != DEFERRED_PRODUCT_STATE
        or receipt.get("preparation_pool_admitted") is not True
        or receipt.get("a2v1_accepted") is not False
        or receipt.get("scientific_training_eligible") is not False
        or receipt.get("panel_admission_authorized") is not False
        or receipt.get("temporal_panel_frozen") is not False
        or receipt.get("model_training_authorized") is not False
        or receipt.get("sealed_test_access_authorized") is not False
        or receipt.get("s78_included") is not False
        or receipt.get("s79_s80_touched") is not False
        or receipt.get("sealed_aperture_photometry_opened") is not False
        or receipt.get("sealed_shards_written") is not False
        or receipt.get("sector_wide_hdf5_identity_cadence_quality_qa_performed")
        is not True
    ):
        raise FM0ContractError("preparation-pool admission contract drifted")
    nonsealed = _integer(
        receipt.get("n_nonsealed_preparation_rows"),
        label="n_nonsealed_preparation_rows",
    )
    if (
        _integer(receipt.get("n_six_view_shards"), label="n_six_view_shards")
        != nonsealed
    ):
        raise FM0ContractError("preparation-pool shard/source count drifted")
    _integer(receipt.get("n_source_rows"), label="n_source_rows")
    _integer(
        receipt.get("n_sealed_identity_rows"),
        label="n_sealed_identity_rows",
        allow_zero=True,
    )
    sector_receipts = tuple(
        _mapping(value, label="admission sector receipt")
        for value in _sequence(
            receipt.get("sector_receipts"), label="admission sector receipts"
        )
    )
    if tuple(int(value.get("sector", -1)) for value in sector_receipts) != (
        TEMPORAL_PANEL_SECTORS
    ):
        raise FM0ContractError(
            "admission sector receipts must be chronological S66--S77"
        )
    return sector_receipts


def _panel_row(
    row: Mapping[str, str],
    *,
    cohort: str,
    binding: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "panel_schema_version": TEMPORAL_PANEL_SCHEMA_VERSION,
        "observation_key": row["observation_key"],
        "product_instance_id": row["product_instance_id"],
        "physical_source_id": row["physical_source_id"],
        "gaia_dr3_source_id": row["gaia_dr3_source_id"],
        "tic_id": row["tic_id"],
        "sector": row["sector"],
        "leakage_component_id": row["leakage_component_id"],
        "source_partition": row["source_partition"],
        "temporal_cohort": cohort,
        "sector_release_root": binding["sector_release_root"],
        "sector_receipt_sha256": binding["sector_receipt_sha256"],
        "sector_manifest_sha256": binding["sector_manifest_sha256"],
        "relative_path": row["relative_path"],
        "sha256": row["sha256"],
        "input_source_sha256": row["input_source_sha256"],
        "source_row_sha256": row["source_row_sha256"],
        "n_cadences": row["n_cadences"],
        "n_segments": row["n_segments"],
        "view_present_json": row["view_present_json"],
        "input_adapter": row["input_adapter"],
        "mission_quality_provider": row["mission_quality_provider"],
        "mission_quality_reference_manifest_sha256": row[
            "mission_quality_reference_manifest_sha256"
        ],
        "hdf5_quality_receipt_sha256": row["hdf5_quality_receipt_sha256"],
        "full_visit_shard": "true",
        "model_context_length_bound": "false",
        "scientific_training_eligible": "false",
        "development_evaluation_eligible": "true",
    }


def _make_tree_read_only(root: Path) -> None:
    for path in root.rglob("*"):
        path.chmod(0o555 if path.is_dir() else 0o444)
    root.chmod(0o555)


def _output_leaf_path(path: str | Path) -> Path:
    raw = Path(path).expanduser()
    if raw.name in {"", ".", ".."}:
        raise FM0ContractError("temporal-panel output must be a named directory")
    return raw.absolute()


def freeze_temporal_panel(
    *,
    config_path: str | Path,
    config_sha256: str,
    admission_receipt_path: str | Path,
    admission_receipt_sha256: str,
    baseline_manifest_path: str | Path,
    baseline_manifest_sha256: str,
    producer_git_sha: str,
    output_dir: str | Path,
) -> TemporalPanelFreeze:
    """Freeze an immutable identity-only S66--S77 development panel."""

    producer = str(producer_git_sha).strip().lower()
    if _GIT_SHA.fullmatch(producer) is None:
        raise FM0ContractError("producer_git_sha must be a full lowercase Git SHA")
    config, config_source, config_hash, policy_baseline_hash = _load_policy(
        config_path, expected_sha256=config_sha256
    )
    runtime_baseline_hash = _digest(
        baseline_manifest_sha256, label="runtime baseline manifest SHA-256"
    )
    if runtime_baseline_hash != policy_baseline_hash:
        raise FM0ContractError("runtime baseline manifest does not match frozen policy")
    (
        baseline_source,
        baseline_hash,
        baseline_components,
        n_baseline_rows,
    ) = _load_baseline_components(
        baseline_manifest_path, expected_sha256=runtime_baseline_hash
    )

    admission_source, admission_hash = _verify_file(
        admission_receipt_path,
        admission_receipt_sha256,
        label="S66--S77 admission receipt",
        require_read_only=True,
    )
    admission = _read_json(admission_source, label="S66--S77 admission receipt")
    sector_receipts = _validate_admission(admission)

    seen_observations: set[str] = set()
    seen_products: set[str] = set()
    panel_rows: list[dict[str, Any]] = []
    bindings: list[dict[str, Any]] = []
    for sector, admission_sector in zip(
        TEMPORAL_PANEL_SECTORS, sector_receipts, strict=True
    ):
        manifest_rows, binding = _load_sector_manifest(
            admission_sector=admission_sector,
            sector=sector,
            seen_observations=seen_observations,
            seen_products=seen_products,
        )
        bindings.append(binding)
        for row in manifest_rows:
            if row["source_partition"] == TRAIN_PARTITION:
                continue
            if row["source_partition"] != DEVELOPMENT_PARTITION:
                raise FM0ContractError("non-development row reached panel selection")
            component = row["leakage_component_id"]
            baseline_partition = baseline_components.get(component)
            if (
                baseline_partition is not None
                and baseline_partition != DEVELOPMENT_PARTITION
            ):
                raise FM0ContractError(
                    "later development component crosses its frozen baseline partition"
                )
            cohort = "repeated" if baseline_partition is not None else "new"
            panel_rows.append(_panel_row(row, cohort=cohort, binding=binding))

    if not panel_rows:
        raise FM0ContractError("temporal panel has no development rows")
    panel_rows.sort(key=lambda row: (int(row["sector"]), row["observation_key"]))
    panel_observations = {str(row["observation_key"]) for row in panel_rows}
    if len(panel_observations) != len(panel_rows):
        raise FM0ContractError("temporal panel contains duplicate observations")
    assert_no_search_columns(TEMPORAL_PANEL_FIELDS, context="temporal panel")

    total_manifest_rows = sum(int(row["n_manifest_rows"]) for row in bindings)
    total_train_excluded = sum(int(row["n_train_rows_excluded"]) for row in bindings)
    admission_nonsealed = int(admission["n_nonsealed_preparation_rows"])
    if total_manifest_rows != admission_nonsealed:
        raise FM0ContractError("panel source manifests do not close to admission total")
    if total_train_excluded + len(panel_rows) != total_manifest_rows:
        raise FM0ContractError("panel development/train partition closure failed")

    repeated_rows = [row for row in panel_rows if row["temporal_cohort"] == "repeated"]
    new_rows = [row for row in panel_rows if row["temporal_cohort"] == "new"]
    if not repeated_rows or not new_rows:
        raise FM0ContractError(
            "temporal panel requires nonempty repeated and new development cohorts"
        )
    repeated_components = {str(row["leakage_component_id"]) for row in repeated_rows}
    new_components = {str(row["leakage_component_id"]) for row in new_rows}
    if repeated_components & new_components:
        raise FM0ContractError(
            "one component was classified into both temporal cohorts"
        )

    panel_payload = _csv_bytes(panel_rows, TEMPORAL_PANEL_FIELDS)
    bindings_payload = _csv_bytes(bindings, TEMPORAL_PANEL_SECTOR_BINDING_FIELDS)
    panel_hash = hashlib.sha256(panel_payload).hexdigest()
    bindings_hash = hashlib.sha256(bindings_payload).hexdigest()
    receipt = {
        "schema_version": TEMPORAL_PANEL_RECEIPT_SCHEMA_VERSION,
        "ready_state": TEMPORAL_PANEL_READY_STATE,
        "passed": True,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "producer_git_sha": producer,
        "campaign_id": str(config.get("campaign_id", "")).strip(),
        "config_path": str(config_source),
        "config_sha256": config_hash,
        "admission_receipt_path": str(admission_source),
        "admission_receipt_sha256": admission_hash,
        "baseline_manifest_path": str(baseline_source),
        "baseline_manifest_sha256": baseline_hash,
        "baseline_manifest_sectors": list(BASELINE_SECTORS),
        "panel_sectors": list(TEMPORAL_PANEL_SECTORS),
        "excluded_sectors": list(EXCLUDED_SECTORS),
        "n_baseline_manifest_rows": n_baseline_rows,
        "n_baseline_components": len(baseline_components),
        "n_admission_source_rows": int(admission["n_source_rows"]),
        "n_admission_nonsealed_rows": admission_nonsealed,
        "n_admission_sealed_identity_rows": int(admission["n_sealed_identity_rows"]),
        "n_panel_rows": len(panel_rows),
        "n_panel_components": len(repeated_components | new_components),
        "n_repeated_rows": len(repeated_rows),
        "n_repeated_components": len(repeated_components),
        "n_new_rows": len(new_rows),
        "n_new_components": len(new_components),
        "n_train_rows_excluded": total_train_excluded,
        "n_sealed_rows_emitted": 0,
        "panel_path": "panel.csv",
        "panel_sha256": panel_hash,
        "sector_bindings_path": "sector_bindings.csv",
        "sector_bindings_sha256": bindings_hash,
        "label_blind": True,
        "identity_only": True,
        "shard_payloads_opened": False,
        "temporal_panel_frozen": True,
        "development_evaluation_eligible": True,
        "scientific_training_eligible": False,
        "model_training_authorized": False,
        "sealed_test_access_authorized": False,
        "event_retention_authorized": False,
        "formal_model_gate_authorized": False,
        "production_model_claim": False,
        "prospective_test_claim": False,
        "claim_limit": str(config.get("claim_limit", "")).strip(),
    }
    receipt_payload = (
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    receipt_hash = hashlib.sha256(receipt_payload).hexdigest()

    final = _output_leaf_path(output_dir)
    partial = final.with_name(final.name + ".partial")
    if final.exists() or final.is_symlink() or partial.exists() or partial.is_symlink():
        raise FM0ContractError(f"refusing to overwrite temporal panel: {final}")
    final.parent.mkdir(parents=True, exist_ok=True)
    try:
        partial.mkdir()
        publish_immutable(partial / "panel.csv", panel_payload)
        publish_immutable(partial / "sector_bindings.csv", bindings_payload)
        publish_immutable(partial / "receipt.json", receipt_payload)
        publish_immutable(partial / "READY", receipt_hash.encode("ascii") + b"\n")
        _make_tree_read_only(partial)
        os.replace(partial, final)
    except Exception:
        if partial.exists() and not partial.is_symlink():
            for path in (partial, *partial.rglob("*")):
                try:
                    path.chmod(0o755 if path.is_dir() else 0o644)
                except OSError:
                    pass
            shutil.rmtree(partial, ignore_errors=True)
        raise
    return TemporalPanelFreeze(
        output_dir=final,
        panel_path=final / "panel.csv",
        sector_bindings_path=final / "sector_bindings.csv",
        receipt_path=final / "receipt.json",
        ready_path=final / "READY",
        receipt_sha256=receipt_hash,
        receipt=receipt,
    )


__all__ = [
    "BASELINE_SECTORS",
    "DEVELOPMENT_PARTITION",
    "TEMPORAL_PANEL_CAMPAIGN_ID",
    "TEMPORAL_PANEL_CONFIG_SCHEMA_VERSION",
    "TEMPORAL_PANEL_FIELDS",
    "TEMPORAL_PANEL_FROZEN_POLICY_SHA256",
    "TEMPORAL_PANEL_READY_STATE",
    "TEMPORAL_PANEL_RECEIPT_SCHEMA_VERSION",
    "TEMPORAL_PANEL_SCHEMA_VERSION",
    "TEMPORAL_PANEL_SECTORS",
    "TEMPORAL_PANEL_SECTOR_BINDING_FIELDS",
    "TEMPORAL_PANEL_SECTOR_BINDING_SCHEMA_VERSION",
    "TemporalPanelFreeze",
    "freeze_temporal_panel",
]
