"""Fail-closed admission of the deferred S66--S77 FM preparation pool.

This is deliberately separate from :mod:`temporal_admission`, whose frozen
v1 policy requires accepted A2v1 products beginning at Sector 65.  The v2
policy proves only that every deferred ORCD sector in the exact S66--S77 pool
has four mutually bound preparation receipts.  It cannot convert those
checkpoints into accepted Stage-1 products, freeze a temporal panel, or make
their light curves scientifically training eligible.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path
from typing import Any

import yaml

from twirl.lightcurves.mission_quality_reference import (
    MISSION_QUALITY_REFERENCE_SCHEMA_VERSION,
    MissionQualityReference,
    load_mission_quality_reference,
)

from .hdf5_quality_admission import (
    HDF5_QUALITY_READY_STATE,
    HDF5_QUALITY_RECEIPT_SCHEMA_VERSION,
)
from .input_release import FLUX_VIEW_NAMES
from .later_hdf5_adapter import later_input_source_sha256
from .later_sector_release import (
    _load_hdf5_quality_receipt,
    _load_inventory_bundle,
    validate_later_sector_release,
)
from .later_source_inventory import (
    SUMMARY_SCHEMA_VERSION,
    load_later_source_rows,
)
from .registry import FM0ContractError, publish_immutable, sha256_file

CONFIG_SCHEMA_VERSION = "twirl_fm0_2_later_sector_admission_config_v2"
PRODUCTION_CAMPAIGN_ID = "twirl_fm0_2_later_sector_admission_v2"
# Filled from the separately reviewed immutable config.  It is intentionally
# independent of the frozen v1 hash in temporal_admission.py.
FROZEN_POLICY_SHA256 = (
    "10536c89d2656604c97f6289cb03fc9f5aca459aa5340303324314e4905dfd76"
)
EXCLUSION_LEDGER_SCHEMA_VERSION = "twirl_fm0_later_sector_exclusion_ledger_v1"
EXCLUSION_DECISION_ID = "twirl_fm0_later_sector_exclusions_v1"
POOL_RECEIPT_SCHEMA_VERSION = "twirl_fm0_later_sector_preparation_pool_receipt_v2"
SIX_VIEW_RECEIPT_SCHEMA_VERSION = "twirl_fm0_later_six_view_sector_release_v1"
SIX_VIEW_READY_STATE = "FM_SIX_VIEW_DEFERRED_READY"
DEFERRED_PRODUCT_STATE = "ORCD_COMPLETE_DEFERRED"
PREPARATION_SECTORS = tuple(range(66, 78))
EXCLUDED_SECTORS = (65,)
RECEIPT_ROLES = (
    "mission_quality_reference",
    "hdf5_quality",
    "source_inventory",
    "six_view_release",
)

_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_GIT_SHA = re.compile(r"^[0-9a-f]{40}$")
_LEAKAGE_COMPONENT = re.compile(r"^leakage_[0-9a-f]{64}$")
_ALLOWED_PARTITIONS = frozenset({"poc_train", "poc_development"})
_SEALED_PARTITIONS = frozenset({"poc_sealed_test"})


@dataclass(frozen=True)
class ReceiptBinding:
    """One exact runtime path/digest binding."""

    path: Path
    sha256: str


@dataclass(frozen=True)
class SectorPreparationBindings:
    """The four exact preparation receipts required for one sector."""

    sector: int
    mission_quality_reference: ReceiptBinding
    hdf5_quality: ReceiptBinding
    source_inventory: ReceiptBinding
    six_view_release: ReceiptBinding


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


def _positive_int(value: Any, *, label: str, allow_zero: bool = False) -> int:
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


def _read_json(path: Path, *, label: str) -> Mapping[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FM0ContractError(f"invalid {label}: {path}") from exc
    return _mapping(payload, label=label)


def _verify_file(binding: ReceiptBinding, *, label: str) -> Path:
    expected = _digest(binding.sha256, label=f"{label} SHA-256")
    source = Path(binding.path).expanduser().resolve()
    if binding.path.is_symlink() or not source.is_file() or source.stat().st_size <= 0:
        raise FM0ContractError(f"missing materialized {label}: {source}")
    observed = sha256_file(source)
    if observed != expected:
        raise FM0ContractError(
            f"{label} hash mismatch: expected {expected}, observed {observed}"
        )
    return source


def _verify_declared_file(
    *, declaring_path: Path, raw_path: Any, raw_sha256: Any, label: str
) -> tuple[Path, str]:
    text = str(raw_path).strip()
    if not text:
        raise FM0ContractError(f"{label} path is missing")
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = declaring_path.parent / path
    digest = _digest(raw_sha256, label=f"{label} SHA-256")
    return _verify_file(ReceiptBinding(path=path, sha256=digest), label=label), digest


def _expected_orbits(sector: int) -> tuple[int, int]:
    return 2 * sector + 7, 2 * sector + 8


def _load_policy(
    path: str | Path, *, expected_sha256: str
) -> tuple[Mapping[str, Any], Path, str]:
    source = _verify_file(
        ReceiptBinding(Path(path), expected_sha256), label="later-sector v2 config"
    )
    try:
        payload = yaml.safe_load(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        raise FM0ContractError(f"invalid later-sector v2 config: {source}") from exc
    config = _mapping(payload, label="later-sector v2 config")
    observed = _digest(expected_sha256, label="later-sector v2 config SHA-256")
    if config.get("schema_version") != CONFIG_SCHEMA_VERSION:
        raise FM0ContractError("later-sector v2 config schema mismatch")
    if (
        config.get("campaign_id") == PRODUCTION_CAMPAIGN_ID
        and observed != FROZEN_POLICY_SHA256
    ):
        raise FM0ContractError("production later-sector v2 policy hash drifted")
    _validate_policy(config)
    return config, source, observed


def _validate_policy(config: Mapping[str, Any]) -> None:
    scope = _mapping(config.get("scope"), label="v2 scope")
    exclusion = _mapping(config.get("exclusion_ledger"), label="v2 exclusion")
    selection = _mapping(config.get("selection"), label="v2 selection")
    evidence = _mapping(config.get("evidence"), label="v2 evidence")
    boundary = _mapping(config.get("product_boundary"), label="v2 product boundary")

    if scope.get("policy_only") is not True or scope.get(
        "preparation_pool_admission_only"
    ) is not True:
        raise FM0ContractError("v2 is not restricted to preparation admission")
    false_scope = (
        "results_recorded",
        "panel_frozen",
        "model_training_authorized",
        "sealed_test_access_authorized",
        "event_retention_authorized",
        "formal_model_gate_authorized",
        "production_model_claim",
        "foundation_model_claim",
        "prospective_test_claim",
    )
    if any(scope.get(field) is not False for field in false_scope):
        raise FM0ContractError("v2 scope improperly authorizes a model or panel claim")

    if (
        exclusion.get("required_schema_version")
        != EXCLUSION_LEDGER_SCHEMA_VERSION
        or exclusion.get("required_decision_id") != EXCLUSION_DECISION_ID
        or tuple(int(value) for value in exclusion.get("excluded_sectors", ()))
        != EXCLUDED_SECTORS
        or exclusion.get("exclusion_must_precede_model_evaluation") is not True
        or exclusion.get("repaired_sector_retroactive_insertion") != "forbidden"
    ):
        raise FM0ContractError("v2 exclusion-ledger policy drifted")
    _digest(exclusion.get("required_sha256"), label="exclusion-ledger policy SHA-256")

    sectors = tuple(int(value) for value in selection.get("ordered_sectors", ()))
    if sectors != PREPARATION_SECTORS:
        raise FM0ContractError("v2 preparation pool must be exactly S66--S77")
    if any(right != left + 1 for left, right in pairwise(sectors)):
        raise FM0ContractError("v2 preparation pool is not contiguous")
    if (
        int(selection.get("first_sector", -1)) != sectors[0]
        or int(selection.get("final_sector", -1)) != sectors[-1]
        or selection.get("exact_sector_set_required") is not True
        or selection.get("chronological_order_required") is not True
        or selection.get("require_contiguous_numeric_sectors") is not True
        or selection.get("skipped_sector_policy") != "fail"
        or selection.get("excluded_sector_policy") != "fail"
        or selection.get("unmapped_sector_policy") != "fail"
        or selection.get("expected_orbits_must_be_explicit") is not True
        or selection.get("s78_included") is not False
        or selection.get("s79_s80_touched") is not False
    ):
        raise FM0ContractError("v2 exact-sector selection policy drifted")
    orbit_map = _mapping(
        selection.get("expected_orbits_by_sector"), label="v2 orbit map"
    )
    try:
        normalized_orbits = {
            int(raw_sector): tuple(int(value) for value in raw_orbits)
            for raw_sector, raw_orbits in orbit_map.items()
        }
    except (TypeError, ValueError) as exc:
        raise FM0ContractError("v2 orbit map is invalid") from exc
    if normalized_orbits != {
        sector: _expected_orbits(sector) for sector in PREPARATION_SECTORS
    }:
        raise FM0ContractError("v2 explicit orbit map drifted")

    required_receipts = _mapping(
        evidence.get("required_receipts_per_sector"), label="v2 receipt policy"
    )
    if tuple(required_receipts) != RECEIPT_ROLES:
        raise FM0ContractError("v2 required receipt roles or order drifted")
    expected_receipts = {
        "mission_quality_reference": (
            "schema_version",
            MISSION_QUALITY_REFERENCE_SCHEMA_VERSION,
        ),
        "hdf5_quality": ("schema_version", HDF5_QUALITY_RECEIPT_SCHEMA_VERSION),
        "source_inventory": ("schema_version", SUMMARY_SCHEMA_VERSION),
        "six_view_release": ("schema_version", SIX_VIEW_RECEIPT_SCHEMA_VERSION),
    }
    for role, (field, expected) in expected_receipts.items():
        if required_receipts[role].get(field) != expected:
            raise FM0ContractError(f"v2 {role} schema policy drifted")
    if (
        required_receipts["hdf5_quality"].get("quality_state")
        != HDF5_QUALITY_READY_STATE
        or required_receipts["six_view_release"].get("release_state")
        != SIX_VIEW_READY_STATE
        or tuple(required_receipts["six_view_release"].get("flux_view_names", ()))
        != tuple(FLUX_VIEW_NAMES)
    ):
        raise FM0ContractError("v2 receipt state/view policy drifted")
    true_evidence = (
        "fail_closed",
        "exact_path_and_sha256_supplied_at_runtime",
        "distinct_receipt_paths_and_sha256_required",
        "every_receipt_must_pass",
        "cross_receipt_sha256_bindings_required",
        "model_outcome_blind_required",
        "sealed_aperture_photometry_or_derived_shard_access_forbidden",
        "label_blind_sector_wide_hdf5_quality_qa_allowed",
    )
    if any(evidence.get(field) is not True for field in true_evidence) or evidence.get(
        "missing_or_unbound_receipt_policy"
    ) != "fail":
        raise FM0ContractError("v2 evidence policy is not fail-closed")

    if (
        boundary.get("required_source_product_state") != DEFERRED_PRODUCT_STATE
        or boundary.get("preparation_pool_state") != DEFERRED_PRODUCT_STATE
        or boundary.get("preparation_pool_admitted") is not True
    ):
        raise FM0ContractError("v2 deferred preparation state drifted")
    false_boundary = (
        "a2v1_accepted",
        "scientific_training_eligible",
        "panel_admission_authorized",
        "temporal_panel_frozen",
        "model_training_authorized",
        "deferred_checkpoint_is_stage1_acceptance",
        "deferred_checkpoint_is_scientific_release",
    )
    if any(boundary.get(field) is not False for field in false_boundary) or boundary.get(
        "separate_future_acceptance_contract_required"
    ) is not True:
        raise FM0ContractError("v2 product/training boundary drifted")


def _verify_exclusion_ledger(
    *,
    config: Mapping[str, Any],
    ledger_path: str | Path,
    expected_sha256: str,
) -> tuple[Path, str]:
    policy = _mapping(config.get("exclusion_ledger"), label="exclusion policy")
    digest = _digest(expected_sha256, label="runtime exclusion-ledger SHA-256")
    if digest != policy.get("required_sha256"):
        raise FM0ContractError("runtime exclusion-ledger hash differs from v2 policy")
    path = _verify_file(
        ReceiptBinding(Path(ledger_path), digest), label="S65 exclusion ledger"
    )
    try:
        payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        raise FM0ContractError(f"invalid S65 exclusion ledger: {path}") from exc
    ledger = _mapping(payload, label="S65 exclusion ledger")
    scope = _mapping(ledger.get("scope"), label="S65 exclusion scope")
    rows = _sequence(ledger.get("excluded_sectors"), label="excluded sectors")
    if (
        ledger.get("schema_version") != EXCLUSION_LEDGER_SCHEMA_VERSION
        or ledger.get("decision_id") != EXCLUSION_DECISION_ID
        or ledger.get("status") != "frozen_pre_model_evaluation"
        or len(rows) != 1
        or scope.get("selection_is_label_blind") is not True
        or scope.get("selection_precedes_later_sector_embedding_or_model_score_access")
        is not True
        or scope.get("model_training_authorized") is not False
        or scope.get("temporal_panel_frozen") is not False
        or scope.get("sealed_test_access_authorized") is not False
        or scope.get("stage1_product_state_changed") is not False
    ):
        raise FM0ContractError("S65 exclusion ledger identity/scope drifted")
    row = _mapping(rows[0], label="S65 exclusion row")
    if (
        int(row.get("sector", -1)) != 65
        or row.get("decision") != "exclude"
        or row.get("reason_code") != "incomplete_mission_quality_authority"
        or row.get("mission_quality_provider") != "spoc"
        or int(row.get("missing_cadence_count", -1)) != 1462
        or row.get("exclusion_is_model_outcome_based") is not False
    ):
        raise FM0ContractError("S65 exclusion decision drifted")
    return path, digest


def _verify_mission_quality(
    *, sector: int, binding: ReceiptBinding
) -> tuple[Path, MissionQualityReference, dict[str, Any]]:
    path = _verify_file(binding, label=f"S{sector} mission-quality reference")
    if path.name != "manifest.json":
        raise FM0ContractError(
            f"S{sector} mission-quality binding must name manifest.json"
        )
    orbits = _expected_orbits(sector)
    provider = "spoc" if sector < 67 else "tica"
    reference = load_mission_quality_reference(
        reference_dir=path.parent,
        sector=sector,
        expected_orbits=orbits,
    )
    if (
        reference.sector != sector
        or reference.expected_orbits != orbits
        or reference.provider != provider
        or reference.manifest_path.resolve() != path
        or reference.manifest_sha256
        != _digest(binding.sha256, label=f"S{sector} mission-quality SHA-256")
    ):
        raise FM0ContractError(
            f"S{sector} fully loaded mission-quality reference drifted"
        )
    reference.assert_unchanged()
    return path, reference, {
        "schema_version": MISSION_QUALITY_REFERENCE_SCHEMA_VERSION,
        "sha256": binding.sha256,
        "mission_quality_provider": provider,
        "table_sha256": reference.table_sha256,
        "path": str(path),
    }


def _verify_hdf5_quality(
    *,
    sector: int,
    binding: ReceiptBinding,
    inventory_summary: Mapping[str, Any],
    reference: MissionQualityReference,
) -> tuple[Path, Mapping[str, Any], dict[str, Any]]:
    path = _verify_file(binding, label=f"S{sector} HDF5-quality receipt")
    receipt, loaded_path, loaded_sha = _load_hdf5_quality_receipt(
        path=path,
        sector=sector,
        inventory_summary=inventory_summary,
        reference=reference,
        expected_receipt_sha256=binding.sha256,
    )
    if (
        loaded_path != path
        or loaded_sha
        != _digest(binding.sha256, label=f"S{sector} HDF5-quality SHA-256")
    ):
        raise FM0ContractError(f"S{sector} HDF5-quality exact binding drifted")
    products = _positive_int(
        receipt.get("n_hdf5_products"), label=f"S{sector} HDF5 product count"
    )
    cadences = _positive_int(
        receipt.get("n_cadences_checked"), label=f"S{sector} HDF5 cadence count"
    )
    reference.assert_unchanged()
    return path, receipt, {
        "schema_version": HDF5_QUALITY_RECEIPT_SCHEMA_VERSION,
        "quality_state": HDF5_QUALITY_READY_STATE,
        "sha256": loaded_sha,
        "n_hdf5_products": products,
        "n_cadences_checked": cadences,
        "path": str(path),
    }


def _verify_source_inventory(
    *, sector: int, binding: ReceiptBinding
) -> tuple[
    Path,
    Mapping[str, Any],
    tuple[dict[str, str], ...],
    dict[str, Any],
]:
    path = _verify_file(binding, label=f"S{sector} source-inventory summary")
    if path.name != "summary.json":
        raise FM0ContractError(
            f"S{sector} source-inventory binding must name summary.json"
        )
    (
        summary,
        loaded_summary_path,
        loaded_summary_sha,
        all_rows,
        sources_path,
        sources_sha,
    ) = _load_inventory_bundle(
        source_inventory_dir=path.parent,
        sector=sector,
        expected_summary_sha256=binding.sha256,
    )
    if (
        loaded_summary_path != path
        or loaded_summary_sha
        != _digest(binding.sha256, label=f"S{sector} source-inventory SHA-256")
    ):
        raise FM0ContractError(f"S{sector} source-inventory exact binding drifted")
    selected_rows = load_later_source_rows(
        path.parent,
        expected_summary_sha256=loaded_summary_sha,
        allowed_source_partitions=tuple(sorted(_ALLOWED_PARTITIONS)),
    )
    expected_selected = tuple(
        row for row in all_rows if row["source_partition"] in _ALLOWED_PARTITIONS
    )
    if not selected_rows:
        raise FM0ContractError(f"S{sector} has no nonsealed preparation rows")
    if selected_rows != expected_selected:
        raise FM0ContractError(
            f"S{sector} strict source-row identity/evidence closure drifted"
        )
    n_source_rows = len(all_rows)
    n_sealed = n_source_rows - len(selected_rows)
    if n_sealed != sum(
        int(summary["source_rows_by_partition"].get(partition, 0))
        for partition in _SEALED_PARTITIONS
    ):
        raise FM0ContractError(f"S{sector} sealed source-row closure drifted")
    return path, summary, selected_rows, {
        "schema_version": SUMMARY_SCHEMA_VERSION,
        "sha256": loaded_summary_sha,
        "sources_sha256": sources_sha,
        "n_source_rows": n_source_rows,
        "n_usable_source_rows": len(selected_rows),
        "n_sealed_identity_rows": n_sealed,
        "path": str(path),
        "sources_path": str(sources_path),
    }


def _binding_matches(
    *, receipt_path: Path, receipt: Mapping[str, Any], prefix: str, expected: ReceiptBinding
) -> bool:
    text = str(receipt.get(f"{prefix}_path", "")).strip()
    if not text:
        return False
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = receipt_path.parent / path
    return path.resolve() == expected.path.expanduser().resolve() and _digest(
        receipt.get(f"{prefix}_sha256"), label=f"six-view {prefix} SHA-256"
    ) == _digest(expected.sha256, label=f"expected {prefix} SHA-256")


def _stable_id(prefix: str, payload: Any) -> str:
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(encoded).hexdigest()}"


def _verify_six_view(
    *,
    sector: int,
    binding: ReceiptBinding,
    mission_quality: ReceiptBinding,
    hdf5_quality: ReceiptBinding,
    source_inventory: ReceiptBinding,
    expected_usable_rows: int,
    reference: MissionQualityReference,
    inventory_summary: Mapping[str, Any],
    selected_source_rows: Sequence[Mapping[str, str]],
) -> tuple[Path, Mapping[str, Any], dict[str, Any]]:
    path = _verify_file(binding, label=f"S{sector} six-view release receipt")
    loaded_path, receipt, bundle = validate_later_sector_release(
        path.parent,
        expected_receipt_sha256=binding.sha256,
        require_read_only=True,
        verify_shard_payloads=False,
    )
    observations = int(bundle["n_observations"])
    shards = int(bundle["n_shards"])
    if (
        loaded_path != path
        or int(bundle["sector"]) != sector
        or bundle["mission_quality_provider"] != reference.provider
        or observations != expected_usable_rows
        or shards != observations
    ):
        raise FM0ContractError(f"S{sector} six-view validated summary drifted")
    expected_bindings = (
        ("mission_quality_reference_manifest", mission_quality),
        ("hdf5_quality_receipt", hdf5_quality),
        ("source_inventory_summary", source_inventory),
    )
    if any(
        not _binding_matches(
            receipt_path=path, receipt=receipt, prefix=prefix, expected=expected
        )
        for prefix, expected in expected_bindings
    ):
        raise FM0ContractError(f"S{sector} six-view cross-receipt binding drifted")

    inventory_root = Path(source_inventory.path).expanduser().resolve().parent
    inventory_sources = inventory_root / "sources.csv"
    if (
        not _binding_matches(
            receipt_path=path,
            receipt=receipt,
            prefix="source_inventory_sources",
            expected=ReceiptBinding(
                inventory_sources,
                _digest(
                    inventory_summary.get("sources_sha256"),
                    label=f"S{sector} inventory sources SHA-256",
                ),
            ),
        )
        or not _binding_matches(
            receipt_path=path,
            receipt=receipt,
            prefix="corpus_selection",
            expected=ReceiptBinding(
                Path(str(inventory_summary.get("identity_binding_path", ""))),
                _digest(
                    inventory_summary.get("identity_binding_sha256"),
                    label=f"S{sector} corpus selection SHA-256",
                ),
            ),
        )
        or not _binding_matches(
            receipt_path=path,
            receipt=receipt,
            prefix="mission_quality_reference_table",
            expected=ReceiptBinding(reference.table_path, reference.table_sha256),
        )
    ):
        raise FM0ContractError(f"S{sector} six-view upstream bundle binding drifted")
    manifest_rows = list(bundle["manifest_rows"])
    source_rows = list(bundle["source_rows"])
    expected_source_rows = sorted(
        (dict(row) for row in selected_source_rows),
        key=lambda row: (int(row["gaia_dr3_source_id"]), int(row["tic_id"])),
    )
    if (
        source_rows != expected_source_rows
        or _positive_int(
            receipt.get("n_sealed_source_rows_excluded"),
            label="six-view sealed source count",
            allow_zero=True,
        )
        != int(inventory_summary.get("n_source_rows", -1)) - observations
    ):
        raise FM0ContractError(f"S{sector} six-view source closure drifted")

    source_by_hash = {row["source_row_sha256"]: row for row in source_rows}
    for row in manifest_rows:
        observation_key = str(row["observation_key"])
        source = source_by_hash.get(row["source_row_sha256"])
        if source is None:
            raise FM0ContractError(f"S{sector} six-view manifest/source closure drifted")
        physical_source_id = f"gaia_dr3:{source['gaia_dr3_source_id']}"
        expected_observation_key = _stable_id(
            "observation",
            {"physical_source_id": physical_source_id, "sector": sector},
        )
        expected_input_source_sha = later_input_source_sha256(
            source_row_sha256=source["source_row_sha256"],
            mission_quality_reference_sha256=reference.manifest_sha256,
            hdf5_quality_receipt_sha256=hdf5_quality.sha256,
        )
        expected_product_id = _stable_id(
            "product",
            {
                "observation_key": observation_key,
                "tic_id": source["tic_id"],
                "a2v1_product_version": "A2v1_ORCD_DEFERRED",
                "source_sha256": expected_input_source_sha,
            },
        )
        if (
            observation_key != expected_observation_key
            or row["product_instance_id"] != expected_product_id
            or row["input_source_sha256"] != expected_input_source_sha
        ):
            raise FM0ContractError(f"S{sector} six-view derived identity drifted")
    reference.assert_unchanged()
    return path, receipt, {
        "schema_version": SIX_VIEW_RECEIPT_SCHEMA_VERSION,
        "release_state": SIX_VIEW_READY_STATE,
        "sha256": binding.sha256,
        "n_observations": observations,
        "n_shards": shards,
        "n_cadences": int(bundle["n_cadences"]),
        "manifest_sha256": receipt["manifest_sha256"],
        "source_manifest_sha256": receipt["source_manifest_sha256"],
        "visit_timing_sha256": receipt["visit_timing_sha256"],
        "path": str(path),
    }


def construct_preparation_pool_receipt(
    *,
    config_path: str | Path,
    expected_config_sha256: str,
    exclusion_ledger_path: str | Path,
    expected_exclusion_ledger_sha256: str,
    ordered_sector_bindings: Sequence[SectorPreparationBindings],
    producer_git_sha: str,
) -> dict[str, Any]:
    """Verify and admit exactly the deferred S66--S77 preparation pool."""

    if _GIT_SHA.fullmatch(str(producer_git_sha)) is None:
        raise FM0ContractError("producer_git_sha must be a full lowercase Git SHA")
    config, policy_path, policy_sha = _load_policy(
        config_path, expected_sha256=expected_config_sha256
    )
    ledger_path, ledger_sha = _verify_exclusion_ledger(
        config=config,
        ledger_path=exclusion_ledger_path,
        expected_sha256=expected_exclusion_ledger_sha256,
    )
    bindings = tuple(ordered_sector_bindings)
    sectors = tuple(int(value.sector) for value in bindings)
    if sectors != PREPARATION_SECTORS:
        raise FM0ContractError(
            "runtime sector receipts must cover exactly chronological S66--S77"
        )

    seen_paths: set[Path] = set()
    seen_hashes: set[str] = set()
    for item in bindings:
        for role in RECEIPT_ROLES:
            binding = getattr(item, role)
            path = Path(binding.path).expanduser().resolve()
            digest = _digest(binding.sha256, label=f"S{item.sector} {role} SHA-256")
            if path in seen_paths or digest in seen_hashes:
                raise FM0ContractError("preparation receipt paths/hashes must be distinct")
            seen_paths.add(path)
            seen_hashes.add(digest)

    sector_receipts: list[dict[str, Any]] = []
    total_source_rows = 0
    total_usable_rows = 0
    total_sealed_rows = 0
    total_shards = 0
    for item in bindings:
        sector = int(item.sector)
        _mission_path, reference, mission_summary = _verify_mission_quality(
            sector=sector, binding=item.mission_quality_reference
        )
        (
            _source_path,
            inventory_summary,
            selected_source_rows,
            source_summary,
        ) = _verify_source_inventory(
            sector=sector, binding=item.source_inventory
        )
        _hdf5_path, _hdf5, hdf5_summary = _verify_hdf5_quality(
            sector=sector,
            binding=item.hdf5_quality,
            inventory_summary=inventory_summary,
            reference=reference,
        )
        _six_path, _six, six_summary = _verify_six_view(
            sector=sector,
            binding=item.six_view_release,
            mission_quality=item.mission_quality_reference,
            hdf5_quality=item.hdf5_quality,
            source_inventory=item.source_inventory,
            expected_usable_rows=int(source_summary["n_usable_source_rows"]),
            reference=reference,
            inventory_summary=inventory_summary,
            selected_source_rows=selected_source_rows,
        )
        total_source_rows += int(source_summary["n_source_rows"])
        total_usable_rows += int(source_summary["n_usable_source_rows"])
        total_sealed_rows += int(source_summary["n_sealed_identity_rows"])
        total_shards += int(six_summary["n_shards"])
        sector_receipts.append(
            {
                "sector": sector,
                "expected_orbits": list(_expected_orbits(sector)),
                "product_state": DEFERRED_PRODUCT_STATE,
                "preparation_admitted": True,
                "a2v1_accepted": False,
                "scientific_training_eligible": False,
                "evidence": {
                    "mission_quality_reference": mission_summary,
                    "hdf5_quality": hdf5_summary,
                    "source_inventory": source_summary,
                    "six_view_release": six_summary,
                },
            }
        )

    if total_shards != total_usable_rows:
        raise FM0ContractError("pool shard total differs from nonsealed source rows")
    return {
        "schema_version": POOL_RECEIPT_SCHEMA_VERSION,
        "campaign_id": PRODUCTION_CAMPAIGN_ID,
        "producer_git_sha": producer_git_sha,
        "policy_path": str(policy_path),
        "policy_sha256": policy_sha,
        "exclusion_ledger_path": str(ledger_path),
        "exclusion_ledger_sha256": ledger_sha,
        "excluded_sectors": list(EXCLUDED_SECTORS),
        "preparation_pool_sectors": list(PREPARATION_SECTORS),
        "product_state": DEFERRED_PRODUCT_STATE,
        "preparation_pool_admitted": True,
        "a2v1_accepted": False,
        "scientific_training_eligible": False,
        "panel_admission_authorized": False,
        "temporal_panel_frozen": False,
        "model_training_authorized": False,
        "sealed_test_access_authorized": False,
        "formal_model_gate_authorized": False,
        "production_model_claim": False,
        "foundation_model_claim": False,
        "prospective_test_claim": False,
        "s78_included": False,
        "s79_s80_touched": False,
        "n_source_rows": total_source_rows,
        "n_nonsealed_preparation_rows": total_usable_rows,
        "n_sealed_identity_rows": total_sealed_rows,
        "n_six_view_shards": total_shards,
        "sector_wide_hdf5_identity_cadence_quality_qa_performed": True,
        "sealed_aperture_photometry_opened": False,
        "sealed_shards_written": False,
        "sector_receipts": sector_receipts,
        "claim_limit": str(config.get("claim_limit", "")).strip(),
    }


def write_preparation_pool_receipt(
    *, output_path: str | Path, receipt: Mapping[str, Any]
) -> tuple[Path, str]:
    """Publish one immutable pool receipt and return its exact digest."""

    if (
        receipt.get("schema_version") != POOL_RECEIPT_SCHEMA_VERSION
        or receipt.get("preparation_pool_admitted") is not True
        or receipt.get("a2v1_accepted") is not False
        or receipt.get("scientific_training_eligible") is not False
        or receipt.get("model_training_authorized") is not False
    ):
        raise FM0ContractError("refusing to publish an invalid v2 pool receipt")
    payload = (
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    path = Path(output_path).expanduser().resolve()
    publish_immutable(path, payload)
    return path, hashlib.sha256(payload).hexdigest()


def admit_preparation_pool(
    *,
    config_path: str | Path,
    expected_config_sha256: str,
    exclusion_ledger_path: str | Path,
    expected_exclusion_ledger_sha256: str,
    ordered_sector_bindings: Sequence[SectorPreparationBindings],
    producer_git_sha: str,
    output_path: str | Path,
) -> tuple[dict[str, Any], str]:
    """Verify S66--S77 and publish the deferred preparation-pool receipt."""

    receipt = construct_preparation_pool_receipt(
        config_path=config_path,
        expected_config_sha256=expected_config_sha256,
        exclusion_ledger_path=exclusion_ledger_path,
        expected_exclusion_ledger_sha256=expected_exclusion_ledger_sha256,
        ordered_sector_bindings=ordered_sector_bindings,
        producer_git_sha=producer_git_sha,
    )
    _path, digest = write_preparation_pool_receipt(
        output_path=output_path, receipt=receipt
    )
    return receipt, digest


__all__ = [
    "CONFIG_SCHEMA_VERSION",
    "DEFERRED_PRODUCT_STATE",
    "EXCLUDED_SECTORS",
    "FROZEN_POLICY_SHA256",
    "POOL_RECEIPT_SCHEMA_VERSION",
    "PREPARATION_SECTORS",
    "PRODUCTION_CAMPAIGN_ID",
    "RECEIPT_ROLES",
    "SIX_VIEW_READY_STATE",
    "SIX_VIEW_RECEIPT_SCHEMA_VERSION",
    "ReceiptBinding",
    "SectorPreparationBindings",
    "admit_preparation_pool",
    "construct_preparation_pool_receipt",
    "write_preparation_pool_receipt",
]
