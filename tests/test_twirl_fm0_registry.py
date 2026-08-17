from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from twirl.models.fm0.registry import (
    FM0ContractError,
    DIAGNOSTIC_ADMISSION_RECEIPT_SCHEMA_VERSION,
    DIAGNOSTIC_GATES,
    SPLIT_SALT,
    assert_no_search_columns,
    build_alias_registry,
    build_observation_registry,
    deterministic_source_partition,
    load_frozen_contract,
    write_registry_release,
)


def _diagnostic_receipt(
    path: Path,
    *,
    sector: int = 66,
    passed: bool = True,
    failed_gate: str | None = None,
    schema: str = DIAGNOSTIC_ADMISSION_RECEIPT_SCHEMA_VERSION,
) -> tuple[str, str]:
    gates = {}
    for gate in DIAGNOSTIC_GATES:
        evidence_path = path.parent / f"{path.stem}-{gate}.json"
        gate_passed = gate != failed_gate
        evidence_path.write_text(
            json.dumps(
                {
                    "passed": gate_passed,
                    "sector": sector,
                    "schema_version": "twirl_fm0_1_diagnostic_evidence_v1",
                },
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        gates[gate] = {
            "passed": gate_passed,
            "path": evidence_path.name,
            "sha256": hashlib.sha256(evidence_path.read_bytes()).hexdigest(),
        }
    payload = {
        "receipt_schema_version": schema,
        "sector": sector,
        "passed": passed,
        "evidence_gates": gates,
    }
    path.write_text(json.dumps(payload, sort_keys=True) + "\n", encoding="utf-8")
    return str(path), hashlib.sha256(path.read_bytes()).hexdigest()


def test_frozen_contract_hashes_are_bound() -> None:
    contract = load_frozen_contract()
    assert contract.config["campaign_id"] == "twirl_fm0_1_s56_s67_poc"
    assert contract.freeze_receipt["scientific_contract_status"] == "frozen"


def test_connected_alias_components_are_deterministic_and_quarantined() -> None:
    rows = [
        {"gaia_dr3_source_id": "100", "tic_id": "11"},
        {"gaia_dr3_source_id": "100", "tic_id": "12"},
        {"gaia_dr3_source_id": "200", "tic_id": "21"},
        {"gaia_dr3_source_id": "300", "tic_id": "21"},
        {"gaia_dr3_source_id": "400", "tic_id": ""},
    ]
    registry = build_alias_registry(rows)
    reversed_registry = build_alias_registry(reversed(rows))
    assert registry == reversed_registry
    alias = {(row["gaia_dr3_source_id"], row["tic_id"]): row for row in registry.aliases}
    assert alias[("100", "11")]["leakage_component_id"] == alias[("100", "12")][
        "leakage_component_id"
    ]
    assert not alias[("100", "11")]["quarantined"]
    assert alias[("200", "21")]["quarantined"]
    assert alias[("300", "21")]["source_partition"] == "quarantine"
    assert len(registry.components) == 3


def test_source_split_uses_full_digest_big_endian_modulo() -> None:
    component = "leakage_example"
    expected = int.from_bytes(
        hashlib.sha256(f"{SPLIT_SALT}:{component}".encode()).digest(), "big"
    ) % 10000
    partition, bucket = deterministic_source_partition(component)
    assert bucket == expected
    assert partition == (
        "poc_train"
        if expected < 8000
        else "poc_development"
        if expected < 9000
        else "poc_sealed_test"
    )


def test_observation_keys_and_diagnostic_admission_gates(tmp_path: Path) -> None:
    aliases = build_alias_registry(
        [{"gaia_dr3_source_id": "123456789012345678", "tic_id": "42"}]
    )
    accepted = {
        "gaia_dr3_source_id": "123456789012345678",
        "tic_id": "42",
        "sector": 56,
        "a2v1_product_version": "A2v1",
        "source_sha256": "a" * 64,
        "product_state": "A2V1_ACCEPTED",
    }
    first = build_observation_registry([accepted], aliases)
    second = build_observation_registry([dict(accepted)], aliases)
    assert first == second
    assert first[0]["observation_key"].startswith("observation_")
    assert first[0]["product_instance_id"].startswith("product_")

    diagnostic = dict(accepted, sector=66, product_state="ORCD_COMPLETE_DEFERRED")
    with pytest.raises(FM0ContractError, match="admission_receipt_path"):
        build_observation_registry([diagnostic], aliases)
    diagnostic["immutable_cell_receipts_passed"] = True
    with pytest.raises(FM0ContractError, match="inline diagnostic gate"):
        build_observation_registry([diagnostic], aliases)
    diagnostic.pop("immutable_cell_receipts_passed")

    receipt_path, receipt_hash = _diagnostic_receipt(tmp_path / "admission.json")
    diagnostic["diagnostic_admission_receipt_path"] = receipt_path
    diagnostic["diagnostic_admission_receipt_sha256"] = receipt_hash
    admitted = build_observation_registry([diagnostic], aliases)
    assert admitted[0]["product_state"] == "ORCD_COMPLETE_DEFERRED"
    assert admitted[0]["diagnostic_admission_receipt_sha256"] == receipt_hash

    failed_path, failed_hash = _diagnostic_receipt(
        tmp_path / "failed.json", failed_gate="external_quality_join"
    )
    failed = dict(
        diagnostic,
        diagnostic_admission_receipt_path=failed_path,
        diagnostic_admission_receipt_sha256=failed_hash,
    )
    with pytest.raises(FM0ContractError, match="failed evidence gates"):
        build_observation_registry([failed], aliases)

    not_passed_path, not_passed_hash = _diagnostic_receipt(
        tmp_path / "not-passed.json", passed=False
    )
    not_passed = dict(
        diagnostic,
        diagnostic_admission_receipt_path=not_passed_path,
        diagnostic_admission_receipt_sha256=not_passed_hash,
    )
    with pytest.raises(FM0ContractError, match="is not passed"):
        build_observation_registry([not_passed], aliases)

    bad_schema_path, bad_schema_hash = _diagnostic_receipt(
        tmp_path / "bad-schema.json", schema="wrong"
    )
    bad_schema = dict(
        diagnostic,
        diagnostic_admission_receipt_path=bad_schema_path,
        diagnostic_admission_receipt_sha256=bad_schema_hash,
    )
    with pytest.raises(FM0ContractError, match="schema mismatch"):
        build_observation_registry([bad_schema], aliases)

    wrong_path, wrong_hash = _diagnostic_receipt(tmp_path / "wrong-sector.json", sector=67)
    wrong = dict(
        diagnostic,
        diagnostic_admission_receipt_path=wrong_path,
        diagnostic_admission_receipt_sha256=wrong_hash,
    )
    with pytest.raises(FM0ContractError, match="does not match"):
        build_observation_registry([wrong], aliases)


def test_strict_no_bls_gate_and_immutable_registry(tmp_path) -> None:
    with pytest.raises(FM0ContractError, match="forbidden"):
        assert_no_search_columns(["raw_flux", "bls_period"])
    for search_field in (
        "period",
        "epoch",
        "depth",
        "duration",
        "sde",
        "teacher_score",
    ):
        with pytest.raises(FM0ContractError, match="forbidden"):
            assert_no_search_columns(["observation_key", search_field])
    aliases = build_alias_registry(
        [{"gaia_dr3_source_id": "100", "tic_id": "10"}]
    )
    summary = write_registry_release(tmp_path, aliases)
    assert not summary["certifies_full_campaign"]
    assert write_registry_release(tmp_path, aliases) == summary
    (tmp_path / "aliases.csv").write_text("different\n", encoding="utf-8")
    with pytest.raises(FM0ContractError, match="refusing to replace"):
        write_registry_release(tmp_path, aliases)
