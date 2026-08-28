from __future__ import annotations

import csv
import hashlib
import json
import re
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from twirl.models.fm0 import later_sector_admission_v2 as admission
from twirl.models.fm0.later_source_inventory import SOURCE_FIELDS
from twirl.models.fm0.registry import FM0ContractError
from twirl.models.fm0.temporal_admission import (
    FROZEN_POLICY_SHA256 as V1_FROZEN_POLICY_SHA256,
)

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs/models/twirl_fm0_2_later_sector_admission_v2.yaml"
LEDGER = ROOT / "configs/models/twirl_fm0_later_sector_exclusions_v1.yaml"
PRODUCER = "a" * 40


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, payload: object) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return _sha(path)


def _binding(path: Path) -> admission.ReceiptBinding:
    return admission.ReceiptBinding(path=path, sha256=_sha(path))


def _source_row(*, sector: int, partition: str, suffix: int) -> dict[str, object]:
    row: dict[str, object] = {field: "" for field in SOURCE_FIELDS}
    gaia = sector * 10_000 + suffix
    tic = sector * 1_000 + suffix
    component = hashlib.sha256(f"{gaia}:{tic}".encode("ascii")).hexdigest()
    row.update(
        {
            "source_row_schema_version": "twirl_fm0_later_retained_source_row_v1",
            "sector": sector,
            "gaia_dr3_source_id": gaia,
            "tic_id": tic,
            "leakage_component_id": f"leakage_{component}",
            "source_partition": partition,
            "identity_source": "checksum_bound_corpus_selection_v1",
            "target_authority_sha256": "b" * 64,
            "identity_ambiguous": "false",
            "product_state": admission.DEFERRED_PRODUCT_STATE,
            "n_hdf5_products": 2,
            "orbits_json": json.dumps([2 * sector + 7, 2 * sector + 8]),
            "hdf5_paths_json": "[]",
            "hdf5_sha256_json": "[]",
            "hdf5_evidence_json": "[]",
            "retained_source_sha256": "c" * 64,
            "source_row_sha256": hashlib.sha256(
                f"source:{sector}:{suffix}".encode("ascii")
            ).hexdigest(),
            "source_receipt_path": f"/source/s{sector}.json",
            "source_receipt_sha256": "d" * 64,
            "model_outcome_blind": "true",
            "panel_admission_authorized": "false",
        }
    )
    return row


def _make_sector(tmp_path: Path, sector: int) -> admission.SectorPreparationBindings:
    root = tmp_path / f"s{sector:04d}"
    root.mkdir(parents=True)
    orbits = [2 * sector + 7, 2 * sector + 8]

    quality_root = root / "quality"
    quality_root.mkdir()
    table = quality_root / "cadence_quality.csv"
    table.write_text(f"sector,cadence\n{sector},{sector * 100}\n", encoding="utf-8")
    mission = {
        "schema_version": "twirl_fm0_mission_quality_reference_v2",
        "producer_git_sha": PRODUCER,
        "sector": sector,
        "expected_orbits": orbits,
        "mission_quality_provider": "spoc" if sector == 66 else "tica",
        "hdf5_quality_join_required": True,
        "six_view_shards_verified": False,
        "panel_admission_authorized": False,
        "table_file": table.name,
        "table_sha256": _sha(table),
        "n_rows": 1,
    }
    mission_path = quality_root / "manifest.json"
    _write_json(mission_path, mission)
    (quality_root / "READY").write_text(PRODUCER + "\n", encoding="utf-8")

    upstream = root / "upstream"
    source_receipt = upstream / "source.json"
    quality_transfer = upstream / "quality-transfer.json"
    _write_json(source_receipt, {"sector": sector, "kind": "source"})
    _write_json(quality_transfer, {"sector": sector, "kind": "quality"})
    hdf5_path = root / "hdf5-quality.json"
    hdf5 = {
        "schema_version": "twirl_fm0_hdf5_quality_admission_v1",
        "sector": sector,
        "expected_orbits": orbits,
        "quality_state": "FM_HDF5_QUALITY_READY",
        "passed": True,
        "n_hdf5_products": sector + 10,
        "n_hdf5_opened": sector + 10,
        "n_unreadable_hdf5": 0,
        "n_cadences_checked": sector * 1000,
        "hdf5_openability_verified": True,
        "internal_cadence_quality_verified": True,
        "external_cadence_quality_verified": True,
        "six_view_shards_verified": False,
        "panel_admission_authorized": False,
        "source_receipt_path": str(source_receipt),
        "source_receipt_sha256": _sha(source_receipt),
        "quality_transfer_manifest_path": str(quality_transfer),
        "quality_transfer_manifest_sha256": _sha(quality_transfer),
    }
    _write_json(hdf5_path, hdf5)

    inventory_root = root / "inventory"
    inventory_root.mkdir()
    source_rows = [
        _source_row(sector=sector, partition="poc_train", suffix=1),
        _source_row(sector=sector, partition="poc_development", suffix=2),
        _source_row(sector=sector, partition="poc_sealed_test", suffix=3),
    ]
    sources_path = inventory_root / "sources.csv"
    with sources_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SOURCE_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(source_rows)
    aliases_path = inventory_root / "aliases.csv"
    aliases_path.write_text(f"sector,gaia\n{sector},{sector}\n", encoding="utf-8")
    source_summary = {
        "summary_schema_version": "twirl_fm0_later_retained_source_summary_v1",
        "sector": sector,
        "n_source_rows": len(source_rows),
        "source_rows_by_partition": {
            "poc_train": 1,
            "poc_development": 1,
            "poc_sealed_test": 1,
        },
        "sources_sha256": _sha(sources_path),
        "aliases_sha256": _sha(aliases_path),
        "hdf5_content_opened": False,
        "sealed_hdf5_content_opened": False,
        "model_outcome_blind": True,
        "six_view_shards_verified": False,
        "panel_admission_authorized": False,
    }
    source_summary_path = inventory_root / "summary.json"
    source_summary_sha = _write_json(source_summary_path, source_summary)
    (inventory_root / "READY").write_text(source_summary_sha + "\n", encoding="utf-8")

    six_path = root / "six-view.json"
    six_view = {
        "schema_version": admission.SIX_VIEW_RECEIPT_SCHEMA_VERSION,
        "sector": sector,
        "expected_orbits": orbits,
        "release_state": admission.SIX_VIEW_READY_STATE,
        "product_state": admission.DEFERRED_PRODUCT_STATE,
        "passed": True,
        "flux_view_names": [
            "raw_relative_1x1",
            "raw_relative_3x3",
            "adp_1x1",
            "adp_3x3",
            "adp015_1x1",
            "adp015_3x3",
        ],
        "n_observations": 2,
        "n_shards": 2,
        "six_view_shards_verified": True,
        "scientific_training_eligible": False,
        "panel_admission_authorized": False,
        "sealed_hdf5_content_opened": False,
        "sealed_shards_written": False,
        "model_outcome_blind": True,
        "mission_quality_reference_manifest_path": str(mission_path),
        "mission_quality_reference_manifest_sha256": _sha(mission_path),
        "hdf5_quality_receipt_path": str(hdf5_path),
        "hdf5_quality_receipt_sha256": _sha(hdf5_path),
        "source_inventory_summary_path": str(source_summary_path),
        "source_inventory_summary_sha256": _sha(source_summary_path),
    }
    _write_json(six_path, six_view)
    return admission.SectorPreparationBindings(
        sector=sector,
        mission_quality_reference=_binding(mission_path),
        hdf5_quality=_binding(hdf5_path),
        source_inventory=_binding(source_summary_path),
        six_view_release=_binding(six_path),
    )


@pytest.fixture()
def sector_bindings(tmp_path: Path) -> list[admission.SectorPreparationBindings]:
    return [_make_sector(tmp_path, sector) for sector in admission.PREPARATION_SECTORS]


def _construct(
    bindings: list[admission.SectorPreparationBindings],
) -> dict[str, object]:
    return admission.construct_preparation_pool_receipt(
        config_path=CONFIG,
        expected_config_sha256=admission.FROZEN_POLICY_SHA256,
        exclusion_ledger_path=LEDGER,
        expected_exclusion_ledger_sha256=_sha(LEDGER),
        ordered_sector_bindings=bindings,
        producer_git_sha=PRODUCER,
    )


def _patch_lightweight_validators(monkeypatch: pytest.MonkeyPatch) -> None:
    """Keep pool-policy tests small; bundle integrity has dedicated tests."""

    def mission(*, sector: int, binding: admission.ReceiptBinding):
        admission._verify_file(binding, label=f"S{sector} mission fixture")
        provider = "spoc" if sector == 66 else "tica"
        reference = SimpleNamespace(
            provider=provider,
            manifest_sha256=binding.sha256,
            table_path=binding.path.parent / "cadence_quality.csv",
            table_sha256=_sha(binding.path.parent / "cadence_quality.csv"),
            assert_unchanged=lambda: None,
        )
        return binding.path.resolve(), reference, {
            "schema_version": "twirl_fm0_mission_quality_reference_v2",
            "sha256": binding.sha256,
            "mission_quality_provider": provider,
            "path": str(binding.path.resolve()),
        }

    def source(*, sector: int, binding: admission.ReceiptBinding):
        admission._verify_file(binding, label=f"S{sector} source fixture")
        return binding.path.resolve(), {"n_source_rows": 3}, ({}, {}), {
            "schema_version": "twirl_fm0_later_retained_source_summary_v1",
            "sha256": binding.sha256,
            "n_source_rows": 3,
            "n_usable_source_rows": 2,
            "n_sealed_identity_rows": 1,
            "path": str(binding.path.resolve()),
        }

    def hdf5(
        *,
        sector: int,
        binding: admission.ReceiptBinding,
        inventory_summary,
        reference,
    ):
        del inventory_summary, reference
        admission._verify_file(binding, label=f"S{sector} HDF5 fixture")
        return binding.path.resolve(), {}, {
            "schema_version": "twirl_fm0_hdf5_quality_admission_v1",
            "sha256": binding.sha256,
            "n_hdf5_products": sector + 10,
            "n_cadences_checked": sector * 1000,
            "path": str(binding.path.resolve()),
        }

    def six(
        *,
        sector: int,
        binding: admission.ReceiptBinding,
        mission_quality: admission.ReceiptBinding,
        hdf5_quality: admission.ReceiptBinding,
        source_inventory: admission.ReceiptBinding,
        expected_usable_rows: int,
        **_kwargs,
    ):
        path = admission._verify_file(binding, label=f"S{sector} six-view fixture")
        receipt = json.loads(path.read_text(encoding="utf-8"))
        if (
            receipt.get("product_state") != admission.DEFERRED_PRODUCT_STATE
            or receipt.get("scientific_training_eligible") is not False
            or receipt.get("sealed_hdf5_content_opened") is not False
        ):
            raise FM0ContractError("six-view receipt did not pass")
        for prefix, expected in (
            ("mission_quality_reference_manifest", mission_quality),
            ("hdf5_quality_receipt", hdf5_quality),
            ("source_inventory_summary", source_inventory),
        ):
            if not admission._binding_matches(
                receipt_path=path,
                receipt=receipt,
                prefix=prefix,
                expected=expected,
            ):
                raise FM0ContractError("six-view cross-receipt binding drifted")
        return path, receipt, {
            "schema_version": admission.SIX_VIEW_RECEIPT_SCHEMA_VERSION,
            "release_state": admission.SIX_VIEW_READY_STATE,
            "sha256": binding.sha256,
            "n_observations": expected_usable_rows,
            "n_shards": expected_usable_rows,
            "path": str(path),
        }

    monkeypatch.setattr(admission, "_verify_mission_quality", mission)
    monkeypatch.setattr(admission, "_verify_source_inventory", source)
    monkeypatch.setattr(admission, "_verify_hdf5_quality", hdf5)
    monkeypatch.setattr(admission, "_verify_six_view", six)


def _rewrite_json_binding(
    binding: admission.ReceiptBinding,
    update: dict[str, object],
) -> admission.ReceiptBinding:
    payload = json.loads(binding.path.read_text(encoding="utf-8"))
    payload.update(update)
    _write_json(binding.path, payload)
    return _binding(binding.path)


def test_v2_admits_only_deferred_preparation_and_keeps_all_claims_false(
    sector_bindings: list[admission.SectorPreparationBindings],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _patch_lightweight_validators(monkeypatch)
    receipt = _construct(sector_bindings)
    assert receipt["preparation_pool_sectors"] == list(range(66, 78))
    assert receipt["excluded_sectors"] == [65]
    assert receipt["product_state"] == "ORCD_COMPLETE_DEFERRED"
    assert receipt["preparation_pool_admitted"] is True
    assert receipt["n_source_rows"] == 36
    assert receipt["n_nonsealed_preparation_rows"] == 24
    assert receipt["n_sealed_identity_rows"] == 12
    assert receipt["n_six_view_shards"] == 24
    for field in (
        "a2v1_accepted",
        "scientific_training_eligible",
        "panel_admission_authorized",
        "temporal_panel_frozen",
        "model_training_authorized",
        "sealed_test_access_authorized",
        "formal_model_gate_authorized",
        "production_model_claim",
        "foundation_model_claim",
        "prospective_test_claim",
        "s78_included",
        "s79_s80_touched",
        "sealed_aperture_photometry_opened",
        "sealed_shards_written",
    ):
        assert receipt[field] is False
    assert (
        receipt["sector_wide_hdf5_identity_cadence_quality_qa_performed"] is True
    )
    output = tmp_path / "pool-receipt.json"
    path, digest = admission.write_preparation_pool_receipt(
        output_path=output, receipt=receipt
    )
    assert path == output.resolve()
    assert digest == _sha(output)
    changed_receipt = dict(receipt, n_source_rows=37)
    with pytest.raises(FM0ContractError, match="refusing to replace"):
        admission.write_preparation_pool_receipt(
            output_path=output, receipt=changed_receipt
        )


@pytest.mark.parametrize("mutation", ("missing", "reversed", "s65"))
def test_runtime_pool_must_be_exact_chronological_s66_s77(
    sector_bindings: list[admission.SectorPreparationBindings], mutation: str
) -> None:
    if mutation == "missing":
        candidate = sector_bindings[:-1]
    elif mutation == "reversed":
        candidate = list(reversed(sector_bindings))
    else:
        candidate = [replace(sector_bindings[0], sector=65), *sector_bindings[1:]]
    with pytest.raises(FM0ContractError, match="exactly chronological S66--S77"):
        _construct(candidate)


def test_exclusion_ledger_hash_is_a_frozen_runtime_binding(
    sector_bindings: list[admission.SectorPreparationBindings], tmp_path: Path
) -> None:
    altered = tmp_path / "ledger.yaml"
    altered.write_bytes(LEDGER.read_bytes() + b"\n")
    with pytest.raises(FM0ContractError, match="differs from v2 policy"):
        admission.construct_preparation_pool_receipt(
            config_path=CONFIG,
            expected_config_sha256=admission.FROZEN_POLICY_SHA256,
            exclusion_ledger_path=altered,
            expected_exclusion_ledger_sha256=_sha(altered),
            ordered_sector_bindings=sector_bindings,
            producer_git_sha=PRODUCER,
        )


def test_six_view_cannot_promote_deferred_sector_to_a2v1_accepted(
    sector_bindings: list[admission.SectorPreparationBindings],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _patch_lightweight_validators(monkeypatch)
    first = sector_bindings[0]
    changed = _rewrite_json_binding(
        first.six_view_release,
        {"product_state": "A2V1_ACCEPTED", "scientific_training_eligible": True},
    )
    sector_bindings[0] = replace(first, six_view_release=changed)
    with pytest.raises(FM0ContractError, match="six-view receipt did not pass"):
        _construct(sector_bindings)


def test_six_view_must_bind_the_exact_other_three_receipts(
    sector_bindings: list[admission.SectorPreparationBindings],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _patch_lightweight_validators(monkeypatch)
    first = sector_bindings[0]
    changed = _rewrite_json_binding(
        first.six_view_release,
        {"hdf5_quality_receipt_sha256": "f" * 64},
    )
    sector_bindings[0] = replace(first, six_view_release=changed)
    with pytest.raises(FM0ContractError, match="cross-receipt binding drifted"):
        _construct(sector_bindings)


def test_sealed_content_access_fails_closed(
    sector_bindings: list[admission.SectorPreparationBindings],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _patch_lightweight_validators(monkeypatch)
    first = sector_bindings[0]
    changed = _rewrite_json_binding(
        first.six_view_release, {"sealed_hdf5_content_opened": True}
    )
    sector_bindings[0] = replace(first, six_view_release=changed)
    with pytest.raises(FM0ContractError, match="six-view receipt did not pass"):
        _construct(sector_bindings)


def test_receipt_bytes_are_bound_not_only_their_paths(
    sector_bindings: list[admission.SectorPreparationBindings],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _patch_lightweight_validators(monkeypatch)
    sector_bindings[0].hdf5_quality.path.write_text("{}\n", encoding="utf-8")
    with pytest.raises(FM0ContractError, match="hash mismatch"):
        _construct(sector_bindings)


def test_production_config_is_frozen_and_v1_policy_is_untouched() -> None:
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    assert _sha(CONFIG) == admission.FROZEN_POLICY_SHA256
    assert config["selection"]["ordered_sectors"] == list(range(66, 78))
    assert config["exclusion_ledger"]["required_sha256"] == _sha(LEDGER)
    assert config["product_boundary"]["a2v1_accepted"] is False
    assert config["product_boundary"]["scientific_training_eligible"] is False
    assert config["selection"]["s79_s80_touched"] is False
    assert V1_FROZEN_POLICY_SHA256 == (
        "527a7b4d9f9c452f02576eea7f155abe65ad8439057317e8bc10c80b0ed93da3"
    )


def test_mission_quality_gate_uses_the_full_reference_loader(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = tmp_path / "quality"
    manifest = root / "manifest.json"
    table = root / "cadence_quality.csv"
    _write_json(manifest, {"fixture": True})
    table.write_text("cadence\n1\n", encoding="utf-8")
    binding = _binding(manifest)
    unchanged_calls = 0

    def assert_unchanged() -> None:
        nonlocal unchanged_calls
        unchanged_calls += 1

    reference = SimpleNamespace(
        sector=66,
        expected_orbits=(139, 140),
        provider="spoc",
        manifest_path=manifest.resolve(),
        manifest_sha256=binding.sha256,
        table_path=table.resolve(),
        table_sha256=_sha(table),
        assert_unchanged=assert_unchanged,
    )
    observed: dict[str, object] = {}

    def load_reference(**kwargs):
        observed.update(kwargs)
        return reference

    monkeypatch.setattr(admission, "load_mission_quality_reference", load_reference)
    _path, loaded, summary = admission._verify_mission_quality(
        sector=66, binding=binding
    )
    assert loaded is reference
    assert observed == {
        "reference_dir": root.resolve(),
        "sector": 66,
        "expected_orbits": (139, 140),
    }
    assert summary["mission_quality_provider"] == "spoc"
    assert unchanged_calls == 1


def test_source_inventory_gate_cross_checks_both_strict_loaders_and_bound_sha(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = tmp_path / "inventory"
    summary_path = root / "summary.json"
    sources_path = root / "sources.csv"
    summary = {
        "n_source_rows": 2,
        "source_rows_by_partition": {
            "poc_train": 1,
            "poc_development": 0,
            "poc_sealed_test": 1,
        },
    }
    _write_json(summary_path, summary)
    sources_path.write_text("fixture\n", encoding="utf-8")
    binding = _binding(summary_path)
    all_rows = [
        {
            "source_partition": "poc_train",
            "source_row_sha256": "1" * 64,
            "hdf5_evidence_json": "[]",
        },
        {
            "source_partition": "poc_sealed_test",
            "source_row_sha256": "2" * 64,
            "hdf5_evidence_json": "[]",
        },
    ]
    calls: dict[str, object] = {}

    def inventory_bundle(**kwargs):
        calls["bundle"] = kwargs
        return (
            summary,
            summary_path.resolve(),
            binding.sha256,
            all_rows,
            sources_path.resolve(),
            _sha(sources_path),
        )

    def source_rows(*args, **kwargs):
        calls["rows"] = (args, kwargs)
        drifted = dict(all_rows[0], hdf5_evidence_json='[{"drift":true}]')
        return (drifted,)

    monkeypatch.setattr(admission, "_load_inventory_bundle", inventory_bundle)
    monkeypatch.setattr(admission, "load_later_source_rows", source_rows)
    with pytest.raises(FM0ContractError, match="identity/evidence closure drifted"):
        admission._verify_source_inventory(sector=67, binding=binding)
    assert calls["bundle"] == {
        "source_inventory_dir": root.resolve(),
        "sector": 67,
        "expected_summary_sha256": binding.sha256,
    }
    assert calls["rows"] == (
        (root.resolve(),),
        {
            "expected_summary_sha256": binding.sha256,
            "allowed_source_partitions": ("poc_development", "poc_train"),
        },
    )


def test_hdf5_gate_forwards_exact_sha_to_strict_upstream_validator(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "hdf5-quality.json"
    _write_json(path, {"fixture": True})
    binding = _binding(path)
    receipt = {"n_hdf5_products": 4, "n_cadences_checked": 320}
    unchanged_calls = 0

    def assert_unchanged() -> None:
        nonlocal unchanged_calls
        unchanged_calls += 1

    reference = SimpleNamespace(assert_unchanged=assert_unchanged)
    inventory_summary = {"n_hdf5_products_declared": 4}
    calls: dict[str, object] = {}

    def load_hdf5(**kwargs):
        calls.update(kwargs)
        return receipt, path.resolve(), binding.sha256

    monkeypatch.setattr(admission, "_load_hdf5_quality_receipt", load_hdf5)
    _path, loaded, summary = admission._verify_hdf5_quality(
        sector=67,
        binding=binding,
        inventory_summary=inventory_summary,
        reference=reference,
    )
    assert loaded is receipt
    assert calls == {
        "path": path.resolve(),
        "sector": 67,
        "inventory_summary": inventory_summary,
        "reference": reference,
        "expected_receipt_sha256": binding.sha256,
    }
    assert summary["n_cadences_checked"] == 320
    assert unchanged_calls == 1


@pytest.mark.parametrize(
    "validator_error",
    (
        "six-view release root closure drifted; missing shard",
        "later six-view shard hash drifted",
        "later six-view shard closure drifted; extra shard",
        "later six-view receipt contract drifted; wrong provider",
        "six-view release is not read-only",
    ),
)
def test_six_view_gate_propagates_every_full_bundle_validation_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    validator_error: str,
) -> None:
    receipt_path = tmp_path / "release" / "receipt.json"
    _write_json(receipt_path, {"fixture": True})
    binding = _binding(receipt_path)
    observed: dict[str, object] = {}

    def validate(*args, **kwargs):
        observed["args"] = args
        observed["kwargs"] = kwargs
        raise FM0ContractError(validator_error)

    monkeypatch.setattr(admission, "validate_later_sector_release", validate)
    with pytest.raises(FM0ContractError, match=re.escape(validator_error)):
        admission._verify_six_view(
            sector=67,
            binding=binding,
            mission_quality=binding,
            hdf5_quality=binding,
            source_inventory=binding,
            expected_usable_rows=1,
            reference=SimpleNamespace(provider="tica"),
            inventory_summary={},
            selected_source_rows=(),
        )
    assert observed == {
        "args": (receipt_path.parent.resolve(),),
        "kwargs": {
            "expected_receipt_sha256": binding.sha256,
            "require_read_only": True,
        },
    }
