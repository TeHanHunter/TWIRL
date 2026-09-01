from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest
import yaml

from twirl.models.fm0 import temporal_panel as panel
from twirl.models.fm0.input_release import (
    INPUT_RELEASE_SCHEMA_VERSION,
    MANIFEST_COLUMNS,
)
from twirl.models.fm0.later_sector_admission_v2 import POOL_RECEIPT_SCHEMA_VERSION
from twirl.models.fm0.later_sector_release import (
    LATER_SIX_VIEW_MANIFEST_FIELDS,
    LATER_SIX_VIEW_MANIFEST_SCHEMA_VERSION,
    LATER_SIX_VIEW_READY_STATE,
    LATER_SIX_VIEW_RECEIPT_SCHEMA_VERSION,
)
from twirl.models.fm0.registry import FM0ContractError, sha256_file

REPO_ROOT = Path(__file__).resolve().parents[1]
POLICY_PATH = REPO_ROOT / "configs/models/twirl_fm0_2_s66_s77_temporal_panel_v1.yaml"
GIT_SHA = "6" * 40
REPEATED_COMPONENT = "leakage_" + "a" * 64


def _hash_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _write_csv(
    path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, value: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _baseline_row(*, partition: str) -> dict[str, object]:
    return {
        "input_release_schema_version": INPUT_RELEASE_SCHEMA_VERSION,
        "observation_key": "baseline-observation",
        "product_instance_id": "baseline-product",
        "source_sha256": _hash_text("baseline-source"),
        "leakage_component_id": REPEATED_COMPONENT,
        "source_partition": partition,
        "product_state": "A2V1_ACCEPTED",
        "relative_path": "shards/baseline.npz",
        "sha256": _hash_text("baseline-shard"),
        "input_source_sha256": _hash_text("baseline-input"),
        "n_cadences": 100,
        "n_segments": 1,
        "view_present_json": "[1,1,1,1,1,1]",
        "host_visit_offset_cadences": 0,
        "host_visit_gap_cadences": "",
        "host_visit_overlaps_previous": "false",
        "input_adapter": "fixture",
        "scientific_training_eligible": "true",
    }


def _manifest_row(
    *, sector: int, name: str, component: str, partition: str
) -> dict[str, object]:
    gaia = str(10_000_000 + sector * 10 + {"repeat": 1, "new": 2, "train": 3}[name])
    return {
        "manifest_schema_version": LATER_SIX_VIEW_MANIFEST_SCHEMA_VERSION,
        "observation_key": f"observation-s{sector}-{name}",
        "product_instance_id": f"product-s{sector}-{name}",
        "physical_source_id": f"gaia_dr3:{gaia}",
        "gaia_dr3_source_id": gaia,
        "tic_id": str(20_000_000 + sector * 10 + len(name)),
        "sector": sector,
        "leakage_component_id": component,
        "source_partition": partition,
        "product_state": "ORCD_COMPLETE_DEFERRED",
        "relative_path": f"shards/observation-s{sector}-{name}.npz",
        "sha256": _hash_text(f"shard-s{sector}-{name}"),
        "input_source_sha256": _hash_text(f"input-s{sector}-{name}"),
        "source_row_sha256": _hash_text(f"row-s{sector}-{name}"),
        "n_cadences": 100,
        "n_segments": 1,
        "view_present_json": "[1,1,1,1,1,1]",
        "input_adapter": "twirl_fm0_later_hdf5_adapter_v1",
        "mission_quality_provider": "spoc" if sector == 66 else "tica",
        "mission_quality_reference_manifest_sha256": _hash_text(f"quality-s{sector}"),
        "hdf5_quality_receipt_sha256": _hash_text(f"hdf5-s{sector}"),
        "full_visit_shard": "true",
        "model_context_length_bound": "false",
        "scientific_training_eligible": "false",
        "panel_admission_authorized": "false",
    }


def _make_fixture(
    tmp_path: Path,
    *,
    baseline_partition: str = "poc_development",
    baseline_component: str = REPEATED_COMPONENT,
    duplicate_observation: bool = False,
    sealed_release_row: bool = False,
) -> dict[str, object]:
    baseline = tmp_path / "baseline.csv"
    _write_csv(
        baseline,
        MANIFEST_COLUMNS,
        [
            {
                **_baseline_row(partition=baseline_partition),
                "leakage_component_id": baseline_component,
            }
        ],
    )
    baseline_sha = sha256_file(baseline)

    config = yaml.safe_load(POLICY_PATH.read_text(encoding="utf-8"))
    config["campaign_id"] = "fixture_temporal_panel"
    config["selection"]["baseline_manifest_sha256"] = baseline_sha
    config_path = tmp_path / "policy.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    sector_receipts: list[dict[str, object]] = []
    for sector in panel.TEMPORAL_PANEL_SECTORS:
        root = tmp_path / f"s{sector}"
        new_component = "leakage_" + _hash_text(f"new-s{sector}")
        train_component = "leakage_" + _hash_text(f"train-s{sector}")
        train_partition = "poc_sealed_test" if sealed_release_row else "poc_train"
        rows = [
            _manifest_row(
                sector=sector,
                name="repeat",
                component=REPEATED_COMPONENT,
                partition="poc_development",
            ),
            _manifest_row(
                sector=sector,
                name="new",
                component=new_component,
                partition="poc_development",
            ),
            _manifest_row(
                sector=sector,
                name="train",
                component=train_component,
                partition=train_partition,
            ),
        ]
        if duplicate_observation and sector == 67:
            rows[0]["observation_key"] = "observation-s66-repeat"
        manifest = root / "manifest.csv"
        _write_csv(manifest, LATER_SIX_VIEW_MANIFEST_FIELDS, rows)
        manifest_sha = sha256_file(manifest)
        source_manifest_sha = _hash_text(f"source-manifest-s{sector}")
        visit_timing_sha = _hash_text(f"visit-timing-s{sector}")
        receipt = {
            "schema_version": LATER_SIX_VIEW_RECEIPT_SCHEMA_VERSION,
            "release_state": LATER_SIX_VIEW_READY_STATE,
            "passed": True,
            "sector": sector,
            "expected_orbits": [2 * sector + 7, 2 * sector + 8],
            "product_state": "ORCD_COMPLETE_DEFERRED",
            "n_observations": len(rows),
            "n_shards": len(rows),
            "n_cadences": 300,
            "n_segments": 3,
            "source_rows_by_partition": {
                "poc_development": 2,
                train_partition: 1,
            },
            "manifest_sha256": manifest_sha,
            "source_manifest_sha256": source_manifest_sha,
            "visit_timing_sha256": visit_timing_sha,
            "full_visit_shards": True,
            "model_context_length_bound": False,
            "sealed_hdf5_content_opened": False,
            "sealed_shards_written": False,
            "model_outcome_blind": True,
            "six_view_shards_verified": True,
            "scientific_training_eligible": False,
            "panel_admission_authorized": False,
        }
        receipt_path = root / "receipt.json"
        _write_json(receipt_path, receipt)
        receipt_sha = sha256_file(receipt_path)
        ready = root / "READY"
        ready.write_text(receipt_sha + "\n", encoding="utf-8")
        manifest.chmod(0o444)
        receipt_path.chmod(0o444)
        ready.chmod(0o444)
        root.chmod(0o555)
        sector_receipts.append(
            {
                "sector": sector,
                "expected_orbits": [2 * sector + 7, 2 * sector + 8],
                "product_state": "ORCD_COMPLETE_DEFERRED",
                "preparation_admitted": True,
                "a2v1_accepted": False,
                "scientific_training_eligible": False,
                "evidence": {
                    "six_view_release": {
                        "schema_version": LATER_SIX_VIEW_RECEIPT_SCHEMA_VERSION,
                        "release_state": LATER_SIX_VIEW_READY_STATE,
                        "sha256": receipt_sha,
                        "n_observations": len(rows),
                        "n_shards": len(rows),
                        "n_cadences": 300,
                        "manifest_sha256": manifest_sha,
                        "source_manifest_sha256": source_manifest_sha,
                        "visit_timing_sha256": visit_timing_sha,
                        "path": str(receipt_path.resolve()),
                    }
                },
            }
        )

    admission = {
        "schema_version": POOL_RECEIPT_SCHEMA_VERSION,
        "preparation_pool_sectors": list(panel.TEMPORAL_PANEL_SECTORS),
        "excluded_sectors": [65],
        "product_state": "ORCD_COMPLETE_DEFERRED",
        "preparation_pool_admitted": True,
        "a2v1_accepted": False,
        "scientific_training_eligible": False,
        "panel_admission_authorized": False,
        "temporal_panel_frozen": False,
        "model_training_authorized": False,
        "sealed_test_access_authorized": False,
        "s78_included": False,
        "s79_s80_touched": False,
        "sealed_aperture_photometry_opened": False,
        "sealed_shards_written": False,
        "sector_wide_hdf5_identity_cadence_quality_qa_performed": True,
        "n_source_rows": 48,
        "n_nonsealed_preparation_rows": 36,
        "n_sealed_identity_rows": 12,
        "n_six_view_shards": 36,
        "sector_receipts": sector_receipts,
    }
    admission_path = tmp_path / "admission.json"
    _write_json(admission_path, admission)
    admission_path.chmod(0o444)
    return {
        "config_path": config_path,
        "config_sha256": sha256_file(config_path),
        "admission_receipt_path": admission_path,
        "admission_receipt_sha256": sha256_file(admission_path),
        "baseline_manifest_path": baseline,
        "baseline_manifest_sha256": baseline_sha,
        "producer_git_sha": GIT_SHA,
        "output_dir": tmp_path / "panel",
    }


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def test_freeze_is_identity_only_and_classifies_repeated_and_new(
    tmp_path: Path,
) -> None:
    args = _make_fixture(tmp_path)
    frozen = panel.freeze_temporal_panel(**args)

    rows = _read_rows(frozen.panel_path)
    assert len(rows) == 24
    assert {row["source_partition"] for row in rows} == {"poc_development"}
    assert {row["scientific_training_eligible"] for row in rows} == {"false"}
    assert {row["development_evaluation_eligible"] for row in rows} == {"true"}
    assert sum(row["temporal_cohort"] == "repeated" for row in rows) == 12
    assert sum(row["temporal_cohort"] == "new" for row in rows) == 12
    assert not any(
        (Path(row["sector_release_root"]) / row["relative_path"]).exists()
        for row in rows
    )

    receipt = frozen.receipt
    assert receipt["n_repeated_components"] == 1
    assert receipt["n_new_components"] == 12
    assert receipt["n_train_rows_excluded"] == 12
    assert receipt["n_sealed_rows_emitted"] == 0
    assert receipt["shard_payloads_opened"] is False
    assert receipt["scientific_training_eligible"] is False
    assert receipt["development_evaluation_eligible"] is True
    assert (
        frozen.ready_path.read_text(encoding="utf-8").strip() == frozen.receipt_sha256
    )
    assert all(
        not (path.stat().st_mode & 0o222) for path in frozen.output_dir.rglob("*")
    )
    with pytest.raises(FM0ContractError, match="refusing to overwrite"):
        panel.freeze_temporal_panel(**args)


def test_sector_manifest_hash_drift_fails_closed(tmp_path: Path) -> None:
    args = _make_fixture(tmp_path)
    manifest = tmp_path / "s70/manifest.csv"
    (tmp_path / "s70").chmod(0o755)
    manifest.chmod(0o644)
    manifest.write_text(manifest.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    manifest.chmod(0o444)
    (tmp_path / "s70").chmod(0o555)
    with pytest.raises(FM0ContractError, match="manifest hash mismatch"):
        panel.freeze_temporal_panel(**args)


def test_duplicate_later_observation_fails_closed(tmp_path: Path) -> None:
    args = _make_fixture(tmp_path, duplicate_observation=True)
    with pytest.raises(FM0ContractError, match="identities are not unique"):
        panel.freeze_temporal_panel(**args)


def test_sealed_release_row_fails_closed(tmp_path: Path) -> None:
    args = _make_fixture(tmp_path, sealed_release_row=True)
    with pytest.raises(FM0ContractError, match="manifest identity drifted"):
        panel.freeze_temporal_panel(**args)


@pytest.mark.parametrize("index,forbidden_sector", ((0, 65), (-1, 78)))
def test_s65_or_later_sector_admission_drift_fails_closed(
    tmp_path: Path, index: int, forbidden_sector: int
) -> None:
    args = _make_fixture(tmp_path)
    admission_path = Path(args["admission_receipt_path"])
    admission_path.chmod(0o644)
    admission = json.loads(admission_path.read_text(encoding="utf-8"))
    admission["preparation_pool_sectors"][index] = forbidden_sector
    _write_json(admission_path, admission)
    admission_path.chmod(0o444)
    args["admission_receipt_sha256"] = sha256_file(admission_path)
    with pytest.raises(FM0ContractError, match="admission contract drifted"):
        panel.freeze_temporal_panel(**args)


def test_baseline_partition_drift_fails_closed(tmp_path: Path) -> None:
    args = _make_fixture(tmp_path, baseline_partition="poc_train")
    with pytest.raises(FM0ContractError, match="crosses its frozen baseline partition"):
        panel.freeze_temporal_panel(**args)


def test_empty_temporal_cohort_fails_before_publication(tmp_path: Path) -> None:
    args = _make_fixture(
        tmp_path,
        baseline_component="leakage_" + _hash_text("not-in-later-sectors"),
    )
    with pytest.raises(FM0ContractError, match="nonempty repeated and new"):
        panel.freeze_temporal_panel(**args)
    assert not Path(args["output_dir"]).exists()
