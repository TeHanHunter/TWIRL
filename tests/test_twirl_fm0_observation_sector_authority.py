from __future__ import annotations

import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "scripts"
    / "stage5_validation"
    / "build_twirl_fm0_observation_sector_authority.py"
)
SPEC = importlib.util.spec_from_file_location("fm0_sector_authority", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


FIELDS = (
    "input_release_schema_version",
    "observation_key",
    "product_instance_id",
    "source_sha256",
    "leakage_component_id",
    "source_partition",
    "product_state",
    "relative_path",
    "sha256",
    "input_source_sha256",
    "n_cadences",
    "n_segments",
    "view_present_json",
    "host_visit_offset_cadences",
    "host_visit_gap_cadences",
    "host_visit_overlaps_previous",
    "input_adapter",
    "scientific_training_eligible",
)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _row(key: str, partition: str) -> dict[str, str]:
    return {
        "input_release_schema_version": "twirl_fm0_1_input_release_v1",
        "observation_key": key,
        "product_instance_id": f"product_{key}",
        "source_sha256": "a" * 64,
        "leakage_component_id": f"leakage_{key}",
        "source_partition": partition,
        "product_state": "A2V1_ACCEPTED",
        "relative_path": f"shards/{key}.npz",
        "sha256": "b" * 64,
        "input_source_sha256": "c" * 64,
        "n_cadences": "2048",
        "n_segments": "1",
        "view_present_json": "[1,1,1,1,1,1]",
        "host_visit_offset_cadences": "0.0",
        "host_visit_gap_cadences": "0.0",
        "host_visit_overlaps_previous": "False",
        "input_adapter": "a2v1_hdf5_quality_aware_v1",
        "scientific_training_eligible": "True",
    }


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _fixture(tmp_path: Path):
    s56 = tmp_path / "s56.csv"
    s57 = tmp_path / "s57.csv"
    merged = tmp_path / "merged.csv"
    summary = tmp_path / "merge_summary.json"
    train = _row("observation_train", "poc_train")
    development = _row("observation_development", "poc_development")
    sealed = _row("observation_sealed", "poc_sealed_test")
    _write_csv(s56, [train, sealed])
    _write_csv(s57, [development])
    _write_csv(merged, [sealed, development, train])
    summary.write_text(
        json.dumps(
            {
                "manifest_sha256": _sha(merged),
                "sectors": [56, 57],
                "sector_bindings": {
                    "56": {"manifest_sha256": _sha(s56), "n_observations": 2},
                    "57": {"manifest_sha256": _sha(s57), "n_observations": 1},
                },
            },
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return merged, summary, [(56, s56), (57, s57)]


def test_builds_exact_nonsealed_four_column_authority(tmp_path: Path) -> None:
    merged, summary, sectors = _fixture(tmp_path)
    payload, receipt = MODULE.build_authority(
        merged_manifest=merged,
        merge_summary=summary,
        merge_summary_sha256=_sha(summary),
        sector_manifests=sectors,
    )
    rows = list(csv.DictReader(io.StringIO(payload.decode("utf-8"))))
    assert list(rows[0]) == list(MODULE.OUTPUT_FIELDS)
    assert [row["observation_key"] for row in rows] == [
        "observation_development",
        "observation_train",
    ]
    assert {row["sector"] for row in rows} == {"56", "57"}
    assert receipt["partition_counts"]["poc_sealed_test"] == 1
    assert receipt["sealed_rows_emitted"] == 0


def test_rejects_merge_or_identity_drift(tmp_path: Path) -> None:
    merged, summary, sectors = _fixture(tmp_path)
    with pytest.raises(ValueError, match="merge-summary SHA"):
        MODULE.build_authority(
            merged_manifest=merged,
            merge_summary=summary,
            merge_summary_sha256="0" * 64,
            sector_manifests=sectors,
        )
    rows = list(csv.DictReader(merged.read_text(encoding="utf-8").splitlines()))
    rows[0]["leakage_component_id"] = "leakage_tampered"
    _write_csv(merged, rows)
    data = json.loads(summary.read_text(encoding="utf-8"))
    data["manifest_sha256"] = _sha(merged)
    summary.write_text(json.dumps(data, sort_keys=True) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="merged/per-sector mismatch"):
        MODULE.build_authority(
            merged_manifest=merged,
            merge_summary=summary,
            merge_summary_sha256=_sha(summary),
            sector_manifests=sectors,
        )
