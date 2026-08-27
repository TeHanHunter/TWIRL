from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from twirl.models.fm0.orcd_source_admission import (
    SOURCE_READY_STATE,
    verify_archive_source,
    verify_retained_sector_source,
)
from twirl.models.fm0.registry import FM0ContractError


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True) + "\n")


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_archive_source_is_ready_but_not_panel_admitted(tmp_path: Path) -> None:
    sector = 65
    archive = tmp_path / "s0065_A2v1_raw_hdf5.tar"
    archive.write_bytes(b"test-archive")
    digest = _sha(archive)
    (tmp_path / f"{archive.name}.sha256").write_text(f"{digest}  {archive.name}\n")
    (tmp_path / "READY").write_text(digest + "\n")
    _write_json(
        tmp_path / "READY_ORCD",
        {"tag": "s0065", "task_id": "task-1", "verified_at": "now"},
    )
    _write_json(
        tmp_path / "summary.json",
        {
            "schema_version": "twirl_fm0_1_a2v1_sector_tar_v1",
            "sector": sector,
            "orbits": [137, 138],
            "archive": archive.name,
            "archive_bytes": archive.stat().st_size,
            "archive_sha256": digest,
            "accepted_sources_modified": False,
            "n_hdf5_files": 12,
        },
    )
    receipt = verify_archive_source(
        archive_dir=tmp_path, sector=sector, expected_orbits=(137, 138)
    )
    assert receipt["source_state"] == SOURCE_READY_STATE
    assert receipt["product_state"] == "A2V1_ACCEPTED"
    assert receipt["n_hdf5_products_declared"] == 12
    assert receipt["hdf5_openability_verified"] is False
    assert receipt["panel_admission_authorized"] is False


def _retained_sector(tmp_path: Path, *, sector: int = 66) -> Path:
    root = tmp_path / f"s{sector:04d}"
    for orbit in (139, 140):
        for camera in range(1, 5):
            for ccd in range(1, 5):
                cell = f"s{sector:04d}_o{orbit}_cam{camera}_ccd{ccd}"
                cell_root = root / "outputs" / cell
                complete_path = cell_root / "complete" / "slurm-1.json"
                complete = {
                    "schema": "twirl-a2v1-orcd-cell-complete-v1",
                    "state": "complete",
                    "cell": cell,
                    "attempt_id": "slurm-1",
                    "input_manifest_sha256": "a" * 64,
                    "output_manifest_sha256": "b" * 64,
                    "environment_manifest_sha256": "c" * 64,
                }
                _write_json(complete_path, complete)
                retained_root = cell_root / "retained" / "slurm-1-v2"
                retained_root.mkdir(parents=True)
                retention = {
                    "schema": "twirl-a2v1-orcd-retention-v1",
                    "ok": True,
                    "cell": cell,
                    "attempt_id": "slurm-1",
                    "pdo_return_deferred": True,
                    "pdo_sector_accepted": False,
                    "retained_root": str(retained_root),
                    "retained_root_validation": {
                        "ok": True,
                        "checks": {
                            "completion_marker_gate_reused": True,
                            "exact_epsf_count_nonzero": True,
                            "hdf5_count_nonzero_tic_authority": True,
                            "output_manifest_paths": True,
                        },
                        "completion_json_sha256": _sha(complete_path),
                        "input_manifest_sha256": "a" * 64,
                        "output_manifest_sha256": "b" * 64,
                        "environment_manifest_sha256": "c" * 64,
                        "outputs": {
                            "hdf5": {
                                "validated_tics": 3,
                                "missing_expected": [],
                                "extra_unrequested": [],
                            }
                        },
                    },
                }
                _write_json(
                    cell_root / "retained" / "slurm-1-v2.retention.json",
                    retention,
                )
    return root


def test_retained_sector_requires_all_32_receipt_bound_cells(tmp_path: Path) -> None:
    root = _retained_sector(tmp_path)
    historical = next(root.glob("outputs/*/retained/*.retention.json"))
    payload = json.loads(historical.read_text())
    label = payload["cell"]
    payload["cell"] = {
        "sector": 66,
        "orbit": int(label.split("_o", 1)[1].split("_", 1)[0]),
        "camera": int(label.split("_cam", 1)[1].split("_", 1)[0]),
        "ccd": int(label.rsplit("_ccd", 1)[1]),
        "label": label,
    }
    validation = payload["retained_root_validation"]
    validation["hashes"] = {
        key: validation.pop(key)
        for key in (
            "completion_json_sha256",
            "input_manifest_sha256",
            "output_manifest_sha256",
            "environment_manifest_sha256",
        )
    }
    historical.unlink()
    _write_json(historical.parent / "legacy" / "retention.json", payload)
    modern = next(root.glob("outputs/*/retained/*.retention.json"))
    _write_json(
        modern.parent / "older-attempt" / "retention.json",
        json.loads(modern.read_text()),
    )
    receipt = verify_retained_sector_source(
        sector_root=root, sector=66, expected_orbits=(139, 140)
    )
    assert receipt["source_state"] == SOURCE_READY_STATE
    assert receipt["product_state"] == "ORCD_COMPLETE_DEFERRED"
    assert receipt["pdo_sector_accepted"] is False
    assert receipt["n_cells"] == 32
    assert receipt["n_hdf5_products_declared"] == 96
    assert receipt["panel_admission_authorized"] is False


def test_retained_sector_fails_on_missing_cell(tmp_path: Path) -> None:
    root = _retained_sector(tmp_path)
    missing = root / "outputs" / "s0066_o139_cam1_ccd1"
    for path in sorted(missing.rglob("*"), reverse=True):
        if path.is_file():
            path.unlink()
        else:
            path.rmdir()
    missing.rmdir()
    with pytest.raises(FM0ContractError, match="retained cells differ"):
        verify_retained_sector_source(
            sector_root=root, sector=66, expected_orbits=(139, 140)
        )


def test_retained_sector_fails_on_unaccepted_deferred_state(tmp_path: Path) -> None:
    root = _retained_sector(tmp_path)
    path = next(root.glob("outputs/*/retained/*.retention.json"))
    payload = json.loads(path.read_text())
    payload["pdo_sector_accepted"] = True
    _write_json(path, payload)
    with pytest.raises(FM0ContractError, match="identity drifted"):
        verify_retained_sector_source(
            sector_root=root, sector=66, expected_orbits=(139, 140)
        )
