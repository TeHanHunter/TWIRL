from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from twirl.models.fm0.mission_quality_admission import (
    MISSION_QUALITY_READY_STATE,
    MISSION_QUALITY_RECEIPT_SCHEMA_VERSION,
)
from twirl.models.fm0.mission_quality_transfer import (
    MISSION_QUALITY_TRANSFER_SCHEMA_VERSION,
    package_mission_quality_transfer,
)
from twirl.models.fm0.registry import FM0ContractError


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _receipt(tmp_path: Path, *, sector: int = 67) -> Path:
    quat = tmp_path / f"orbit-{2 * sector + 7}/ffi/run/cam1_quat.txt"
    qflag = tmp_path / f"orbit-{2 * sector + 7}/ffi/run/cam1ccd1_qflag.txt"
    mission = tmp_path / f"ticaffiflag_s{sector}_cam1_ccd1.txt"
    for path, content in (
        (quat, "cadence,flag\n100,0\n"),
        (qflag, "100, 0\n"),
        (mission, "100, 4\n"),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
    receipt = {
        "schema_version": MISSION_QUALITY_RECEIPT_SCHEMA_VERSION,
        "sector": sector,
        "quality_state": MISSION_QUALITY_READY_STATE,
        "passed": True,
        "panel_admission_authorized": False,
        "source_bindings": [
            {
                "role": "qlp_quaternion",
                "orbit": 2 * sector + 7,
                "camera": 1,
                "path": str(quat),
                "sha256": _sha(quat),
            },
            {
                "role": "qlp_detector_qflag",
                "orbit": 2 * sector + 7,
                "camera": 1,
                "ccd": 1,
                "path": str(qflag),
                "sha256": _sha(qflag),
            },
            {
                "role": "tica_mission_quality",
                "camera": 1,
                "ccd": 1,
                "path": str(mission),
                "sha256": _sha(mission),
            },
        ],
    }
    path = tmp_path / f"s{sector}.receipt.json"
    path.write_text(json.dumps(receipt) + "\n", encoding="utf-8")
    return path


def test_package_mission_quality_transfer_copies_bound_sources(tmp_path: Path) -> None:
    receipt = _receipt(tmp_path)
    output = tmp_path / "release"
    manifest = package_mission_quality_transfer(
        receipt_paths=[receipt],
        output_dir=output,
        producer_git_sha="a" * 40,
    )
    assert manifest["schema_version"] == MISSION_QUALITY_TRANSFER_SCHEMA_VERSION
    assert manifest["sectors"] == [67]
    assert manifest["n_source_files"] == 3
    assert (output / "receipts/s0067.mission_quality.json").is_file()
    assert (output / "sources/s0067/mission/ticaffiflag_s67_cam1_ccd1.txt").is_file()
    assert (output / "sources/s0067/qlp/orbit-141/ffi/run/cam1_quat.txt").is_file()


def test_package_mission_quality_transfer_rejects_source_drift(
    tmp_path: Path,
) -> None:
    receipt = _receipt(tmp_path)
    payload = json.loads(receipt.read_text(encoding="utf-8"))
    Path(payload["source_bindings"][0]["path"]).write_text(
        "cadence,flag\n101,0\n", encoding="utf-8"
    )
    with pytest.raises(FM0ContractError, match="hash-mismatched"):
        package_mission_quality_transfer(
            receipt_paths=[receipt],
            output_dir=tmp_path / "release",
            producer_git_sha="a" * 40,
        )
    assert not (tmp_path / "release.partial").exists()


def test_package_mission_quality_transfer_refuses_overwrite(tmp_path: Path) -> None:
    receipt = _receipt(tmp_path)
    output = tmp_path / "release"
    output.mkdir()
    with pytest.raises(FM0ContractError, match="refusing to overwrite"):
        package_mission_quality_transfer(
            receipt_paths=[receipt],
            output_dir=output,
            producer_git_sha="a" * 40,
        )
