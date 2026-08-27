from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from twirl.models.fm0.mission_quality_admission import (
    MISSION_QUALITY_READY_STATE,
    mission_quality_type,
    verify_mission_quality_sources,
    write_mission_quality_receipt,
)
from twirl.models.fm0.registry import FM0ContractError


def _write_quat(path: Path, cadences: tuple[int, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("cadence,flag\n" + "".join(f"{value},0\n" for value in cadences))


def _write_flags(
    path: Path, cadences: tuple[int, ...], *, nonzero_value: int = 1
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "".join(
            f"{value}, {nonzero_value if index % 2 else 0}\n"
            for index, value in enumerate(cadences)
        )
    )


def _quality_tree(tmp_path: Path, *, sector: int) -> Path:
    orbits = (2 * sector + 7, 2 * sector + 8)
    provider = mission_quality_type(sector)
    for camera in range(1, 5):
        cadence_by_orbit = {
            orbits[0]: (1000 + camera * 10, 1001 + camera * 10),
            orbits[1]: (2000 + camera * 10, 2001 + camera * 10),
        }
        for orbit, cadences in cadence_by_orbit.items():
            run = tmp_path / f"orbit-{orbit}/ffi/run"
            _write_quat(run / f"cam{camera}_quat.txt", cadences)
            for ccd in range(1, 5):
                _write_flags(run / f"cam{camera}ccd{ccd}_qflag.txt", cadences)
        all_cadences = cadence_by_orbit[orbits[0]] + cadence_by_orbit[orbits[1]]
        prefix = "spocffiflag" if provider == "spoc" else "ticaffiflag"
        directory = "spocflags" if provider == "spoc" else "ticaflags"
        for ccd in range(1, 5):
            _write_flags(
                tmp_path / directory / f"{prefix}_s{sector}_cam{camera}_ccd{ccd}.txt",
                all_cadences,
                nonzero_value=4 if provider == "tica" else 1,
            )
    if provider == "tica":
        directory = tmp_path / "ticaflags"
        detectors = []
        for camera in range(1, 5):
            for ccd in range(1, 5):
                path = directory / f"ticaffiflag_s{sector}_cam{camera}_ccd{ccd}.txt"
                detectors.append(
                    {
                        "camera": camera,
                        "ccd": ccd,
                        "path": path.name,
                        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                        "n_rows": len(path.read_text(encoding="utf-8").splitlines()),
                        "n_nonzero": 2,
                    }
                )
        producer = "a" * 40
        summary = {
            "schema_version": "twirl_fm0_tica_quality_materialization_v1",
            "producer_git_sha": producer,
            "sector": sector,
            "mission_quality_type": "tica",
            "source": "qlp.lctools.bin.hlsp.query_ticaflags",
            "qlp_version": "0.14.6",
            "qlp_source_path": "/qlp/lctools/bin/hlsp.py",
            "qlp_source_sha256": "b" * 64,
            "n_detectors": 16,
            "n_rows": 64,
            "n_nonzero": 32,
            "cadence_coverage_verified": False,
            "detectors": detectors,
        }
        (directory / "summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        (directory / "READY").write_text(producer + "\n", encoding="utf-8")
    return tmp_path


@pytest.mark.parametrize(
    ("sector", "provider"), ((66, "spoc"), (67, "tica"))
)
def test_quality_gate_uses_uniform_sector_provider_policy(
    tmp_path: Path, sector: int, provider: str
) -> None:
    root = _quality_tree(tmp_path, sector=sector)
    receipt = verify_mission_quality_sources(
        sector=sector,
        expected_orbits=(2 * sector + 7, 2 * sector + 8),
        qlp_root=root,
    )
    assert receipt["mission_quality_type"] == provider
    assert receipt["quality_state"] == MISSION_QUALITY_READY_STATE
    assert receipt["n_qflag_files"] == 32
    assert receipt["n_mission_quality_files"] == 16
    assert receipt["n_mission_quality_rows_excluded_by_quat"] == 0
    assert receipt["hdf5_quality_join_verified"] is False
    assert receipt["panel_admission_authorized"] is False


def test_quality_gate_rejects_missing_mission_cadence(tmp_path: Path) -> None:
    root = _quality_tree(tmp_path, sector=66)
    path = root / "spocflags/spocffiflag_s66_cam4_ccd4.txt"
    path.write_text("2040, 0\n2041, 0\n")
    with pytest.raises(FM0ContractError, match="missing quaternion cadences"):
        verify_mission_quality_sources(
            sector=66,
            expected_orbits=(139, 140),
            qlp_root=root,
        )


def test_quality_gate_records_mission_rows_excluded_by_quat(
    tmp_path: Path,
) -> None:
    root = _quality_tree(tmp_path, sector=67)
    path = root / "ticaflags/ticaffiflag_s67_cam1_ccd1.txt"
    with path.open("a", encoding="utf-8") as handle:
        handle.write("9998, 0\n9999, 4\n")
    summary_path = root / "ticaflags/summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    entry = summary["detectors"][0]
    entry["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    entry["n_rows"] = 6
    entry["n_nonzero"] = 3
    summary["n_rows"] = 66
    summary["n_nonzero"] = 33
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    receipt = verify_mission_quality_sources(
        sector=67,
        expected_orbits=(141, 142),
        qlp_root=root,
    )
    assert receipt["n_mission_quality_rows_excluded_by_quat"] == 2
    assert receipt["n_mission_quality_rows_retained"] == 64
    assert receipt["n_nonzero_mission_quality_rows_retained"] == 32
    assert receipt["mission_quality_authority_exclusions"]["by_detector"][
        "cam1_ccd1"
    ] == {
        "n_rows": 2,
        "rows": [
            {"cadence": 9998, "mission_quality": 0},
            {"cadence": 9999, "mission_quality": 4},
        ],
    }
    assert len(receipt["mission_quality_authority_exclusions_sha256"]) == 64


def test_quality_gate_rejects_extra_qlp_qflag_cadence(tmp_path: Path) -> None:
    root = _quality_tree(tmp_path, sector=67)
    path = root / "orbit-141/ffi/run/cam1ccd1_qflag.txt"
    with path.open("a", encoding="utf-8") as handle:
        handle.write("9999, 0\n")
    with pytest.raises(FM0ContractError, match="QLP qflag cadence coverage differs"):
        verify_mission_quality_sources(
            sector=67,
            expected_orbits=(141, 142),
            qlp_root=root,
        )


def test_quality_gate_rejects_tica_bits_outside_current_qlp_mask(
    tmp_path: Path,
) -> None:
    root = _quality_tree(tmp_path, sector=67)
    path = root / "ticaflags/ticaffiflag_s67_cam1_ccd1.txt"
    text = path.read_text(encoding="utf-8").replace("1010, 0", "1010, 1")
    path.write_text(text, encoding="utf-8")
    summary_path = root / "ticaflags/summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["detectors"][0]["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    with pytest.raises(FM0ContractError, match="three-bit authority"):
        verify_mission_quality_sources(
            sector=67,
            expected_orbits=(141, 142),
            qlp_root=root,
        )


def test_quality_gate_rejects_missing_tica_materialization_summary(
    tmp_path: Path,
) -> None:
    root = _quality_tree(tmp_path, sector=67)
    (root / "ticaflags/summary.json").unlink()
    with pytest.raises(FM0ContractError, match="lacks summary.json"):
        verify_mission_quality_sources(
            sector=67,
            expected_orbits=(141, 142),
            qlp_root=root,
        )


def test_quality_gate_rejects_wrong_orbits_and_overwrite(tmp_path: Path) -> None:
    root = _quality_tree(tmp_path, sector=67)
    with pytest.raises(FM0ContractError, match="expected orbits"):
        verify_mission_quality_sources(
            sector=67,
            expected_orbits=(140, 141),
            qlp_root=root,
        )
    receipt = verify_mission_quality_sources(
        sector=67,
        expected_orbits=(141, 142),
        qlp_root=root,
    )
    output = tmp_path / "receipt.json"
    write_mission_quality_receipt(receipt, output)
    assert json.loads(output.read_text())["passed"] is True
    with pytest.raises(FM0ContractError, match="refusing to overwrite"):
        write_mission_quality_receipt(receipt, output)
