from __future__ import annotations

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


def _write_flags(path: Path, cadences: tuple[int, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{value}, {index % 2}\n" for index, value in enumerate(cadences)))


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
            )
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
    assert receipt["hdf5_quality_join_verified"] is False
    assert receipt["panel_admission_authorized"] is False


def test_quality_gate_rejects_missing_mission_cadence(tmp_path: Path) -> None:
    root = _quality_tree(tmp_path, sector=66)
    path = root / "spocflags/spocffiflag_s66_cam4_ccd4.txt"
    path.write_text("2040, 0\n2041, 0\n")
    with pytest.raises(FM0ContractError, match="cadence coverage differs"):
        verify_mission_quality_sources(
            sector=66,
            expected_orbits=(139, 140),
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
