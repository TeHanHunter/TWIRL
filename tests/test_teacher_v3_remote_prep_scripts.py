from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script(relative: str, name: str):
    path = REPO_ROOT / relative
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_generic_orbit_contract_preserves_locked_s56() -> None:
    spoc = _load_script(
        "scripts/stage1_lightcurves/build_a2v1_spoc_quality_table.py",
        "build_a2v1_spoc_quality_table_test",
    )
    cadence = _load_script(
        "scripts/stage1_lightcurves/build_a2v1_cadence_reference.py",
        "build_a2v1_cadence_reference_test",
    )
    for module in (spoc, cadence):
        assert module._expected_orbits(sector=56, supplied=None) == (119, 120)
        assert module._expected_orbits(
            sector=57,
            supplied=[122, 121],
        ) == (121, 122)
        with pytest.raises(ValueError, match="required"):
            module._expected_orbits(sector=57, supplied=None)
        with pytest.raises(ValueError, match="locked"):
            module._expected_orbits(sector=56, supplied=[121, 122])


def test_generic_s57_quality_and_cadence_clis(tmp_path: Path) -> None:
    sector = 57
    orbits = (121, 122)
    quat_specs: list[str] = []
    qflag_specs: list[str] = []
    spoc_specs: list[str] = []
    camera_cadences: dict[int, list[int]] = {}

    for orbit_index, orbit in enumerate(orbits):
        for camera in range(1, 5):
            cadences = [
                1_000_000 + orbit_index * 1000 + camera * 10 + offset
                for offset in range(2)
            ]
            camera_cadences.setdefault(camera, []).extend(cadences)
            quat = tmp_path / f"o{orbit}_cam{camera}_quat.txt"
            quat.write_text(
                "cadence\n" + "".join(f"{value}\n" for value in cadences),
                encoding="utf-8",
            )
            quat_specs.append(f"{orbit},{camera},{quat}")
            for ccd in range(1, 5):
                qflag = tmp_path / f"o{orbit}_cam{camera}_ccd{ccd}_qflag.txt"
                qflag.write_text(
                    "".join(f"{value} 0\n" for value in cadences),
                    encoding="utf-8",
                )
                qflag_specs.append(f"{orbit},{camera},{ccd},{qflag}")

    for camera in range(1, 5):
        for ccd in range(1, 5):
            spoc = (
                tmp_path
                / f"spocffiflag_s{sector}_cam{camera}_ccd{ccd}.txt"
            )
            spoc.write_text(
                "".join(
                    f"{value},0\n" for value in camera_cadences[camera]
                ),
                encoding="utf-8",
            )
            spoc_specs.append(f"{camera},{ccd},{spoc}")

    quality_table = tmp_path / "spoc_quality.csv"
    quality_provenance = tmp_path / "spoc_quality.provenance.json"
    spoc_command = [
        sys.executable,
        str(
            REPO_ROOT
            / "scripts/stage1_lightcurves/build_a2v1_spoc_quality_table.py"
        ),
        "--sector",
        str(sector),
        "--expected-orbit",
        str(orbits[0]),
        "--expected-orbit",
        str(orbits[1]),
        "--output-table",
        str(quality_table),
        "--output-provenance",
        str(quality_provenance),
    ]
    for value in quat_specs:
        spoc_command.extend(["--quat", value])
    for value in spoc_specs:
        spoc_command.extend(["--spoc-flag", value])
    completed = subprocess.run(
        spoc_command,
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr

    cadence_table = tmp_path / "cadence_reference.csv"
    cadence_manifest = tmp_path / "cadence_reference.manifest.json"
    cadence_command = [
        sys.executable,
        str(
            REPO_ROOT
            / "scripts/stage1_lightcurves/build_a2v1_cadence_reference.py"
        ),
        "--sector",
        str(sector),
        "--expected-orbit",
        str(orbits[0]),
        "--expected-orbit",
        str(orbits[1]),
        "--spoc-quality-table",
        str(quality_table),
        "--spoc-quality-provenance",
        str(quality_provenance),
        "--output-table",
        str(cadence_table),
        "--output-manifest",
        str(cadence_manifest),
    ]
    for value in quat_specs:
        cadence_command.extend(["--quat", value])
    for value in qflag_specs:
        cadence_command.extend(["--qlp-qflag", value])
    completed = subprocess.run(
        cadence_command,
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    manifest = json.loads(cadence_manifest.read_text(encoding="utf-8"))
    assert manifest["contract_version"] == "s57_a2v1_cadence_reference_v1"
    assert manifest["orbits"] == [121, 122]
    assert manifest["n_rows"] == 64
    assert manifest["n_qlp_qflag_files_verified"] == 32
    assert manifest["n_spoc_authority_files_verified"] == 16


def test_teacher_v3_shell_wrappers_have_valid_bash_syntax() -> None:
    scripts = (
        "scripts/stage1_lightcurves/run_a2v1_cadence_reference_pdo.sh",
        "scripts/stage5_validation/run_teacher_v3_raw_export_pdo.sh",
        "scripts/orcd/slurm_teacher_v3_native_sector_cpu.sbatch",
        "scripts/orcd/slurm_teacher_v3_native_registry_cpu.sbatch",
        "scripts/orcd/slurm_teacher_v3_metadata_bootstrap_cpu.sbatch",
        "scripts/orcd/slurm_teacher_v3_h200.sbatch",
    )
    for relative in scripts:
        completed = subprocess.run(
            ["bash", "-n", str(REPO_ROOT / relative)],
            cwd=REPO_ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
        assert completed.returncode == 0, f"{relative}: {completed.stderr}"

    metadata_wrapper = (
        REPO_ROOT
        / "scripts/orcd/slurm_teacher_v3_metadata_bootstrap_cpu.sbatch"
    ).read_text(encoding="utf-8")
    assert "--metadata-only-bootstrap" in metadata_wrapper
    assert "--allow-cpu" in metadata_wrapper
    assert "--workers 0" in metadata_wrapper
