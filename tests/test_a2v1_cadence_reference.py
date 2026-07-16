from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest

from twirl.lightcurves.a2v1_cadence_reference import (
    CADENCE_REFERENCE_COLUMNS,
    QuatSource,
    SPOC_QUALITY_DERIVATION_CONTRACT,
    build_cadence_reference,
    file_sha256,
    write_cadence_reference,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


def _write_quat(path: Path, cadences: list[int]) -> None:
    lines = ["q1_mean,flag,cadence"]
    lines.extend(f"0.0,0,{cadence}" for cadence in cadences)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _make_inputs(
    root: Path,
    *,
    cameras: tuple[int, ...],
    detectors: tuple[tuple[int, int], ...],
    include_orbitid: bool = True,
) -> tuple[tuple[QuatSource, ...], Path, Path, dict[tuple[int, int], int]]:
    cadence_orbits: dict[tuple[int, int], int] = {}
    sources: list[QuatSource] = []
    for orbit_index, orbitid in enumerate((119, 120)):
        for camera in cameras:
            cadences = [
                100_000 + 1000 * orbit_index + 10 * camera,
                100_001 + 1000 * orbit_index + 10 * camera,
            ]
            path = root / f"o{orbitid}_cam{camera}_quat.txt"
            _write_quat(path, cadences)
            sources.append(QuatSource(orbitid=orbitid, camera=camera, path=path))
            cadence_orbits.update({(camera, cadence): orbitid for cadence in cadences})

    rows: list[dict[str, int]] = []
    for camera, ccd in detectors:
        for (candidate_camera, cadence), orbitid in sorted(cadence_orbits.items()):
            if candidate_camera != camera:
                continue
            row = {
                "sector": 56,
                "camera": camera,
                "ccd": ccd,
                "cadenceno": cadence,
                "quality": 8 if cadence % 2 else 0,
            }
            if include_orbitid:
                row["orbitid"] = orbitid
            rows.append(row)
    quality_path = root / "s56_spoc_quality.csv"
    quality = pd.DataFrame(rows)
    quality.to_csv(quality_path, index=False)
    flag_sources: list[dict[str, object]] = []
    for camera, ccd in detectors:
        detector = quality.loc[
            (quality["camera"] == camera) & (quality["ccd"] == ccd)
        ]
        flag_path = root / f"spocffiflag_s56_cam{camera}_ccd{ccd}.txt"
        lines = [
            f"{int(row.cadenceno)}.000000,{int(row.quality)}"
            for row in detector.itertuples(index=False)
        ]
        if camera == 3:
            # S56 cam3 has one SPOC cadence intentionally absent from the
            # authoritative orbit-120 quaternion table.
            lines.append("199999.000000,0")
        flag_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        flag_sources.append(
            {
                "camera": camera,
                "ccd": ccd,
                "path": str(flag_path),
                "sha256": file_sha256(flag_path),
            }
        )
    provenance_path = root / "s56_spoc_quality.provenance.json"
    provenance = {
        "contract_version": "s56_spoc_quality_table_provenance_v1",
        "derivation_contract": SPOC_QUALITY_DERIVATION_CONTRACT,
        "sector": 56,
        "quality_authority": "spoc_quality_flags",
        "table_sha256": file_sha256(quality_path),
        "n_rows": len(quality),
        "detectors": [f"cam{camera}_ccd{ccd}" for camera, ccd in detectors],
        "source_flag_files": flag_sources,
    }
    provenance_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return tuple(sources), quality_path, provenance_path, cadence_orbits


def _refresh_table_binding(provenance_path: Path, quality_path: Path) -> None:
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    provenance["table_sha256"] = file_sha256(quality_path)
    provenance["n_rows"] = len(pd.read_csv(quality_path))
    provenance_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def test_build_and_atomically_write_hash_bound_reference(tmp_path: Path) -> None:
    detectors = ((1, 1), (1, 2), (2, 3))
    sources, quality_path, provenance_path, cadence_orbits = _make_inputs(
        tmp_path, cameras=(1, 2), detectors=detectors
    )
    table_path = tmp_path / "s56_cadence_reference.csv"
    manifest_path = tmp_path / "s56_cadence_reference.json"

    manifest = write_cadence_reference(
        sector=56,
        quat_sources=sources,
        spoc_quality_table=quality_path,
        spoc_quality_provenance=provenance_path,
        output_table=table_path,
        output_manifest=manifest_path,
        expected_orbits=(119, 120),
        expected_detectors=detectors,
    )

    table = pd.read_csv(table_path)
    assert tuple(table.columns) == CADENCE_REFERENCE_COLUMNS
    assert len(table) == 12
    assert {
        (int(row.camera), int(row.cadenceno)): int(row.orbitid)
        for row in table.itertuples(index=False)
    } == cadence_orbits
    assert manifest["contract_version"] == "s56_a2v1_cadence_reference_v1"
    assert manifest["cadence_authority"] == "qlp_cam_quat"
    assert manifest["quality_authority"] == "spoc_quality_flags"
    assert manifest["table_sha256"] == file_sha256(table_path)
    assert manifest["n_rows"] == len(table)
    assert manifest["detectors"] == ["cam1_ccd1", "cam1_ccd2", "cam2_ccd3"]
    assert manifest["orbits"] == [119, 120]
    assert manifest["n_spoc_authority_files_verified"] == len(detectors)
    assert manifest["spoc_quality_provenance_sha256"] == file_sha256(
        provenance_path
    )
    assert len(manifest["source_file_sha256"]) == (
        len(sources) + len(detectors) + 2
    )
    assert json.loads(manifest_path.read_text(encoding="utf-8")) == manifest
    assert not list(tmp_path.glob(".s56_cadence_reference.*"))

    with pytest.raises(FileExistsError, match="output already exists"):
        write_cadence_reference(
            sector=56,
            quat_sources=sources,
            spoc_quality_table=quality_path,
            spoc_quality_provenance=provenance_path,
            output_table=table_path,
            output_manifest=manifest_path,
            expected_orbits=(119, 120),
            expected_detectors=detectors,
        )


def test_builder_rejects_nonexact_spoc_cadence_join(tmp_path: Path) -> None:
    detectors = ((1, 1),)
    sources, quality_path, provenance_path, _ = _make_inputs(
        tmp_path, cameras=(1,), detectors=detectors, include_orbitid=False
    )
    quality = pd.read_csv(quality_path)
    quality.iloc[:-1].to_csv(quality_path, index=False)
    _refresh_table_binding(provenance_path, quality_path)

    with pytest.raises(ValueError, match="join is not exact"):
        build_cadence_reference(
            sector=56,
            quat_sources=sources,
            spoc_quality_table=quality_path,
            spoc_quality_provenance=provenance_path,
            expected_orbits=(119, 120),
            expected_detectors=detectors,
        )


def test_builder_rejects_missing_spoc_provenance(tmp_path: Path) -> None:
    detectors = ((1, 1),)
    sources, quality_path, _, _ = _make_inputs(
        tmp_path, cameras=(1,), detectors=detectors
    )

    with pytest.raises(FileNotFoundError):
        build_cadence_reference(
            sector=56,
            quat_sources=sources,
            spoc_quality_table=quality_path,
            spoc_quality_provenance=tmp_path / "missing.json",
            expected_orbits=(119, 120),
            expected_detectors=detectors,
        )


def test_builder_rejects_forged_spoc_source_hash(tmp_path: Path) -> None:
    detectors = ((1, 1),)
    sources, quality_path, provenance_path, _ = _make_inputs(
        tmp_path, cameras=(1,), detectors=detectors
    )
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    provenance["source_flag_files"][0]["sha256"] = "f" * 64
    provenance_path.write_text(json.dumps(provenance), encoding="utf-8")

    with pytest.raises(ValueError, match="authority file hash mismatch"):
        build_cadence_reference(
            sector=56,
            quat_sources=sources,
            spoc_quality_table=quality_path,
            spoc_quality_provenance=provenance_path,
            expected_orbits=(119, 120),
            expected_detectors=detectors,
        )


def test_builder_rejects_incomplete_spoc_source_inventory(tmp_path: Path) -> None:
    detectors = ((1, 1), (1, 2))
    sources, quality_path, provenance_path, _ = _make_inputs(
        tmp_path, cameras=(1,), detectors=detectors
    )
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    provenance["source_flag_files"] = provenance["source_flag_files"][:-1]
    provenance_path.write_text(json.dumps(provenance), encoding="utf-8")

    with pytest.raises(ValueError, match="source inventory is incomplete"):
        build_cadence_reference(
            sector=56,
            quat_sources=sources,
            spoc_quality_table=quality_path,
            spoc_quality_provenance=provenance_path,
            expected_orbits=(119, 120),
            expected_detectors=detectors,
        )


def test_builder_rejects_rehashed_table_that_disagrees_with_raw_spoc(
    tmp_path: Path,
) -> None:
    detectors = ((1, 1),)
    sources, quality_path, provenance_path, _ = _make_inputs(
        tmp_path, cameras=(1,), detectors=detectors
    )
    quality = pd.read_csv(quality_path)
    quality.loc[0, "quality"] = int(quality.loc[0, "quality"]) ^ 8
    quality.to_csv(quality_path, index=False)
    _refresh_table_binding(provenance_path, quality_path)

    with pytest.raises(ValueError, match="does not match its original authority"):
        build_cadence_reference(
            sector=56,
            quat_sources=sources,
            spoc_quality_table=quality_path,
            spoc_quality_provenance=provenance_path,
            expected_orbits=(119, 120),
            expected_detectors=detectors,
        )


def test_builder_rejects_duplicate_orbit_assignment(tmp_path: Path) -> None:
    detectors = ((1, 1),)
    sources, quality_path, provenance_path, _ = _make_inputs(
        tmp_path, cameras=(1,), detectors=detectors
    )
    _write_quat(sources[1].path, [100_010, 101_011])

    with pytest.raises(ValueError, match="multiple orbits"):
        build_cadence_reference(
            sector=56,
            quat_sources=sources,
            spoc_quality_table=quality_path,
            spoc_quality_provenance=provenance_path,
            expected_orbits=(119, 120),
            expected_detectors=detectors,
        )


def test_s56_cli_builds_full_detector_reference(tmp_path: Path) -> None:
    detectors = tuple(
        (camera, ccd) for camera in range(1, 5) for ccd in range(1, 5)
    )
    sources, quality_path, provenance_path, _ = _make_inputs(
        tmp_path, cameras=(1, 2, 3, 4), detectors=detectors
    )
    table_path = tmp_path / "reference.csv"
    manifest_path = tmp_path / "reference.json"
    command = [
        sys.executable,
        str(
            REPO_ROOT
            / "scripts"
            / "stage1_lightcurves"
            / "build_a2v1_cadence_reference.py"
        ),
        "--spoc-quality-table",
        str(quality_path),
        "--spoc-quality-provenance",
        str(provenance_path),
        "--output-table",
        str(table_path),
        "--output-manifest",
        str(manifest_path),
    ]
    for source in sources:
        command.extend(
            ["--quat", f"{source.orbitid},{source.camera},{source.path}"]
        )

    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["n_rows"] == 64
    assert len(manifest["detectors"]) == 16
    assert manifest["n_spoc_rows_excluded_by_quat"] == 4
    assert file_sha256(table_path) == manifest["table_sha256"]
