from __future__ import annotations

import json
from pathlib import Path
import stat
import subprocess
import sys

import pandas as pd
import pytest

import twirl.lightcurves.a2v1_cadence_reference as cadence_reference
from twirl.lightcurves.a2v1_cadence_reference import (
    CADENCE_REFERENCE_COLUMNS,
    QLP_QUALITY_EXTERNAL_BIT,
    SPOC_QUALITY_TABLE_COLUMNS,
    QlpQflagSource,
    QuatSource,
    SPOC_QUALITY_DERIVATION_CONTRACT,
    SpocFlagInput,
    build_cadence_reference,
    file_sha256,
    load_spoc_quality_provenance,
    write_cadence_reference,
    write_spoc_quality_table,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


def _write_quat(path: Path, cadences: list[int]) -> None:
    lines = ["q1_mean,flag,cadence"]
    lines.extend(f"0.0,0,{cadence}" for cadence in cadences)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _make_qlp_qflags(
    root: Path,
    *,
    quat_sources: tuple[QuatSource, ...],
    detectors: tuple[tuple[int, int], ...],
) -> tuple[QlpQflagSource, ...]:
    sources: list[QlpQflagSource] = []
    for quat in quat_sources:
        cadences = cadence_reference.read_qlp_quat_cadences(quat.path)
        for camera, ccd in detectors:
            if camera != quat.camera:
                continue
            path = (
                root
                / f"o{quat.orbitid}_cam{quat.camera}_ccd{ccd}_qflag.txt"
            )
            rows = [
                f"{int(cadence)} {1 if int(cadence) % 2 else 0}"
                for cadence in cadences
            ]
            path.write_text("\n".join(rows) + "\n", encoding="utf-8")
            sources.append(
                QlpQflagSource(
                    orbitid=quat.orbitid,
                    camera=quat.camera,
                    ccd=ccd,
                    path=path,
                )
            )
    return tuple(sources)


def _make_raw_authorities(
    root: Path,
    *,
    cameras: tuple[int, ...],
    detectors: tuple[tuple[int, int], ...],
) -> tuple[
    tuple[QuatSource, ...],
    tuple[SpocFlagInput, ...],
    dict[tuple[int, int], int],
]:
    cadence_orbits: dict[tuple[int, int], int] = {}
    quat_sources: list[QuatSource] = []
    for orbit_index, orbitid in enumerate((119, 120)):
        for camera in cameras:
            cadences = [
                100_000 + 1000 * orbit_index + 10 * camera,
                100_001 + 1000 * orbit_index + 10 * camera,
            ]
            quat_path = root / f"o{orbitid}_cam{camera}_quat.txt"
            _write_quat(quat_path, cadences)
            quat_sources.append(
                QuatSource(orbitid=orbitid, camera=camera, path=quat_path)
            )
            cadence_orbits.update(
                {(camera, cadence): orbitid for cadence in cadences}
            )

    spoc_sources: list[SpocFlagInput] = []
    for camera, ccd in detectors:
        rows = [
            f"{cadence}.000000,{8 if cadence % 2 else 0}"
            for candidate_camera, cadence in sorted(cadence_orbits)
            if candidate_camera == camera
        ]
        if camera == 3:
            rows.append(f"{199_900 + ccd}.000000,{ccd}")
        path = root / f"spocffiflag_s56_cam{camera}_ccd{ccd}.txt"
        path.write_text("\n".join(rows) + "\n", encoding="utf-8")
        spoc_sources.append(SpocFlagInput(camera=camera, ccd=ccd, path=path))
    return tuple(quat_sources), tuple(spoc_sources), cadence_orbits


def _make_inputs(
    root: Path,
    *,
    cameras: tuple[int, ...],
    detectors: tuple[tuple[int, int], ...],
    include_orbitid: bool = True,
) -> tuple[
    tuple[QuatSource, ...],
    tuple[QlpQflagSource, ...],
    Path,
    Path,
    dict[tuple[int, int], int],
]:
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
    qflags = _make_qlp_qflags(
        root,
        quat_sources=tuple(sources),
        detectors=detectors,
    )
    return tuple(sources), qflags, quality_path, provenance_path, cadence_orbits


def _refresh_table_binding(provenance_path: Path, quality_path: Path) -> None:
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    provenance["table_sha256"] = file_sha256(quality_path)
    provenance["n_rows"] = len(pd.read_csv(quality_path))
    provenance_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def test_spoc_builder_writes_loader_compatible_atomic_pair(tmp_path: Path) -> None:
    detectors = tuple(
        (camera, ccd) for camera in range(1, 5) for ccd in range(1, 5)
    )
    quats, flags, cadence_orbits = _make_raw_authorities(
        tmp_path, cameras=(1, 2, 3, 4), detectors=detectors
    )
    table_path = tmp_path / "s56_spoc_quality.csv"
    provenance_path = tmp_path / "s56_spoc_quality.provenance.json"

    manifest = write_spoc_quality_table(
        sector=56,
        quat_sources=quats,
        spoc_flag_sources=flags,
        output_table=table_path,
        output_provenance=provenance_path,
        expected_orbits=(119, 120),
        expected_detectors=detectors,
    )

    table = pd.read_csv(table_path)
    assert tuple(table.columns) == SPOC_QUALITY_TABLE_COLUMNS
    assert len(table) == 64
    assert manifest["contract_version"] == "s56_spoc_quality_table_provenance_v1"
    assert manifest["derivation_contract"] == SPOC_QUALITY_DERIVATION_CONTRACT
    assert manifest["table_sha256"] == file_sha256(table_path)
    assert manifest["n_rows"] == len(table)
    assert manifest["n_rows_raw"] == 68
    assert manifest["n_rows_retained"] == 64
    assert manifest["n_rows_excluded_by_quat"] == 4
    assert manifest["n_quaternion_files"] == 8
    assert manifest["n_spoc_flag_files"] == 16
    assert len(manifest["source_file_sha256"]) == 24
    assert manifest["exclusions"]["n_rows"] == 4
    assert manifest["exclusions"]["by_detector"]["cam3_ccd2"]["rows"] == [
        {"cadenceno": 199902, "quality": 2}
    ]
    assert {
        (int(row.camera), int(row.cadenceno)): int(row.orbitid)
        for row in table.itertuples(index=False)
    } == cadence_orbits

    payload, verified_sources = load_spoc_quality_provenance(
        provenance_path,
        spoc_quality_table=table_path,
        sector=56,
        expected_detectors=detectors,
    )
    assert payload == manifest
    assert len(verified_sources) == 16
    qflags = _make_qlp_qflags(
        tmp_path,
        quat_sources=quats,
        detectors=detectors,
    )
    rebuilt = build_cadence_reference(
        sector=56,
        quat_sources=quats,
        qlp_qflag_sources=qflags,
        spoc_quality_table=table_path,
        spoc_quality_provenance=provenance_path,
        expected_orbits=(119, 120),
        expected_detectors=detectors,
    )
    assert tuple(rebuilt.columns) == CADENCE_REFERENCE_COLUMNS
    pd.testing.assert_series_equal(
        rebuilt["spoc_quality"], table["quality"], check_names=False
    )
    assert not list(tmp_path.glob(".s56_spoc_quality.*"))

    with pytest.raises(FileExistsError, match="output already exists"):
        write_spoc_quality_table(
            sector=56,
            quat_sources=quats,
            spoc_flag_sources=flags,
            output_table=table_path,
            output_provenance=provenance_path,
            expected_orbits=(119, 120),
            expected_detectors=detectors,
        )


def test_spoc_builder_fails_closed_on_missing_quaternion_cadence(
    tmp_path: Path,
) -> None:
    detectors = ((1, 1),)
    quats, flags, _ = _make_raw_authorities(
        tmp_path, cameras=(1,), detectors=detectors
    )
    raw_rows = flags[0].path.read_text(encoding="utf-8").splitlines()
    flags[0].path.write_text("\n".join(raw_rows[1:]) + "\n", encoding="utf-8")
    table_path = tmp_path / "quality.csv"
    provenance_path = tmp_path / "quality.json"

    with pytest.raises(ValueError, match="missing quaternion cadences"):
        write_spoc_quality_table(
            sector=56,
            quat_sources=quats,
            spoc_flag_sources=flags,
            output_table=table_path,
            output_provenance=provenance_path,
            expected_orbits=(119, 120),
            expected_detectors=detectors,
        )
    assert not table_path.exists()
    assert not provenance_path.exists()


def test_spoc_builder_does_not_publish_when_manifest_write_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    detectors = ((1, 1),)
    quats, flags, _ = _make_raw_authorities(
        tmp_path, cameras=(1,), detectors=detectors
    )
    table_path = tmp_path / "quality.csv"
    provenance_path = tmp_path / "quality.json"

    def fail_write_json(payload: dict[str, object], path: Path) -> None:
        raise OSError("synthetic manifest failure")

    monkeypatch.setattr(cadence_reference, "_write_json", fail_write_json)
    with pytest.raises(OSError, match="synthetic manifest failure"):
        write_spoc_quality_table(
            sector=56,
            quat_sources=quats,
            spoc_flag_sources=flags,
            output_table=table_path,
            output_provenance=provenance_path,
            expected_orbits=(119, 120),
            expected_detectors=detectors,
        )
    assert not table_path.exists()
    assert not provenance_path.exists()
    assert not list(tmp_path.glob(".quality.*"))


def test_atomic_pair_publish_restores_existing_outputs_on_second_replace_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_table = tmp_path / "quality.csv"
    output_manifest = tmp_path / "quality.json"
    output_table.write_text("old table\n", encoding="utf-8")
    output_manifest.write_text("old manifest\n", encoding="utf-8")
    table_temporary = tmp_path / ".new-table.csv"
    manifest_temporary = tmp_path / ".new-manifest.json"
    table_temporary.write_text("new table\n", encoding="utf-8")
    manifest_temporary.write_text("new manifest\n", encoding="utf-8")
    real_replace = cadence_reference.os.replace

    def fail_second_replace(source: Path, destination: Path) -> None:
        if Path(source) == manifest_temporary and Path(destination) == output_manifest:
            raise OSError("synthetic second replace failure")
        real_replace(source, destination)

    monkeypatch.setattr(cadence_reference.os, "replace", fail_second_replace)
    with pytest.raises(OSError, match="synthetic second replace failure"):
        cadence_reference._publish_pair_with_rollback(
            table_temporary=table_temporary,
            manifest_temporary=manifest_temporary,
            output_table=output_table,
            output_manifest=output_manifest,
        )
    assert output_table.read_text(encoding="utf-8") == "old table\n"
    assert output_manifest.read_text(encoding="utf-8") == "old manifest\n"


def test_build_and_publish_hash_bound_reference_with_stable_permissions(
    tmp_path: Path,
) -> None:
    detectors = ((1, 1), (1, 2), (2, 3))
    sources, qflags, quality_path, provenance_path, cadence_orbits = _make_inputs(
        tmp_path, cameras=(1, 2), detectors=detectors
    )
    table_path = tmp_path / "s56_cadence_reference.csv"
    manifest_path = tmp_path / "s56_cadence_reference.json"

    manifest = write_cadence_reference(
        sector=56,
        quat_sources=sources,
        qlp_qflag_sources=qflags,
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
    assert manifest["quality_authority"] == "spoc_and_qlp_quality_flags"
    assert manifest["quality_composition"] == {
        "external_quality": "spoc_quality | (qlp_quality << 30)",
        "qlp_quality_raw_values": [0, 1],
        "qlp_quality_external_bit": QLP_QUALITY_EXTERNAL_BIT,
    }
    assert manifest["table_sha256"] == file_sha256(table_path)
    assert manifest["n_rows"] == len(table)
    assert manifest["detectors"] == ["cam1_ccd1", "cam1_ccd2", "cam2_ccd3"]
    assert manifest["orbits"] == [119, 120]
    assert manifest["n_spoc_authority_files_verified"] == len(detectors)
    assert manifest["n_qlp_qflag_files_verified"] == 2 * len(detectors)
    assert manifest["spoc_quality_provenance_sha256"] == file_sha256(
        provenance_path
    )
    assert len(manifest["source_file_sha256"]) == (
        len(sources) + len(qflags) + len(detectors) + 2
    )
    expected_external = pd.Series(
        table["spoc_quality"].to_numpy(dtype="int64")
        | (
            table["qlp_quality"].to_numpy(dtype="int64")
            << QLP_QUALITY_EXTERNAL_BIT
        )
    )
    pd.testing.assert_series_equal(
        table["external_quality"], expected_external, check_names=False
    )
    assert json.loads(manifest_path.read_text(encoding="utf-8")) == manifest
    assert stat.S_IMODE(table_path.stat().st_mode) == 0o640
    assert stat.S_IMODE(manifest_path.stat().st_mode) == 0o640
    assert not list(tmp_path.glob(".s56_cadence_reference.*"))

    with pytest.raises(FileExistsError, match="output already exists"):
        write_cadence_reference(
            sector=56,
            quat_sources=sources,
            qlp_qflag_sources=qflags,
            spoc_quality_table=quality_path,
            spoc_quality_provenance=provenance_path,
            output_table=table_path,
            output_manifest=manifest_path,
            expected_orbits=(119, 120),
            expected_detectors=detectors,
        )


def test_builder_rejects_nonexact_spoc_cadence_join(tmp_path: Path) -> None:
    detectors = ((1, 1),)
    sources, qflags, quality_path, provenance_path, _ = _make_inputs(
        tmp_path, cameras=(1,), detectors=detectors, include_orbitid=False
    )
    quality = pd.read_csv(quality_path)
    quality.iloc[:-1].to_csv(quality_path, index=False)
    _refresh_table_binding(provenance_path, quality_path)

    with pytest.raises(ValueError, match="join is not exact"):
        build_cadence_reference(
            sector=56,
            quat_sources=sources,
            qlp_qflag_sources=qflags,
            spoc_quality_table=quality_path,
            spoc_quality_provenance=provenance_path,
            expected_orbits=(119, 120),
            expected_detectors=detectors,
        )


def test_builder_rejects_missing_spoc_provenance(tmp_path: Path) -> None:
    detectors = ((1, 1),)
    sources, qflags, quality_path, _, _ = _make_inputs(
        tmp_path, cameras=(1,), detectors=detectors
    )

    with pytest.raises(FileNotFoundError):
        build_cadence_reference(
            sector=56,
            quat_sources=sources,
            qlp_qflag_sources=qflags,
            spoc_quality_table=quality_path,
            spoc_quality_provenance=tmp_path / "missing.json",
            expected_orbits=(119, 120),
            expected_detectors=detectors,
        )


def test_builder_rejects_forged_spoc_source_hash(tmp_path: Path) -> None:
    detectors = ((1, 1),)
    sources, qflags, quality_path, provenance_path, _ = _make_inputs(
        tmp_path, cameras=(1,), detectors=detectors
    )
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    provenance["source_flag_files"][0]["sha256"] = "f" * 64
    provenance_path.write_text(json.dumps(provenance), encoding="utf-8")

    with pytest.raises(ValueError, match="authority file hash mismatch"):
        build_cadence_reference(
            sector=56,
            quat_sources=sources,
            qlp_qflag_sources=qflags,
            spoc_quality_table=quality_path,
            spoc_quality_provenance=provenance_path,
            expected_orbits=(119, 120),
            expected_detectors=detectors,
        )


def test_builder_rejects_incomplete_spoc_source_inventory(tmp_path: Path) -> None:
    detectors = ((1, 1), (1, 2))
    sources, qflags, quality_path, provenance_path, _ = _make_inputs(
        tmp_path, cameras=(1,), detectors=detectors
    )
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    provenance["source_flag_files"] = provenance["source_flag_files"][:-1]
    provenance_path.write_text(json.dumps(provenance), encoding="utf-8")

    with pytest.raises(ValueError, match="source inventory is incomplete"):
        build_cadence_reference(
            sector=56,
            quat_sources=sources,
            qlp_qflag_sources=qflags,
            spoc_quality_table=quality_path,
            spoc_quality_provenance=provenance_path,
            expected_orbits=(119, 120),
            expected_detectors=detectors,
        )


def test_builder_rejects_rehashed_table_that_disagrees_with_raw_spoc(
    tmp_path: Path,
) -> None:
    detectors = ((1, 1),)
    sources, qflags, quality_path, provenance_path, _ = _make_inputs(
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
            qlp_qflag_sources=qflags,
            spoc_quality_table=quality_path,
            spoc_quality_provenance=provenance_path,
            expected_orbits=(119, 120),
            expected_detectors=detectors,
        )


def test_builder_rejects_duplicate_orbit_assignment(tmp_path: Path) -> None:
    detectors = ((1, 1),)
    sources, qflags, quality_path, provenance_path, _ = _make_inputs(
        tmp_path, cameras=(1,), detectors=detectors
    )
    _write_quat(sources[1].path, [100_010, 101_011])

    with pytest.raises(ValueError, match="multiple orbits"):
        build_cadence_reference(
            sector=56,
            quat_sources=sources,
            qlp_qflag_sources=qflags,
            spoc_quality_table=quality_path,
            spoc_quality_provenance=provenance_path,
            expected_orbits=(119, 120),
            expected_detectors=detectors,
        )


@pytest.mark.parametrize(
    ("replacement", "match"),
    [
        ("100010 0\n", "coverage is not exact"),
        ("100010 0\n100010 1\n", "duplicate cadence"),
        ("100010 2\n100011 0\n", "raw 0 or 1"),
        ("100010,0\n100011,1\n", "expected whitespace-delimited"),
    ],
)
def test_builder_rejects_bad_qlp_qflag_coverage_or_format(
    tmp_path: Path, replacement: str, match: str
) -> None:
    detectors = ((1, 1),)
    sources, qflags, quality_path, provenance_path, _ = _make_inputs(
        tmp_path, cameras=(1,), detectors=detectors
    )
    qflags[0].path.write_text(replacement, encoding="utf-8")

    with pytest.raises(ValueError, match=match):
        build_cadence_reference(
            sector=56,
            quat_sources=sources,
            qlp_qflag_sources=qflags,
            spoc_quality_table=quality_path,
            spoc_quality_provenance=provenance_path,
            expected_orbits=(119, 120),
            expected_detectors=detectors,
        )


def test_builder_rejects_incomplete_qlp_qflag_inventory(tmp_path: Path) -> None:
    detectors = ((1, 1),)
    sources, qflags, quality_path, provenance_path, _ = _make_inputs(
        tmp_path, cameras=(1,), detectors=detectors
    )

    with pytest.raises(ValueError, match="QLP qflag source inventory mismatch"):
        build_cadence_reference(
            sector=56,
            quat_sources=sources,
            qlp_qflag_sources=qflags[:-1],
            spoc_quality_table=quality_path,
            spoc_quality_provenance=provenance_path,
            expected_orbits=(119, 120),
            expected_detectors=detectors,
        )


def test_writer_rejects_qlp_qflag_source_change_during_build(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    detectors = ((1, 1),)
    sources, qflags, quality_path, provenance_path, _ = _make_inputs(
        tmp_path, cameras=(1,), detectors=detectors
    )
    table_path = tmp_path / "reference.csv"
    manifest_path = tmp_path / "reference.json"
    real_builder = cadence_reference.build_cadence_reference

    def mutate_after_build(**kwargs: object) -> pd.DataFrame:
        frame = real_builder(**kwargs)
        qflags[0].path.write_text(
            qflags[0].path.read_text(encoding="utf-8") + "\n",
            encoding="utf-8",
        )
        return frame

    monkeypatch.setattr(
        cadence_reference, "build_cadence_reference", mutate_after_build
    )
    with pytest.raises(RuntimeError, match="authority input changed"):
        write_cadence_reference(
            sector=56,
            quat_sources=sources,
            qlp_qflag_sources=qflags,
            spoc_quality_table=quality_path,
            spoc_quality_provenance=provenance_path,
            output_table=table_path,
            output_manifest=manifest_path,
            expected_orbits=(119, 120),
            expected_detectors=detectors,
        )
    assert not table_path.exists()
    assert not manifest_path.exists()


def test_s56_spoc_quality_cli_builds_full_authority_pair(tmp_path: Path) -> None:
    detectors = tuple(
        (camera, ccd) for camera in range(1, 5) for ccd in range(1, 5)
    )
    quats, flags, _ = _make_raw_authorities(
        tmp_path, cameras=(1, 2, 3, 4), detectors=detectors
    )
    table_path = tmp_path / "spoc_quality.csv"
    provenance_path = tmp_path / "spoc_quality.json"
    command = [
        sys.executable,
        str(
            REPO_ROOT
            / "scripts"
            / "stage1_lightcurves"
            / "build_s56_spoc_quality_table.py"
        ),
        "--output-table",
        str(table_path),
        "--output-provenance",
        str(provenance_path),
    ]
    for source in quats:
        command.extend(
            ["--quat", f"{source.orbitid},{source.camera},{source.path}"]
        )
    for source in flags:
        command.extend(
            ["--spoc-flag", f"{source.camera},{source.ccd},{source.path}"]
        )

    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr
    manifest = json.loads(provenance_path.read_text(encoding="utf-8"))
    assert manifest["n_rows"] == 64
    assert manifest["n_spoc_flag_files"] == 16
    assert manifest["n_quaternion_files"] == 8
    assert manifest["n_rows_excluded_by_quat"] == 4
    assert file_sha256(table_path) == manifest["table_sha256"]


def test_s56_cli_builds_full_detector_reference(tmp_path: Path) -> None:
    detectors = tuple(
        (camera, ccd) for camera in range(1, 5) for ccd in range(1, 5)
    )
    sources, qflags, quality_path, provenance_path, _ = _make_inputs(
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
    for source in qflags:
        command.extend(
            [
                "--qlp-qflag",
                f"{source.orbitid},{source.camera},{source.ccd},{source.path}",
            ]
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
    assert manifest["n_qlp_qflag_files_verified"] == 32
    assert manifest["quality_authority"] == "spoc_and_qlp_quality_flags"
    assert manifest["n_spoc_rows_excluded_by_quat"] == 4
    assert file_sha256(table_path) == manifest["table_sha256"]
