from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from twirl.lightcurves import mission_quality_reference as quality_reference
from twirl.lightcurves.mission_quality_reference import (
    MISSION_QUALITY_REFERENCE_COLUMNS,
    MISSION_QUALITY_REFERENCE_SCHEMA_VERSION,
    build_mission_quality_reference,
    load_mission_quality_reference,
)
from twirl.models.fm0.mission_quality_admission import (
    mission_quality_type,
    verify_mission_quality_sources,
    write_mission_quality_receipt,
)
from twirl.models.fm0.mission_quality_transfer import (
    package_mission_quality_transfer,
)
from twirl.models.fm0.registry import FM0ContractError


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_quat(path: Path, cadences: tuple[int, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "cadence,flag\n" + "".join(f"{cadence},0\n" for cadence in cadences),
        encoding="utf-8",
    )


def _write_flags(path: Path, cadences: tuple[int, ...], *, nonzero_value: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "".join(
            f"{cadence}, {nonzero_value if index % 2 else 0}\n"
            for index, cadence in enumerate(cadences)
        ),
        encoding="utf-8",
    )


def _quality_tree(
    root: Path, *, sector: int, exclusion: tuple[int, int] | None = None
) -> Path:
    orbits = (2 * sector + 7, 2 * sector + 8)
    provider = mission_quality_type(sector)
    prefix = "spocffiflag" if provider == "spoc" else "ticaffiflag"
    directory_name = "spocflags" if provider == "spoc" else "ticaflags"
    for camera in range(1, 5):
        first = (10_000 + camera * 100, 10_001 + camera * 100)
        second = (20_000 + camera * 100, 20_001 + camera * 100)
        for orbit, cadences in zip(orbits, (first, second), strict=True):
            run = root / f"orbit-{orbit}/ffi/run"
            _write_quat(run / f"cam{camera}_quat.txt", cadences)
            for ccd in range(1, 5):
                _write_flags(
                    run / f"cam{camera}ccd{ccd}_qflag.txt",
                    cadences,
                    nonzero_value=1,
                )
        for ccd in range(1, 5):
            path = (
                root / directory_name / f"{prefix}_s{sector}_cam{camera}_ccd{ccd}.txt"
            )
            _write_flags(
                path,
                first + second,
                nonzero_value=4 if provider == "tica" else 2,
            )
            if exclusion is not None and (camera, ccd) == (1, 1):
                cadence, value = exclusion
                with path.open("a", encoding="utf-8") as handle:
                    handle.write(f"{cadence}, {value}\n")

    if provider == "tica":
        directory = root / "ticaflags"
        detectors = []
        n_rows = 0
        n_nonzero = 0
        for camera in range(1, 5):
            for ccd in range(1, 5):
                path = directory / f"ticaffiflag_s{sector}_cam{camera}_ccd{ccd}.txt"
                rows = [
                    line
                    for line in path.read_text(encoding="utf-8").splitlines()
                    if line
                ]
                nonzero = sum(int(line.split(",")[1]) != 0 for line in rows)
                detectors.append(
                    {
                        "camera": camera,
                        "ccd": ccd,
                        "path": path.name,
                        "sha256": _sha(path),
                        "n_rows": len(rows),
                        "n_nonzero": nonzero,
                    }
                )
                n_rows += len(rows)
                n_nonzero += nonzero
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
            "n_rows": n_rows,
            "n_nonzero": n_nonzero,
            "cadence_coverage_verified": False,
            "detectors": detectors,
        }
        (directory / "summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        (directory / "READY").write_text(producer + "\n", encoding="utf-8")
    return root


def _transfer(
    tmp_path: Path, *, sector: int, exclusion: tuple[int, int] | None = None
) -> Path:
    authority_root = _quality_tree(
        tmp_path / "authorities", sector=sector, exclusion=exclusion
    )
    receipt = verify_mission_quality_sources(
        sector=sector,
        expected_orbits=(2 * sector + 7, 2 * sector + 8),
        qlp_root=authority_root,
    )
    receipt_path = tmp_path / f"s{sector}.mission-quality.json"
    write_mission_quality_receipt(receipt, receipt_path)
    transfer = tmp_path / "transfer"
    package_mission_quality_transfer(
        receipt_paths=[receipt_path],
        output_dir=transfer,
        producer_git_sha="c" * 40,
    )
    return transfer


@pytest.mark.parametrize(("sector", "provider"), ((66, "spoc"), (67, "tica")))
def test_build_and_load_provider_neutral_reference(
    tmp_path: Path, sector: int, provider: str
) -> None:
    transfer = _transfer(tmp_path, sector=sector)
    output = tmp_path / "reference"
    manifest = build_mission_quality_reference(
        quality_transfer_root=transfer,
        expected_quality_transfer_manifest_sha256=_sha(transfer / "manifest.json"),
        sector=sector,
        output_dir=output,
        producer_git_sha="d" * 40,
    )

    assert manifest["schema_version"] == MISSION_QUALITY_REFERENCE_SCHEMA_VERSION
    assert manifest["mission_quality_provider"] == provider
    assert manifest["table_columns"] == list(MISSION_QUALITY_REFERENCE_COLUMNS)
    assert manifest["n_rows"] == 64
    assert manifest["hdf5_quality_join_required"] is True
    assert manifest["six_view_shards_verified"] is False
    assert manifest["panel_admission_authorized"] is False
    assert all(
        not path.is_symlink() and not path.stat().st_mode & 0o222
        for path in (output, *output.rglob("*"))
    )

    reference = load_mission_quality_reference(
        reference_dir=output,
        sector=sector,
        expected_orbits=(2 * sector + 7, 2 * sector + 8),
    )
    result = reference.lookup(
        orbit=2 * sector + 7,
        camera=1,
        ccd=1,
        cadence=np.asarray((10_100, 10_101)),
    )
    expected_mission = 4 if provider == "tica" else 2
    assert result.mission_quality.tolist() == [0, expected_mission]
    assert result.qlp_quality.tolist() == [0, 1]
    assert result.authority_excluded.tolist() == [False, False]
    assert result.external_quality.tolist() == [
        0,
        expected_mission | (1 << 30),
    ]
    assert result.external_bad.tolist() == [False, True]
    reference.assert_unchanged()
    with pytest.raises(FM0ContractError, match="refusing to overwrite"):
        build_mission_quality_reference(
            quality_transfer_root=transfer,
            expected_quality_transfer_manifest_sha256=_sha(
                transfer / "manifest.json"
            ),
            sector=sector,
            output_dir=output,
            producer_git_sha="d" * 40,
        )


def test_mission_only_rows_remain_explicit_authority_exclusions(
    tmp_path: Path,
) -> None:
    transfer = _transfer(tmp_path, sector=67, exclusion=(99_999, 4))
    output = tmp_path / "reference"
    build_mission_quality_reference(
        quality_transfer_root=transfer,
        expected_quality_transfer_manifest_sha256=_sha(transfer / "manifest.json"),
        sector=67,
        output_dir=output,
        producer_git_sha="d" * 40,
    )
    reference = load_mission_quality_reference(reference_dir=output, sector=67)
    result = reference.lookup(
        orbit=141,
        camera=1,
        ccd=1,
        cadence=np.asarray((10_100, 99_999)),
    )
    assert result.mission_quality.tolist() == [0, 4]
    assert result.qlp_quality.tolist() == [0, 0]
    assert result.authority_excluded.tolist() == [False, True]
    assert result.external_bad.tolist() == [False, True]
    assert reference.provenance["n_authority_exclusions"] == 1
    with pytest.raises(FM0ContractError, match="unexplained cadence"):
        reference.lookup(
            orbit=141,
            camera=1,
            ccd=1,
            cadence=np.asarray((88_888,)),
        )


def test_builder_requires_controller_authorized_transfer_manifest(tmp_path: Path) -> None:
    transfer = _transfer(tmp_path, sector=66)
    with pytest.raises(FM0ContractError, match="transfer manifest hash drifted"):
        build_mission_quality_reference(
            quality_transfer_root=transfer,
            expected_quality_transfer_manifest_sha256="f" * 64,
            sector=66,
            output_dir=tmp_path / "reference",
            producer_git_sha="d" * 40,
        )


def test_builder_cleans_owned_published_tree_and_allows_retry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    transfer = _transfer(tmp_path, sector=66)
    output = tmp_path / "reference"
    sibling = tmp_path / "reference-unrelated"
    sibling.write_text("keep\n", encoding="utf-8")
    original = quality_reference._assert_output_tree_read_only

    def fail_after_publish(root: Path) -> None:
        original(root)
        if root == output:
            raise KeyboardInterrupt("simulated post-publish interruption")

    monkeypatch.setattr(
        quality_reference, "_assert_output_tree_read_only", fail_after_publish
    )
    with pytest.raises(KeyboardInterrupt, match="post-publish interruption"):
        build_mission_quality_reference(
            quality_transfer_root=transfer,
            expected_quality_transfer_manifest_sha256=_sha(
                transfer / "manifest.json"
            ),
            sector=66,
            output_dir=output,
            producer_git_sha="d" * 40,
        )
    assert not output.exists()
    assert not output.with_name(output.name + ".partial").exists()
    assert sibling.read_text(encoding="utf-8") == "keep\n"

    monkeypatch.setattr(
        quality_reference, "_assert_output_tree_read_only", original
    )
    build_mission_quality_reference(
        quality_transfer_root=transfer,
        expected_quality_transfer_manifest_sha256=_sha(transfer / "manifest.json"),
        sector=66,
        output_dir=output,
        producer_git_sha="d" * 40,
    )
    assert output.is_dir()


def test_builder_refuses_broken_output_symlink_without_touching_target(
    tmp_path: Path,
) -> None:
    transfer = _transfer(tmp_path, sector=66)
    target = tmp_path / "outside-target"
    output = tmp_path / "reference"
    output.symlink_to(target, target_is_directory=True)

    with pytest.raises(FM0ContractError, match="refusing to overwrite"):
        build_mission_quality_reference(
            quality_transfer_root=transfer,
            expected_quality_transfer_manifest_sha256=_sha(
                transfer / "manifest.json"
            ),
            sector=66,
            output_dir=output,
            producer_git_sha="d" * 40,
        )
    assert output.is_symlink()
    assert not target.exists()


def test_builder_rejects_transfer_role_drift(tmp_path: Path) -> None:
    transfer = _transfer(tmp_path, sector=67)
    manifest_path = transfer / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    mission_row = next(
        row for row in manifest["sources"] if row["role"] == "tica_mission_quality"
    )
    mission_row["role"] = "spoc_mission_quality"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    with pytest.raises(FM0ContractError, match="role/path/hash drifted"):
        build_mission_quality_reference(
            quality_transfer_root=transfer,
            expected_quality_transfer_manifest_sha256=_sha(
                transfer / "manifest.json"
            ),
            sector=67,
            output_dir=tmp_path / "reference",
            producer_git_sha="d" * 40,
        )


def test_builder_rechecks_qflag_coverage_after_fully_rebound_transfer(
    tmp_path: Path,
) -> None:
    transfer = _transfer(tmp_path, sector=67)
    manifest_path = transfer / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    source_row = next(
        row
        for row in manifest["sources"]
        if row["role"] == "qlp_detector_qflag"
        and "orbit-141" in row["staged_path"]
        and row["staged_path"].endswith("cam1ccd1_qflag.txt")
    )
    source_path = transfer / source_row["staged_path"]
    with source_path.open("a", encoding="utf-8") as handle:
        handle.write("99999, 0\n")
    source_row["sha256"] = _sha(source_path)
    source_row["bytes"] = source_path.stat().st_size

    receipt_row = manifest["receipts"][0]
    receipt_path = transfer / receipt_row["staged_path"]
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    binding = next(
        row
        for row in receipt["source_bindings"]
        if row["role"] == "qlp_detector_qflag"
        and row["orbit"] == 141
        and row["camera"] == 1
        and row["ccd"] == 1
    )
    binding["sha256"] = source_row["sha256"]
    binding["n_rows"] += 1
    receipt_path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    receipt_row["sha256"] = _sha(receipt_path)
    receipt_row["bytes"] = receipt_path.stat().st_size
    manifest["n_bytes"] = sum(
        row["bytes"] for row in manifest["receipts"] + manifest["sources"]
    )
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    with pytest.raises(FM0ContractError, match="QLP qflag cadence coverage differs"):
        build_mission_quality_reference(
            quality_transfer_root=transfer,
            expected_quality_transfer_manifest_sha256=_sha(
                transfer / "manifest.json"
            ),
            sector=67,
            output_dir=tmp_path / "reference",
            producer_git_sha="d" * 40,
        )


def test_loader_rejects_rehashed_table_semantic_drift(tmp_path: Path) -> None:
    transfer = _transfer(tmp_path, sector=67)
    output = tmp_path / "reference"
    build_mission_quality_reference(
        quality_transfer_root=transfer,
        expected_quality_transfer_manifest_sha256=_sha(transfer / "manifest.json"),
        sector=67,
        output_dir=output,
        producer_git_sha="d" * 40,
    )
    table_path = output / "cadence_quality.csv"
    output.chmod(0o750)
    table_path.chmod(0o640)
    frame = pd.read_csv(table_path)
    index = int(frame.index[frame["mission_quality"] == 4][0])
    frame.loc[index, "mission_quality"] = 32
    frame.to_csv(table_path, index=False, lineterminator="\n")
    manifest_path = output / "manifest.json"
    manifest_path.chmod(0o640)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["table_sha256"] = _sha(table_path)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    with pytest.raises(FM0ContractError, match="disagrees with transferred"):
        load_mission_quality_reference(reference_dir=output, sector=67)


def test_loaded_reference_detects_bound_source_drift(tmp_path: Path) -> None:
    transfer = _transfer(tmp_path, sector=66)
    output = tmp_path / "reference"
    build_mission_quality_reference(
        quality_transfer_root=transfer,
        expected_quality_transfer_manifest_sha256=_sha(transfer / "manifest.json"),
        sector=66,
        output_dir=output,
        producer_git_sha="d" * 40,
    )
    reference = load_mission_quality_reference(reference_dir=output, sector=66)
    source = next(iter(reference.source_file_sha256))
    source.write_text(source.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="input changed"):
        reference.assert_unchanged()
