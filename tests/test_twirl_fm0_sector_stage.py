from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path
import tarfile

from twirl.models.fm0.corpus import CORPUS_SELECTION_FIELDS, CORPUS_SELECTION_SCHEMA_VERSION
from twirl.models.fm0.sector_stage import merge_sector_inventories, stage_sector_archive


def _csv(path: Path, fields, rows) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_sector_stage_extracts_only_selected_files_and_merges(tmp_path):
    selection = tmp_path / "selection.csv"
    paths = [
        "orbit-119/ffi/cam1/ccd1/LC/10.h5",
        "orbit-120/ffi/cam1/ccd1/LC/10.h5",
    ]
    _csv(
        selection,
        CORPUS_SELECTION_FIELDS,
        [
            {
                "schema_version": CORPUS_SELECTION_SCHEMA_VERSION,
                "gaia_dr3_source_id": "100",
                "tic_id": "10",
                "sector": 56,
                "camera": 1,
                "ccd": 1,
                "leakage_component_id": "leakage_test",
                "source_partition": "poc_train",
                "orbits_json": "[119,120]",
                "hdf5_paths_json": json.dumps(paths, separators=(",", ":")),
            }
        ],
    )
    selection_sha = hashlib.sha256(selection.read_bytes()).hexdigest()
    archive_dir = tmp_path / "archive"
    archive_dir.mkdir()
    archive = archive_dir / "s0056_A2v1_raw_hdf5.tar"
    with tarfile.open(archive, "w") as handle:
        for index, name in enumerate(paths):
            payload = f"hdf5-{index}".encode()
            info = tarfile.TarInfo(name)
            info.size = len(payload)
            handle.addfile(info, io.BytesIO(payload))
        extra = b"unused"
        info = tarfile.TarInfo("orbit-119/ffi/cam1/ccd1/LC/99.h5")
        info.size = len(extra)
        handle.addfile(info, io.BytesIO(extra))
    archive_sha = hashlib.sha256(archive.read_bytes()).hexdigest()
    (archive_dir / "s0056_A2v1_raw_hdf5.tar.sha256").write_text(
        f"{archive_sha}  s0056_A2v1_raw_hdf5.tar\n"
    )
    (archive_dir / "READY_ORCD").write_text("{}\n")
    quality_table = tmp_path / "cadence_reference.csv"
    quality_manifest = tmp_path / "cadence_reference.manifest.json"
    quality_table.write_text("sector,cadenceno\n56,1\n")
    quality_manifest.write_text("{}\n")
    sha = "1" * 40
    sector_root = tmp_path / "sectors"
    summary = stage_sector_archive(
        sector=56,
        selection_path=selection,
        selection_sha256=selection_sha,
        archive_dir=archive_dir,
        quality_table=quality_table,
        quality_manifest=quality_manifest,
        output_root=sector_root,
        producer_git_sha=sha,
        workers=2,
    )
    assert summary["n_hdf5_files"] == 2
    assert not (sector_root / "s0056" / "orbit-119/ffi/cam1/ccd1/LC/99.h5").exists()
    merged = merge_sector_inventories(
        sector_root=sector_root,
        sectors=(56,),
        out_dir=tmp_path / "merged",
        producer_git_sha=sha,
    )
    assert merged["n_observations"] == 1
    assert merged["sector_bindings"]["56"]["producer_git_sha"] == sha


def test_sector_merge_permits_only_explicit_reviewed_producers(tmp_path):
    sector_root = tmp_path / "sectors"
    sector_dir = sector_root / "s0056"
    sector_dir.mkdir(parents=True)
    inventory = sector_dir / "source_inventory.csv"
    from twirl.models.fm0.a2v1_adapter import A2V1_HDF5_SOURCE_INVENTORY_FIELDS

    _csv(
        inventory,
        A2V1_HDF5_SOURCE_INVENTORY_FIELDS,
        [
            {
                "gaia_dr3_source_id": "100",
                "tic_id": "10",
                "sector": 56,
                "a2v1_product_version": "A2v1",
                "product_state": "A2V1_ACCEPTED",
                "diagnostic_admission_receipt_path": "",
                "diagnostic_admission_receipt_sha256": "",
                "hdf5_paths_json": "[]",
                "hdf5_sha256_json": "[]",
                "quality_table_path": "/authority/table",
                "quality_table_sha256": "a" * 64,
                "quality_manifest_path": "/authority/manifest",
                "quality_manifest_sha256": "b" * 64,
            }
        ],
    )
    sector_sha = "1" * 40
    summary = {
        "schema_version": "twirl_fm0_1_sector_stage_v1",
        "producer_git_sha": sector_sha,
        "sector": 56,
        "source_inventory_sha256": hashlib.sha256(inventory.read_bytes()).hexdigest(),
    }
    (sector_dir / "summary.json").write_text(json.dumps(summary))
    (sector_dir / "READY").write_text(sector_sha + "\n")
    merged = merge_sector_inventories(
        sector_root=sector_root,
        sectors=(56,),
        out_dir=tmp_path / "merged",
        producer_git_sha="2" * 40,
        allowed_sector_producer_git_shas=(sector_sha,),
    )
    assert merged["allowed_sector_producer_git_shas"] == [sector_sha]
