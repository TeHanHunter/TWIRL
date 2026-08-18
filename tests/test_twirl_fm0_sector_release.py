from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

from twirl.models.fm0.a2v1_adapter import A2V1_HDF5_MANIFEST_FIELDS
from twirl.models.fm0.input_release import MANIFEST_COLUMNS, VISIT_TIMING_COLUMNS
from twirl.models.fm0.registry import build_alias_registry, build_observation_registry, write_registry_release
from twirl.models.fm0.sector_release import merge_sector_input_releases


def _csv(path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _stage(
    *,
    root: Path,
    sector: int,
    observation: dict[str, object],
    registry_sha: str,
    source_sha: str,
    start: float,
    end: float,
    producer_sha: str,
) -> None:
    directory = root / f"s{sector:04d}"
    (directory / "shards").mkdir(parents=True)
    shard_name = f"shards/{observation['observation_key']}.npz"
    shard = directory / shard_name
    shard.write_bytes(f"sector-{sector}".encode())
    manifest = {
        "input_release_schema_version": "twirl_fm0_1_input_release_v1",
        "observation_key": observation["observation_key"],
        "product_instance_id": observation["product_instance_id"],
        "source_sha256": observation["source_sha256"],
        "leakage_component_id": observation["leakage_component_id"],
        "source_partition": observation["source_partition"],
        "product_state": "A2V1_ACCEPTED",
        "relative_path": shard_name,
        "sha256": _sha(shard),
        "input_source_sha256": observation["source_sha256"],
        "n_cadences": 80,
        "n_segments": 1,
        "view_present_json": "[1,1,1,1,1,1]",
        "host_visit_offset_cadences": "0.0",
        "host_visit_gap_cadences": "0.0",
        "host_visit_overlaps_previous": False,
        "input_adapter": "a2v1_hdf5_quality_aware_v1",
        "scientific_training_eligible": True,
    }
    _csv(directory / "manifest.csv", MANIFEST_COLUMNS, [manifest])
    timing = {
        "observation_key": observation["observation_key"],
        "physical_source_id": observation["physical_source_id"],
        "absolute_visit_start": repr(start),
        "absolute_visit_end": repr(end),
    }
    _csv(directory / "visit_timing.csv", VISIT_TIMING_COLUMNS, [timing])
    source = {
        "observation_key": observation["observation_key"],
        "product_instance_id": observation["product_instance_id"],
        "source_sha256": observation["source_sha256"],
        "hdf5_paths_json": "[]",
        "hdf5_sha256_json": "[]",
        "quality_table_path": "/quality/table",
        "quality_table_sha256": "a" * 64,
        "quality_manifest_path": "/quality/manifest",
        "quality_manifest_sha256": "b" * 64,
    }
    _csv(directory / "a2v1_hdf5_manifest.csv", A2V1_HDF5_MANIFEST_FIELDS, [source])
    (directory / "summary.json").write_text("{}\n")
    sector_summary = {
        "schema_version": "twirl_fm0_1_sector_input_release_v1",
        "producer_git_sha": producer_sha,
        "sector": sector,
        "registry_observations_sha256": registry_sha,
        "a2v1_hdf5_manifest_sha256": _sha(directory / "a2v1_hdf5_manifest.csv"),
        "manifest_sha256": _sha(directory / "manifest.csv"),
        "summary_sha256": _sha(directory / "summary.json"),
        "visit_timing_sha256": _sha(directory / "visit_timing.csv"),
        "n_observations": 1,
        "n_cadences": 80,
        "input_adapter": "a2v1_hdf5_quality_aware_v1",
        "scientific_training_eligible": True,
    }
    (directory / "sector_summary.json").write_text(
        json.dumps(sector_summary, sort_keys=True) + "\n"
    )
    (directory / "READY").write_text(producer_sha + "\n")


def test_sector_input_merge_recomputes_cross_sector_host_timing(tmp_path: Path) -> None:
    aliases = build_alias_registry([{"gaia_dr3_source_id": "100", "tic_id": "10"}])
    observations = build_observation_registry(
        [
            {
                "gaia_dr3_source_id": "100",
                "tic_id": "10",
                "sector": sector,
                "a2v1_product_version": "A2v1",
                "source_sha256": digest * 64,
                "product_state": "A2V1_ACCEPTED",
            }
            for sector, digest in ((56, "a"), (57, "b"))
        ],
        aliases,
    )
    registry = tmp_path / "registry"
    write_registry_release(registry, aliases, observations)
    producer = "1" * 40
    stages = tmp_path / "stages"
    observations_by_sector = {int(row["sector"]): row for row in observations}
    _stage(
        root=stages,
        sector=56,
        observation=observations_by_sector[56],
        registry_sha=_sha(registry / "observations.csv"),
        source_sha="a" * 64,
        start=100.0,
        end=110.0,
        producer_sha=producer,
    )
    _stage(
        root=stages,
        sector=57,
        observation=observations_by_sector[57],
        registry_sha=_sha(registry / "observations.csv"),
        source_sha="b" * 64,
        start=130.0,
        end=140.0,
        producer_sha=producer,
    )
    merged = merge_sector_input_releases(
        sector_root=stages,
        sectors=(56, 57),
        registry_dir=registry,
        out_dir=tmp_path / "merged",
        producer_git_sha=producer,
    )
    assert merged["n_observations"] == 2
    rows = list(csv.DictReader((tmp_path / "merged" / "manifest.csv").open(newline="")))
    rows_by_key = {row["observation_key"]: row for row in rows}
    assert float(rows_by_key[observations_by_sector[56]["observation_key"]]["host_visit_offset_cadences"]) == 0.0
    assert float(rows_by_key[observations_by_sector[57]["observation_key"]]["host_visit_offset_cadences"]) == 30.0 * 86400.0 / 200.0
    assert float(rows_by_key[observations_by_sector[57]["observation_key"]]["host_visit_gap_cadences"]) == 20.0 * 86400.0 / 200.0
    assert (tmp_path / "merged" / "shards" / f"{observations_by_sector[56]['observation_key']}.npz").is_file()
