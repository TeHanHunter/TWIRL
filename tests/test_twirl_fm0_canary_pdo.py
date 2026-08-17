from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys


SCRIPT = Path(__file__).parents[1] / "scripts/stage5_validation/prepare_twirl_fm0_s56_canary_pdo.py"


def test_pdo_canary_stage_hashes_sources_without_modifying_them(tmp_path: Path) -> None:
    accepted = tmp_path / "accepted"
    paths = []
    for orbit in (119, 120):
        path = accepted / f"orbit-{orbit}/ffi/cam4/ccd1/LC/42.h5"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f"orbit-{orbit}".encode())
        paths.append(path)
    selection = tmp_path / "selection.csv"
    fields = (
        "schema_version", "gaia_dr3_source_id", "tic_id", "sector", "camera", "ccd",
        "leakage_component_id", "source_partition", "rank_sha256", "is_benchmark",
        "orbits_json", "hdf5_paths_json",
    )
    with selection.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerow(
            {
                "schema_version": "twirl_fm0_1_s56_canary_selection_v1",
                "gaia_dr3_source_id": "123", "tic_id": "42", "sector": 56,
                "camera": 4, "ccd": 1, "leakage_component_id": "leakage_x",
                "source_partition": "poc_train", "rank_sha256": "a" * 64,
                "is_benchmark": True, "orbits_json": "[119,120]",
                "hdf5_paths_json": json.dumps([str(path) for path in paths]),
            }
        )
    selection_sha = hashlib.sha256(selection.read_bytes()).hexdigest()
    staging = tmp_path / "staging"
    completed = subprocess.run(
        [
            sys.executable, str(SCRIPT), "--selection", str(selection),
            "--selection-sha256", selection_sha, "--accepted-root", str(accepted),
            "--staging-dir", str(staging), "--orcd-source-root", "/orcd/fm0/source",
            "--quality-table", Path("/orcd/quality.csv"),
            "--quality-table-sha256", "b" * 64,
            "--quality-manifest", Path("/orcd/quality.json"),
            "--quality-manifest-sha256", "c" * 64,
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert json.loads(completed.stdout)["n_hdf5_files"] == 2
    assert paths[0].read_bytes() == b"orbit-119"
    inventory = list(csv.DictReader((staging / "source_inventory.csv").open()))
    assert len(inventory) == 1
    assert inventory[0]["product_state"] == "A2V1_ACCEPTED"
    assert json.loads(inventory[0]["hdf5_paths_json"])[0].startswith("/orcd/fm0/source/")
