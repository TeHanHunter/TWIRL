#!/usr/bin/env python3
"""Hash a frozen S56 selection and prepare its bounded PDO transfer manifest.

This utility reads accepted A2v1 HDF5 files but never modifies them.  It writes
only below the caller-supplied user-owned staging directory.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
from typing import Any


SELECTION_FIELDS = (
    "schema_version",
    "gaia_dr3_source_id",
    "tic_id",
    "sector",
    "camera",
    "ccd",
    "leakage_component_id",
    "source_partition",
    "rank_sha256",
    "is_benchmark",
    "orbits_json",
    "hdf5_paths_json",
)
INVENTORY_FIELDS = (
    "gaia_dr3_source_id",
    "tic_id",
    "sector",
    "a2v1_product_version",
    "product_state",
    "diagnostic_admission_receipt_path",
    "diagnostic_admission_receipt_sha256",
    "hdf5_paths_json",
    "hdf5_sha256_json",
    "quality_table_path",
    "quality_table_sha256",
    "quality_manifest_path",
    "quality_manifest_sha256",
)
SHA_RE = re.compile(r"^[0-9a-f]{64}$")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write(path: Path, payload: bytes) -> None:
    if path.exists():
        raise RuntimeError(f"refusing to overwrite {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    if temporary.exists():
        raise RuntimeError(f"stale temporary output: {temporary}")
    with temporary.open("xb") as handle:
        handle.write(payload)
        handle.flush()
    temporary.replace(path)


def _csv_bytes(rows: list[dict[str, Any]], fields: tuple[str, ...]) -> bytes:
    from io import StringIO

    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(fields), lineterminator="\n")
    writer.writeheader()
    writer.writerows({field: row[field] for field in fields} for row in rows)
    return stream.getvalue().encode("utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", required=True, type=Path)
    parser.add_argument("--selection-sha256", required=True)
    parser.add_argument("--accepted-root", required=True, type=Path)
    parser.add_argument("--staging-dir", required=True, type=Path)
    parser.add_argument("--orcd-source-root", required=True, type=Path)
    parser.add_argument("--quality-table", required=True, type=Path)
    parser.add_argument("--quality-table-sha256", required=True)
    parser.add_argument("--quality-manifest", required=True, type=Path)
    parser.add_argument("--quality-manifest-sha256", required=True)
    args = parser.parse_args()

    accepted_root = args.accepted_root.resolve()
    staging = args.staging_dir.resolve()
    if staging == accepted_root or accepted_root in staging.parents:
        raise RuntimeError("staging directory must not be inside the accepted A2v1 tree")
    for value, name in (
        (args.selection_sha256, "selection SHA-256"),
        (args.quality_table_sha256, "quality-table SHA-256"),
        (args.quality_manifest_sha256, "quality-manifest SHA-256"),
    ):
        if not SHA_RE.fullmatch(value):
            raise RuntimeError(f"invalid {name}")
    if sha256_file(args.selection) != args.selection_sha256:
        raise RuntimeError("selection checksum mismatch")

    with args.selection.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != SELECTION_FIELDS:
            raise RuntimeError("selection columns do not match the frozen canary schema")
        selected = list(reader)
    if not selected:
        raise RuntimeError("empty canary selection")

    inventory: list[dict[str, Any]] = []
    relative_files: list[str] = []
    checksum_rows: list[dict[str, str]] = []
    total_bytes = 0
    seen_sources: set[tuple[str, str]] = set()
    for row in selected:
        if row["schema_version"] != "twirl_fm0_1_s56_canary_selection_v1":
            raise RuntimeError("selection schema version mismatch")
        if int(row["sector"]) != 56:
            raise RuntimeError("canary contains a non-S56 row")
        identity = (row["gaia_dr3_source_id"], row["tic_id"])
        if identity in seen_sources:
            raise RuntimeError(f"duplicate canary source: {identity}")
        seen_sources.add(identity)
        source_paths = [Path(value).resolve() for value in json.loads(row["hdf5_paths_json"])]
        if len(source_paths) != 2:
            raise RuntimeError(f"canary source does not have two orbit files: {identity}")
        hashes: list[str] = []
        remote_paths: list[str] = []
        for source in source_paths:
            if accepted_root not in source.parents or not source.is_file():
                raise RuntimeError(f"missing or unsafe accepted source: {source}")
            relative = source.relative_to(accepted_root)
            if relative.parts[0] not in {"orbit-119", "orbit-120"}:
                raise RuntimeError(f"unexpected S56 orbit path: {relative}")
            observed_sha = sha256_file(source)
            size = source.stat().st_size
            total_bytes += size
            relative_text = str(relative)
            relative_files.append(relative_text)
            hashes.append(observed_sha)
            remote_paths.append(str(args.orcd_source_root / relative))
            checksum_rows.append(
                {"relative_path": relative_text, "size_bytes": str(size), "sha256": observed_sha}
            )
        inventory.append(
            {
                "gaia_dr3_source_id": row["gaia_dr3_source_id"],
                "tic_id": row["tic_id"],
                "sector": "56",
                "a2v1_product_version": "A2v1",
                "product_state": "A2V1_ACCEPTED",
                "diagnostic_admission_receipt_path": "",
                "diagnostic_admission_receipt_sha256": "",
                "hdf5_paths_json": json.dumps(remote_paths, separators=(",", ":")),
                "hdf5_sha256_json": json.dumps(hashes, separators=(",", ":")),
                "quality_table_path": str(args.quality_table),
                "quality_table_sha256": args.quality_table_sha256,
                "quality_manifest_path": str(args.quality_manifest),
                "quality_manifest_sha256": args.quality_manifest_sha256,
            }
        )

    if len(relative_files) != len(set(relative_files)):
        raise RuntimeError("the canary selection reuses an HDF5 path")
    checksum_rows.sort(key=lambda item: item["relative_path"])
    relative_files.sort()
    inventory.sort(key=lambda item: (int(item["gaia_dr3_source_id"]), int(item["tic_id"])))
    inventory_bytes = _csv_bytes(inventory, INVENTORY_FIELDS)
    checksums_bytes = _csv_bytes(
        checksum_rows, ("relative_path", "size_bytes", "sha256")
    )
    files_bytes = ("\n".join(relative_files) + "\n").encode("utf-8")
    summary = {
        "schema_version": "twirl_fm0_1_s56_canary_pdo_stage_v1",
        "selection_sha256": args.selection_sha256,
        "n_sources": len(inventory),
        "n_hdf5_files": len(relative_files),
        "total_hdf5_bytes": total_bytes,
        "outputs_sha256": {
            "source_inventory.csv": hashlib.sha256(inventory_bytes).hexdigest(),
            "hdf5_checksums.csv": hashlib.sha256(checksums_bytes).hexdigest(),
            "files.txt": hashlib.sha256(files_bytes).hexdigest(),
        },
        "accepted_sources_modified": False,
        "claim_limit": "bounded input canary staging; not a model result",
    }
    atomic_write(staging / "source_inventory.csv", inventory_bytes)
    atomic_write(staging / "hdf5_checksums.csv", checksums_bytes)
    atomic_write(staging / "files.txt", files_bytes)
    atomic_write(
        staging / "stage_summary.json",
        (json.dumps(summary, indent=2, sort_keys=True) + "\n").encode("utf-8"),
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
