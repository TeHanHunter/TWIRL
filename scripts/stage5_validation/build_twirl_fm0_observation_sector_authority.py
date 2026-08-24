#!/usr/bin/env python3
"""Build the sealed-safe FM0 evaluator observation-to-sector authority."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
from typing import Any


OUTPUT_FIELDS = (
    "observation_key",
    "leakage_component_id",
    "source_partition",
    "sector",
)
IDENTITY_FIELDS = (
    "product_instance_id",
    "source_sha256",
    "leakage_component_id",
    "source_partition",
    "product_state",
    "sha256",
    "input_source_sha256",
    "n_cadences",
    "n_segments",
    "view_present_json",
    "input_adapter",
    "scientific_training_eligible",
)
ALLOWED_OUTPUT_PARTITIONS = {"poc_train", "poc_development"}
SEALED_PARTITION = "poc_sealed_test"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")
        required = {"observation_key", *IDENTITY_FIELDS}
        if not required.issubset(reader.fieldnames):
            raise ValueError(f"CSV schema is incomplete: {path}")
        return [dict(row) for row in reader]


def write_atomic(path: Path, payload: bytes) -> str:
    if path.exists() or path.with_name(path.name + ".sha256").exists():
        raise FileExistsError(f"refusing to overwrite authority artifact: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    temporary.write_bytes(payload)
    os.replace(temporary, path)
    digest = hashlib.sha256(payload).hexdigest()
    sidecar = path.with_name(path.name + ".sha256")
    sidecar.write_text(f"{digest}  {path.name}\n", encoding="utf-8")
    return digest


def parse_sector_manifest(value: str) -> tuple[int, Path]:
    raw_sector, separator, raw_path = value.partition("=")
    if separator != "=" or not raw_sector.isdigit() or not raw_path:
        raise argparse.ArgumentTypeError("sector manifest must be SECTOR=PATH")
    sector = int(raw_sector)
    if sector < 56:
        raise argparse.ArgumentTypeError("FM0 sector authority requires Sector >= 56")
    return sector, Path(raw_path)


def build_authority(
    *,
    merged_manifest: Path,
    merge_summary: Path,
    merge_summary_sha256: str,
    sector_manifests: list[tuple[int, Path]],
) -> tuple[bytes, dict[str, Any]]:
    observed_merge_summary_sha = sha256_file(merge_summary)
    if observed_merge_summary_sha != merge_summary_sha256:
        raise ValueError("merge-summary SHA-256 mismatch")
    summary = json.loads(merge_summary.read_text(encoding="utf-8"))
    if not isinstance(summary, dict):
        raise ValueError("merge summary is not a mapping")
    merged_manifest_sha = sha256_file(merged_manifest)
    if summary.get("manifest_sha256") != merged_manifest_sha:
        raise ValueError("merge summary does not bind the merged manifest")
    sectors = summary.get("sectors")
    bindings = summary.get("sector_bindings")
    if not isinstance(sectors, list) or not isinstance(bindings, dict):
        raise ValueError("merge summary has no sector bindings")
    supplied = {sector: path.resolve(strict=True) for sector, path in sector_manifests}
    if len(supplied) != len(sector_manifests) or sorted(supplied) != sectors:
        raise ValueError("sector manifests do not exactly cover the merged release")

    by_observation: dict[str, tuple[int, dict[str, str]]] = {}
    sector_receipts: dict[str, dict[str, Any]] = {}
    for sector in sectors:
        path = supplied[sector]
        digest = sha256_file(path)
        binding = bindings.get(str(sector))
        if not isinstance(binding, dict) or binding.get("manifest_sha256") != digest:
            raise ValueError(f"per-sector manifest binding differs for S{sector}")
        rows = read_rows(path)
        if binding.get("n_observations") != len(rows):
            raise ValueError(f"per-sector row count differs for S{sector}")
        for row in rows:
            key = row["observation_key"]
            if not key or key in by_observation:
                raise ValueError(f"duplicate/empty observation key: {key!r}")
            by_observation[key] = (sector, row)
        sector_receipts[str(sector)] = {
            "path": str(path),
            "sha256": digest,
            "n_observations": len(rows),
        }

    merged_rows = read_rows(merged_manifest)
    if len(merged_rows) != len(by_observation):
        raise ValueError("per-sector and merged manifest row counts differ")
    output_rows: list[dict[str, str | int]] = []
    partition_counts: dict[str, int] = {}
    seen: set[str] = set()
    for row in merged_rows:
        key = row["observation_key"]
        if key in seen or key not in by_observation:
            raise ValueError(f"merged observation is duplicated or unbound: {key}")
        seen.add(key)
        sector, sector_row = by_observation[key]
        for field in IDENTITY_FIELDS:
            if row[field] != sector_row[field]:
                raise ValueError(f"merged/per-sector mismatch: {key}.{field}")
        partition = row["source_partition"]
        partition_counts[partition] = partition_counts.get(partition, 0) + 1
        if partition in ALLOWED_OUTPUT_PARTITIONS:
            output_rows.append(
                {
                    "observation_key": key,
                    "leakage_component_id": row["leakage_component_id"],
                    "source_partition": partition,
                    "sector": sector,
                }
            )
        elif partition != SEALED_PARTITION:
            raise ValueError(f"unexpected source partition: {partition}")
    if seen != set(by_observation):
        raise ValueError("per-sector manifests contain rows absent from merged manifest")
    output_rows.sort(key=lambda row: str(row["observation_key"]))
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=OUTPUT_FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(output_rows)
    payload = buffer.getvalue().encode("utf-8")
    receipt = {
        "schema_version": "twirl_fm0_2_observation_sector_authority_receipt_v1",
        "passed": True,
        "construction": "per_sector_manifest_to_merged_manifest_one_to_one_join_v1",
        "merged_manifest": {
            "path": str(merged_manifest.resolve(strict=True)),
            "sha256": merged_manifest_sha,
            "n_observations": len(merged_rows),
        },
        "merge_summary": {
            "path": str(merge_summary.resolve(strict=True)),
            "sha256": observed_merge_summary_sha,
        },
        "sector_manifests": sector_receipts,
        "output_columns": list(OUTPUT_FIELDS),
        "output_rows": len(output_rows),
        "partition_counts": partition_counts,
        "sealed_rows_emitted": 0,
        "model_visible": False,
        "claim_limit": (
            "Evaluator-only sector authority; no model input, sealed metric, "
            "production use, or foundation-model claim."
        ),
    }
    return payload, receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--merged-manifest", type=Path, required=True)
    parser.add_argument("--merge-summary", type=Path, required=True)
    parser.add_argument("--merge-summary-sha256", required=True)
    parser.add_argument(
        "--sector-manifest",
        action="append",
        type=parse_sector_manifest,
        required=True,
    )
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-receipt", type=Path, required=True)
    args = parser.parse_args()
    if len(args.merge_summary_sha256) != 64:
        raise ValueError("merge-summary SHA-256 must be 64 hexadecimal characters")
    csv_payload, receipt = build_authority(
        merged_manifest=args.merged_manifest.resolve(strict=True),
        merge_summary=args.merge_summary.resolve(strict=True),
        merge_summary_sha256=args.merge_summary_sha256,
        sector_manifests=args.sector_manifest,
    )
    csv_sha = write_atomic(args.output_csv, csv_payload)
    receipt["output_authority"] = {
        "path": str(args.output_csv.resolve(strict=True)),
        "sha256": csv_sha,
    }
    receipt_payload = (
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    write_atomic(args.output_receipt, receipt_payload)
    print(args.output_receipt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
