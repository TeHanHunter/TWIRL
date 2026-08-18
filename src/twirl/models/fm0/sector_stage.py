"""Sector-at-a-time expansion of checksum-bound A2v1 archives for FM0.1."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
from typing import Any, Mapping, Sequence

from .a2v1_adapter import A2V1_HDF5_SOURCE_INVENTORY_FIELDS
from .corpus import CORPUS_SELECTION_FIELDS, CORPUS_SELECTION_SCHEMA_VERSION, expected_orbits
from .registry import FM0ContractError, publish_immutable, sha256_file


SECTOR_STAGE_SCHEMA_VERSION = "twirl_fm0_1_sector_stage_v1"
SECTOR_MERGE_SCHEMA_VERSION = "twirl_fm0_1_sector_inventory_merge_v1"
SHA_RE = re.compile(r"^[0-9a-f]{64}$")


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if not rows:
        raise FM0ContractError(f"CSV is empty: {path}")
    return rows


def _csv_bytes(rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> bytes:
    from io import StringIO

    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(fields), lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row[field] for field in fields})
    return stream.getvalue().encode("utf-8")


def _declared_archive_hash(path: Path) -> str:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) != 1:
        raise FM0ContractError("sector archive SHA sidecar must contain one line")
    digest = lines[0].split()[0].strip().lower()
    if not SHA_RE.fullmatch(digest):
        raise FM0ContractError("sector archive sidecar has an invalid SHA-256")
    return digest


def _safe_relative_hdf5(value: str, *, sector: int) -> Path:
    path = Path(value)
    if path.is_absolute() or ".." in path.parts or path.suffix != ".h5":
        raise FM0ContractError(f"unsafe selected HDF5 path: {value!r}")
    if len(path.parts) != 6:
        raise FM0ContractError(f"selected HDF5 path has unexpected layout: {value}")
    allowed_orbits = {f"orbit-{orbit}" for orbit in expected_orbits(sector)}
    if path.parts[0] not in allowed_orbits or path.parts[1] != "ffi":
        raise FM0ContractError(f"selected HDF5 path has the wrong sector orbit: {value}")
    if not re.fullmatch(r"cam[1-4]", path.parts[2]):
        raise FM0ContractError(f"selected HDF5 path has invalid camera: {value}")
    if not re.fullmatch(r"ccd[1-4]", path.parts[3]) or path.parts[4] != "LC":
        raise FM0ContractError(f"selected HDF5 path has invalid CCD/layout: {value}")
    if not path.stem.isdigit() or int(path.stem) <= 0:
        raise FM0ContractError(f"selected HDF5 path has invalid TIC name: {value}")
    return path


def _hash_pair(item: tuple[str, Path]) -> tuple[str, int, str]:
    relative, path = item
    return relative, path.stat().st_size, sha256_file(path)


def stage_sector_archive(
    *,
    sector: int,
    selection_path: str | Path,
    selection_sha256: str,
    archive_dir: str | Path,
    quality_table: str | Path,
    quality_manifest: str | Path,
    output_root: str | Path,
    producer_git_sha: str,
    workers: int = 4,
) -> dict[str, Any]:
    """Extract selected members and publish one scientific source inventory."""

    sector = int(sector)
    if sector not in range(56, 66):
        raise FM0ContractError("accepted-sector staging is bounded to S56-S65")
    if not re.fullmatch(r"[0-9a-f]{40}", producer_git_sha):
        raise FM0ContractError("producer_git_sha must be full lowercase hex")
    if workers <= 0:
        raise FM0ContractError("workers must be positive")
    selection = Path(selection_path).resolve()
    if sha256_file(selection) != selection_sha256:
        raise FM0ContractError("corpus selection hash mismatch")
    rows = _read_csv(selection)
    if tuple(rows[0]) != CORPUS_SELECTION_FIELDS:
        raise FM0ContractError("corpus selection columns differ from the contract")
    selected = [row for row in rows if int(row["sector"]) == sector]
    if not selected:
        raise FM0ContractError(f"selection contains no S{sector} observations")

    tag = f"s{sector:04d}"
    archive_root = Path(archive_dir).resolve()
    archive = archive_root / f"{tag}_A2v1_raw_hdf5.tar"
    sidecar = archive.with_suffix(archive.suffix + ".sha256")
    if not (archive.is_file() and sidecar.is_file() and (archive_root / "READY_ORCD").is_file()):
        raise FM0ContractError(f"S{sector} archive is not ORCD-ready")
    archive_sha = _declared_archive_hash(sidecar)
    if sha256_file(archive) != archive_sha:
        raise FM0ContractError(f"S{sector} archive hash mismatch")

    quality_table_path = Path(quality_table).resolve()
    quality_manifest_path = Path(quality_manifest).resolve()
    if not quality_table_path.is_file() or not quality_manifest_path.is_file():
        raise FM0ContractError(f"S{sector} cadence authority is missing")
    quality_table_sha = sha256_file(quality_table_path)
    quality_manifest_sha = sha256_file(quality_manifest_path)

    final = Path(output_root).resolve() / tag
    partial = final.with_name(final.name + ".partial")
    if final.exists() or partial.exists():
        raise FM0ContractError(f"refusing to overwrite sector stage: {final}")
    partial.mkdir(parents=True)
    requested: list[Path] = []
    seen_visits: set[tuple[str, str, int]] = set()
    for row in selected:
        if row["schema_version"] != CORPUS_SELECTION_SCHEMA_VERSION:
            raise FM0ContractError("corpus selection schema drift")
        identity = (row["gaia_dr3_source_id"], row["tic_id"], sector)
        if identity in seen_visits:
            raise FM0ContractError(f"duplicate selected visit: {identity}")
        seen_visits.add(identity)
        paths = json.loads(row["hdf5_paths_json"])
        if not isinstance(paths, list) or len(paths) != 2:
            raise FM0ContractError(f"selected visit lacks two HDF5 paths: {identity}")
        requested.extend(_safe_relative_hdf5(str(value), sector=sector) for value in paths)
    if len(requested) != len(set(requested)):
        raise FM0ContractError("selected visits reuse an HDF5 member")
    requested.sort(key=str)
    files_path = partial / "selected_hdf5_files.txt"
    files_path.write_text("".join(f"{path}\n" for path in requested), encoding="utf-8")
    subprocess.run(
        ["tar", "-xf", str(archive), "-C", str(partial), "--files-from", str(files_path)],
        check=True,
    )
    extracted = [(str(path), partial / path) for path in requested]
    if any(not path.is_file() for _, path in extracted):
        raise FM0ContractError("sector tar omitted one or more selected HDF5 members")
    with ThreadPoolExecutor(max_workers=workers) as pool:
        hashed = dict(
            (relative, (size, digest))
            for relative, size, digest in pool.map(_hash_pair, extracted)
        )

    inventory: list[dict[str, Any]] = []
    total_bytes = 0
    for row in selected:
        relative_paths = [str(_safe_relative_hdf5(value, sector=sector)) for value in json.loads(row["hdf5_paths_json"])]
        hashes = [hashed[value][1] for value in relative_paths]
        total_bytes += sum(hashed[value][0] for value in relative_paths)
        inventory.append(
            {
                "gaia_dr3_source_id": row["gaia_dr3_source_id"],
                "tic_id": row["tic_id"],
                "sector": sector,
                "a2v1_product_version": "A2v1",
                "product_state": "A2V1_ACCEPTED",
                "diagnostic_admission_receipt_path": "",
                "diagnostic_admission_receipt_sha256": "",
                "hdf5_paths_json": json.dumps(
                    [str(final / value) for value in relative_paths], separators=(",", ":")
                ),
                "hdf5_sha256_json": json.dumps(hashes, separators=(",", ":")),
                "quality_table_path": str(quality_table_path),
                "quality_table_sha256": quality_table_sha,
                "quality_manifest_path": str(quality_manifest_path),
                "quality_manifest_sha256": quality_manifest_sha,
            }
        )
    inventory.sort(key=lambda row: int(row["gaia_dr3_source_id"]))
    inventory_payload = _csv_bytes(inventory, A2V1_HDF5_SOURCE_INVENTORY_FIELDS)
    publish_immutable(partial / "source_inventory.csv", inventory_payload)
    summary = {
        "schema_version": SECTOR_STAGE_SCHEMA_VERSION,
        "producer_git_sha": producer_git_sha,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sector": sector,
        "archive_sha256": archive_sha,
        "selection_sha256": selection_sha256,
        "quality_table_sha256": quality_table_sha,
        "quality_manifest_sha256": quality_manifest_sha,
        "n_observations": len(inventory),
        "n_hdf5_files": len(requested),
        "total_hdf5_bytes": total_bytes,
        "source_inventory_sha256": hashlib.sha256(inventory_payload).hexdigest(),
        "accepted_sources_modified": False,
        "claim_limit": "source staging only; not a model result",
    }
    publish_immutable(
        partial / "summary.json",
        (json.dumps(summary, indent=2, sort_keys=True) + "\n").encode("utf-8"),
    )
    publish_immutable(partial / "READY", (producer_git_sha + "\n").encode("utf-8"))
    os.replace(partial, final)
    return summary


def merge_sector_inventories(
    *,
    sector_root: str | Path,
    sectors: Sequence[int],
    out_dir: str | Path,
    producer_git_sha: str,
) -> dict[str, Any]:
    """Merge completed sector inventories without rereading HDF5 payloads."""

    root = Path(sector_root).resolve()
    output = Path(out_dir)
    rows: list[dict[str, str]] = []
    sector_bindings: dict[str, Any] = {}
    identities: set[tuple[str, str, int]] = set()
    for sector in sorted({int(value) for value in sectors}):
        sector_dir = root / f"s{sector:04d}"
        inventory = sector_dir / "source_inventory.csv"
        summary_path = sector_dir / "summary.json"
        ready = sector_dir / "READY"
        if not (inventory.is_file() and summary_path.is_file() and ready.is_file()):
            raise FM0ContractError(f"S{sector} sector stage is incomplete")
        if ready.read_text(encoding="utf-8").strip() != producer_git_sha:
            raise FM0ContractError(f"S{sector} sector stage has the wrong code revision")
        sector_rows = _read_csv(inventory)
        if tuple(sector_rows[0]) != A2V1_HDF5_SOURCE_INVENTORY_FIELDS:
            raise FM0ContractError(f"S{sector} source inventory columns drifted")
        for row in sector_rows:
            identity = (row["gaia_dr3_source_id"], row["tic_id"], int(row["sector"]))
            if identity in identities:
                raise FM0ContractError(f"duplicate merged source identity: {identity}")
            identities.add(identity)
            rows.append(row)
        sector_bindings[str(sector)] = {
            "source_inventory_sha256": sha256_file(inventory),
            "summary_sha256": sha256_file(summary_path),
            "n_observations": len(sector_rows),
        }
    rows.sort(key=lambda row: (int(row["sector"]), int(row["gaia_dr3_source_id"])))
    payload = _csv_bytes(rows, A2V1_HDF5_SOURCE_INVENTORY_FIELDS)
    summary = {
        "schema_version": SECTOR_MERGE_SCHEMA_VERSION,
        "producer_git_sha": producer_git_sha,
        "sectors": sorted({int(value) for value in sectors}),
        "n_observations": len(rows),
        "source_inventory_sha256": hashlib.sha256(payload).hexdigest(),
        "sector_bindings": sector_bindings,
        "claim_limit": "source-inventory merge only; not a model result",
    }
    publish_immutable(output / "source_inventory.csv", payload)
    publish_immutable(
        output / "summary.json",
        (json.dumps(summary, indent=2, sort_keys=True) + "\n").encode("utf-8"),
    )
    return summary


__all__ = [
    "SECTOR_MERGE_SCHEMA_VERSION",
    "SECTOR_STAGE_SCHEMA_VERSION",
    "merge_sector_inventories",
    "stage_sector_archive",
]
