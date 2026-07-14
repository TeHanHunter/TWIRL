#!/usr/bin/env python3
"""Verify transferred sector A2v1 recovery inputs before ORCD compute."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import h5py
import pandas as pd

from twirl.injections.a2v1_recovery import compare_adp_compact_products


EXPECTED_BY_SECTOR = {
    56: {
        "fits_count": 31_450,
        "fits_bytes": 27_080_968_320,
        "raw_bytes": 7_176_576_282,
        "adp_bytes": 4_379_183_054,
        "raw_targets": 31_446,
        "adp_targets": 31_450,
        "teacher_rows": 2_159,
        "teacher_unique_tics": 2_044,
        "adp_without_raw": 4,
        "raw_without_adp": 0,
    },
    57: {
        "fits_count": 27_213,
        "fits_bytes": 21_552_696_000,
        "raw_bytes": 5_818_577_929,
        "adp_bytes": 3_586_289_254,
        "raw_targets": 27_212,
        "adp_targets": 27_213,
        "teacher_rows": 2_159,
        "teacher_unique_tics": 2_044,
        "adp_without_raw": 1,
        "raw_without_adp": 0,
    },
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _write_fits_checksum_manifest(paths: list[Path], root: Path, out_csv: Path) -> str:
    tree_digest = hashlib.sha256()
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("relative_path", "bytes", "sha256"))
        writer.writeheader()
        for index, path in enumerate(paths, start=1):
            relative = path.relative_to(root).as_posix()
            sha256 = _sha256(path)
            size = path.stat().st_size
            writer.writerow(
                {"relative_path": relative, "bytes": size, "sha256": sha256}
            )
            tree_digest.update(relative.encode("utf-8"))
            tree_digest.update(b"\0")
            tree_digest.update(str(size).encode("ascii"))
            tree_digest.update(b"\0")
            tree_digest.update(sha256.encode("ascii"))
            tree_digest.update(b"\n")
            if index % 500 == 0:
                print(f"[fits-sha256] {index:,}/{len(paths):,}", flush=True)
    return tree_digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector", type=int, default=56)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--hash-large-files", action="store_true")
    parser.add_argument("--fits-checksum-manifest", type=Path)
    args = parser.parse_args()

    if args.sector not in EXPECTED_BY_SECTOR:
        raise SystemExit(f"no locked input expectations for Sector {args.sector}")
    expected = EXPECTED_BY_SECTOR[args.sector]
    root = args.input_root
    sector_tag = f"s{args.sector}"
    sector_padded = f"s{args.sector:04d}"
    fits_root = root / f"hlsp_{sector_padded}_A2v1"
    raw_path = root / f"{sector_tag}_A2v1_tglc_raw_sources.h5"
    adp_reference_path = root / f"{sector_tag}_A2v1_adp_pair.h5"
    adp_active_path = root / f"{sector_tag}_A2v1_adp_pair_rebuilt.h5"
    teacher_path = root / "human_vetting_training_table_adjudicated.csv"
    validation_path = root / f"{sector_tag}_A2v1_validation_full_schema.json"
    fits = sorted(fits_root.rglob("*.fits"))
    fits_bytes = sum(path.stat().st_size for path in fits)
    teacher = pd.read_csv(teacher_path, low_memory=False)
    validation = json.loads(validation_path.read_text())
    with h5py.File(raw_path, "r") as raw, h5py.File(adp_reference_path, "r") as adp:
        raw_keys = set(raw["targets"].keys())
        adp_keys = set(adp["targets"].keys())
        observed = {
            "fits_count": len(fits),
            "fits_bytes": fits_bytes,
            "raw_bytes": raw_path.stat().st_size,
            "adp_bytes": adp_reference_path.stat().st_size,
            "active_adp_bytes": adp_active_path.stat().st_size,
            "raw_targets": len(raw_keys),
            "adp_targets": len(adp_keys),
            "teacher_rows": len(teacher),
            "teacher_unique_tics": int(
                pd.to_numeric(teacher["tic"], errors="coerce").nunique()
            ),
            "adp_without_raw": len(adp_keys - raw_keys),
            "raw_without_adp": len(raw_keys - adp_keys),
            "raw_contract_version": str(raw.attrs.get("contract_version", "")),
            "adp_sector": int(adp.attrs.get("sector", -1)),
            "adp_flux_columns": str(adp.attrs.get("flux_columns", "")),
        }
    failures = [
        f"{key}={observed[key]!r}; expected {value!r}"
        for key, value in expected.items()
        if observed[key] != value
    ]
    if observed["adp_sector"] != args.sector:
        failures.append("ADP HDF5 sector does not match the requested sector")
    accepted_raw_contracts = {
        "s56_tglc_raw_pair_v1",
        f"s{args.sector}_tglc_raw_pair_v1",
        "a2v1_tglc_raw_pair_v1",
    }
    if observed["raw_contract_version"] not in accepted_raw_contracts:
        failures.append("raw HDF5 has the wrong contract version")
    if (
        "DET_FLUX_ADP_SML" not in observed["adp_flux_columns"]
        or "DET_FLUX_ADP" not in observed["adp_flux_columns"]
    ):
        failures.append("ADP HDF5 does not declare the two required apertures")
    if not bool(validation.get("ok", False)):
        failures.append("source A2v1 full-schema validation is not ok")
    compact_parity = compare_adp_compact_products(
        adp_reference_path,
        adp_active_path,
    )
    if not compact_parity["passed"]:
        failures.extend(compact_parity["failures"])
    hashes = {}
    fits_checksum_manifest = None
    if args.hash_large_files:
        for path in (
            raw_path,
            adp_reference_path,
            adp_active_path,
            teacher_path,
            validation_path,
        ):
            hashes[path.name] = _sha256(path)
        fits_checksum_manifest = (
            args.fits_checksum_manifest
            or root / f"hlsp_{sector_padded}_A2v1_sha256_manifest.csv"
        )
        hashes[f"hlsp_{sector_padded}_A2v1_tree"] = _write_fits_checksum_manifest(
            fits,
            fits_root,
            fits_checksum_manifest,
        )
    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sector": args.sector,
        "input_root": str(root.resolve()),
        "expected": expected,
        "observed": observed,
        "hashes": hashes,
        "fits_checksum_manifest": (
            str(fits_checksum_manifest.resolve()) if fits_checksum_manifest else None
        ),
        "compact_rebuild_parity": compact_parity,
        "source_validation_ok": bool(validation.get("ok", False)),
        "passed": not failures,
        "failures": failures,
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0 if payload["passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
