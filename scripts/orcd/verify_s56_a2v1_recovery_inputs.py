#!/usr/bin/env python3
"""Verify the transferred S56 A2v1 recovery inputs before ORCD compute."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import h5py
import pandas as pd

from twirl.injections.a2v1_recovery import compare_adp_compact_products


EXPECTED = {
    "fits_count": 31_450,
    "fits_bytes": 27_080_968_320,
    "raw_bytes": 7_176_576_282,
    "adp_bytes": 4_379_183_054,
    "raw_targets": 31_446,
    "adp_targets": 31_450,
    "teacher_rows": 2_159,
    "teacher_unique_tics": 2_044,
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--hash-large-files", action="store_true")
    args = parser.parse_args()

    root = args.input_root
    fits_root = root / "hlsp_s0056_A2v1"
    raw_path = root / "s56_A2v1_tglc_raw_sources.h5"
    adp_reference_path = root / "s56_A2v1_adp_pair.h5"
    adp_active_path = root / "s56_A2v1_adp_pair_rebuilt.h5"
    teacher_path = root / "human_vetting_training_table_adjudicated.csv"
    validation_path = root / "s56_A2v1_validation_full_schema.json"
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
            "adp_flux_columns": str(adp.attrs.get("flux_columns", "")),
        }
    failures = [
        f"{key}={observed[key]!r}; expected {value!r}"
        for key, value in EXPECTED.items()
        if observed[key] != value
    ]
    if observed["adp_without_raw"] != 4 or observed["raw_without_adp"] != 0:
        failures.append(
            "raw/ADP intersection does not have the expected four FITS-only targets"
        )
    if observed["raw_contract_version"] != "s56_tglc_raw_pair_v1":
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
    if args.hash_large_files:
        for path in (
            raw_path,
            adp_reference_path,
            adp_active_path,
            teacher_path,
            validation_path,
        ):
            hashes[path.name] = _sha256(path)
    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_root": str(root.resolve()),
        "expected": EXPECTED,
        "observed": observed,
        "hashes": hashes,
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
