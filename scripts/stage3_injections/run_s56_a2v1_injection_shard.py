#!/usr/bin/env python3
"""Generate one restartable fresh S56 A2v1 injection shard."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import pandas as pd

from twirl.injections.a2v1_recovery import (
    EPOCH_QUALITY_POLICY_CONTRACT,
    fresh_injection_contract,
    load_recovery_config,
    write_fresh_injection_shard,
)
from twirl.lightcurves.a2v1_cadence_reference import file_sha256


def _complete(
    path: Path,
    expected: int,
    selection_mode: str,
    contract: str,
    cadence_reference_table_sha256: str,
    cadence_reference_manifest_sha256: str,
) -> bool:
    if not path.exists():
        return False
    try:
        with h5py.File(path, "r") as h5:
            return (
                str(h5.attrs.get("contract_version", "")) == contract
                and len(h5["injections"]) == expected
                and str(h5.attrs.get("selection_mode", "shard")) == selection_mode
                and str(h5.attrs.get("epoch_quality_policy_contract", ""))
                == EPOCH_QUALITY_POLICY_CONTRACT
                and str(h5.attrs.get("cadence_reference_table_sha256", ""))
                == cadence_reference_table_sha256
                and str(h5.attrs.get("cadence_reference_manifest_sha256", ""))
                == cadence_reference_manifest_sha256
            )
    except (KeyError, OSError):
        return False


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--raw-h5", type=Path, required=True)
    parser.add_argument("--adp-h5", type=Path, required=True)
    parser.add_argument("--cadence-reference-table", type=Path, required=True)
    parser.add_argument("--cadence-reference-manifest", type=Path, required=True)
    parser.add_argument("--schedule", type=Path, required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--out-h5", type=Path, required=True)
    parser.add_argument("--limit", type=int)
    parser.add_argument(
        "--selection-mode",
        choices=("shard", "parameter_spanning"),
        default="shard",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    config = load_recovery_config(args.config)
    contract = fresh_injection_contract(config.sector)
    expected = (
        config.rows_per_shard
        if args.limit is None
        else min(config.rows_per_shard, args.limit)
    )
    cadence_reference_table_sha256 = file_sha256(args.cadence_reference_table)
    cadence_reference_manifest_sha256 = file_sha256(
        args.cadence_reference_manifest
    )
    if (
        _complete(
            args.out_h5,
            expected,
            args.selection_mode,
            contract,
            cadence_reference_table_sha256,
            cadence_reference_manifest_sha256,
        )
        and not args.overwrite
    ):
        print(f"[a2v1-injection-shard] complete output exists: {args.out_h5}")
        return 0
    if args.out_h5.exists() and not args.overwrite:
        raise SystemExit(
            f"incomplete output exists; inspect or pass --overwrite: {args.out_h5}"
        )
    schedule = (
        pd.read_parquet(args.schedule)
        if args.schedule.suffix.lower() == ".parquet"
        else pd.read_csv(args.schedule, low_memory=False)
    )
    summary = write_fresh_injection_shard(
        raw_h5=args.raw_h5,
        adp_h5=args.adp_h5,
        cadence_reference_table=args.cadence_reference_table,
        cadence_reference_manifest=args.cadence_reference_manifest,
        schedule=schedule,
        shard_index=args.shard_index,
        config=config,
        out_h5=args.out_h5,
        limit=args.limit,
        selection_mode=args.selection_mode,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
