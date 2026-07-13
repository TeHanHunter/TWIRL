#!/usr/bin/env python3
"""Generate one restartable fresh S56 A2v1 injection shard."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import pandas as pd

from twirl.injections.a2v1_recovery import (
    FRESH_INJECTION_CONTRACT,
    load_recovery_config,
    write_fresh_injection_shard,
)


def _complete(path: Path, expected: int) -> bool:
    if not path.exists():
        return False
    try:
        with h5py.File(path, "r") as h5:
            return (
                str(h5.attrs.get("contract_version", "")) == FRESH_INJECTION_CONTRACT
                and len(h5["injections"]) == expected
            )
    except OSError:
        return False


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--raw-h5", type=Path, required=True)
    parser.add_argument("--adp-h5", type=Path, required=True)
    parser.add_argument("--schedule", type=Path, required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--out-h5", type=Path, required=True)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    config = load_recovery_config(args.config)
    expected = (
        config.rows_per_shard
        if args.limit is None
        else min(config.rows_per_shard, args.limit)
    )
    if _complete(args.out_h5, expected) and not args.overwrite:
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
        schedule=schedule,
        shard_index=args.shard_index,
        config=config,
        out_h5=args.out_h5,
        limit=args.limit,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
