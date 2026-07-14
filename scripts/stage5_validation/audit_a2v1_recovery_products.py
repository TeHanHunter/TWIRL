#!/usr/bin/env python3
"""Audit complete A2v1 injection shards against their source light curves."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from twirl.injections.a2v1_recovery import (
    audit_fresh_injection_shards,
    load_recovery_config,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--schedule", type=Path, required=True)
    parser.add_argument("--shards", type=Path, nargs="+", required=True)
    parser.add_argument("--raw-h5", type=Path, required=True)
    parser.add_argument("--adp-h5", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--relative-tolerance", type=float, default=1.0e-6)
    args = parser.parse_args()

    schedule = (
        pd.read_parquet(args.schedule)
        if args.schedule.suffix.lower() == ".parquet"
        else pd.read_csv(args.schedule, low_memory=False)
    )
    payload = audit_fresh_injection_shards(
        shard_paths=args.shards,
        schedule=schedule,
        raw_h5=args.raw_h5,
        adp_h5=args.adp_h5,
        config=load_recovery_config(args.config),
        relative_tolerance=args.relative_tolerance,
    )
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0 if payload["passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
