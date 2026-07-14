#!/usr/bin/env python3
"""Build and parity-check the fresh S56 A2v1 injection schedule."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from twirl.injections.a2v1_recovery import (
    build_fresh_injection_schedule,
    load_recovery_config,
    run_adp_roundtrip_parity,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--raw-h5", type=Path, required=True)
    parser.add_argument("--adp-h5", type=Path, required=True)
    parser.add_argument("--teacher-table", type=Path, required=True)
    parser.add_argument(
        "--additional-exclusion-table",
        type=Path,
        action="append",
        default=[],
        help="Optional prior evaluation table whose TICs must also be excluded.",
    )
    parser.add_argument(
        "--host-overlap-audit-table",
        type=Path,
        action="append",
        default=[],
        help="Optional table whose TIC overlap is reported but not excluded.",
    )
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    schedule_path = args.out_dir / "injection_schedule.parquet"
    if schedule_path.exists() and not args.overwrite:
        schedule = pd.read_parquet(schedule_path)
        print(f"[a2v1-schedule] reusing {schedule_path} ({len(schedule):,} rows)")
    else:
        schedule, _, summary = build_fresh_injection_schedule(
            raw_h5=args.raw_h5,
            adp_h5=args.adp_h5,
            teacher_table=args.teacher_table,
            config=load_recovery_config(args.config),
            out_dir=args.out_dir,
            additional_exclusion_tables=args.additional_exclusion_table,
            host_overlap_audit_tables=args.host_overlap_audit_table,
        )
        print(json.dumps(summary, indent=2, sort_keys=True))

    _, parity = run_adp_roundtrip_parity(
        raw_h5=args.raw_h5,
        adp_h5=args.adp_h5,
        schedule=schedule,
        config=load_recovery_config(args.config),
        out_dir=args.out_dir,
    )
    print(json.dumps(parity, indent=2, sort_keys=True))
    if not parity["passed"]:
        print("[a2v1-schedule] ADP parity failed; full injection generation is blocked")
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
