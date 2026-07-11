#!/usr/bin/env python3
"""Merge and strictly verify parallel S56 harmonic native-input shards."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.harmonic_export import merge_raw_pair_shards  # noqa: E402
from twirl.vetting.harmonic_dataset import prepare_harmonic_training_rows  # noqa: E402
from twirl.vetting.harmonic_inputs import verify_raw_pair_contract  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shards", type=Path, nargs="+", required=True)
    parser.add_argument("--out-h5", type=Path, required=True)
    parser.add_argument("--training-table", type=Path, required=True)
    args = parser.parse_args()
    merge = merge_raw_pair_shards(shard_paths=args.shards, out_h5=args.out_h5)
    verification = verify_raw_pair_contract(
        args.out_h5,
        require_errors=True,
        require_periodograms=True,
    )
    rows = prepare_harmonic_training_rows(pd.read_csv(args.training_table, low_memory=False))
    injection = rows["native_group_path"].str.startswith("injections/")
    expected = {
        "targets": int(rows.loc[~injection, "native_group_path"].nunique()),
        "injections": int(rows.loc[injection, "native_group_path"].nunique()),
    }
    count_match = merge["counts"] == expected
    summary = {
        "merge": merge,
        "verification": verification,
        "expected_counts": expected,
        "exact_count_match": count_match,
    }
    args.out_h5.with_suffix(".summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    if not verification["passed"] or not count_match:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
