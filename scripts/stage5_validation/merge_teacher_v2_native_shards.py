#!/usr/bin/env python3
"""Merge and verify native-input shards for Teacher-v2 candidate tables."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from twirl.vetting.harmonic_export import merge_raw_pair_shards
from twirl.vetting.harmonic_inputs import native_group_path, verify_raw_pair_contract


def _read(path: Path) -> pd.DataFrame:
    return (
        pd.read_parquet(path)
        if path.suffix.lower() == ".parquet"
        else pd.read_csv(path, low_memory=False)
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shards", type=Path, nargs="+", required=True)
    parser.add_argument("--out-h5", type=Path, required=True)
    parser.add_argument("--candidate-table", type=Path, required=True)
    args = parser.parse_args()

    candidates = _read(args.candidate_table)
    if "native_input_include" in candidates:
        include = candidates["native_input_include"]
        if include.dtype != bool:
            include = include.fillna("").astype(str).str.lower().isin(
                {"1", "true", "t", "yes", "y"}
            )
        candidates = candidates.loc[include].copy()
    if "native_group_path" not in candidates:
        candidates["native_group_path"] = [
            native_group_path(row) for row in candidates.to_dict("records")
        ]
    paths = candidates["native_group_path"].fillna("").astype(str)
    if paths.eq("").any():
        raise ValueError("candidate table contains empty native group paths")
    expected = {
        "targets": int(paths[~paths.str.startswith("injections/")].nunique()),
        "injections": int(paths[paths.str.startswith("injections/")].nunique()),
    }
    merge = merge_raw_pair_shards(shard_paths=args.shards, out_h5=args.out_h5)
    verification = verify_raw_pair_contract(
        args.out_h5,
        require_errors=True,
        require_periodograms=True,
    )
    summary = {
        "merge": merge,
        "verification": verification,
        "expected_counts": expected,
        "exact_count_match": merge["counts"] == expected,
    }
    args.out_h5.with_suffix(".summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if verification["passed"] and summary["exact_count_match"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
