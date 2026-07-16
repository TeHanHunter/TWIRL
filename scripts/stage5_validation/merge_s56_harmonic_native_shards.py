#!/usr/bin/env python3
"""Merge and strictly verify parallel S56 harmonic native-input shards."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.harmonic_export import (  # noqa: E402
    merge_raw_pair_shards,
    read_candidate_table,
)
from twirl.vetting.harmonic_dataset import prepare_harmonic_training_rows  # noqa: E402
from twirl.vetting.harmonic_inputs import (  # noqa: E402
    native_group_path,
    verify_raw_pair_contract,
)


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
    source_rows = read_candidate_table(args.training_table)
    if "native_input_include" in source_rows:
        include = source_rows["native_input_include"]
        if include.dtype != bool:
            include = (
                include.fillna("")
                .astype(str)
                .str.lower()
                .isin({"1", "true", "yes", "y"})
            )
        rows = source_rows.loc[include].copy()
        if "native_group_path" not in rows:
            rows["native_group_path"] = [
                native_group_path(row) for row in rows.to_dict("records")
            ]
    else:
        rows = prepare_harmonic_training_rows(source_rows)
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
