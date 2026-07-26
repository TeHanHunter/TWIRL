#!/usr/bin/env python3
"""Merge and strictly verify one sector's Teacher-v3 native-input shards."""
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
from twirl.vetting.harmonic_inputs import (  # noqa: E402
    native_group_path,
    verify_raw_pair_contract,
)


def _truth(values):
    if values.dtype == bool:
        return values.fillna(False)
    return (
        values.fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        .isin({"1", "true", "t", "yes", "y"})
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector", type=int, required=True)
    parser.add_argument("--shards", type=Path, nargs="+", required=True)
    parser.add_argument("--out-h5", type=Path, required=True)
    parser.add_argument("--training-table", type=Path, required=True)
    args = parser.parse_args()
    if args.sector < 56 or args.sector > 62:
        raise SystemExit("Teacher-v3 native merge is bounded to sectors 56-62")
    merge = merge_raw_pair_shards(
        shard_paths=args.shards,
        out_h5=args.out_h5,
    )
    verification = verify_raw_pair_contract(
        args.out_h5,
        require_errors=True,
        require_periodograms=True,
    )
    source_rows = read_candidate_table(args.training_table)
    if "native_input_include" not in source_rows:
        raise KeyError("Teacher-v3 sector table lacks native_input_include")
    rows = source_rows.loc[_truth(source_rows["native_input_include"])].copy()
    if "native_group_path" not in rows:
        rows["native_group_path"] = [
            native_group_path(row) for row in rows.to_dict("records")
        ]
    injection = rows["native_group_path"].str.startswith("injections/")
    expected = {
        "targets": int(rows.loc[~injection, "native_group_path"].nunique()),
        "injections": int(rows.loc[injection, "native_group_path"].nunique()),
    }
    count_match = merge["counts"] == expected
    summary = {
        "sector": int(args.sector),
        "merge": merge,
        "verification": verification,
        "expected_counts": expected,
        "exact_count_match": count_match,
    }
    args.out_h5.with_suffix(".summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    if not verification["passed"] or not count_match:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
