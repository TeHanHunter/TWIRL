#!/usr/bin/env python3
"""Merge all fresh-recovery candidate and native-input shards."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from twirl.vetting.harmonic_export import merge_raw_pair_shards
from twirl.vetting.harmonic_inputs import verify_raw_pair_contract


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-shards", type=Path, nargs="+", required=True)
    parser.add_argument("--native-shards", type=Path, nargs="+", required=True)
    parser.add_argument("--out-candidates", type=Path, required=True)
    parser.add_argument("--out-native-h5", type=Path, required=True)
    parser.add_argument("--maximum-injections", type=int, default=20_000)
    args = parser.parse_args()

    candidates = pd.concat(
        [pd.read_csv(path, low_memory=False) for path in sorted(args.candidate_shards)],
        ignore_index=True,
    )
    if candidates["review_id"].duplicated().any():
        raise ValueError("candidate shards contain duplicate review_id values")
    n_candidate_injections = int(candidates["injection_id"].nunique())
    if not 0 < n_candidate_injections <= args.maximum_injections:
        raise ValueError(
            f"candidate shards cover {n_candidate_injections:,} injections; "
            f"expected between 1 and {args.maximum_injections:,}"
        )
    args.out_candidates.parent.mkdir(parents=True, exist_ok=True)
    candidates.to_csv(args.out_candidates, index=False)
    merge = merge_raw_pair_shards(
        shard_paths=sorted(args.native_shards),
        out_h5=args.out_native_h5,
    )
    verification = verify_raw_pair_contract(
        args.out_native_h5,
        require_errors=True,
        require_periodograms=True,
    )
    expected_counts = {"targets": 0, "injections": n_candidate_injections}
    passed = verification["passed"] and merge["counts"] == expected_counts
    summary = {
        "n_candidate_rows": int(len(candidates)),
        "n_injections": n_candidate_injections,
        "merge": merge,
        "verification": verification,
        "expected_counts": expected_counts,
        "passed": bool(passed),
    }
    args.out_native_h5.with_suffix(".merge_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if passed else 3


if __name__ == "__main__":
    raise SystemExit(main())
