#!/usr/bin/env python3
"""Build candidate metadata and native Teacher-v1 inputs for one injection shard."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from twirl.vetting.harmonic_export import build_raw_pair_export
from twirl.vetting.harmonic_inputs import verify_raw_pair_contract
from twirl.vetting.injection_teacher_recovery import (
    enrich_injection_candidate_metadata,
    normalize_injection_peak_candidates,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peak-table", type=Path, required=True)
    parser.add_argument("--pair-h5", type=Path, required=True)
    parser.add_argument("--raw-h5", type=Path, required=True)
    parser.add_argument("--adp-h5", type=Path, required=True)
    parser.add_argument("--out-candidates", type=Path, required=True)
    parser.add_argument("--out-native-h5", type=Path, required=True)
    parser.add_argument("--n-periods", type=int, default=4096)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    if (
        args.out_candidates.exists() or args.out_native_h5.exists()
    ) and not args.overwrite:
        raise SystemExit("teacher-shard output exists; inspect or pass --overwrite")
    peaks = (
        pd.read_parquet(args.peak_table)
        if args.peak_table.suffix == ".parquet"
        else pd.read_csv(args.peak_table, low_memory=False)
    )
    candidates = normalize_injection_peak_candidates(peaks, pair_h5=args.pair_h5)
    candidates, metadata_summary = enrich_injection_candidate_metadata(
        candidates,
        pair_h5=args.pair_h5,
    )
    if not metadata_summary["passed"]:
        raise RuntimeError("candidate metadata generation was incomplete")
    args.out_candidates.parent.mkdir(parents=True, exist_ok=True)
    candidates.to_csv(args.out_candidates, index=False)
    build = build_raw_pair_export(
        training_table=args.out_candidates,
        raw_source_h5=args.raw_h5,
        compact_adp_h5=args.adp_h5,
        injection_pair_h5=args.pair_h5,
        out_h5=args.out_native_h5,
        repo_root=Path.cwd(),
        n_periods=args.n_periods,
    )
    verification = verify_raw_pair_contract(
        args.out_native_h5,
        require_errors=True,
        require_periodograms=True,
    )
    summary = {
        "n_candidates": int(len(candidates)),
        "n_injections": int(candidates["injection_id"].nunique()),
        "metadata": metadata_summary,
        "native_build": build,
        "native_verification": verification,
    }
    args.out_native_h5.with_suffix(".teacher_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if verification["passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
