#!/usr/bin/env python3
"""Assemble and verify the S56 harmonic-CNN native input HDF5 on ORCD."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.harmonic_export import build_raw_pair_export  # noqa: E402
from twirl.vetting.harmonic_inputs import (  # noqa: E402
    verify_native_candidate_binding,
    verify_raw_pair_contract,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--training-table", type=Path, required=True)
    parser.add_argument("--training-summary", type=Path)
    parser.add_argument("--raw-source-h5", type=Path, required=True)
    parser.add_argument("--compact-adp-h5", type=Path, required=True)
    parser.add_argument("--injection-pair-h5", type=Path)
    parser.add_argument("--cadence-reference-table", type=Path, required=True)
    parser.add_argument("--cadence-reference-manifest", type=Path, required=True)
    parser.add_argument("--out-h5", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument("--n-periods", type=int, default=4096)
    parser.add_argument("--sector", type=int, default=56)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--n-shards", type=int, default=1)
    args = parser.parse_args()
    build = build_raw_pair_export(
        training_table=args.training_table,
        training_summary=args.training_summary,
        raw_source_h5=args.raw_source_h5,
        compact_adp_h5=args.compact_adp_h5,
        injection_pair_h5=args.injection_pair_h5,
        cadence_reference_table=args.cadence_reference_table,
        cadence_reference_manifest=args.cadence_reference_manifest,
        out_h5=args.out_h5,
        repo_root=args.repo_root,
        sector=args.sector,
        n_periods=args.n_periods,
        shard_index=args.shard_index,
        n_shards=args.n_shards,
    )
    verification = (
        verify_native_candidate_binding(
            args.out_h5,
            candidate_table=args.training_table,
            candidate_summary=args.training_summary,
            require_periodograms=True,
            expected_periodogram_n=args.n_periods,
            expected_shard_index=args.shard_index,
            expected_n_shards=args.n_shards,
        )
        if args.training_summary is not None
        else verify_raw_pair_contract(
            args.out_h5,
            require_errors=True,
            require_periodograms=True,
        )
    )
    summary = {"build": build, "verification": verification}
    summary_path = args.out_h5.with_suffix(".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    if not verification["passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
