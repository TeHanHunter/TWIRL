#!/usr/bin/env python3
"""Run fail-closed A2v1 Tier-1 QA for bounded enrichment evidence."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.lightcurves.a2v1_tier1_qa import run_a2v1_tier1_qa  # noqa: E402


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector", type=int, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--tier0-summary", type=Path, required=True)
    parser.add_argument("--compact-lc", type=Path, required=True)
    parser.add_argument(
        "--cadence-reference",
        type=Path,
        required=True,
        help="Authoritative detector/orbit cadence-and-quality CSV or Parquet.",
    )
    parser.add_argument(
        "--cadence-reference-manifest",
        type=Path,
        required=True,
        help="Hash-bound QLP-quaternion/SPOC-quality provenance manifest.",
    )
    parser.add_argument(
        "--injection-source-parity",
        type=Path,
        required=True,
        help="Locked JSON proving parity between the injection ADP source and compact product.",
    )
    parser.add_argument("--independent-metrics", type=Path, required=True)
    parser.add_argument("--independent-manifest", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--gate-json", type=Path, required=True)
    injection = parser.add_mutually_exclusive_group(required=True)
    injection.add_argument(
        "--injection-metrics",
        type=Path,
        help="Precomputed model-weighted injection-retention CSV or Parquet.",
    )
    injection.add_argument(
        "--injection-shard",
        type=Path,
        action="append",
        dest="injection_shards",
        help="Fresh A2v1 injection HDF5 shard; repeat for multiple shards.",
    )
    parser.add_argument(
        "--injection-manifest",
        type=Path,
        help="Required with --injection-metrics; contains the frozen contract.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.injection_metrics is not None and args.injection_manifest is None:
        raise SystemExit("--injection-metrics requires --injection-manifest")
    if args.injection_shards and args.injection_manifest is not None:
        raise SystemExit("--injection-manifest is only used with precomputed metrics")
    summary = run_a2v1_tier1_qa(
        sector=args.sector,
        config_path=args.config,
        tier0_summary_path=args.tier0_summary,
        compact_lc=args.compact_lc,
        cadence_reference_path=args.cadence_reference,
        cadence_reference_manifest_path=args.cadence_reference_manifest,
        injection_source_parity_path=args.injection_source_parity,
        injection_metrics_path=args.injection_metrics,
        injection_manifest_path=args.injection_manifest,
        injection_shards=tuple(args.injection_shards or ()),
        independent_metrics_path=args.independent_metrics,
        independent_manifest_path=args.independent_manifest,
        out_dir=args.out_dir,
        gate_json=args.gate_json,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
