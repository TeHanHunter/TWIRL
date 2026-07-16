#!/usr/bin/env python3
"""Build the A2v1 Tier-0 integrity/benchmark gate before teacher inference."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.lightcurves.a2v1_qa import run_a2v1_photometric_qa  # noqa: E402


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector", type=int, required=True)
    parser.add_argument("--orbits", type=int, nargs="+", required=True)
    parser.add_argument("--compact-lc", type=Path, required=True)
    parser.add_argument("--schema-summary", type=Path, required=True)
    parser.add_argument("--bls-peaks", type=Path, required=True)
    parser.add_argument("--a2v1-root", type=Path, required=True)
    parser.add_argument("--reference-raw-root", type=Path, required=True)
    parser.add_argument("--reference-qa-summary", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--gate-json", type=Path, required=True)
    parser.add_argument("--sample-size", type=int, default=256)
    parser.add_argument("--seed", type=int, default=56)
    parser.add_argument("--workers", type=int, default=8)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.sector != 56 and args.reference_qa_summary is None:
        raise SystemExit("S57+ QA requires --reference-qa-summary for the passed S56 benchmark")
    summary = run_a2v1_photometric_qa(
        sector=args.sector,
        orbits=args.orbits,
        compact_lc=args.compact_lc,
        schema_summary=args.schema_summary,
        bls_peaks=args.bls_peaks,
        a2v1_root=args.a2v1_root,
        reference_raw_root=args.reference_raw_root,
        reference_qa_summary=args.reference_qa_summary,
        out_dir=args.out_dir,
        gate_json=args.gate_json,
        sample_size=args.sample_size,
        seed=args.seed,
        workers=args.workers,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
