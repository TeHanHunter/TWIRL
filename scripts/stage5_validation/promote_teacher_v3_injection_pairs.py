#!/usr/bin/env python3
"""Promote the two frozen S56 Teacher-v3 injection-pair files into one."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.teacher_v3_injection_pairs import (  # noqa: E402
    TEACHER_V3_S56_EXPECTED_INJECTIONS,
    TEACHER_V3_S56_EXPECTED_SOURCE_COUNTS,
    promote_teacher_v3_injection_pairs,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection-table", type=Path, required=True)
    parser.add_argument(
        "--source-pair-h5",
        type=Path,
        nargs=2,
        metavar=("PAIR_300_H5", "PAIR_156_H5"),
        required=True,
    )
    parser.add_argument(
        "--canonical-injection-h5",
        type=Path,
        required=True,
        help="Stable raw-injection HDF5 containing all selected injection groups.",
    )
    parser.add_argument("--out-h5", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument(
        "--expected-count",
        type=int,
        default=TEACHER_V3_S56_EXPECTED_INJECTIONS,
    )
    parser.add_argument(
        "--expected-source-counts",
        type=int,
        nargs=2,
        default=TEACHER_V3_S56_EXPECTED_SOURCE_COUNTS,
        metavar=("FIRST", "SECOND"),
    )
    parser.add_argument("--expected-sector", type=int, default=56)
    parser.add_argument("--summary-json", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    summary = promote_teacher_v3_injection_pairs(
        selection_table=args.selection_table,
        source_pair_h5s=args.source_pair_h5,
        canonical_injection_h5=args.canonical_injection_h5,
        out_h5=args.out_h5,
        repo_root=args.repo_root,
        expected_count=args.expected_count,
        expected_source_counts=args.expected_source_counts,
        expected_sector=args.expected_sector,
        overwrite=args.overwrite,
    )
    summary_path = (
        args.summary_json
        if args.summary_json is not None
        else args.out_h5.with_suffix(".summary.json")
    )
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
