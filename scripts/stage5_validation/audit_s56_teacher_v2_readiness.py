#!/usr/bin/env python3
"""Gate the five-way teacher and student on real-label support and test metrics."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.teacher_active_learning import teacher_v2_readiness  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--test-metrics", type=Path)
    parser.add_argument("--unresolved-ephemerides", type=int, default=0)
    parser.add_argument("--join-mismatches", type=int, default=0)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    labels = (
        pd.read_parquet(args.labels)
        if args.labels.suffix.lower() == ".parquet"
        else pd.read_csv(args.labels, low_memory=False)
    )
    metrics = json.loads(args.test_metrics.read_text()) if args.test_metrics else None
    summary = teacher_v2_readiness(
        labels,
        test_metrics=metrics,
        unresolved_ephemerides=args.unresolved_ephemerides,
        join_mismatches=args.join_mismatches,
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=True) + "\n"
    )
    lines = [
        "# S56 Teacher v2 / Student Readiness",
        "",
        f"- Teacher outputs: `{summary['teacher_output_count']}`",
        f"- Broad class promoted: `{summary['broad_class_promoted']}`",
        f"- Student ready: `{summary['student_ready']}`",
        "",
        "## Checks",
        "",
    ]
    lines.extend(
        f"- {name}: `{passed}`" for name, passed in summary["student_checks"].items()
    )
    (args.out_dir / "summary.md").write_text("\n".join(lines) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=True))


if __name__ == "__main__":
    main()
