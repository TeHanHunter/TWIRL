#!/usr/bin/env python3
"""Audit S56 adjudication labels and build the deduplicated training table."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.vetting.adjudication_audit import run_adjudication_audit  # noqa: E402


TEACHER_ROOT = REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue_2k"
PILOT_ROOT = REPO_ROOT / "reports/stage5_validation/s56_eb_miner_adp_only/review_queue_eb_priority_pilot100"
DEPRECATED_ROOT = REPO_ROOT / "reports/stage5_validation/s56_eb_miner/review_queue_eb_priority_500"
DEFAULT_ROOT = REPO_ROOT / "reports/stage5_validation/s56_label_adjudication_real343"


def parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    ap.add_argument(
        "--active-training-table",
        type=Path,
        default=TEACHER_ROOT / "human_training_table_adp_only/human_vetting_training_table.csv",
    )
    ap.add_argument("--eb-pilot-queue", type=Path, default=PILOT_ROOT / "review_queue_eb_priority_100.csv")
    ap.add_argument("--eb-pilot-labels", type=Path, default=PILOT_ROOT / "human_labels_vetted.csv")
    ap.add_argument(
        "--deprecated-eb-queue", type=Path, default=DEPRECATED_ROOT / "review_queue_eb_priority_500.csv"
    )
    ap.add_argument("--deprecated-eb-labels", type=Path, default=DEPRECATED_ROOT / "human_labels_vetted.csv")
    ap.add_argument(
        "--current-adp-bls",
        type=Path,
        default=REPO_ROOT / "reports/stage5_validation/s56_adp_real_bls_peaks/real_adp_bls_peaks.parquet",
    )
    ap.add_argument("--require-complete", action="store_true")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    summary = run_adjudication_audit(
        root=args.root,
        active_training_table=args.active_training_table,
        eb_pilot_queue=args.eb_pilot_queue,
        eb_pilot_labels=args.eb_pilot_labels,
        deprecated_eb_queue=args.deprecated_eb_queue,
        deprecated_eb_labels=args.deprecated_eb_labels,
        current_adp_bls=args.current_adp_bls,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, default=str))
    if args.require_complete and not summary["cleanup_complete"]:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
