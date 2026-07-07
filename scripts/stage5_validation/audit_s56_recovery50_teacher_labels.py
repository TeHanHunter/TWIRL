#!/usr/bin/env python3
"""Audit human labels and injected-signal visibility for the S56 recovery50 queue."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.vetting.recovery50_teacher import audit_recovery50_labels, json_default  # noqa: E402


DEFAULT_QUEUE = REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue/review_queue_1k.csv"
DEFAULT_LABELS = DEFAULT_QUEUE.with_name("human_labels_vetted.csv")
DEFAULT_OUT_DIR = DEFAULT_QUEUE.with_name("human_label_audit")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue-csv", type=Path, default=DEFAULT_QUEUE)
    parser.add_argument("--labels-csv", type=Path, default=DEFAULT_LABELS)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    summary = audit_recovery50_labels(args.queue_csv, args.labels_csv, args.out_dir)
    print(json.dumps(summary, indent=2, sort_keys=True, default=json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
