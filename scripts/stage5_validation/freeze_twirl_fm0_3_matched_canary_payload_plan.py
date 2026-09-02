#!/usr/bin/env python3
"""Freeze the payload-screened direct-128 matched FM0.3 crop schedule."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.models.fm0.matched_canary_payload_plan import (
    freeze_matched_canary_payload_plan,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--identity-plan-root", required=True, type=Path)
    parser.add_argument("--identity-plan-receipt-sha256", required=True)
    parser.add_argument("--producer-git-sha", required=True)
    parser.add_argument("--output", required=True, type=Path)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    result = freeze_matched_canary_payload_plan(
        args.output,
        identity_plan_root=args.identity_plan_root,
        identity_plan_receipt_sha256=args.identity_plan_receipt_sha256,
        producer_git_sha=args.producer_git_sha,
    )
    print(
        json.dumps(
            {
                "root": str(result.root),
                "receipt_sha256": result.receipt_sha256,
                "schedule_sha256": result.schedule_sha256,
                "n_schedule_rows": result.receipt["schedule"]["n_rows"],
                "source_payloads_opened": True,
                "events_injected": False,
                "checkpoints_loaded": False,
                "evaluation_executed": False,
                "model_training_executed": False,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
