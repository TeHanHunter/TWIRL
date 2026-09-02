#!/usr/bin/env python3
"""Freeze the identity-only matched FM0.3 canary evaluation schedule."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.models.fm0.matched_canary_plan import (
    freeze_matched_canary_evaluation_plan,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--temporal-panel-root", required=True, type=Path)
    parser.add_argument("--temporal-panel-receipt-sha256", required=True)
    parser.add_argument("--composite-root", required=True, type=Path)
    parser.add_argument("--composite-receipt-sha256", required=True)
    parser.add_argument("--composite-source-bindings-sha256", required=True)
    parser.add_argument("--composite-role-index-sha256", required=True)
    parser.add_argument("--producer-git-sha", required=True)
    parser.add_argument("--output", required=True, type=Path)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    result = freeze_matched_canary_evaluation_plan(
        args.output,
        temporal_panel_root=args.temporal_panel_root,
        temporal_panel_receipt_sha256=args.temporal_panel_receipt_sha256,
        composite_root=args.composite_root,
        composite_receipt_sha256=args.composite_receipt_sha256,
        composite_source_bindings_sha256=args.composite_source_bindings_sha256,
        composite_role_index_sha256=args.composite_role_index_sha256,
        producer_git_sha=args.producer_git_sha,
    )
    print(
        json.dumps(
            {
                "root": str(result.root),
                "receipt_sha256": result.receipt_sha256,
                "schedule_sha256": result.schedule_sha256,
                "n_schedule_rows": result.receipt["schedule"]["n_rows"],
                "split_counts": result.receipt["split_contract"]["counts"],
                "identity_only": True,
                "evaluation_executed": False,
                "model_training_authorized": False,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
