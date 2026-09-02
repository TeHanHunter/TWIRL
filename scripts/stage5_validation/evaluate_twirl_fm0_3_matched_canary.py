#!/usr/bin/env python3
"""Evaluate the frozen matched TWIRL-FM0.3 canary payload schedule."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.models.fm0.matched_canary_evaluation import (
    MODEL_ARM_CONTRACT,
    ModelArm,
    evaluate_matched_canary,
)
from twirl.models.fm0.validation import require_clean_git_revision


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--payload-plan-root", required=True, type=Path)
    parser.add_argument("--payload-plan-receipt-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-git-sha", required=True)
    parser.add_argument(
        "--controls-only",
        action="store_true",
        help="run only the exact-center raw and quality-only CPU controls",
    )
    parser.add_argument(
        "--model-arm",
        action="append",
        default=[],
        nargs=3,
        metavar=("NAME", "RUN_DIR", "CHECKPOINT"),
        help=(
            "add one frozen checkpoint arm; NAME must be one of: "
            + ", ".join(MODEL_ARM_CONTRACT)
        ),
    )
    parser.add_argument("--device", choices=("auto", "cpu", "cuda"), default="auto")
    parser.add_argument("--batch-size", type=int, default=32)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.controls_only and args.model_arm:
        raise ValueError("--controls-only cannot be combined with --model-arm")
    if not args.controls_only and not args.model_arm:
        raise ValueError("select --controls-only or supply frozen --model-arm entries")
    producer = require_clean_git_revision(ROOT, args.expected_git_sha)
    model_arms = tuple(
        ModelArm(name=name, run_dir=Path(run_dir), checkpoint_path=Path(checkpoint))
        for name, run_dir, checkpoint in args.model_arm
    )
    result = evaluate_matched_canary(
        args.output,
        payload_plan_root=args.payload_plan_root,
        payload_plan_receipt_sha256=args.payload_plan_receipt_sha256,
        producer_git_sha=producer,
        model_arms=model_arms,
        device=args.device,
        batch_size=args.batch_size,
    )
    print(
        json.dumps(
            {
                "root": str(result.root),
                "receipt_sha256": result.receipt_sha256,
                "results_sha256": result.results_sha256,
                "scores_sha256": result.scores_sha256,
                "evaluation_scope": result.receipt["gate_summary"]["evaluation_scope"],
                "controls_preflight_passed": result.receipt["gate_summary"][
                    "controls_preflight"
                ]["passed"],
                "feature_arms": result.receipt["feature_arms"],
                "architecture_promotion_authorized": False,
                "foundation_model_claim_authorized": False,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
