#!/usr/bin/env python3
"""Evaluate exact FM0 step-0/step-2000 checkpoints on frozen S66--S77."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.models.fm0.temporal_zero_shot import (
    evaluate_temporal_zero_shot,
)
from twirl.models.fm0.validation import (
    require_clean_git_revision,
    write_json_with_sha256,
)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--step0-checkpoint", type=Path, required=True)
    parser.add_argument("--step2000-checkpoint", type=Path, required=True)
    parser.add_argument("--temporal-panel-dir", type=Path, required=True)
    parser.add_argument("--temporal-panel-receipt-sha256", required=True)
    parser.add_argument("--baseline-manifest", type=Path, required=True)
    parser.add_argument("--baseline-manifest-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--max-repeated-components", type=int, default=256)
    parser.add_argument("--max-new-components", type=int, default=256)
    parser.add_argument("--max-baseline-train-components", type=int, default=512)
    parser.add_argument("--pca-components", type=int, default=32)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--expected-git-sha", required=True)
    return parser


def main() -> int:
    args = _parser().parse_args()
    if (
        args.output.exists()
        or args.output.with_name(args.output.name + ".sha256").exists()
    ):
        raise FileExistsError("refusing to overwrite temporal zero-shot output")
    require_clean_git_revision(ROOT, args.expected_git_sha)
    payload = evaluate_temporal_zero_shot(
        run_dir=args.run_dir,
        step0_checkpoint_path=args.step0_checkpoint,
        step2000_checkpoint_path=args.step2000_checkpoint,
        temporal_panel_dir=args.temporal_panel_dir,
        temporal_panel_receipt_sha256=args.temporal_panel_receipt_sha256,
        baseline_manifest_path=args.baseline_manifest,
        baseline_manifest_sha256=args.baseline_manifest_sha256,
        device=args.device,
        max_repeated_components=args.max_repeated_components,
        max_new_components=args.max_new_components,
        max_baseline_train_components=args.max_baseline_train_components,
        pca_components=args.pca_components,
        batch_size=args.batch_size,
    )
    payload["evaluator_git_sha"] = args.expected_git_sha
    write_json_with_sha256(args.output, payload)
    os.chmod(args.output, 0o444)
    os.chmod(args.output.with_name(args.output.name + ".sha256"), 0o444)
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
