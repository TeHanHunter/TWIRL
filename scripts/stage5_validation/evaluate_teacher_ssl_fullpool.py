#!/usr/bin/env python3
"""Evaluate completed full-pool SSL encoders on sealed development folds."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.vetting.teacher_ssl_fullpool_evaluation import (  # noqa: E402
    run_fullpool_ssl_evaluation,
    run_fullpool_ssl_evaluation_preflight,
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training-table", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--registry-summary", type=Path, required=True)
    parser.add_argument("--completion-release", type=Path, required=True)
    parser.add_argument("--baseline-teacher-v3-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--fine-tune-epochs", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=560064)
    parser.add_argument("--linear-probe-steps", type=int, default=500)
    parser.add_argument("--bootstrap-draws", type=int, default=2000)
    parser.add_argument(
        "--preflight-only",
        action="store_true",
        help="Validate authorities and run one bounded encoder forward pass.",
    )
    parser.add_argument(
        "--allow-cpu",
        action="store_true",
        help="Allow a bounded local test without CUDA; production uses one H200.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help=(
            "Contract-check and resume the known interrupted evaluation root; "
            "completed artifacts are validated before reuse."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.preflight_only:
        if args.resume:
            raise ValueError("--resume cannot be combined with --preflight-only")
        summary = run_fullpool_ssl_evaluation_preflight(
            training_table_path=args.training_table,
            registry_path=args.registry,
            registry_summary_path=args.registry_summary,
            completion_release_path=args.completion_release,
            batch_size=int(args.batch_size),
            workers=0,
            seed=int(args.seed),
            require_cuda=not bool(args.allow_cpu),
        )
        print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
        return 0
    summary = run_fullpool_ssl_evaluation(
        training_table_path=args.training_table,
        registry_path=args.registry,
        registry_summary_path=args.registry_summary,
        completion_release_path=args.completion_release,
        baseline_teacher_v3_root=args.baseline_teacher_v3_root,
        output_dir=args.out_dir,
        fine_tune_epochs=int(args.fine_tune_epochs),
        batch_size=int(args.batch_size),
        workers=int(args.workers),
        seed=int(args.seed),
        linear_probe_steps=int(args.linear_probe_steps),
        bootstrap_draws=int(args.bootstrap_draws),
        require_cuda=not bool(args.allow_cpu),
        resume=bool(args.resume),
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
