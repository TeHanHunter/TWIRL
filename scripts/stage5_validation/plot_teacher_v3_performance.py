#!/usr/bin/env python3
"""Render the frozen Teacher v3 fixed-test performance comparison."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.plotting.teacher_v3_performance import (  # noqa: E402
    build_teacher_v3_performance_figure,
    validate_teacher_v3_performance_inputs,
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Completed Teacher v3 `full` output directory.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        help="Output directory; defaults to RUN_DIR/performance.",
    )
    parser.add_argument(
        "--expected-bootstrap-draws",
        type=int,
        default=2000,
    )
    parser.add_argument(
        "--expected-bootstrap-seed",
        type=int,
        default=560063,
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Validate the completed run and write no figure products.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    run_dir = args.run_dir.expanduser().resolve()
    out_dir = (
        args.out_dir.expanduser().resolve()
        if args.out_dir is not None
        else run_dir / "performance"
    )
    if args.check_only:
        result = validate_teacher_v3_performance_inputs(
            run_dir=run_dir,
            expected_draws=int(args.expected_bootstrap_draws),
            expected_seed=int(args.expected_bootstrap_seed),
        )
    else:
        result = build_teacher_v3_performance_figure(
            run_dir=run_dir,
            out_dir=out_dir,
            expected_draws=int(args.expected_bootstrap_draws),
            expected_seed=int(args.expected_bootstrap_seed),
        )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
