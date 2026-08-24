#!/usr/bin/env python3
"""Build the matched development-only FM0.1.1 versus FM0.1.2 report plot."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.plotting.fm0_development_comparison import (  # noqa: E402
    make_fm0_development_comparison,
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fm0-1-1-log", type=Path, required=True)
    parser.add_argument("--fm0-1-2-log", type=Path, required=True)
    parser.add_argument("--fm0-1-1-representation", type=Path, required=True)
    parser.add_argument("--fm0-1-2-representation", type=Path, required=True)
    parser.add_argument("--fm0-1-1-summary", type=Path)
    parser.add_argument("--fm0-1-2-summary", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--effective-rank-gate", type=float, default=26.0)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    result = make_fm0_development_comparison(
        fm0_1_1_log_path=args.fm0_1_1_log,
        fm0_1_2_log_path=args.fm0_1_2_log,
        fm0_1_1_representation_path=args.fm0_1_1_representation,
        fm0_1_2_representation_path=args.fm0_1_2_representation,
        fm0_1_1_summary_path=args.fm0_1_1_summary,
        fm0_1_2_summary_path=args.fm0_1_2_summary,
        output_dir=args.output_dir,
        effective_rank_gate=args.effective_rank_gate,
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
