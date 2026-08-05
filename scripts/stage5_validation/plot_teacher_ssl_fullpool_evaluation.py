#!/usr/bin/env python3
"""Render UMAP and matched-development figures for full-pool SSL."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.plotting.teacher_ssl_fullpool_evaluation import (  # noqa: E402
    make_fullpool_evaluation_figures,
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embedding", type=Path, required=True)
    parser.add_argument("--performance", type=Path, required=True)
    parser.add_argument("--per-class", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--random-state", type=int, default=560064)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    summary = make_fullpool_evaluation_figures(
        embedding_path=args.embedding,
        performance_path=args.performance,
        per_class_path=args.per_class,
        output_dir=args.out_dir,
        random_state=int(args.random_state),
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
