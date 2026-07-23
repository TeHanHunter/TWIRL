#!/usr/bin/env python3
"""Freeze a deterministic TIC-grouped split registry for Teacher-v1."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.vetting.teacher_split_registry import (
    DEFAULT_OBSERVATION_IDENTITY_COLUMNS,
    write_tic_split_registry,
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--out-registry", type=Path, required=True)
    parser.add_argument("--out-summary", type=Path, required=True)
    parser.add_argument(
        "--seed",
        type=int,
        required=True,
        help="Declared deterministic split seed.",
    )
    parser.add_argument("--stratum-column", default="split_stratum")
    parser.add_argument(
        "--identity-column",
        action="append",
        default=[],
        help=(
            "Observation identity column; repeat to override the default "
            "(sector, tic, candidate_key)."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    identity_columns = tuple(
        args.identity_column or DEFAULT_OBSERVATION_IDENTITY_COLUMNS
    )
    summary = write_tic_split_registry(
        corpus_path=args.corpus,
        registry_path=args.out_registry,
        summary_path=args.out_summary,
        seed=args.seed,
        stratum_column=args.stratum_column,
        identity_columns=identity_columns,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
