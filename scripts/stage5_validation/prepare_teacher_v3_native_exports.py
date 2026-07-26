#!/usr/bin/env python3
"""Prepare split-bound per-sector native-export tables for Teacher v3."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.vetting.teacher_v3_preparation import (
    write_teacher_v3_native_preparation,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ROOT = (
    ROOT
    / "reports/stage5_validation/teacher_v3_s56_s62_a2v1_current_adp"
)
DEFAULT_FROZEN = DEFAULT_ROOT / "frozen"
DEFAULT_NATIVE_ROOT = Path(
    "/orcd/data/mki_aryeh/001/twirl/reports/stage5_validation/"
    "teacher_v3_s56_s62_a2v1_current_adp/native"
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--corpus",
        type=Path,
        default=DEFAULT_FROZEN / "observation_morphology_corpus.csv",
    )
    parser.add_argument(
        "--corpus-summary",
        type=Path,
        default=(
            DEFAULT_FROZEN / "observation_morphology_corpus.summary.json"
        ),
    )
    parser.add_argument(
        "--split-registry",
        type=Path,
        default=DEFAULT_FROZEN / "tic_split_registry.csv",
    )
    parser.add_argument(
        "--split-summary",
        type=Path,
        default=DEFAULT_FROZEN / "tic_split_registry.summary.json",
    )
    parser.add_argument(
        "--native-root",
        type=Path,
        default=DEFAULT_NATIVE_ROOT,
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_ROOT / "native_preparation",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    summary = write_teacher_v3_native_preparation(
        corpus_path=args.corpus,
        corpus_summary_path=args.corpus_summary,
        split_registry_path=args.split_registry,
        split_summary_path=args.split_summary,
        native_root=args.native_root,
        out_dir=args.out_dir,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
