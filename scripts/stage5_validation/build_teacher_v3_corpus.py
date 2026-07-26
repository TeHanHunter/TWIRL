#!/usr/bin/env python3
"""Freeze the seven-sector observation-level corpus for Teacher v3."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.vetting.teacher_v3 import write_teacher_v3_corpus


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_S56 = (
    ROOT
    / "reports/stage5_validation/s56_label_adjudication_real343"
    / "adjudicated_training_table/human_vetting_training_table_adjudicated.csv"
)
DEFAULT_S57_S59 = (
    ROOT
    / "reports/stage5_validation/franklin_s57_s59_label_return_20260721"
    / "accepted_morphology_labels.csv"
)
DEFAULT_S60_S62 = (
    ROOT
    / "reports/stage5_validation/franklin_s60_s62_label_return_20260724"
    / "accepted_morphology_labels.csv"
)
DEFAULT_SIGNAL_FREEZE = (
    ROOT
    / "reports/stage5_validation/s56_s62_morphology_corpus_teacher_v3_v1"
    / "accepted_signal_rereview.csv"
)
DEFAULT_OUT_DIR = (
    ROOT
    / "reports/stage5_validation/teacher_v3_s56_s62_a2v1_current_adp"
    / "frozen"
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--s56-training", type=Path, default=DEFAULT_S56)
    parser.add_argument("--s57-s59-labels", type=Path, default=DEFAULT_S57_S59)
    parser.add_argument("--s60-s62-labels", type=Path, default=DEFAULT_S60_S62)
    parser.add_argument(
        "--signal-freeze",
        type=Path,
        default=DEFAULT_SIGNAL_FREEZE,
    )
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    summary = write_teacher_v3_corpus(
        s56_training_path=args.s56_training,
        s57_s59_labels_path=args.s57_s59_labels,
        s60_s62_labels_path=args.s60_s62_labels,
        signal_freeze_path=args.signal_freeze,
        output_path=args.out_dir / "observation_morphology_corpus.csv",
        summary_path=args.out_dir / "observation_morphology_corpus.summary.json",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
