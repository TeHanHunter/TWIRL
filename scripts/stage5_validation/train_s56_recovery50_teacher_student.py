#!/usr/bin/env python3
"""Train first recovery50 teacher/student smoke models for S56."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.vetting.recovery50_teacher import (  # noqa: E402
    json_default,
    load_feature_table,
    train_teacher_student_smoke,
    write_table,
)

DEFAULT_ROOT = REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue"
DEFAULT_TRAINING_TABLE = DEFAULT_ROOT / "human_training_table/human_vetting_training_table.csv"
DEFAULT_SHAPE_FEATURES = DEFAULT_ROOT / "folded_shape_features/folded_shape_features.csv"
DEFAULT_METRICS = (
    DEFAULT_ROOT / "twirl_vet_metrics_real_fullphase_binmatch.csv",
    DEFAULT_ROOT / "twirl_vet_metrics_injected_fullphase_binmatch.csv",
)
DEFAULT_OUT_DIR = DEFAULT_ROOT / "teacher_smoke"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training-table", type=Path, default=DEFAULT_TRAINING_TABLE)
    parser.add_argument("--shape-features", type=Path, default=DEFAULT_SHAPE_FEATURES)
    parser.add_argument("--metrics-table", type=Path, action="append", default=list(DEFAULT_METRICS))
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--min-class-count", type=int, default=40)
    parser.add_argument("--pseudo-min-confidence", type=float, default=0.98)
    parser.add_argument("--pseudo-min-margin", type=float, default=0.50)
    parser.add_argument("--pseudo-weight", type=float, default=0.25)
    parser.add_argument("--random-state", type=int, default=56)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    feature_table = load_feature_table(
        training_table=args.training_table,
        shape_features=args.shape_features if args.shape_features.exists() else None,
        metrics_tables=tuple(path for path in args.metrics_table if path.exists()),
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    feature_table_path = write_table(feature_table, args.out_dir / "model_feature_table.parquet")
    summary = train_teacher_student_smoke(
        feature_table=feature_table,
        out_dir=args.out_dir,
        min_class_count=args.min_class_count,
        pseudo_min_confidence=args.pseudo_min_confidence,
        pseudo_min_margin=args.pseudo_min_margin,
        pseudo_weight=args.pseudo_weight,
        random_state=args.random_state,
    )
    summary = {**summary, "feature_table": str(feature_table_path)}
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=json_default) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True, default=json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
