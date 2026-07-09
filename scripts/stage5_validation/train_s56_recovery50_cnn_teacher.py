#!/usr/bin/env python3
"""Train the S56 recovery50 CNN teacher on human labels."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.vetting.recovery50_cnn import (  # noqa: E402
    CnnTrainConfig,
    train_recovery50_cnn_teacher,
)
from twirl.vetting.recovery50_teacher import json_default  # noqa: E402

DEFAULT_ROOT = REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue_2k"
DEFAULT_TENSOR_NPZ = DEFAULT_ROOT / "cnn_tensors_adp_only/recovery50_cnn_tensors.npz"
DEFAULT_TENSOR_ROWS = DEFAULT_ROOT / "cnn_tensors_adp_only/recovery50_cnn_tensor_rows.csv"
DEFAULT_TRAINING_TABLE = DEFAULT_ROOT / "human_training_table_adp_only/human_vetting_training_table.csv"
DEFAULT_METRICS: tuple[Path, ...] = ()
DEFAULT_OUT_DIR = DEFAULT_ROOT / "cnn_teacher_adp_only"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tensor-npz", type=Path, default=DEFAULT_TENSOR_NPZ)
    parser.add_argument("--tensor-rows", type=Path, default=DEFAULT_TENSOR_ROWS)
    parser.add_argument("--training-table", type=Path, default=DEFAULT_TRAINING_TABLE)
    parser.add_argument("--metrics-table", type=Path, action="append", default=list(DEFAULT_METRICS))
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--learning-rate", type=float, default=5.0e-4)
    parser.add_argument("--weight-decay", type=float, default=1.0e-4)
    parser.add_argument("--dropout", type=float, default=0.20)
    parser.add_argument("--early-stop-patience", type=int, default=14)
    parser.add_argument("--min-class-count", type=int, default=40)
    parser.add_argument("--validation-fraction", type=float, default=0.20)
    parser.add_argument("--test-fraction", type=float, default=0.20)
    parser.add_argument("--seed", type=int, default=56)
    parser.add_argument("--allow-cpu", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    cfg = CnnTrainConfig(
        epochs=args.epochs,
        batch_size=args.batch_size,
        learning_rate=args.learning_rate,
        weight_decay=args.weight_decay,
        dropout=args.dropout,
        early_stop_patience=args.early_stop_patience,
        min_class_count=args.min_class_count,
        validation_fraction=args.validation_fraction,
        test_fraction=args.test_fraction,
        seed=args.seed,
        require_cuda=not args.allow_cpu,
    )
    summary = train_recovery50_cnn_teacher(
        tensor_npz=args.tensor_npz,
        tensor_rows=args.tensor_rows,
        training_table=args.training_table,
        metrics_tables=tuple(path for path in args.metrics_table if path.exists()),
        out_dir=args.out_dir,
        cfg=cfg,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, default=json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
