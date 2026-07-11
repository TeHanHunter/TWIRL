#!/usr/bin/env python3
"""Train and evaluate the post-adjudication S56 harmonic CNN teacher."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.harmonic_cnn import HarmonicTrainConfig  # noqa: E402
from twirl.vetting.harmonic_training import (  # noqa: E402
    DEFAULT_PROFILES,
    run_harmonic_teacher_training,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--training-table", type=Path, required=True)
    parser.add_argument("--native-h5", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--profiles", default=",".join(DEFAULT_PROFILES))
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--pretrain-epochs", type=int, default=20)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--learning-rate", type=float, default=3.0e-4)
    parser.add_argument("--weight-decay", type=float, default=1.0e-4)
    parser.add_argument("--patience", type=int, default=12)
    parser.add_argument("--seed", type=int, default=56)
    parser.add_argument("--allow-cpu", action="store_true")
    args = parser.parse_args()
    profiles = tuple(value.strip() for value in args.profiles.split(",") if value.strip())
    config = HarmonicTrainConfig(
        epochs=args.epochs,
        batch_size=args.batch_size,
        learning_rate=args.learning_rate,
        weight_decay=args.weight_decay,
        patience=args.patience,
        seed=args.seed,
    )
    summary = run_harmonic_teacher_training(
        training_table=args.training_table,
        native_h5=args.native_h5,
        out_dir=args.out_dir,
        profiles=profiles,
        train_config=config,
        workers=args.workers,
        pretrain_epochs=args.pretrain_epochs,
        require_cuda=not args.allow_cpu,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=True))


if __name__ == "__main__":
    main()
