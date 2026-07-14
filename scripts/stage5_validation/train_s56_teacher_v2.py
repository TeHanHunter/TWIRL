#!/usr/bin/env python3
"""Train one or more S56 Teacher-v2 architecture profiles on ORCD."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from twirl.vetting.teacher_v2_cnn import TeacherV2TrainConfig
from twirl.vetting.teacher_v2_training import (
    DEFAULT_TEACHER_V2_PROFILES,
    run_teacher_v2_training,
)


def _balanced_smoke_rows(rows: pd.DataFrame, *, n_rows: int, seed: int) -> pd.DataFrame:
    development = rows.loc[rows["fixed_split"].eq("development")].copy()
    if n_rows >= len(development):
        return development.reset_index(drop=True)
    pieces: list[pd.DataFrame] = []
    strata = (
        development["cv_fold"].astype(str)
        + "|"
        + development["teacher_v2_role"].fillna("").astype(str)
    )
    quota = max(1, int(n_rows) // max(strata.nunique(), 1))
    for _, group in development.groupby(strata, sort=True):
        pieces.append(group.sample(n=min(quota, len(group)), random_state=seed))
    selected = pd.concat(pieces, ignore_index=False)
    if len(selected) < n_rows:
        remaining = development.drop(index=selected.index)
        selected = pd.concat(
            [
                selected,
                remaining.sample(
                    n=min(n_rows - len(selected), len(remaining)), random_state=seed + 1
                ),
            ]
        )
    # A smoke never opens or scores the holdout. Keep zero holdout rows in its
    # serialized input so accidental holdout access cannot occur.
    return selected.sample(frac=1.0, random_state=seed).reset_index(drop=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training-table", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--profiles", default=",".join(DEFAULT_TEACHER_V2_PROFILES))
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--patience", type=int, default=12)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=56)
    parser.add_argument("--max-development-rows", type=int)
    parser.add_argument("--allow-cpu", action="store_true")
    args = parser.parse_args()

    training_table = args.training_table
    if args.max_development_rows is not None:
        rows = (
            pd.read_parquet(training_table)
            if training_table.suffix.lower() == ".parquet"
            else pd.read_csv(training_table, low_memory=False)
        )
        rows = _balanced_smoke_rows(
            rows, n_rows=int(args.max_development_rows), seed=int(args.seed)
        )
        args.out_dir.mkdir(parents=True, exist_ok=True)
        training_table = args.out_dir / "smoke_training_rows.parquet"
        rows.to_parquet(training_table, compression="zstd", index=False)
    profiles = tuple(value.strip() for value in args.profiles.split(",") if value.strip())
    summary = run_teacher_v2_training(
        training_table=training_table,
        out_dir=args.out_dir,
        profiles=profiles,
        train_config=TeacherV2TrainConfig(
            epochs=int(args.epochs),
            patience=int(args.patience),
            batch_size=int(args.batch_size),
            seed=int(args.seed),
        ),
        workers=int(args.workers),
        require_cuda=not args.allow_cpu,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
