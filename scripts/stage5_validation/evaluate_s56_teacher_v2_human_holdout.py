#!/usr/bin/env python3
"""Evaluate Teacher v2 and v1 on the same locked S56 human rows."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from twirl.vetting.teacher_v2_evaluation import evaluate_locked_human_holdout


def _read(path: Path) -> pd.DataFrame:
    return (
        pd.read_parquet(path)
        if path.suffix.lower() == ".parquet"
        else pd.read_csv(path, low_memory=False)
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--human-role-table", type=Path, required=True)
    parser.add_argument("--teacher-v2-scores", type=Path, required=True)
    parser.add_argument("--teacher-v1-scores", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    rows, summary = evaluate_locked_human_holdout(
        _read(args.human_role_table),
        _read(args.teacher_v2_scores),
        _read(args.teacher_v1_scores),
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows.to_parquet(
        args.out_dir / "same_row_human_holdout_predictions.parquet",
        compression="zstd",
        index=False,
    )
    (args.out_dir / "same_row_human_holdout_metrics.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=True))
    return 0 if summary["acceptance"]["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
