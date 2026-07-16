#!/usr/bin/env python3
"""Build and optionally score the real-label transfer gate onto S56 A2v1."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.teacher_active_learning import (  # noqa: E402
    build_a2v1_label_transfer_table,
    evaluate_a2v1_scored_transfer,
)


def read(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path) if path.suffix.lower() == ".parquet" else pd.read_csv(path, low_memory=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--adp-peaks", type=Path, required=True)
    parser.add_argument("--teacher-scores", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    transfer, compatibility = build_a2v1_label_transfer_table(
        read(args.labels), read(args.adp_peaks)
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    transfer.to_csv(args.out_dir / "a2v1_real_label_transfer_rows.csv", index=False)
    summary: dict = {"compatibility": compatibility}
    if args.teacher_scores is not None:
        evaluation, scored_summary = evaluate_a2v1_scored_transfer(
            transfer, read(args.teacher_scores)
        )
        evaluation.to_csv(args.out_dir / "a2v1_real_label_transfer_predictions.csv", index=False)
        summary["scored_transfer"] = scored_summary
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=True))


if __name__ == "__main__":
    main()
