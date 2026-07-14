#!/usr/bin/env python3
"""Freeze a candidate-score threshold at a maximum real-TIC review load."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import json
from pathlib import Path

import pandas as pd

from twirl.vetting.teacher_v2 import freeze_real_tic_workload_threshold


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--score-column", default="p_preserve")
    parser.add_argument("--max-fraction", type=float, default=0.05)
    parser.add_argument("--model-name", required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    args = parser.parse_args()

    scores = (
        pd.read_parquet(args.scores)
        if args.scores.suffix.lower() == ".parquet"
        else pd.read_csv(args.scores, low_memory=False)
    )
    threshold = freeze_real_tic_workload_threshold(
        scores,
        score_column=args.score_column,
        max_fraction=args.max_fraction,
    )
    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "model_name": args.model_name,
        "workload_threshold": asdict(threshold),
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
