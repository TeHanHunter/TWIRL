#!/usr/bin/env python3
"""Verify a blinded existing-Teacher enrichment batch and rendered sheets."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from twirl.vetting.teacher_v2_active_learning import verify_enrichment_batch


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue", type=Path, required=True)
    parser.add_argument("--overlap", type=Path, required=True)
    parser.add_argument("--hidden", type=Path, required=True)
    parser.add_argument("--sheet-dir", type=Path, default=None)
    parser.add_argument("--out-json", type=Path, required=True)
    args = parser.parse_args()
    result = verify_enrichment_batch(
        pd.read_csv(args.queue),
        pd.read_csv(args.overlap),
        pd.read_parquet(args.hidden),
        sheet_dir=args.sheet_dir,
    )
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
