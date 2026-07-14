#!/usr/bin/env python3
"""Transfer clean real human labels onto current S56 A2v1 ADP ephemerides."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import pandas as pd

from twirl.vetting.teacher_v2 import (
    normalize_real_adp_candidates,
    transfer_human_labels_to_a2v1_candidates,
)


def _read(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path) if path.suffix.lower() == ".parquet" else pd.read_csv(path, low_memory=False)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--human-table", type=Path, required=True)
    parser.add_argument("--real-bls-peaks", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--small-peaks", type=int, default=10)
    args = parser.parse_args()

    human = _read(args.human_table)
    candidates = normalize_real_adp_candidates(
        _read(args.real_bls_peaks), small_peaks_per_tic=args.small_peaks
    )
    transferred, compatibility = transfer_human_labels_to_a2v1_candidates(
        human, candidates
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    candidates.to_parquet(
        args.out_dir / "s56_a2v1_real_candidates_top10.parquet",
        compression="zstd",
        index=False,
    )
    top5 = candidates.loc[
        pd.to_numeric(candidates["rep_peak_rank"], errors="coerce").le(5)
    ].copy()
    top5.to_parquet(
        args.out_dir / "s56_a2v1_real_candidates_top5.parquet",
        compression="zstd",
        index=False,
    )
    compatibility.to_csv(args.out_dir / "human_a2v1_compatibility.csv", index=False)
    transferred.to_csv(args.out_dir / "human_a2v1_transferred_training_rows.csv", index=False)
    label_counts = (
        transferred.get("human_label", pd.Series(dtype=str))
        .fillna("")
        .astype(str)
        .value_counts()
    )
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "n_human_input": int(len(human)),
        "n_real_human_input": int(len(compatibility)),
        "n_candidates": int(len(candidates)),
        "n_top5_candidates": int(len(top5)),
        "n_candidate_tics": int(candidates["tic"].nunique()),
        "n_transferred": int(len(transferred)),
        "transfer_fraction": float(compatibility["a2v1_transfer_ok"].mean()) if len(compatibility) else 0.0,
        "transferred_label_counts": {str(key): int(value) for key, value in label_counts.items()},
        "period_tolerance": 0.02,
        "minimum_window_overlap": 0.5,
        "outputs": {
            "candidates": str(args.out_dir / "s56_a2v1_real_candidates_top10.parquet"),
            "candidates_top5": str(args.out_dir / "s56_a2v1_real_candidates_top5.parquet"),
            "compatibility": str(args.out_dir / "human_a2v1_compatibility.csv"),
            "transferred": str(args.out_dir / "human_a2v1_transferred_training_rows.csv"),
        },
    }
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
