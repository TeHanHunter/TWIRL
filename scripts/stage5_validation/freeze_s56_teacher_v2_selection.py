#!/usr/bin/env python3
"""Freeze Teacher-v2 architecture and compact threshold without opening holdouts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.vetting.teacher_v2_selection import freeze_teacher_v2_selection
from twirl.vetting.teacher_v2_training import DEFAULT_TEACHER_V2_PROFILES


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training-root", type=Path, required=True)
    parser.add_argument("--real-score-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--profiles", default=",".join(DEFAULT_TEACHER_V2_PROFILES))
    parser.add_argument("--max-real-tic-fraction", type=float, default=0.05)
    args = parser.parse_args()

    profiles = tuple(value.strip() for value in args.profiles.split(",") if value.strip())
    prediction_paths = {
        profile: [
            args.training_root / profile / f"fold_{fold}" / "validation_predictions.parquet"
            for fold in range(5)
        ]
        for profile in profiles
    }
    checkpoint_paths = {
        profile: [
            args.training_root / profile / f"fold_{fold}" / "teacher_v2.pt"
            for fold in range(5)
        ]
        for profile in profiles
    }
    real_score_paths = {
        profile: args.real_score_dir / f"{profile}.parquet" for profile in profiles
    }
    ranking, frozen = freeze_teacher_v2_selection(
        profile_prediction_paths=prediction_paths,
        profile_real_score_paths=real_score_paths,
        profile_checkpoint_paths=checkpoint_paths,
        max_real_tic_fraction=args.max_real_tic_fraction,
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    ranking.to_csv(args.out_dir / "development_profile_ranking.csv", index=False)
    (args.out_dir / "frozen_selection.json").write_text(
        json.dumps(frozen, indent=2, sort_keys=True, allow_nan=True) + "\n"
    )
    print(json.dumps(frozen, indent=2, sort_keys=True, allow_nan=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
