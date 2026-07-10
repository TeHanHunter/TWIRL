#!/usr/bin/env python3
"""Build EB-miner tables from the final adjudicated training table only."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.vetting.eb_miner import (  # noqa: E402
    EBTrainingConfig,
    build_candidate_scoring_pool,
    build_eb_miner_training_table,
)


DEFAULT_ADJUDICATION = REPO_ROOT / "reports/stage5_validation/s56_label_adjudication_real343"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--training-table",
        type=Path,
        default=DEFAULT_ADJUDICATION / "adjudicated_training_table/human_vetting_training_table_adjudicated.csv",
    )
    parser.add_argument(
        "--primary-candidates",
        type=Path,
        default=REPO_ROOT / "reports/stage5_validation/s56_adp_real_bls_peaks/real_adp_bls_peaks.parquet",
    )
    parser.add_argument("--out-root", type=Path, default=DEFAULT_ADJUDICATION / "retrained/eb_miner")
    parser.add_argument("--max-uncertain-negatives", type=int, default=200)
    parser.add_argument("--max-negatives-per-label", type=int, default=300)
    args = parser.parse_args(argv)
    config = EBTrainingConfig(
        random_state=56017,
        max_uncertain_negatives=args.max_uncertain_negatives,
        max_negatives_per_label=args.max_negatives_per_label,
        min_positive_rows=5,
    )
    training = build_eb_miner_training_table(
        joined_tables=(args.training_table,),
        queue_label_pairs=(),
        out_dir=args.out_root / "training",
        cfg=config,
    )
    candidates = build_candidate_scoring_pool(
        primary_candidates=args.primary_candidates,
        fallback_candidates=None,
        exclude_tables=(args.training_table,),
        out_dir=args.out_root / "candidate_pool",
        small_peaks_per_tic=3,
        primary_peaks_per_tic=1,
    )
    summary = {"training": training, "candidate_pool": candidates}
    args.out_root.mkdir(parents=True, exist_ok=True)
    (args.out_root / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=str) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
