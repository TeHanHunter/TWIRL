#!/usr/bin/env python3
"""Build S56 EB/PCEB miner training and candidate-scan tables."""
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
from twirl.vetting.recovery50_teacher import json_default  # noqa: E402


DEFAULT_ROOT = REPO_ROOT / "reports/stage5_validation/s56_eb_miner_adp_only"


DEFAULT_JOINED_TABLES = (
    REPO_ROOT
    / "reports/stage5_validation/s56_recovery50_teacher_queue_2k/"
    / "human_training_table_adp_only/human_vetting_training_table.csv",
)
DEFAULT_QUEUE_LABEL_PAIRS: tuple[tuple[Path, Path, str], ...] = ()
DEFAULT_PRIMARY_CANDIDATES = (
    REPO_ROOT
    / "reports/stage5_validation/s56_adp_real_bls_peaks/real_adp_bls_peaks.parquet"
)
DEFAULT_FALLBACK_CANDIDATES = None
DEFAULT_EXCLUDE_TABLES = (
    REPO_ROOT / "reports/stage5_validation/s56_eb_miner_adp_only/training/eb_miner_labeled_audit.csv",
    REPO_ROOT / "reports/stage5_validation/s56_franklin_real5k_handoff/franklin_review_queue_5k_real.csv",
    REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue_next4k/review_queue_4k.csv",
    REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue/review_queue_1k.csv",
    REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue_next1k/review_queue_1k.csv",
    REPO_ROOT / "reports/stage5_validation/s56_mixed_teacher_queue_pdo/review_queue_1k.csv",
    REPO_ROOT
    / "reports/stage5_validation/s56_eb_miner/review_queue_eb_priority_500/"
    / "human_labels_vetted.csv",
)


def _parse_pair(raw: str) -> tuple[Path, Path, str]:
    parts = raw.split(":", 2)
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("queue-label pair must be queue.csv:labels.csv:name")
    return Path(parts[0]), Path(parts[1]), parts[2]


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-root", type=Path, default=DEFAULT_ROOT)
    ap.add_argument("--joined-table", type=Path, action="append", default=list(DEFAULT_JOINED_TABLES))
    ap.add_argument(
        "--queue-label-pair",
        type=_parse_pair,
        action="append",
        default=list(DEFAULT_QUEUE_LABEL_PAIRS),
        help="Label source as queue.csv:labels.csv:name; may be repeated.",
    )
    ap.add_argument("--primary-candidates", type=Path, default=DEFAULT_PRIMARY_CANDIDATES)
    ap.add_argument(
        "--fallback-candidates",
        type=Path,
        default=DEFAULT_FALLBACK_CANDIDATES,
        help="Optional second ADP-only BLS peak table. Canonical search tables are rejected.",
    )
    ap.add_argument("--small-peaks-per-tic", type=int, default=3)
    ap.add_argument("--primary-peaks-per-tic", type=int, default=1)
    ap.add_argument("--exclude-table", type=Path, action="append", default=list(DEFAULT_EXCLUDE_TABLES))
    ap.add_argument("--max-uncertain-negatives", type=int, default=100)
    ap.add_argument("--max-negatives-per-label", type=int, default=200)
    ap.add_argument("--min-positive-rows", type=int, default=5)
    ap.add_argument("--random-state", type=int, default=56017)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    cfg = EBTrainingConfig(
        random_state=args.random_state,
        max_uncertain_negatives=args.max_uncertain_negatives,
        max_negatives_per_label=args.max_negatives_per_label,
        min_positive_rows=args.min_positive_rows,
    )
    training_summary = build_eb_miner_training_table(
        joined_tables=tuple(args.joined_table or ()),
        queue_label_pairs=tuple(args.queue_label_pair or ()),
        out_dir=args.out_root / "training",
        cfg=cfg,
    )
    candidate_summary = build_candidate_scoring_pool(
        primary_candidates=args.primary_candidates,
        fallback_candidates=args.fallback_candidates,
        exclude_tables=tuple(args.exclude_table or ()),
        out_dir=args.out_root / "candidate_pool",
        small_peaks_per_tic=args.small_peaks_per_tic,
        primary_peaks_per_tic=args.primary_peaks_per_tic,
    )
    summary = {
        "training": training_summary,
        "candidate_pool": candidate_summary,
        "outputs": {
            "training_table": training_summary["outputs"]["training_table"],
            "candidate_pool": candidate_summary["outputs"]["candidate_pool"],
        },
    }
    (args.out_root / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=json_default) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True, default=json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
