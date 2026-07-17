#!/usr/bin/env python3
"""Build a blinded, rank-1 Teacher-v1 enrichment queue for Franklin."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.vetting.franklin_multisector import write_franklin_multisector_batch
from twirl.vetting.teacher_v2_active_learning import EnrichmentQuotas


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sector-score",
        nargs=2,
        action="append",
        metavar=("SECTOR", "TEACHER_SCORES"),
        required=True,
    )
    parser.add_argument(
        "--sector-quota",
        nargs=6,
        action="append",
        metavar=("SECTOR", "COMPACT", "EB", "VARIABLE", "DISAGREEMENT", "CONTROL"),
        required=True,
    )
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--exclude", type=Path, nargs="*", default=[])
    parser.add_argument("--batch-index", type=int, default=0)
    parser.add_argument("--double-review-count", type=int, default=0)
    parser.add_argument("--seed", type=int, default=560717)
    args = parser.parse_args()

    score_paths: dict[int, Path] = {}
    for sector_text, path_text in args.sector_score:
        sector = int(sector_text)
        if sector in score_paths:
            raise ValueError(f"duplicate score path for Sector {sector}")
        score_paths[sector] = Path(path_text)

    quotas: dict[int, EnrichmentQuotas] = {}
    for values in args.sector_quota:
        sector = int(values[0])
        if sector in quotas:
            raise ValueError(f"duplicate quota for Sector {sector}")
        counts = [int(value) for value in values[1:]]
        quotas[sector] = EnrichmentQuotas(
            compact_transit=counts[0],
            eclipse_contact=counts[1],
            smooth_variable=counts[2],
            model_disagreement=counts[3],
            stratified_control=counts[4],
        )
    if set(score_paths) != set(quotas):
        raise ValueError(
            f"score sectors {sorted(score_paths)} differ from quota sectors {sorted(quotas)}"
        )

    summary = write_franklin_multisector_batch(
        sector_score_paths=score_paths,
        sector_quotas=quotas,
        out_dir=args.out_dir,
        exclude_paths=args.exclude,
        batch_index=args.batch_index,
        double_review_count=args.double_review_count,
        seed=args.seed,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
