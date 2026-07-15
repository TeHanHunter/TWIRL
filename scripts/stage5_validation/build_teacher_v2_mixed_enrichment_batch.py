#!/usr/bin/env python3
"""Build one blinded, sector-balanced batch from frozen Teacher-v2 scores."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.vetting.teacher_v2_active_learning import (
    EnrichmentQuotas,
    write_mixed_existing_teacher_enrichment_batch,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sector-score-set",
        nargs=3,
        action="append",
        metavar=("SECTOR", "COMPACT_SCORES", "MORPHOLOGY_SCORES"),
        required=True,
    )
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--batch-index", type=int, default=0)
    parser.add_argument("--exclude", type=Path, nargs="*", default=[])
    parser.add_argument("--compact-count-per-sector", type=int, default=200)
    parser.add_argument("--eb-count-per-sector", type=int, default=150)
    parser.add_argument("--variable-count-per-sector", type=int, default=50)
    parser.add_argument("--disagreement-count-per-sector", type=int, default=50)
    parser.add_argument("--control-count-per-sector", type=int, default=50)
    parser.add_argument("--double-review-count", type=int, default=100)
    args = parser.parse_args()

    score_paths: dict[int, tuple[Path, Path]] = {}
    for sector_text, compact_text, morphology_text in args.sector_score_set:
        sector = int(sector_text)
        if sector in score_paths:
            raise ValueError(f"duplicate score set for Sector {sector}")
        score_paths[sector] = (Path(compact_text), Path(morphology_text))
    summary = write_mixed_existing_teacher_enrichment_batch(
        sector_score_paths=score_paths,
        out_dir=args.out_dir,
        batch_index=args.batch_index,
        exclude_paths=args.exclude,
        per_sector_quotas=EnrichmentQuotas(
            compact_transit=args.compact_count_per_sector,
            eclipse_contact=args.eb_count_per_sector,
            smooth_variable=args.variable_count_per_sector,
            model_disagreement=args.disagreement_count_per_sector,
            stratified_control=args.control_count_per_sector,
        ),
        double_review_count=args.double_review_count,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
