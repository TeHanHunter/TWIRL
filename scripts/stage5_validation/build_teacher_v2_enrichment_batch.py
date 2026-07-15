#!/usr/bin/env python3
"""Build one blinded enrichment batch from existing Teacher-v2 models."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.vetting.teacher_v2_active_learning import (
    EnrichmentQuotas,
    write_existing_teacher_enrichment_batch,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compact-scores", type=Path, required=True)
    parser.add_argument("--morphology-scores", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--sector", type=int, required=True)
    parser.add_argument("--batch-index", type=int, default=0)
    parser.add_argument("--exclude", type=Path, nargs="*", default=[])
    parser.add_argument("--compact-count", type=int, default=400)
    parser.add_argument("--eb-count", type=int, default=300)
    parser.add_argument("--variable-count", type=int, default=100)
    parser.add_argument("--disagreement-count", type=int, default=100)
    parser.add_argument("--control-count", type=int, default=100)
    args = parser.parse_args()
    summary = write_existing_teacher_enrichment_batch(
        compact_scores_path=args.compact_scores,
        morphology_scores_path=args.morphology_scores,
        out_dir=args.out_dir,
        sector=args.sector,
        batch_index=args.batch_index,
        exclude_paths=args.exclude,
        quotas=EnrichmentQuotas(
            compact_transit=args.compact_count,
            eclipse_contact=args.eb_count,
            smooth_variable=args.variable_count,
            model_disagreement=args.disagreement_count,
            stratified_control=args.control_count,
        ),
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
