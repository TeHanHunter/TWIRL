#!/usr/bin/env python3
"""Freeze the identity-only S66--S77 TWIRL-FM development panel."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.models.fm0.temporal_panel import freeze_temporal_panel


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--config-sha256", required=True)
    parser.add_argument("--admission-receipt", required=True, type=Path)
    parser.add_argument("--admission-receipt-sha256", required=True)
    parser.add_argument("--baseline-manifest", required=True, type=Path)
    parser.add_argument("--baseline-manifest-sha256", required=True)
    parser.add_argument("--producer-git-sha", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    frozen = freeze_temporal_panel(
        config_path=args.config,
        config_sha256=args.config_sha256,
        admission_receipt_path=args.admission_receipt,
        admission_receipt_sha256=args.admission_receipt_sha256,
        baseline_manifest_path=args.baseline_manifest,
        baseline_manifest_sha256=args.baseline_manifest_sha256,
        producer_git_sha=args.producer_git_sha,
        output_dir=args.output_dir,
    )
    print(
        json.dumps(
            {
                "output_dir": str(frozen.output_dir),
                "receipt_path": str(frozen.receipt_path),
                "receipt_sha256": frozen.receipt_sha256,
                "n_panel_rows": frozen.receipt["n_panel_rows"],
                "n_repeated_rows": frozen.receipt["n_repeated_rows"],
                "n_new_rows": frozen.receipt["n_new_rows"],
                "scientific_training_eligible": False,
                "development_evaluation_eligible": True,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
