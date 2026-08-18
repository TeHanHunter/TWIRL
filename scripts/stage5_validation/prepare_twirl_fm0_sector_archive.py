#!/usr/bin/env python3
"""Expand one transferred A2v1 sector tar into an FM0.1 source inventory."""
from __future__ import annotations

import argparse
import json

from twirl.models.fm0.sector_stage import stage_sector_archive


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector", required=True, type=int)
    parser.add_argument("--selection", required=True)
    parser.add_argument("--selection-sha256", required=True)
    parser.add_argument("--archive-dir", required=True)
    parser.add_argument("--quality-table", required=True)
    parser.add_argument("--quality-manifest", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--producer-git-sha", required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    summary = stage_sector_archive(
        sector=args.sector,
        selection_path=args.selection,
        selection_sha256=args.selection_sha256,
        archive_dir=args.archive_dir,
        quality_table=args.quality_table,
        quality_manifest=args.quality_manifest,
        output_root=args.output_root,
        producer_git_sha=args.producer_git_sha,
        workers=args.workers,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
