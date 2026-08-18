#!/usr/bin/env python3
"""Merge completed FM0.1 sector source inventories."""
from __future__ import annotations

import argparse
import json

from twirl.models.fm0.sector_stage import merge_sector_inventories


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector-root", required=True)
    parser.add_argument("--from-sector", type=int, default=56)
    parser.add_argument("--through-sector", type=int, default=64)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--producer-git-sha", required=True)
    parser.add_argument(
        "--allowed-sector-producer-git-sha",
        action="append",
        default=[],
        help="Repeat for each reviewed sector-stage commit; defaults to the merge commit.",
    )
    args = parser.parse_args()
    if args.through_sector < args.from_sector:
        raise SystemExit("invalid sector range")
    summary = merge_sector_inventories(
        sector_root=args.sector_root,
        sectors=range(args.from_sector, args.through_sector + 1),
        out_dir=args.out_dir,
        producer_git_sha=args.producer_git_sha,
        allowed_sector_producer_git_shas=(
            args.allowed_sector_producer_git_sha or [args.producer_git_sha]
        ),
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
