#!/usr/bin/env python3
"""Merge immutable FM0.1 sector input sub-releases into one full release."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.models.fm0.registry import load_frozen_contract
from twirl.models.fm0.sector_release import merge_sector_input_releases


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector-root", required=True, type=Path)
    parser.add_argument("--registry-dir", required=True, type=Path)
    parser.add_argument("--from-sector", type=int, default=56)
    parser.add_argument("--through-sector", type=int, default=64)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--producer-git-sha", required=True)
    parser.add_argument(
        "--allowed-sector-producer-git-sha",
        action="append",
        default=[],
        help="Repeat for each reviewed sector-input commit; defaults to the merge commit.",
    )
    parser.add_argument("--design", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--freeze-receipt", type=Path)
    args = parser.parse_args()
    if args.through_sector < args.from_sector:
        raise SystemExit("invalid sector range")
    contract = load_frozen_contract(
        design_path=args.design,
        config_path=args.config,
        freeze_receipt_path=args.freeze_receipt,
    )
    summary = merge_sector_input_releases(
        sector_root=args.sector_root,
        sectors=range(args.from_sector, args.through_sector + 1),
        registry_dir=args.registry_dir,
        out_dir=args.out_dir,
        producer_git_sha=args.producer_git_sha,
        allowed_sector_producer_git_shas=(
            args.allowed_sector_producer_git_sha or [args.producer_git_sha]
        ),
        contract=contract,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
