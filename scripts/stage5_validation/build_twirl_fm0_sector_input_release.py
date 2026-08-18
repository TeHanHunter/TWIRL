#!/usr/bin/env python3
"""Build one checksum-bound, BLS-free FM0.1 sector input sub-release."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.models.fm0.registry import load_frozen_contract
from twirl.models.fm0.sector_release import build_sector_input_release


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector", required=True, type=int)
    parser.add_argument("--registry-dir", required=True, type=Path)
    parser.add_argument("--a2v1-hdf5-manifest", required=True, type=Path)
    parser.add_argument("--out-root", required=True, type=Path)
    parser.add_argument("--producer-git-sha", required=True)
    parser.add_argument("--design", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--freeze-receipt", type=Path)
    args = parser.parse_args()
    contract = load_frozen_contract(
        design_path=args.design,
        config_path=args.config,
        freeze_receipt_path=args.freeze_receipt,
    )
    summary = build_sector_input_release(
        sector=args.sector,
        registry_dir=args.registry_dir,
        hdf5_manifest_path=args.a2v1_hdf5_manifest,
        output_root=args.out_root,
        producer_git_sha=args.producer_git_sha,
        contract=contract,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
