#!/usr/bin/env python3
"""Build one deferred, provider-neutral TWIRL-FM later-sector six-view bundle."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.models.fm0.later_sector_release import build_later_sector_release
from twirl.models.fm0.registry import FM0ContractError


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector", required=True, type=int)
    parser.add_argument("--source-inventory-dir", required=True, type=Path)
    parser.add_argument("--source-inventory-summary-sha256", required=True)
    parser.add_argument("--mission-quality-reference-dir", required=True, type=Path)
    parser.add_argument("--mission-quality-reference-manifest-sha256", required=True)
    parser.add_argument("--hdf5-quality-receipt", required=True, type=Path)
    parser.add_argument("--hdf5-quality-receipt-sha256", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--producer-git-sha", required=True)
    args = parser.parse_args()
    receipt = build_later_sector_release(
        sector=args.sector,
        source_inventory_dir=args.source_inventory_dir,
        expected_source_inventory_summary_sha256=(
            args.source_inventory_summary_sha256
        ),
        mission_quality_reference_dir=args.mission_quality_reference_dir,
        expected_mission_quality_reference_manifest_sha256=(
            args.mission_quality_reference_manifest_sha256
        ),
        hdf5_quality_receipt_path=args.hdf5_quality_receipt,
        expected_hdf5_quality_receipt_sha256=args.hdf5_quality_receipt_sha256,
        output_dir=args.output_dir,
        producer_git_sha=args.producer_git_sha,
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FM0ContractError, OSError, ValueError) as exc:
        print(f"[fm0-later-six-view] ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
