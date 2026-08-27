#!/usr/bin/env python3
"""Package passing FM mission-quality receipts and sources for ORCD transfer."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.models.fm0.mission_quality_transfer import (
    package_mission_quality_transfer,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", action="append", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--producer-git-sha", required=True)
    args = parser.parse_args()
    manifest = package_mission_quality_transfer(
        receipt_paths=args.receipt,
        output_dir=args.output_dir,
        producer_git_sha=args.producer_git_sha,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
