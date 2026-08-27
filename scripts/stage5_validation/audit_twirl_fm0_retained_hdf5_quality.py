#!/usr/bin/env python3
"""Audit every retained ORCD HDF5 against transferred cadence authorities."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.models.fm0.hdf5_quality_admission import (
    audit_retained_sector_hdf5_quality,
    write_hdf5_quality_receipt,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector", required=True, type=int)
    parser.add_argument("--expected-orbit", action="append", required=True, type=int)
    parser.add_argument("--source-receipt", required=True, type=Path)
    parser.add_argument("--quality-transfer-root", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--producer-git-sha", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    receipt = audit_retained_sector_hdf5_quality(
        source_receipt_path=args.source_receipt,
        quality_transfer_root=args.quality_transfer_root,
        sector=args.sector,
        expected_orbits=args.expected_orbit,
        workers=args.workers,
        producer_git_sha=args.producer_git_sha,
    )
    path = write_hdf5_quality_receipt(receipt, args.output)
    print(json.dumps({"output": str(path), **receipt}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
