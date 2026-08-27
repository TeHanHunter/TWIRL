#!/usr/bin/env python3
"""Build a fail-closed ORCD FM light-curve source-readiness receipt."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.models.fm0.orcd_source_admission import (
    verify_archive_source,
    verify_retained_sector_source,
    write_source_receipt,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector", required=True, type=int)
    parser.add_argument("--expected-orbit", required=True, type=int, action="append")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--archive-dir", type=Path)
    source.add_argument("--retained-sector-root", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    kwargs = {"sector": args.sector, "expected_orbits": args.expected_orbit}
    if args.archive_dir is not None:
        receipt = verify_archive_source(archive_dir=args.archive_dir, **kwargs)
    else:
        receipt = verify_retained_sector_source(
            sector_root=args.retained_sector_root, **kwargs
        )
    path = write_source_receipt(receipt, args.output)
    print(json.dumps({"output": str(path), **receipt}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
