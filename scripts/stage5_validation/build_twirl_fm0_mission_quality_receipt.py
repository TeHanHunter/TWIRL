#!/usr/bin/env python3
"""Build one fail-closed FM mission-quality source receipt."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.models.fm0.mission_quality_admission import (
    verify_mission_quality_sources,
    write_mission_quality_receipt,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector", required=True, type=int)
    parser.add_argument("--expected-orbit", required=True, type=int, action="append")
    parser.add_argument("--qlp-root", required=True, type=Path)
    parser.add_argument(
        "--mission-flag-root",
        type=Path,
        help="Override the default spocflags (<S67) or ticaflags (S67+) root.",
    )
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    receipt = verify_mission_quality_sources(
        sector=args.sector,
        expected_orbits=args.expected_orbit,
        qlp_root=args.qlp_root,
        mission_flag_root=args.mission_flag_root,
    )
    path = write_mission_quality_receipt(receipt, args.output)
    print(json.dumps({"output": str(path), **receipt}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
