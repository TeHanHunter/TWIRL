#!/usr/bin/env python3
"""Build the hash-bound S56 SPOC quality table used by Tier-1 cadence QA."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.lightcurves.a2v1_cadence_reference import (  # noqa: E402
    S56_EXPECTED_DETECTORS,
    S56_EXPECTED_ORBITS,
    parse_quat_spec,
    parse_spoc_flag_spec,
    write_spoc_quality_table,
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--quat",
        action="append",
        required=True,
        metavar="ORBIT,CAMERA,PATH",
        help=(
            "Explicit QLP camC_quat.txt authority; repeat for all eight "
            "S56 orbit/camera combinations."
        ),
    )
    parser.add_argument(
        "--spoc-flag",
        action="append",
        required=True,
        metavar="CAMERA,CCD,PATH",
        help=(
            "Original detector spocffiflag_s56_camC_ccdD.txt; repeat for "
            "all 16 S56 detectors."
        ),
    )
    parser.add_argument("--output-table", type=Path, required=True)
    parser.add_argument("--output-provenance", type=Path, required=True)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Atomically replace an existing output table/provenance pair.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    try:
        quat_sources = tuple(parse_quat_spec(value) for value in args.quat)
        spoc_sources = tuple(
            parse_spoc_flag_spec(value) for value in args.spoc_flag
        )
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    manifest = write_spoc_quality_table(
        sector=56,
        quat_sources=quat_sources,
        spoc_flag_sources=spoc_sources,
        output_table=args.output_table,
        output_provenance=args.output_provenance,
        expected_orbits=S56_EXPECTED_ORBITS,
        expected_detectors=S56_EXPECTED_DETECTORS,
        overwrite=args.overwrite,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
