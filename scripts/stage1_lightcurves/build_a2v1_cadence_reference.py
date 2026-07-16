#!/usr/bin/env python3
"""Build the hash-bound S56 A2v1 cadence/quality reference for Tier-1 QA."""
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
    write_cadence_reference,
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sector",
        type=int,
        default=56,
        help="TESS sector (the locked command-line contract currently supports S56).",
    )
    parser.add_argument(
        "--quat",
        action="append",
        required=True,
        metavar="ORBIT,CAMERA,PATH",
        help=(
            "QLP camC_quat.txt with explicit orbit and camera; repeat for all "
            "eight S56 orbit/camera combinations."
        ),
    )
    parser.add_argument(
        "--spoc-quality-table",
        type=Path,
        required=True,
        help=(
            "Precomputed CSV/Parquet with sector,camera,ccd,cadenceno,quality "
            "and optional orbitid."
        ),
    )
    parser.add_argument(
        "--spoc-quality-provenance",
        type=Path,
        required=True,
        help=(
            "Mandatory JSON sidecar binding the table to the derivation contract "
            "and all 16 original detector-level SPOC flag files."
        ),
    )
    parser.add_argument("--output-table", type=Path, required=True)
    parser.add_argument("--output-manifest", type=Path, required=True)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Atomically replace existing table/manifest outputs.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.sector != 56:
        raise SystemExit(
            "this locked CLI builds the S56 SPOC-quality reference only; "
            "use the module API for a predeclared later-sector contract"
        )
    try:
        quat_sources = tuple(parse_quat_spec(value) for value in args.quat)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    manifest = write_cadence_reference(
        sector=args.sector,
        quat_sources=quat_sources,
        spoc_quality_table=args.spoc_quality_table,
        spoc_quality_provenance=args.spoc_quality_provenance,
        output_table=args.output_table,
        output_manifest=args.output_manifest,
        expected_orbits=S56_EXPECTED_ORBITS,
        expected_detectors=S56_EXPECTED_DETECTORS,
        overwrite=args.overwrite,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
