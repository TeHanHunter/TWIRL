#!/usr/bin/env python3
"""Build one hash-bound A2v1 cadence/quality reference."""
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
    parse_qlp_qflag_spec,
    parse_quat_spec,
    write_cadence_reference,
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sector",
        type=int,
        default=56,
        help="TESS sector (default: the locked S56 case).",
    )
    parser.add_argument(
        "--expected-orbit",
        action="append",
        type=int,
        dest="expected_orbits",
        help=(
            "Expected absolute TESS orbit; repeat once per sector orbit. "
            "S56 defaults to the locked 119/120 pair. Later sectors require "
            "this option explicitly."
        ),
    )
    parser.add_argument(
        "--quat",
        action="append",
        required=True,
        metavar="ORBIT,CAMERA,PATH",
        help=(
            "QLP camC_quat.txt with explicit orbit and camera; repeat for all "
            "expected orbit/camera combinations."
        ),
    )
    parser.add_argument(
        "--qlp-qflag",
        action="append",
        required=True,
        metavar="ORBIT,CAMERA,CCD,PATH",
        help=(
            "Headerless QLP cadence/qflag authority with explicit orbit and "
            "detector; repeat for every expected orbit/camera/CCD combination."
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


def _expected_orbits(*, sector: int, supplied: list[int] | None) -> tuple[int, ...]:
    if supplied:
        orbits = tuple(sorted(set(int(value) for value in supplied)))
        if len(orbits) != len(supplied):
            raise ValueError("--expected-orbit values must be unique")
    elif int(sector) == 56:
        orbits = S56_EXPECTED_ORBITS
    else:
        raise ValueError(
            "--expected-orbit is required for sectors other than the locked S56 case"
        )
    if not orbits or any(value <= 0 for value in orbits):
        raise ValueError("expected orbit IDs must be positive")
    if int(sector) == 56 and orbits != S56_EXPECTED_ORBITS:
        raise ValueError("S56 expected orbits are locked to 119 and 120")
    return orbits


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.sector < 56:
        raise SystemExit("A2v1 teacher preparation requires sector >= 56")
    try:
        expected_orbits = _expected_orbits(
            sector=args.sector,
            supplied=args.expected_orbits,
        )
        quat_sources = tuple(parse_quat_spec(value) for value in args.quat)
        qlp_qflag_sources = tuple(
            parse_qlp_qflag_spec(value) for value in args.qlp_qflag
        )
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    manifest = write_cadence_reference(
        sector=args.sector,
        quat_sources=quat_sources,
        qlp_qflag_sources=qlp_qflag_sources,
        spoc_quality_table=args.spoc_quality_table,
        spoc_quality_provenance=args.spoc_quality_provenance,
        output_table=args.output_table,
        output_manifest=args.output_manifest,
        expected_orbits=expected_orbits,
        expected_detectors=S56_EXPECTED_DETECTORS,
        overwrite=args.overwrite,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
