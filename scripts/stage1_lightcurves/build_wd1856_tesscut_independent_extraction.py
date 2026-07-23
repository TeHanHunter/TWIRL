#!/usr/bin/env python3
"""Build the independent S56 WD 1856 reference from official TESSCut pixels."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.lightcurves.tesscut_independent_extraction import (  # noqa: E402
    DEFAULT_MAX_TIME_DELTA_SECONDS,
    DEFAULT_POSITION_TOLERANCE_ARCSEC,
    DEFAULT_ROLLING_WINDOW_CADENCES,
    build_wd1856_tesscut_independent_extraction,
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tesscut-fits",
        type=Path,
        required=True,
        help="Raw official odd-square MAST TESSCut FITS for WD 1856 in S56.",
    )
    parser.add_argument(
        "--compact-lc",
        type=Path,
        required=True,
        help=(
            "S56 A2v1 compact HDF5 used only to map timestamps to WD 1856 "
            "CADENCENO values; its flux datasets are never read."
        ),
    )
    parser.add_argument(
        "--cadence-reference-table",
        type=Path,
        required=True,
        help="Authoritative S56 external-quality cadence table.",
    )
    parser.add_argument(
        "--cadence-reference-manifest",
        type=Path,
        required=True,
        help="Hash-bound manifest for --cadence-reference-table.",
    )
    parser.add_argument("--output-fits", type=Path, required=True)
    parser.add_argument("--manifest-json", type=Path, required=True)
    parser.add_argument(
        "--source-url",
        required=True,
        help="Exact HTTPS STScI/MAST URL from which the TESSCut FITS was obtained.",
    )
    parser.add_argument(
        "--rolling-window-cadences",
        type=int,
        default=DEFAULT_ROLLING_WINDOW_CADENCES,
        help="Odd centered quality-zero rolling-median window (default: 433).",
    )
    parser.add_argument(
        "--position-tolerance-arcsec",
        type=float,
        default=DEFAULT_POSITION_TOLERANCE_ARCSEC,
        help="RA/Dec identity tolerance; values above 1 arcsec are rejected.",
    )
    parser.add_argument(
        "--max-time-delta-seconds",
        type=float,
        default=DEFAULT_MAX_TIME_DELTA_SECONDS,
        help="Nearest timestamp tolerance; values above 60 seconds are rejected.",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    manifest = build_wd1856_tesscut_independent_extraction(
        tesscut_fits=args.tesscut_fits,
        compact_lc=args.compact_lc,
        cadence_reference_table=args.cadence_reference_table,
        cadence_reference_manifest=args.cadence_reference_manifest,
        output_fits=args.output_fits,
        manifest_json=args.manifest_json,
        source_url=args.source_url,
        rolling_window_cadences=args.rolling_window_cadences,
        position_tolerance_arcsec=args.position_tolerance_arcsec,
        max_time_delta_seconds=args.max_time_delta_seconds,
        overwrite=args.overwrite,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
