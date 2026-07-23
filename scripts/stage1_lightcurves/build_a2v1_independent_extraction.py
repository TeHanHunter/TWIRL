#!/usr/bin/env python3
"""Build hash-bound S56 WD 1856 independent-extraction Tier-1 evidence."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.lightcurves.a2v1_independent_extraction import (  # noqa: E402
    IndependentExtractionProvenance,
    build_wd1856_independent_metrics,
    parse_reference_flux_column_mappings,
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--compact-lc",
        type=Path,
        required=True,
        help="Exact S56 A2v1 ADP-pair compact HDF5 product.",
    )
    parser.add_argument(
        "--reference-table",
        type=Path,
        required=True,
        help="External WD 1856 cadence table (CSV, Parquet, ECSV, or FITS).",
    )
    parser.add_argument(
        "--reference-product",
        type=Path,
        required=True,
        help=(
            "Original external extraction product to hash and identity-check; "
            "may equal --reference-table."
        ),
    )
    parser.add_argument(
        "--cadence-reference-table",
        type=Path,
        required=True,
        help="Authoritative S56 detector-cadence external-quality CSV/Parquet.",
    )
    parser.add_argument(
        "--cadence-reference-manifest",
        type=Path,
        required=True,
        help="Hash-bound provenance JSON for --cadence-reference-table.",
    )
    parser.add_argument(
        "--reference-time-system",
        choices=("BJD", "BTJD"),
        required=True,
        help="Time system of the external table; BTJD means BJD-2457000.",
    )
    parser.add_argument(
        "--reference-flux-column",
        action="append",
        required=True,
        metavar="CURRENT=REFERENCE",
        help=(
            "Map one active A2v1 aperture to its external reference flux column; "
            "repeat exactly once for each active aperture."
        ),
    )
    parser.add_argument("--reference-cadence-column", default="CADENCENO")
    parser.add_argument("--reference-time-column", default="TIME")
    parser.add_argument("--reference-quality-column", default="QUALITY")
    parser.add_argument(
        "--reference-tic-column",
        default="TICID",
        help="Required identity column for non-FITS reference inputs.",
    )
    parser.add_argument(
        "--reference-sector-column",
        default="SECTOR",
        help="Required identity column for non-FITS reference inputs.",
    )
    parser.add_argument("--current-repository", required=True)
    parser.add_argument("--current-revision", required=True)
    parser.add_argument("--reference-extractor-family", required=True)
    parser.add_argument("--reference-repository", required=True)
    parser.add_argument("--reference-revision", required=True)
    parser.add_argument(
        "--reference-pixel-source",
        required=True,
        help="Explicit detector/pixel/calibration source used by the external extraction.",
    )
    parser.add_argument("--reference-product-aperture", required=True)
    parser.add_argument(
        "--reference-target-id",
        required=True,
        help="Explicit target declaration; must identify TIC 267574918.",
    )
    parser.add_argument(
        "--independence-basis",
        required=True,
        help="Auditable explanation of why the external extraction is independent.",
    )
    parser.add_argument("--metrics-csv", type=Path, required=True)
    parser.add_argument("--manifest-json", type=Path, required=True)
    parser.add_argument("--transit-duration-min", type=float, default=8.0)
    parser.add_argument("--oot-exclusion-min", type=float, default=30.0)
    parser.add_argument("--max-time-delta-seconds", type=float, default=5.0)
    parser.add_argument("--min-common-cadences", type=int, default=100)
    parser.add_argument("--min-in-event-cadences", type=int, default=4)
    parser.add_argument("--min-out-of-event-cadences", type=int, default=100)
    parser.add_argument("--bls-period-relative-half-width", type=float, default=0.01)
    parser.add_argument("--bls-n-periods", type=int, default=4_001)
    parser.add_argument(
        "--bls-duration-min",
        type=float,
        action="append",
        dest="bls_durations_min",
        help="Bounded-BLS duration in minutes; repeat (default: 6, 8, 10).",
    )
    parser.add_argument("--bls-oversample", type=int, default=20)
    parser.add_argument("--min-bls-depth-snr", type=float, default=5.0)
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    try:
        reference_flux_columns = parse_reference_flux_column_mappings(
            args.reference_flux_column
        )
    except (TypeError, ValueError) as exc:
        parser.error(str(exc))
    provenance = IndependentExtractionProvenance(
        current_repository=args.current_repository,
        current_revision=args.current_revision,
        reference_extractor_family=args.reference_extractor_family,
        reference_repository=args.reference_repository,
        reference_revision=args.reference_revision,
        reference_pixel_source=args.reference_pixel_source,
        reference_product_aperture=args.reference_product_aperture,
        reference_target_id=args.reference_target_id,
        independence_basis=args.independence_basis,
    )
    _, manifest = build_wd1856_independent_metrics(
        compact_lc=args.compact_lc,
        reference_table=args.reference_table,
        reference_product=args.reference_product,
        cadence_reference_table=args.cadence_reference_table,
        cadence_reference_manifest=args.cadence_reference_manifest,
        metrics_csv=args.metrics_csv,
        manifest_json=args.manifest_json,
        provenance=provenance,
        reference_time_system=args.reference_time_system,
        reference_flux_columns=reference_flux_columns,
        reference_cadence_column=args.reference_cadence_column,
        reference_time_column=args.reference_time_column,
        reference_quality_column=args.reference_quality_column,
        reference_tic_column=args.reference_tic_column,
        reference_sector_column=args.reference_sector_column,
        transit_duration_min=args.transit_duration_min,
        oot_exclusion_min=args.oot_exclusion_min,
        max_time_delta_seconds=args.max_time_delta_seconds,
        min_common_cadences=args.min_common_cadences,
        min_in_event_cadences=args.min_in_event_cadences,
        min_out_of_event_cadences=args.min_out_of_event_cadences,
        bls_period_relative_half_width=args.bls_period_relative_half_width,
        bls_n_periods=args.bls_n_periods,
        bls_durations_min=tuple(args.bls_durations_min or (6.0, 8.0, 10.0)),
        bls_oversample=args.bls_oversample,
        min_bls_depth_snr=args.min_bls_depth_snr,
        overwrite=args.overwrite,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
