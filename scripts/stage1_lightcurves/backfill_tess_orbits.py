#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

from astropy.table import Table

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.catalogs import (
    attach_tess_observation_columns,
    build_detector_summary_table,
    build_observation_export_table,
    build_sector_orbit_lookup,
    build_sector_summary_table,
    expand_sector_observation_table_to_orbits,
    write_detector_target_tables,
)


DEFAULT_INPUT_CATALOG = Path(
    "data_local/catalogs/twirl_master_catalog/twirl_wd_master_catalog_v0_tesscoverage.fits"
)
DEFAULT_INPUT_OBSERVATIONS = Path(
    "data_local/catalogs/twirl_master_catalog/twirl_wd_tess_observations_v0.fits"
)
DEFAULT_OUTPUT_CATALOG = DEFAULT_INPUT_CATALOG
DEFAULT_OUTPUT_OBSERVATIONS = DEFAULT_INPUT_OBSERVATIONS
DEFAULT_OUTPUT_MANIFEST = Path(
    "data_local/catalogs/twirl_master_catalog/twirl_wd_master_catalog_v0_tesscoverage_manifest.json"
)
DEFAULT_OUTPUT_DETECTOR_SUMMARY = Path(
    "data_local/catalogs/twirl_master_catalog/twirl_wd_tess_detector_summary_v0.csv"
)
DEFAULT_OUTPUT_SECTOR_SUMMARY = Path(
    "data_local/catalogs/twirl_master_catalog/twirl_wd_tess_sector_summary_v0.csv"
)
DEFAULT_OUTPUT_DETECTOR_DIR = Path(
    "data_local/catalogs/twirl_master_catalog/tess_detector_target_tables_v0"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Backfill orbit-aware TESS coverage products from an existing sector-level TWIRL "
            "coverage catalog and observation table, without rerunning tess-point geometry."
        )
    )
    parser.add_argument(
        "--input-catalog",
        type=Path,
        default=DEFAULT_INPUT_CATALOG,
        help="Input coverage-enriched TWIRL master catalog FITS file.",
    )
    parser.add_argument(
        "--input-observations",
        type=Path,
        default=DEFAULT_INPUT_OBSERVATIONS,
        help="Input sector-level observation table FITS file.",
    )
    parser.add_argument(
        "--output-catalog",
        type=Path,
        default=DEFAULT_OUTPUT_CATALOG,
        help="Output FITS master catalog with orbit-aware TESS coverage columns.",
    )
    parser.add_argument(
        "--output-observations",
        type=Path,
        default=DEFAULT_OUTPUT_OBSERVATIONS,
        help="Output FITS observation table with one row per source_id/orbit/sector/camera/ccd hit.",
    )
    parser.add_argument(
        "--output-manifest",
        type=Path,
        default=DEFAULT_OUTPUT_MANIFEST,
        help="Output JSON manifest for the orbit-backfill run.",
    )
    parser.add_argument(
        "--output-detector-summary",
        type=Path,
        default=DEFAULT_OUTPUT_DETECTOR_SUMMARY,
        help="Output CSV with one row per orbit/sector/camera/ccd footprint.",
    )
    parser.add_argument(
        "--output-sector-summary",
        type=Path,
        default=DEFAULT_OUTPUT_SECTOR_SUMMARY,
        help="Output CSV with one row per sector.",
    )
    parser.add_argument(
        "--detector-table-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DETECTOR_DIR,
        help="Directory for per-orbit/camera/ccd TWIRL target tables.",
    )
    parser.add_argument(
        "--export-detector-tables",
        action="store_true",
        help="Write per-orbit/camera/ccd ECSV target tables.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files.",
    )
    return parser.parse_args()


def progress(message: str) -> None:
    print(message, flush=True)


def ensure_writable(paths: list[Path], overwrite: bool) -> None:
    for path in paths:
        if path.exists() and not overwrite:
            raise FileExistsError(f"{path} already exists. Use --overwrite to replace it.")
        path.parent.mkdir(parents=True, exist_ok=True)


def load_table(path: Path) -> Table:
    table = Table.read(path, unit_parse_strict="silent")
    table.convert_bytestring_to_unicode()
    return table


def write_csv_table(path: Path, table: Table, overwrite: bool) -> None:
    table.write(path, format="ascii.csv", overwrite=overwrite)


def clear_existing_detector_tables(output_directory: Path) -> int:
    if not output_directory.exists():
        return 0

    removed = 0
    for path in output_directory.glob("*_twirl_targets.ecsv"):
        path.unlink()
        removed += 1
    return removed


def main() -> None:
    args = parse_args()

    output_paths = [
        args.output_catalog,
        args.output_observations,
        args.output_manifest,
        args.output_detector_summary,
        args.output_sector_summary,
    ]
    ensure_writable(output_paths, overwrite=args.overwrite)

    progress(f"[orbit-backfill] loading catalog from {args.input_catalog}")
    catalog = load_table(args.input_catalog)
    progress(f"[orbit-backfill] loaded catalog rows={len(catalog)}")

    progress(f"[orbit-backfill] loading observations from {args.input_observations}")
    observation_table = load_table(args.input_observations)
    progress(f"[orbit-backfill] loaded observation rows={len(observation_table)}")

    sector_max = int(max(observation_table["sector"])) if len(observation_table) else -1
    sector_orbit_lookup = build_sector_orbit_lookup(sector_max) if sector_max >= 1 else {}
    expanded_observation_table = expand_sector_observation_table_to_orbits(observation_table)
    progress(
        "[orbit-backfill] "
        f"expanded to orbit-aware observation rows={len(expanded_observation_table)}"
    )

    updated_catalog = attach_tess_observation_columns(catalog, expanded_observation_table)
    observation_export_table = build_observation_export_table(updated_catalog, expanded_observation_table)
    detector_summary = build_detector_summary_table(observation_export_table)
    sector_summary = build_sector_summary_table(observation_export_table)

    updated_catalog.write(args.output_catalog, overwrite=args.overwrite)
    observation_export_table.write(args.output_observations, overwrite=args.overwrite)
    write_csv_table(args.output_detector_summary, detector_summary, overwrite=args.overwrite)
    write_csv_table(args.output_sector_summary, sector_summary, overwrite=args.overwrite)

    removed_detector_tables = 0
    detector_table_paths: list[Path] = []
    if args.export_detector_tables:
        args.detector_table_dir.mkdir(parents=True, exist_ok=True)
        if args.overwrite:
            removed_detector_tables = clear_existing_detector_tables(args.detector_table_dir)
        detector_table_paths = write_detector_target_tables(
            observation_export_table=observation_export_table,
            output_directory=args.detector_table_dir,
            overwrite=args.overwrite,
        )
        progress(
            "[orbit-backfill] "
            f"wrote {len(detector_table_paths)} orbit/camera/ccd target tables to {args.detector_table_dir}"
        )

    manifest = {
        "built_at_utc": datetime.now(timezone.utc).isoformat(),
        "operation": "backfill_tess_orbits",
        "input_catalog_path": str(args.input_catalog.resolve()),
        "input_observations_path": str(args.input_observations.resolve()),
        "output_catalog_path": str(args.output_catalog.resolve()),
        "output_observations_path": str(args.output_observations.resolve()),
        "output_detector_summary_path": str(args.output_detector_summary.resolve()),
        "output_sector_summary_path": str(args.output_sector_summary.resolve()),
        "detector_table_dir": str(args.detector_table_dir.resolve()) if args.export_detector_tables else None,
        "removed_detector_table_count": removed_detector_tables,
        "detector_table_count": len(detector_table_paths),
        "n_catalog_rows": int(len(updated_catalog)),
        "n_input_observation_rows": int(len(observation_table)),
        "n_output_observation_rows": int(len(observation_export_table)),
        "max_sector": sector_max,
        "max_orbit": int(max(observation_export_table["orbit"])) if len(observation_export_table) else -1,
        "n_unique_orbits": int(len(set(int(x) for x in observation_export_table["orbit"]))) if len(observation_export_table) else 0,
        "nonstandard_orbits_per_sector": {
            str(sector): len(orbits) for sector, orbits in sector_orbit_lookup.items() if len(orbits) != 2
        },
        "notes": [
            "This run backfills orbit-aware products from an existing sector-level observation table.",
            "It does not rerun tess-point geometry.",
            "Each sector hit is expanded into one row per orbit using the built-in sector-to-orbit mapping.",
            "Detector target tables are grouped by orbit/camera/ccd.",
        ],
    }

    with args.output_manifest.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")

    progress(f"[orbit-backfill] wrote coverage catalog: {args.output_catalog}")
    progress(f"[orbit-backfill] wrote observation table: {args.output_observations}")
    progress(f"[orbit-backfill] wrote detector summary: {args.output_detector_summary}")
    progress(f"[orbit-backfill] wrote sector summary: {args.output_sector_summary}")
    progress(f"[orbit-backfill] wrote manifest: {args.output_manifest}")


if __name__ == "__main__":
    main()
