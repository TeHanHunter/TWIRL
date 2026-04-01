#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

from astropy.io import fits
from astropy.table import Table

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.catalogs import (
    FIRST_200S_SECTOR,
    SectorCoverageConfig,
    attach_tess_observation_columns,
    build_detector_summary_table,
    build_observation_export_table,
    build_sector_summary_table,
    compute_tess_observation_table,
    write_detector_target_tables,
)


DEFAULT_INPUT_CATALOG = Path("data_local/catalogs/twirl_master_catalog/twirl_wd_master_catalog_v0.fits")
DEFAULT_OUTPUT_CATALOG = Path(
    "data_local/catalogs/twirl_master_catalog/twirl_wd_master_catalog_v0_tesscoverage.fits"
)
DEFAULT_OUTPUT_MANIFEST = Path(
    "data_local/catalogs/twirl_master_catalog/twirl_wd_master_catalog_v0_tesscoverage_manifest.json"
)
DEFAULT_OUTPUT_OBSERVATIONS = Path(
    "data_local/catalogs/twirl_master_catalog/twirl_wd_tess_observations_v0.fits"
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
            "Attach TESS Sector >= 56 coverage metadata to the TWIRL master catalog and export "
            "TWIRL-side detector target tables for planning and filtering."
        )
    )
    parser.add_argument(
        "--input-catalog",
        type=Path,
        default=DEFAULT_INPUT_CATALOG,
        help="Input TWIRL master catalog FITS file.",
    )
    parser.add_argument(
        "--output-catalog",
        type=Path,
        default=DEFAULT_OUTPUT_CATALOG,
        help="Output FITS master catalog with attached TESS coverage JSON and summary columns.",
    )
    parser.add_argument(
        "--output-manifest",
        type=Path,
        default=DEFAULT_OUTPUT_MANIFEST,
        help="Output JSON manifest for the coverage-mapping run.",
    )
    parser.add_argument(
        "--output-observations",
        type=Path,
        default=DEFAULT_OUTPUT_OBSERVATIONS,
        help="Output FITS table with one row per source_id/sector/camera/ccd hit.",
    )
    parser.add_argument(
        "--output-detector-summary",
        type=Path,
        default=DEFAULT_OUTPUT_DETECTOR_SUMMARY,
        help="Output CSV with one row per sector/camera/ccd footprint.",
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
        help="Directory for per-sector/camera/ccd TWIRL target tables.",
    )
    parser.add_argument(
        "--min-sector",
        type=int,
        default=FIRST_200S_SECTOR,
        help="Minimum TESS sector to include. Default is the first 200 s sector.",
    )
    parser.add_argument(
        "--max-sector",
        type=int,
        default=None,
        help="Optional maximum TESS sector to include.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=10_000,
        help="Number of catalog rows to process per sector chunk.",
    )
    parser.add_argument(
        "--max-targets",
        type=int,
        default=None,
        help="Optional cap on the number of input targets, useful for smoke tests.",
    )
    parser.add_argument(
        "--export-detector-tables",
        action="store_true",
        help=(
            "Write per-sector/camera/ccd ECSV target tables. "
            "These are TWIRL planning/filtering tables, not replacements for full-field MKI catalogs."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files.",
    )
    return parser.parse_args()


def ensure_writable(paths: list[Path], overwrite: bool) -> None:
    for path in paths:
        if path.exists() and not overwrite:
            raise FileExistsError(f"{path} already exists. Use --overwrite to replace it.")
        path.parent.mkdir(parents=True, exist_ok=True)


def load_catalog(path: Path, max_targets: int | None = None) -> Table:
    if max_targets is None:
        table = Table.read(path, hdu=1, unit_parse_strict="silent")
        table.convert_bytestring_to_unicode()
        return table

    # For smoke tests, avoid loading the full FITS table before truncating it.
    with fits.open(path, memmap=True) as hdul:
        table = Table(hdul[1].data[:max_targets])
    table.convert_bytestring_to_unicode()
    return table


def write_csv_table(path: Path, table: Table, overwrite: bool) -> None:
    table.write(path, format="ascii.csv", overwrite=overwrite)


def progress(message: str) -> None:
    print(message, flush=True)


def main() -> None:
    args = parse_args()
    output_paths = [
        args.output_catalog,
        args.output_manifest,
        args.output_observations,
        args.output_detector_summary,
        args.output_sector_summary,
    ]
    ensure_writable(output_paths, overwrite=args.overwrite)
    if args.export_detector_tables:
        args.detector_table_dir.mkdir(parents=True, exist_ok=True)

    if args.max_targets is None:
        progress(f"[coverage] loading full catalog from {args.input_catalog}")
    else:
        progress(
            f"[coverage] loading up to {args.max_targets} catalog rows from {args.input_catalog}"
        )

    table = load_catalog(args.input_catalog, max_targets=args.max_targets)
    progress(f"[coverage] loaded catalog rows={len(table)} from {args.input_catalog}")

    config = SectorCoverageConfig(
        min_sector=args.min_sector,
        max_sector=args.max_sector,
        chunk_size=args.chunk_size,
    )
    observation_table, coverage_metadata = compute_tess_observation_table(
        table=table,
        config=config,
        progress=progress,
    )
    progress(
        f"[coverage] observation rows={len(observation_table)} across {coverage_metadata['n_sectors']} sectors"
    )

    catalog_with_coverage = attach_tess_observation_columns(table, observation_table)
    observation_export_table = build_observation_export_table(catalog_with_coverage, observation_table)
    detector_summary = build_detector_summary_table(observation_export_table)
    sector_summary = build_sector_summary_table(observation_export_table)

    catalog_with_coverage.write(args.output_catalog, overwrite=args.overwrite)
    observation_export_table.write(args.output_observations, overwrite=args.overwrite)
    write_csv_table(args.output_detector_summary, detector_summary, overwrite=args.overwrite)
    write_csv_table(args.output_sector_summary, sector_summary, overwrite=args.overwrite)

    detector_table_paths: list[Path] = []
    if args.export_detector_tables:
        detector_table_paths = write_detector_target_tables(
            observation_export_table=observation_export_table,
            output_directory=args.detector_table_dir,
            overwrite=args.overwrite,
        )
        progress(
            f"[coverage] wrote {len(detector_table_paths)} detector target tables to {args.detector_table_dir}"
        )

    manifest = {
        "built_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_catalog_path": str(args.input_catalog.resolve()),
        "output_catalog_path": str(args.output_catalog.resolve()),
        "output_observations_path": str(args.output_observations.resolve()),
        "output_detector_summary_path": str(args.output_detector_summary.resolve()),
        "output_sector_summary_path": str(args.output_sector_summary.resolve()),
        "detector_table_dir": str(args.detector_table_dir.resolve()) if args.export_detector_tables else None,
        "detector_table_count": len(detector_table_paths),
        "max_targets": args.max_targets,
        "coverage_metadata": coverage_metadata,
        "catalog_rows_with_coverage": int(catalog_with_coverage["has_tess_200s_coverage"].sum()),
        "n_unique_tic_observation_rows": int(
            (observation_export_table["tic_id"] > 0).sum() if "tic_id" in observation_export_table.colnames else 0
        ),
        "notes": [
            "Sector coverage uses tess-point detector geometry for all requested sectors.",
            "Observed and future sectors are stored together; later extraction can filter by an observed-sector cutoff.",
            "The per-sector/camera/ccd tables are TWIRL target subsets for planning and filtering.",
            "They do not replace the full-field MKI TGLC cached TIC/Gaia catalogs produced by `tglc catalogs`.",
        ],
    }

    with args.output_manifest.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")

    progress(f"[coverage] wrote coverage catalog: {args.output_catalog}")
    progress(f"[coverage] wrote observation table: {args.output_observations}")
    progress(f"[coverage] wrote detector summary: {args.output_detector_summary}")
    progress(f"[coverage] wrote sector summary: {args.output_sector_summary}")
    progress(f"[coverage] wrote manifest: {args.output_manifest}")


if __name__ == "__main__":
    main()
