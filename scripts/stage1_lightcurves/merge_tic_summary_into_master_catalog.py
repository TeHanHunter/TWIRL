#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
from astropy.table import Table


DEFAULT_INPUT_CATALOG = Path("data_local/catalogs/twirl_master_catalog/twirl_wd_master_catalog_v0.fits")
DEFAULT_INPUT_SUMMARY = Path("data_local/catalogs/twirl_master_catalog/gaia_dr3_to_tic_summary.csv")
DEFAULT_OUTPUT_CATALOG = Path(
    "data_local/catalogs/twirl_master_catalog/twirl_wd_master_catalog_v0_ticmatched.fits"
)
DEFAULT_OUTPUT_MANIFEST = Path(
    "data_local/catalogs/twirl_master_catalog/twirl_wd_master_catalog_v0_ticmatched_manifest.json"
)
DEFAULT_OUTPUT_REPORT = Path(
    "data_local/catalogs/twirl_master_catalog/gaia_dr3_to_tic_merge_report.csv"
)

VALID_STATUSES = {
    "no_dr2_match",
    "dr2_no_tic",
    "unique_tic",
    "multiple_tic",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge the Gaia DR3 to TIC summary CSV back into the TWIRL master catalog. "
            "Only unique TIC matches are written into tic_id by default."
        )
    )
    parser.add_argument(
        "--input-catalog",
        type=Path,
        default=DEFAULT_INPUT_CATALOG,
        help="Input TWIRL master catalog FITS file.",
    )
    parser.add_argument(
        "--input-summary",
        type=Path,
        default=DEFAULT_INPUT_SUMMARY,
        help="Input one-row-per-Gaia summary CSV from export_gaia_dr3_tic_matches.py.",
    )
    parser.add_argument(
        "--output-catalog",
        type=Path,
        default=DEFAULT_OUTPUT_CATALOG,
        help="Output FITS file with merged TIC metadata.",
    )
    parser.add_argument(
        "--output-manifest",
        type=Path,
        default=DEFAULT_OUTPUT_MANIFEST,
        help="Output JSON manifest summarizing the merge.",
    )
    parser.add_argument(
        "--output-report",
        type=Path,
        default=DEFAULT_OUTPUT_REPORT,
        help="Output CSV report containing per-status counts after merge.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files.",
    )
    return parser.parse_args()


def ensure_writable(paths, overwrite: bool) -> None:
    for path in paths:
        if path.exists() and not overwrite:
            raise FileExistsError(f"{path} already exists. Use --overwrite to replace it.")
        path.parent.mkdir(parents=True, exist_ok=True)


def load_summary_rows(path: Path) -> Dict[int, Dict[str, str]]:
    summary_by_gaia_id: Dict[int, Dict[str, str]] = {}
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            gaia_id = int(row["gaia_dr3_source_id"])
            summary_by_gaia_id[gaia_id] = row
    if not summary_by_gaia_id:
        raise ValueError(f"No summary rows found in {path}")
    return summary_by_gaia_id


def load_catalog(path: Path) -> Table:
    table = Table.read(path, hdu=1, unit_parse_strict="silent")
    table.convert_bytestring_to_unicode()
    return table


def build_index(source_ids: np.ndarray) -> Dict[int, int]:
    return {int(source_id): idx for idx, source_id in enumerate(source_ids)}


def merge_summary_into_catalog(
    table: Table,
    summary_by_gaia_id: Dict[int, Dict[str, str]],
) -> Tuple[Table, Dict[str, int]]:
    source_ids = np.asarray(table["source_id"], dtype=np.int64)
    index_by_source_id = build_index(source_ids)

    tic_id = np.asarray(table["tic_id"], dtype=np.int64).copy()
    tic_match_status = np.asarray(table["tic_match_status"]).astype("U32")

    counts = {
        "catalog_rows": int(len(table)),
        "summary_rows": int(len(summary_by_gaia_id)),
        "merged_rows": 0,
        "missing_summary_rows": 0,
        "status_no_dr2_match": 0,
        "status_dr2_no_tic": 0,
        "status_unique_tic": 0,
        "status_multiple_tic": 0,
        "status_unknown": 0,
        "tic_id_filled_unique": 0,
        "catalog_rows_without_summary": 0,
    }

    for gaia_id, row in summary_by_gaia_id.items():
        idx = index_by_source_id.get(gaia_id)
        if idx is None:
            counts["missing_summary_rows"] += 1
            continue

        status = row["tic_match_status"].strip()
        if status not in VALID_STATUSES:
            status = "unknown"
            counts["status_unknown"] += 1
        else:
            counts[f"status_{status}"] += 1

        tic_match_status[idx] = status
        tic_id[idx] = -1

        if status == "unique_tic" and row["unique_tic_id"].strip():
            tic_id[idx] = int(row["unique_tic_id"])
            counts["tic_id_filled_unique"] += 1

        counts["merged_rows"] += 1

    counts["catalog_rows_without_summary"] = counts["catalog_rows"] - counts["merged_rows"]

    table["tic_id"] = tic_id
    table["tic_match_status"] = tic_match_status
    table.meta["TICMTSRC"] = "gaia_dr3_to_tic_summary.csv"
    table.meta["TICMTIME"] = datetime.now(timezone.utc).isoformat()

    return table, counts


def write_report(path: Path, counts: Dict[str, int]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        for key in sorted(counts):
            writer.writerow([key, counts[key]])


def main() -> None:
    args = parse_args()
    ensure_writable(
        [args.output_catalog, args.output_manifest, args.output_report],
        overwrite=args.overwrite,
    )

    summary_by_gaia_id = load_summary_rows(args.input_summary)
    table = load_catalog(args.input_catalog)
    merged_table, counts = merge_summary_into_catalog(table, summary_by_gaia_id)

    merged_table.write(args.output_catalog, overwrite=args.overwrite)
    write_report(args.output_report, counts)

    manifest = {
        "built_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_catalog_path": str(args.input_catalog.resolve()),
        "input_summary_path": str(args.input_summary.resolve()),
        "output_catalog_path": str(args.output_catalog.resolve()),
        "output_report_path": str(args.output_report.resolve()),
        "counts": counts,
        "merge_policy": [
            "tic_match_status is populated from the summary CSV for all matched Gaia DR3 rows.",
            "tic_id is filled only for rows with tic_match_status == unique_tic.",
            "ambiguous or unmatched rows keep tic_id = -1.",
            "tic_match_sep_arcsec is left unchanged because the exporter does not yet compute separations.",
        ],
    }
    with args.output_manifest.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(f"[merge] wrote merged catalog: {args.output_catalog}")
    print(f"[merge] wrote merge manifest: {args.output_manifest}")
    print(f"[merge] wrote merge report: {args.output_report}")
    print(
        "[merge] summary: "
        f"unique_tic={counts['status_unique_tic']}, "
        f"multiple_tic={counts['status_multiple_tic']}, "
        f"dr2_no_tic={counts['status_dr2_no_tic']}, "
        f"no_dr2_match={counts['status_no_dr2_match']}, "
        f"filled_tic_id={counts['tic_id_filled_unique']}"
    )


if __name__ == "__main__":
    main()
