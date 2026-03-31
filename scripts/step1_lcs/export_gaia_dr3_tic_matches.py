#!/usr/bin/env python

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits


DEFAULT_INPUT_CATALOG = Path("data_local/catalogs/twirl_master_catalog/twirl_wd_master_catalog_v0.fits")
DEFAULT_OUTPUT_DIR = Path("data_local/catalogs/twirl_master_catalog")

MATCH_FIELDS = [
    "id",
    "ra",
    "dec",
    "tmag",
    "pmra",
    "pmdec",
    "jmag",
    "kmag",
    "vmag",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Export a read-only Gaia DR3 to TIC match table for the Gaia source IDs in a local "
            "TWIRL catalog using the PDO pyticdb databases."
        )
    )
    parser.add_argument(
        "--input-catalog",
        type=Path,
        default=DEFAULT_INPUT_CATALOG,
        help="Path to a local FITS catalog containing a Gaia DR3 source_id column.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for the exported CSV and JSON products.",
    )
    parser.add_argument(
        "--input-hdu",
        type=str,
        default="auto",
        help='FITS HDU to read. Use "auto" to try Joined and then extension 1.',
    )
    parser.add_argument(
        "--gaia-id-column",
        type=str,
        default="source_id",
        help="Column name containing Gaia DR3 source IDs.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=5000,
        help="Number of Gaia DR3 IDs per database query chunk.",
    )
    parser.add_argument(
        "--tic-db-name",
        type=str,
        default="tic_82",
        help="pyticdb database name for the TIC database.",
    )
    parser.add_argument(
        "--dr2-to-dr3-table",
        type=str,
        default="dr2_to_dr3",
        help="Table name for the Gaia DR2 to DR3 mapping table inside the TIC database.",
    )
    parser.add_argument(
        "--tic-table",
        type=str,
        default="ticentries",
        help="Table name for the TIC entries table inside the TIC database.",
    )
    parser.add_argument(
        "--matches-name",
        type=str,
        default="gaia_dr3_to_tic_matches.csv",
        help="Filename for the one-row-per-match CSV export.",
    )
    parser.add_argument(
        "--summary-name",
        type=str,
        default="gaia_dr3_to_tic_summary.csv",
        help="Filename for the one-row-per-Gaia summary CSV export.",
    )
    parser.add_argument(
        "--manifest-name",
        type=str,
        default="gaia_dr3_to_tic_export_manifest.json",
        help="Filename for the JSON sidecar manifest.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files.",
    )
    return parser.parse_args()


def load_gaia_source_ids(path: Path, column_name: str, hdu_name: str) -> np.ndarray:
    with fits.open(path, memmap=True) as hdul:
        if hdu_name == "auto":
            if "Joined" in hdul:
                data = hdul["Joined"].data
            else:
                data = hdul[1].data
        else:
            try:
                hdu = int(hdu_name)
            except ValueError:
                hdu = hdu_name
            data = hdul[hdu].data

        if column_name not in data.columns.names:
            available = ", ".join(data.columns.names[:20])
            raise KeyError(
                f"Column {column_name!r} not found in {path}. "
                f"First available columns: {available}"
            )

        gaia_ids = np.asarray(data[column_name], dtype=np.int64)

    if gaia_ids.size == 0:
        raise ValueError(f"No Gaia source IDs found in {path}")

    return gaia_ids


def chunked(values: np.ndarray, chunk_size: int):
    for start in range(0, len(values), chunk_size):
        yield values[start : start + chunk_size]


class ReflectedDatabase:
    """Tiny wrapper around pyticdb's reflected SQLAlchemy metadata."""

    def __init__(self, database_name: str):
        try:
            from pyticdb import Databases
        except ImportError as exc:
            raise SystemExit(
                "This script requires pyticdb and must be run in the PDO/TGLC environment."
            ) from exc

        self.database_name = database_name
        self.meta, self.Session = Databases[database_name]

    def table(self, table_name: str):
        return self.meta.tables[table_name]

    def execute(self, query):
        with self.Session() as db:
            return list(db.execute(query).fetchall())


def fetch_dr2_to_dr3_rows(
    tic_db: ReflectedDatabase,
    table_name: str,
    dr3_ids: np.ndarray,
) -> list[Any]:
    import sqlalchemy as sa

    table = tic_db.table(table_name)
    query = sa.select(table.c.dr2_source_id, table.c.dr3_source_id).where(
        table.c.dr3_source_id.in_([int(value) for value in dr3_ids])
    )
    return tic_db.execute(query)


def fetch_tic_rows(
    tic_db: ReflectedDatabase,
    table_name: str,
    dr2_ids: list[int],
) -> list[Any]:
    import sqlalchemy as sa

    if not dr2_ids:
        return []

    table = tic_db.table(table_name)
    columns = [getattr(table.c, field) for field in MATCH_FIELDS]
    gaia_dr2_source_id = sa.cast(table.c.gaia, sa.BigInteger).label("gaia_dr2_source_id")
    query = sa.select(gaia_dr2_source_id, *columns).where(
        table.c.gaia.in_([str(int(value)) for value in dr2_ids])
    )
    return tic_db.execute(query)


def ensure_output_paths(paths: list[Path], overwrite: bool) -> None:
    for path in paths:
        if path.exists() and not overwrite:
            raise FileExistsError(
                f"{path} already exists. Use --overwrite to replace existing outputs."
            )
        path.parent.mkdir(parents=True, exist_ok=True)


def status_from_counts(n_dr2_matches: int, n_tic_matches: int) -> str:
    if n_dr2_matches == 0:
        return "no_dr2_match"
    if n_tic_matches == 0:
        return "dr2_no_tic"
    if n_tic_matches == 1:
        return "unique_tic"
    return "multiple_tic"


def export_matches(
    input_catalog: Path,
    output_dir: Path,
    input_hdu: str,
    gaia_id_column: str,
    chunk_size: int,
    tic_db_name: str,
    dr2_to_dr3_table: str,
    tic_table: str,
    matches_name: str,
    summary_name: str,
    manifest_name: str,
    overwrite: bool,
) -> dict[str, Any]:
    gaia_ids = load_gaia_source_ids(input_catalog, gaia_id_column, input_hdu)

    matches_path = output_dir / matches_name
    summary_path = output_dir / summary_name
    manifest_path = output_dir / manifest_name
    ensure_output_paths([matches_path, summary_path, manifest_path], overwrite=overwrite)

    tic_db = ReflectedDatabase(tic_db_name)

    total_match_rows = 0
    total_unique_matches = 0
    total_multiple_matches = 0
    total_no_dr2_matches = 0
    total_dr2_no_tic_matches = 0

    with (
        matches_path.open("w", newline="", encoding="utf-8") as matches_handle,
        summary_path.open("w", newline="", encoding="utf-8") as summary_handle,
    ):
        matches_writer = csv.DictWriter(
            matches_handle,
            fieldnames=[
                "gaia_dr3_source_id",
                "gaia_dr2_source_id",
                "tic_id",
                "tic_match_rank_for_gaia_dr3",
                "tic_ra_deg",
                "tic_dec_deg",
                "tic_tmag",
                "tic_pmra_mas_per_year",
                "tic_pmdec_mas_per_year",
                "tic_jmag",
                "tic_kmag",
                "tic_vmag",
            ],
        )
        summary_writer = csv.DictWriter(
            summary_handle,
            fieldnames=[
                "gaia_dr3_source_id",
                "tic_match_status",
                "n_dr2_matches",
                "n_tic_matches",
                "unique_tic_id",
                "unique_tic_tmag",
                "dr2_source_ids_json",
                "tic_ids_json",
            ],
        )
        matches_writer.writeheader()
        summary_writer.writeheader()

        for gaia_chunk in chunked(gaia_ids, chunk_size):
            dr2_rows = fetch_dr2_to_dr3_rows(
                tic_db=tic_db,
                table_name=dr2_to_dr3_table,
                dr3_ids=gaia_chunk,
            )

            dr3_to_dr2: dict[int, list[int]] = {}
            for dr2_source_id, dr3_source_id in dr2_rows:
                dr3_to_dr2.setdefault(int(dr3_source_id), []).append(int(dr2_source_id))

            unique_dr2_ids = sorted(
                {dr2_id for dr2_ids in dr3_to_dr2.values() for dr2_id in dr2_ids}
            )
            tic_rows = fetch_tic_rows(tic_db=tic_db, table_name=tic_table, dr2_ids=unique_dr2_ids)

            dr2_to_tic_rows: dict[int, list[Any]] = {}
            for tic_row in tic_rows:
                dr2_to_tic_rows.setdefault(int(tic_row.gaia_dr2_source_id), []).append(tic_row)

            for gaia_dr3_id in gaia_chunk:
                gaia_dr3_id = int(gaia_dr3_id)
                dr2_ids = sorted(dr3_to_dr2.get(gaia_dr3_id, []))
                matched_tic_rows = []
                for dr2_id in dr2_ids:
                    matched_tic_rows.extend(dr2_to_tic_rows.get(dr2_id, []))
                matched_tic_rows = sorted(matched_tic_rows, key=lambda row: int(row.id))

                tic_ids = [int(row.id) for row in matched_tic_rows]
                status = status_from_counts(len(dr2_ids), len(matched_tic_rows))

                if status == "no_dr2_match":
                    total_no_dr2_matches += 1
                elif status == "dr2_no_tic":
                    total_dr2_no_tic_matches += 1
                elif status == "unique_tic":
                    total_unique_matches += 1
                else:
                    total_multiple_matches += 1

                summary_writer.writerow(
                    {
                        "gaia_dr3_source_id": gaia_dr3_id,
                        "tic_match_status": status,
                        "n_dr2_matches": len(dr2_ids),
                        "n_tic_matches": len(matched_tic_rows),
                        "unique_tic_id": tic_ids[0] if len(tic_ids) == 1 else "",
                        "unique_tic_tmag": matched_tic_rows[0].tmag if len(tic_ids) == 1 else "",
                        "dr2_source_ids_json": json.dumps(dr2_ids),
                        "tic_ids_json": json.dumps(tic_ids),
                    }
                )

                for rank, tic_row in enumerate(matched_tic_rows, start=1):
                    matches_writer.writerow(
                        {
                            "gaia_dr3_source_id": gaia_dr3_id,
                            "gaia_dr2_source_id": int(tic_row.gaia_dr2_source_id),
                            "tic_id": int(tic_row.id),
                            "tic_match_rank_for_gaia_dr3": rank,
                            "tic_ra_deg": tic_row.ra,
                            "tic_dec_deg": tic_row.dec,
                            "tic_tmag": tic_row.tmag,
                            "tic_pmra_mas_per_year": tic_row.pmra,
                            "tic_pmdec_mas_per_year": tic_row.pmdec,
                            "tic_jmag": tic_row.jmag,
                            "tic_kmag": tic_row.kmag,
                            "tic_vmag": tic_row.vmag,
                        }
                    )
                    total_match_rows += 1

    manifest = {
        "built_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_catalog_path": str(input_catalog.resolve()),
        "input_hdu": input_hdu,
        "gaia_id_column": gaia_id_column,
        "tic_db_name": tic_db_name,
        "dr2_to_dr3_table": dr2_to_dr3_table,
        "tic_table": tic_table,
        "chunk_size": chunk_size,
        "row_count_input": int(len(gaia_ids)),
        "row_count_match_rows": int(total_match_rows),
        "count_no_dr2_match": int(total_no_dr2_matches),
        "count_dr2_no_tic": int(total_dr2_no_tic_matches),
        "count_unique_tic": int(total_unique_matches),
        "count_multiple_tic": int(total_multiple_matches),
        "matches_csv_path": str(matches_path.resolve()),
        "summary_csv_path": str(summary_path.resolve()),
        "notes": [
            "This script is read-only with respect to the PDO databases; it only issues SELECT queries.",
            "Matches are exported using the TGLC-style TIC gaia (Gaia DR2) plus dr2_to_dr3 bridge.",
            "The summary CSV keeps one row per Gaia DR3 source ID; the matches CSV keeps one row per TIC match.",
        ],
    }

    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")

    return manifest


def main() -> None:
    args = parse_args()
    manifest = export_matches(
        input_catalog=args.input_catalog,
        output_dir=args.output_dir,
        input_hdu=args.input_hdu,
        gaia_id_column=args.gaia_id_column,
        chunk_size=args.chunk_size,
        tic_db_name=args.tic_db_name,
        dr2_to_dr3_table=args.dr2_to_dr3_table,
        tic_table=args.tic_table,
        matches_name=args.matches_name,
        summary_name=args.summary_name,
        manifest_name=args.manifest_name,
        overwrite=args.overwrite,
    )

    print(f"Wrote matches CSV: {manifest['matches_csv_path']}")
    print(f"Wrote summary CSV: {manifest['summary_csv_path']}")
    print(
        "Summary counts: "
        f"unique_tic={manifest['count_unique_tic']}, "
        f"multiple_tic={manifest['count_multiple_tic']}, "
        f"dr2_no_tic={manifest['count_dr2_no_tic']}, "
        f"no_dr2_match={manifest['count_no_dr2_match']}"
    )


if __name__ == "__main__":
    main()
