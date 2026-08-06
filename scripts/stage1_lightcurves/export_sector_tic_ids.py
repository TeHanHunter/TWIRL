#!/usr/bin/env python3
"""Export one-row-per-TIC Sector target lists from a TWIRL observation table.

The resulting simple CSVs are suitable for collaborators who need to query
MAST by TIC ID.  They describe the positive-TIC targets requested for a
completed A2v1 sector, rather than applying an additional white-dwarf
probability cut after production.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from collections.abc import Iterable, Sequence
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.table import Table


DEFAULT_OBSERVATIONS = Path(
    "data_local/catalogs/twirl_master_catalog/twirl_wd_tess_observations_v0.fits"
)
DEFAULT_OUTPUT_DIR = Path("reports/stage1_lightcurves/sector_tic_ids")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--observations",
        type=Path,
        default=DEFAULT_OBSERVATIONS,
        help="TWIRL observation-table FITS input.",
    )
    parser.add_argument(
        "--sectors",
        type=int,
        nargs="+",
        required=True,
        help="Completed sectors to export (for example: 60 61 62 63).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Destination directory for CSVs and their manifest.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(
    path: Path,
    rows: Iterable[Sequence[object]],
    header: Sequence[str],
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    sectors = sorted(set(args.sectors))
    if any(sector <= 0 for sector in sectors):
        raise ValueError("Sectors must be positive integers")
    if not args.observations.is_file():
        raise FileNotFoundError(args.observations)

    table = Table.read(args.observations, memmap=True)
    required = {
        "source_id",
        "Pwd",
        "is_highconf_wd",
        "tic_id",
        "sector",
        "orbit",
    }
    missing = sorted(required - set(table.colnames))
    if missing:
        raise ValueError(f"Observation table is missing required columns: {missing}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    requested_rows: list[tuple[int, int]] = []
    context_rows: list[tuple[int, int, int, float, int, str]] = []
    summary: list[dict[str, object]] = []

    for sector in sectors:
        sector_mask = np.asarray(table["sector"], dtype=int) == sector
        positive_tic = np.asarray(table["tic_id"], dtype=np.int64) > 0
        rows = table[sector_mask & positive_tic]
        if len(rows) == 0:
            raise ValueError(f"No positive-TIC observation rows found for Sector {sector}")

        # A target can appear in both sector orbits.  Its target properties are
        # invariant; fail loudly if a malformed observation table disagrees.
        by_tic: dict[int, tuple[int, float, int, list[int]]] = {}
        for row in rows:
            tic = int(row["tic_id"])
            source_id = int(row["source_id"])
            pwd = float(row["Pwd"])
            highconf = int(bool(row["is_highconf_wd"]))
            orbit = int(row["orbit"])
            prior = by_tic.get(tic)
            if prior is None:
                by_tic[tic] = (source_id, pwd, highconf, [orbit])
            else:
                prior_source, prior_pwd, prior_highconf, orbits = prior
                if (source_id, prior_highconf) != (prior_source, highconf) or not np.isclose(
                    pwd, prior_pwd, rtol=0.0, atol=1e-12, equal_nan=True
                ):
                    raise ValueError(
                        f"Sector {sector} TIC {tic} has inconsistent target metadata across visits"
                    )
                if orbit not in orbits:
                    orbits.append(orbit)

        ordered = sorted(by_tic.items())
        tic_rows = [(tic,) for tic, _value in ordered]
        sector_label = f"s{sector:04d}"
        tic_path = args.output_dir / f"tic_ids_{sector_label}_A2v1_targets.csv"
        write_csv(tic_path, tic_rows, ("tic_id",))

        sector_context = []
        for tic, (source_id, pwd, highconf, orbits) in ordered:
            orbit_text = ";".join(str(orbit) for orbit in sorted(orbits))
            sector_context.append((sector, tic, source_id, pwd, highconf, orbit_text))
        context_path = args.output_dir / f"tic_ids_{sector_label}_A2v1_target_context.csv"
        write_csv(
            context_path,
            sector_context,
            ("sector", "tic_id", "gaia_dr3_source_id", "pwd", "is_highconf_wd", "orbits"),
        )
        requested_rows.extend((sector, tic) for tic, _value in ordered)
        context_rows.extend(sector_context)
        orbit_counts = Counter(
            orbit
            for _tic, (_source, _pwd, _highconf, orbits) in ordered
            for orbit in orbits
        )
        n_highconf = sum(
            highconf for _tic, (_source, _pwd, highconf, _orbits) in ordered
        )
        summary.append(
            {
                "sector": sector,
                "n_unique_tics": len(ordered),
                "n_highconf_pwd_gt_0p75": n_highconf,
                "n_lower_pwd_candidates": len(ordered) - n_highconf,
                "orbits": sorted(orbit_counts),
                "n_tics_by_orbit": {
                    str(orbit): orbit_counts[orbit] for orbit in sorted(orbit_counts)
                },
                "tic_csv": tic_path.name,
                "context_csv": context_path.name,
                "tic_csv_sha256": sha256(tic_path),
                "context_csv_sha256": sha256(context_path),
            }
        )

    sector_range = f"s{sectors[0]:04d}_s{sectors[-1]:04d}"
    combined_tic_path = args.output_dir / f"tic_ids_{sector_range}_A2v1_targets.csv"
    write_csv(combined_tic_path, requested_rows, ("sector", "tic_id"))
    combined_context_path = args.output_dir / f"tic_ids_{sector_range}_A2v1_target_context.csv"
    write_csv(
        combined_context_path,
        context_rows,
        ("sector", "tic_id", "gaia_dr3_source_id", "pwd", "is_highconf_wd", "orbits"),
    )
    manifest = {
        "schema_version": "twirl_sector_tic_export_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "selection": (
            "All unique positive TIC IDs in the supplied TWIRL observation table for each sector; "
            "no post-production Pwd cut is applied."
        ),
        "input_observation_table": str(args.observations),
        "input_observation_table_sha256": sha256(args.observations),
        "sectors": summary,
        "combined_tic_csv": combined_tic_path.name,
        "combined_tic_csv_sha256": sha256(combined_tic_path),
        "combined_context_csv": combined_context_path.name,
        "combined_context_csv_sha256": sha256(combined_context_path),
        "n_combined_sector_tic_rows": len(requested_rows),
        "n_unique_tics_across_sectors": len({tic for _sector, tic in requested_rows}),
    }
    manifest_path = args.output_dir / f"tic_ids_{sector_range}_A2v1_targets_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote {len(summary)} sector lists to {args.output_dir}")
    for item in summary:
        print(
            f"S{item['sector']:02d}: {item['n_unique_tics']} TICs "
            f"({item['n_highconf_pwd_gt_0p75']} Pwd > 0.75)"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
