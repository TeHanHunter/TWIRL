#!/usr/bin/env python3
"""Export the exact TIC union represented by completed TGLC HDF5 files."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector", type=int, required=True)
    parser.add_argument(
        "--orbit-root",
        type=Path,
        action="append",
        required=True,
        help="Orbit ffi root to scan; supply once per sector orbit.",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def discover(root: Path) -> tuple[set[int], int, int]:
    if not root.is_dir():
        raise FileNotFoundError(root)
    ids: set[int] = set()
    n_hdf5_paths = 0
    n_duplicate_tic_paths = 0
    # The TGLC production layout is ``ffi/camC/ccdD/LC/<tic>.h5``.  Using
    # this constrained glob avoids traversing the much larger source/ePSF
    # trees when the inventory runs on PDO network storage.
    for path in root.glob("cam*/ccd*/LC/*.h5"):
        if path.stat().st_size == 0:
            continue
        try:
            tic = int(path.stem)
        except ValueError:
            continue
        if tic <= 0:
            raise ValueError(f"Non-positive TIC filename: {path}")
        n_hdf5_paths += 1
        if tic in ids:
            n_duplicate_tic_paths += 1
            continue
        ids.add(tic)
    if not ids:
        raise ValueError(f"No nonzero TIC HDF5 files found in {root}")
    return ids, n_hdf5_paths, n_duplicate_tic_paths


def main() -> int:
    args = parse_args()
    if args.sector <= 0:
        raise ValueError("Sector must be positive")
    roots = [root.resolve() for root in args.orbit_root]
    if len(set(roots)) != len(roots):
        raise ValueError("Each --orbit-root must be unique")

    per_orbit = {str(root): discover(root) for root in roots}
    tic_ids = sorted(set().union(*(ids for ids, _paths, _duplicates in per_orbit.values())))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    orbit_label = "_".join(root.parent.name.removeprefix("orbit-") for root in roots)
    stem = f"tic_ids_s{args.sector:04d}_lightcurves_orbits{orbit_label}"
    csv_path = args.output_dir / f"{stem}.csv"
    txt_path = args.output_dir / f"{stem}.txt"
    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("tic_id",))
        writer.writerows((tic,) for tic in tic_ids)
    txt_path.write_text("".join(f"{tic}\n" for tic in tic_ids))
    metadata = {
        "sector": args.sector,
        "source_roots": list(per_orbit),
        "n_tics": len(tic_ids),
        "n_hdf5_paths_by_orbit_root": {
            root: n_paths for root, (_ids, n_paths, _duplicates) in per_orbit.items()
        },
        "n_unique_tics_by_orbit_root": {
            root: len(ids) for root, (ids, _paths, _duplicates) in per_orbit.items()
        },
        "n_duplicate_tic_paths_by_orbit_root": {
            root: n_duplicates
            for root, (_ids, _paths, n_duplicates) in per_orbit.items()
        },
        "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "csv": str(csv_path),
        "csv_sha256": sha256(csv_path),
        "txt": str(txt_path),
        "txt_sha256": sha256(txt_path),
        "note": (
            "Union of positive TIC IDs with nonzero TGLC HDF5 light curves. "
            "Duplicate TIC paths within an orbit are retained as audit counts but listed once."
        ),
    }
    metadata_path = args.output_dir / f"{stem}_meta.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Wrote {len(tic_ids)} TIC IDs for S{args.sector:02d} to {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
