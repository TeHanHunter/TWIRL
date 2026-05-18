#!/usr/bin/env python3
"""Drop TGLC-only cadences from TGLC h5 light curves so they align with a QLP
camera quaternion file.

Motivating case (2026-04-15): Sector 56 cam3/orbit-120. TGLC wrote `6137`
cadences per light curve while `cam3_quat.txt` covers `6136` cadences. Every
cam3 target then failed `qlp lctools detrend` in `quatspline.fit_trend` with a
length mismatch. The quat file is authoritative for spacecraft pointing, so
the right fix is to remove the TGLC cadence(s) not present in the quat file
(typically one extra frame at the orbit boundary) from the h5s.

Usage on PDO (inside the QLP env, but h5py is fine under the TGLC env too):

    python scripts/stage1_lightcurves/fix_tglc_quat_cadence_mismatch.py \\
        --lc-dir /pdo/users/tehan/tglc-deep-catalogs/orbit-120/ffi/cam3/ccd1/LC \\
        --lc-dir /pdo/users/tehan/tglc-deep-catalogs/orbit-120/ffi/cam3/ccd2/LC \\
        --lc-dir /pdo/users/tehan/tglc-deep-catalogs/orbit-120/ffi/cam3/ccd3/LC \\
        --lc-dir /pdo/users/tehan/tglc-deep-catalogs/orbit-120/ffi/cam3/ccd4/LC \\
        --quat   /pdo/users/tehan/tglc-deep-catalogs/orbit-120/ffi/run/cam3_quat.txt \\
        --dry-run

Drop `--dry-run` to actually mutate the h5 files. The script refuses to run
if it would need to drop more than `--max-drop` cadences per file (default 4),
which catches the case where the quat file and the h5 are from genuinely
different orbits rather than off-by-one at the boundary.
"""

from __future__ import annotations

# Cap BLAS/OpenMP threads before numpy import — required for any
# multiprocessing pool on PDO (see CLAUDE.md).
import os as _os

for _v in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
):
    _os.environ.setdefault(_v, "1")

import argparse
import multiprocessing as mp
import sys
from functools import partial
from pathlib import Path

import h5py
import numpy as np


def load_quat_cadences(quat_path: Path) -> set[int]:
    """Read the `cadence` column of a QLP `camC_quat.txt` file.

    The file is CSV with a header like `q1_mean,...,flag,cadence`. We locate
    the `cadence` column by name rather than by position.
    """
    import csv

    cadences: list[int] = []
    with quat_path.open(newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        try:
            cad_col = header.index("cadence")
        except ValueError:
            raise RuntimeError(
                f"{quat_path}: no 'cadence' column in header {header!r}"
            )
        for row in reader:
            if not row:
                continue
            try:
                cadences.append(int(float(row[cad_col])))
            except (ValueError, IndexError):
                continue
    if not cadences:
        raise RuntimeError(f"no cadences parsed from {quat_path}")
    return set(cadences)


def find_cadence_dataset(f: h5py.File) -> str:
    """Locate the cadence dataset path in a TGLC h5."""
    for candidate in ("LightCurve/Cadence", "Cadence"):
        if candidate in f:
            return candidate
    for key in f:
        group = f[key]
        if isinstance(group, h5py.Group) and "Cadence" in group:
            return f"{key}/Cadence"
    raise KeyError("no Cadence dataset found in h5")


def fix_one_h5(
    h5_path: Path,
    quat_cadences: set[int],
    max_drop: int,
    dry_run: bool,
) -> tuple[str, int]:
    """Returns (status, n_dropped).

    status is one of: 'ok', 'skip_already_aligned', 'skip_cadences_missing',
    'error_too_many_drops'.
    """
    mode = "r" if dry_run else "r+"
    with h5py.File(h5_path, mode) as f:
        cad_key = find_cadence_dataset(f)
        cad = np.asarray(f[cad_key][...], dtype=np.int64)
        n = cad.size

        in_quat = np.fromiter(
            (int(c) in quat_cadences for c in cad), dtype=bool, count=n
        )
        keep_idx = np.flatnonzero(in_quat)
        drop_idx = np.flatnonzero(~in_quat)

        if drop_idx.size == 0:
            return ("skip_already_aligned", 0)
        if drop_idx.size > max_drop:
            return ("error_too_many_drops", int(drop_idx.size))

        # Any equal-length dataset anywhere in the file gets the same trim.
        to_trim: list[str] = []

        def collect(name: str, obj: object) -> None:
            if isinstance(obj, h5py.Dataset) and obj.ndim >= 1 and obj.shape[0] == n:
                to_trim.append(name)

        f.visititems(collect)

        if dry_run:
            return ("ok", int(drop_idx.size))

        for name in to_trim:
            ds = f[name]
            arr = ds[...]
            new_arr = np.take(arr, keep_idx, axis=0)
            attrs = {k: ds.attrs[k] for k in ds.attrs}
            parent_name, _, leaf = name.rpartition("/")
            parent = f[parent_name] if parent_name else f
            # Recreate at the new length. Preserve compression if any.
            compression = ds.compression
            compression_opts = ds.compression_opts
            chunks = ds.chunks
            dtype = ds.dtype
            del parent[leaf]
            new_ds = parent.create_dataset(
                leaf,
                data=new_arr,
                dtype=dtype,
                compression=compression,
                compression_opts=compression_opts,
                chunks=chunks if chunks is None else True,
            )
            for k, v in attrs.items():
                new_ds.attrs[k] = v

        return ("ok", int(drop_idx.size))


def _worker(path: Path, quat_cadences: set[int], max_drop: int, dry_run: bool):
    try:
        status, n_dropped = fix_one_h5(path, quat_cadences, max_drop, dry_run)
        return (str(path), status, n_dropped, None)
    except Exception as e:
        return (str(path), "exception", 0, f"{type(e).__name__}: {e}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--lc-dir",
        action="append",
        required=True,
        type=Path,
        help="directory of TGLC *.h5 files (repeatable)",
    )
    ap.add_argument(
        "--quat", required=True, type=Path, help="path to camC_quat.txt"
    )
    ap.add_argument("--max-drop", type=int, default=4)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--nprocs", type=int, default=1, help="parallel workers (NFS-bound)")
    args = ap.parse_args(argv)

    quat_cadences = load_quat_cadences(args.quat)
    print(f"[quat] {args.quat}: {len(quat_cadences)} cadences", flush=True)

    all_files: list[Path] = []
    for lc_dir in args.lc_dir:
        files = sorted(lc_dir.glob("*.h5"))
        print(f"[dir] {lc_dir}: {len(files)} h5 files", flush=True)
        all_files.extend(files)

    totals: dict[str, int] = {}
    tag = "dry " if args.dry_run else "fix "
    worker = partial(
        _worker,
        quat_cadences=quat_cadences,
        max_drop=args.max_drop,
        dry_run=args.dry_run,
    )

    def _handle(result):
        path_str, status, n_dropped, err = result
        name = Path(path_str).name
        totals[status] = totals.get(status, 0) + 1
        if status == "ok":
            print(f"[{tag}] {name}: dropped {n_dropped} cadence(s)", flush=True)
        elif status == "error_too_many_drops":
            print(
                f"[err ] {name}: would drop {n_dropped} cadences "
                f"(> --max-drop={args.max_drop}); refusing",
                flush=True,
            )
        elif status == "exception":
            print(f"[err ] {name}: {err}", flush=True)

    if args.nprocs <= 1:
        for path in all_files:
            _handle(worker(path))
    else:
        with mp.Pool(processes=args.nprocs) as pool:
            for result in pool.imap_unordered(worker, all_files, chunksize=4):
                _handle(result)

    print(f"[done] {totals}", flush=True)
    return 0 if totals.get("error_too_many_drops", 0) == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
