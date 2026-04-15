#!/usr/bin/env python3
"""Backfill LightCurve/QualityFlag = zeros(len(BJD), uint16) into TGLC h5s.

QLP's `lctools hlsp` reads `LightCurve/QualityFlag` unconditionally from every
per-orbit h5 (see qlp/lctools/bin/hlsp.py:346). TGLC's lightcurves stage does
not write that dataset, so hlsp crashes on the first target.

This script walks the given LC directories and writes a zero-filled uint16
QualityFlag of shape (len(BJD),) into any h5 missing it. Existing datasets
are left alone (idempotent).

Usage (on PDO):
    python scripts/stage1_lcs/backfill_qualityflag.py \\
        /pdo/qlp-data/orbit-119/ffi/cam*/ccd*/LC \\
        /pdo/qlp-data/orbit-120/ffi/cam*/ccd*/LC \\
        --nprocs 16
"""

import argparse
import sys
from multiprocessing import Pool
from pathlib import Path

import h5py
import numpy as np


def backfill_one(path: Path) -> tuple[str, str]:
    try:
        with h5py.File(path, "r+") as f:
            lc = f["LightCurve"]
            if "QualityFlag" in lc:
                return ("skip", str(path))
            n = lc["BJD"].shape[0]
            lc.create_dataset("QualityFlag", data=np.zeros(n, dtype=np.uint16))
            return ("wrote", str(path))
    except Exception as e:
        return (f"error:{e}", str(path))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("lc_dirs", nargs="+", type=Path)
    ap.add_argument("--nprocs", type=int, default=16)
    args = ap.parse_args()

    files: list[Path] = []
    for d in args.lc_dirs:
        files.extend(sorted(d.glob("*.h5")))
    print(f"Found {len(files)} h5 files across {len(args.lc_dirs)} dirs", flush=True)
    if not files:
        return 0

    counts = {"wrote": 0, "skip": 0, "error": 0}
    with Pool(args.nprocs) as pool:
        for i, (status, path) in enumerate(
            pool.imap_unordered(backfill_one, files, chunksize=64), start=1
        ):
            key = "error" if status.startswith("error") else status
            counts[key] += 1
            if status.startswith("error"):
                print(f"ERR {path}: {status}", flush=True)
            if i % 50000 == 0:
                print(f"  progress {i}/{len(files)} {counts}", flush=True)

    print(f"Done. {counts}", flush=True)
    return 0 if counts["error"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
