#!/usr/bin/env python3
"""Normalize one sector's ADP-only BLS peaks for Teacher-v2 inference."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import h5py
import pandas as pd

from twirl.vetting.teacher_v2 import (
    mark_native_input_availability,
    normalize_real_adp_candidates,
)


def _read(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path) if path.suffix.lower() == ".parquet" else pd.read_csv(path, low_memory=False)


def _target_tics(path: Path) -> set[int]:
    with h5py.File(path, "r") as h5:
        if "targets" not in h5:
            raise KeyError(f"raw source has no /targets group: {path}")
        return {int(key) for key in h5["targets"].keys()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--real-bls-peaks", type=Path, required=True)
    parser.add_argument("--raw-source-h5", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--sector", type=int, required=True)
    parser.add_argument("--small-peaks", type=int, default=10)
    args = parser.parse_args()

    candidates = normalize_real_adp_candidates(
        _read(args.real_bls_peaks),
        small_peaks_per_tic=args.small_peaks,
    )
    sectors = sorted(pd.to_numeric(candidates["sector"], errors="raise").astype(int).unique())
    if sectors != [int(args.sector)]:
        raise ValueError(f"candidate sector mismatch: observed={sectors}, expected={args.sector}")
    candidates = mark_native_input_availability(
        candidates,
        available_tics=_target_tics(args.raw_source_h5),
    )
    top5_mask = pd.to_numeric(candidates["rep_peak_rank"], errors="coerce").le(5)
    top5 = candidates.loc[top5_mask & candidates["native_input_include"]].copy()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    top10_path = args.out_dir / "real_candidates_top10.parquet"
    top5_path = args.out_dir / "real_candidates_top5.parquet"
    candidates.to_parquet(top10_path, compression="zstd", index=False)
    top5.to_parquet(top5_path, compression="zstd", index=False)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sector": int(args.sector),
        "raw_source_h5": str(args.raw_source_h5.resolve()),
        "n_candidates": int(len(candidates)),
        "n_top5_candidates": int(len(top5)),
        "n_top5_native_unavailable": int(
            (top5_mask & ~candidates["native_input_include"]).sum()
        ),
        "n_candidate_tics": int(candidates["tic"].nunique()),
        "n_native_available_tics": int(
            candidates.loc[candidates["native_input_include"], "tic"].nunique()
        ),
        "n_native_unavailable_tics": int(
            candidates.loc[~candidates["native_input_include"], "tic"].nunique()
        ),
        "outputs": {
            "top10": str(top10_path),
            "top5": str(top5_path),
        },
    }
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
