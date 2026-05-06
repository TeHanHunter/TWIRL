#!/usr/bin/env python3
"""Render a 3-panel BLS diagnostic plot from a saved periodogram + HLSP.

Standard panels: full SDE periodogram, ±1% zoom on rank-1 peak, phase-folded
light curve with BLS box overlay. Used to visually verify every BLS search.

Example
-------
    python scripts/stage2_search/plot_bls_diagnostic.py \\
        --npz benchmark/wd1856_bls/sector_0056/periodograms/tic_267574918_DET_FLUX.npz \\
        --hlsp benchmark/gpu_production/hlsp_qlp_tess_ffi_s0056-0000000267574918_tess_v01_llc.fits \\
        --candidates benchmark/wd1856_bls/sector_0056/candidates.parquet \\
        --out benchmark/wd1856_bls/wd1856_s56_bls_diagnostic.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

import numpy as np

from twirl.search.diagnostics import plot_bls_diagnostic
from twirl.io.hlsp import read_hlsp


def _peak_from_candidates(candidates_parquet: Path, tic: int, aperture: str,
                          peak_rank: int) -> dict[str, float]:
    import pyarrow.parquet as pq
    df = pq.read_table(candidates_parquet).to_pandas()
    sub = df[(df["tic"] == tic) & (df["aperture"] == aperture)
             & (df["peak_rank"] == peak_rank)]
    if sub.empty:
        raise RuntimeError(
            f"no rank-{peak_rank} peak for TIC {tic} aperture {aperture} "
            f"in {candidates_parquet}"
        )
    r = sub.iloc[0]
    return {
        "period_d": float(r["period_d"]),
        "t0_bjd": float(r["t0_bjd"]),
        "duration_min": float(r["duration_min"]),
        "depth": float(r["depth"]),
        "sde": float(r["sde"]),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--npz", type=Path, required=True,
                    help="Path to periodogram NPZ saved by save_periodogram().")
    ap.add_argument("--hlsp", type=Path, required=True,
                    help="HLSP FITS the periodogram was computed from.")
    ap.add_argument("--candidates", type=Path, default=None,
                    help="Optional candidates.parquet to pull the peak parameters from.")
    ap.add_argument("--peak-rank", type=int, default=1,
                    help="Which peak to render (1=rank-1; only used with --candidates).")
    ap.add_argument("--aperture", type=str, default="DET_FLUX")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--title", type=str, default=None)
    args = ap.parse_args()

    lc = read_hlsp(args.hlsp)
    if lc is None:
        print(f"[plot] could not read HLSP {args.hlsp}", file=sys.stderr)
        return 2

    with np.load(args.npz) as z:
        spectrum = {k: np.asarray(z[k]) for k in z.files}

    if args.candidates:
        peak = _peak_from_candidates(args.candidates, lc.tic, args.aperture,
                                     args.peak_rank)
    else:
        # Fallback: derive rank-1 from the periodogram itself.
        sde = spectrum["sde"]
        period = spectrum["period"]
        idx = int(np.nanargmax(sde))
        from twirl.io.hlsp import BJDREFI
        peak = {
            "period_d": float(period[idx]),
            "t0_bjd": float(spectrum["t0"][idx]) + float(BJDREFI),
            "duration_min": float(spectrum["duration"][idx]) * 1440.0,
            "depth": float(spectrum["depth"][idx]),
            "sde": float(sde[idx]),
        }

    out = plot_bls_diagnostic(spectrum, lc, args.aperture, peak, args.out,
                              title=args.title)
    print(f"[plot] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
