#!/usr/bin/env python3
"""Cross-aperture consolidation of per-sector BLS candidates.

Groups peaks across DET_FLUX_SML / DET_FLUX / DET_FLUX_LAG into clusters
matched in (P, T0). A cluster spanning ≥2 apertures is far more likely a
real transit than a per-aperture systematic.

Usage:
    python scripts/stage2_search/consolidate_candidates.py --sector 56 \\
        --out-dir data_local/stage2/bls_first_pass

Output: <out-dir>/sector_NNNN/consolidated.parquet
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

import pyarrow.parquet as pq

from twirl.search.consolidate import ConsolidateConfig, consolidate_candidates

DEFAULT_OUT_DIR = REPO_ROOT / "data_local/stage2/bls_first_pass"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sector", type=int, required=True)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    ap.add_argument("--out-suffix", type=str, default="",
                    help="Match the sector_run --out-suffix used at search time.")
    ap.add_argument("--sde-min", type=float, default=7.0,
                    help="Drop peaks below this SDE before clustering "
                         "(default 7 — classical BLS detection floor).")
    ap.add_argument("--period-match-rel", type=float, default=0.005,
                    help="Relative period tolerance for clustering (default 0.5%).")
    ap.add_argument("--t0-match-min", type=float, default=5.0,
                    help="Mid-transit phase tolerance in minutes (default 5).")
    args = ap.parse_args()

    sector_dir = args.out_dir / f"sector_{args.sector:04d}{args.out_suffix}"
    candidates = sector_dir / "candidates.parquet"
    if not candidates.exists():
        print(f"[consolidate] FAIL: no candidates parquet at {candidates}",
              file=sys.stderr)
        return 2

    cfg = ConsolidateConfig(
        period_match_rel=float(args.period_match_rel),
        t0_match_min=float(args.t0_match_min),
        sde_min=float(args.sde_min),
    )
    table = consolidate_candidates(candidates, cfg)
    out_path = sector_dir / "consolidated.parquet"
    if table.num_rows == 0:
        print(f"[consolidate] no clusters above SDE>={cfg.sde_min}; "
              f"writing empty table to {out_path}")
        # Write an empty parquet so downstream tooling has a file.
        import pyarrow as pa
        pa.parquet.write_table(pa.table({"tic": []}), out_path)  # minimal
        return 0
    pq.write_table(table, out_path, compression="zstd")
    n_total = table.num_rows
    df = table.to_pandas()
    n_3 = (df["n_apertures_agree"] == 3).sum()
    n_2 = (df["n_apertures_agree"] == 2).sum()
    n_1 = (df["n_apertures_agree"] == 1).sum()
    print(f"[consolidate] sector {args.sector} -> {out_path}")
    print(f"[consolidate] n_clusters={n_total}  3-aperture={n_3}  "
          f"2-aperture={n_2}  1-aperture={n_1}")
    if n_3 + n_2:
        top = df[df["n_apertures_agree"] >= 2].head(10)
        print("[consolidate] top multi-aperture clusters:")
        for _, r in top.iterrows():
            print(f"    TIC {int(r['tic']):>11}  Tmag={r['tmag']:5.2f}  "
                  f"P={r['period_d']:.6f} d  T0={r['t0_bjd']:.4f}  "
                  f"SDE_max={r['sde_max']:.1f}  "
                  f"agree={r['apertures_agree']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
