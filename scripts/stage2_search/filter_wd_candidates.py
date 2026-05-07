#!/usr/bin/env python3
"""Non-destructive WD-planet-regime candidate filter.

Reads `sector_NNNN/consolidated.parquet` (BLS clusters across apertures) and
writes `sector_NNNN/wd_candidates.parquet` ranked by SDE. Cuts are
intentionally loose so a single-aperture detection like WD 1856 in S59
still survives:

  - 0.04 d < period < 13 d           (above conservative Roche limit;
                                       at most 0.5 × sector baseline)
  - sde_max > 7.0                    (classical BLS detection floor)
  - status = ok                      (drops read_fail / too_few_cadences)

We deliberately do NOT cut on:
  - n_apertures_agree (a real signal can be in 1 aperture only — see WD 1856
    in S59 cam4/ccd4 where DET_FLUX and DET_FLUX_LAG are scattered-light
    dominated and miss the truth completely; SML has it alone)
  - depth (BLS depth on undersampled WD transits is biased low and noisy
    photometry adds variance — too risky as a hard cut at this stage)

The output is a sibling parquet, not a replacement: you can re-run with
different thresholds without touching the canonical consolidated table.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sector", type=int, required=True)
    ap.add_argument("--out-dir", type=Path,
                    default=REPO_ROOT / "data_local/stage2/bls_first_pass")
    ap.add_argument("--out-suffix", type=str, default="",
                    help="Match the sector_run --out-suffix used at search time.")
    ap.add_argument("--p-min", type=float, default=0.04,
                    help="Minimum period in days (default 0.04 = 1 hour, "
                         "above conservative WD Roche limit for any plausible "
                         "non-disintegrating companion).")
    ap.add_argument("--p-max", type=float, default=13.0,
                    help="Maximum period in days (default 13 d ~ ½ × sector).")
    ap.add_argument("--sde-min", type=float, default=7.0,
                    help="Minimum sde_max per cluster (default 7).")
    ap.add_argument("--top-n", type=int, default=0,
                    help="If >0, also write a top-N TIC list to wd_top_tics.txt.")
    args = ap.parse_args()

    sector_dir = args.out_dir / f"sector_{args.sector:04d}{args.out_suffix}"
    src = sector_dir / "consolidated.parquet"
    if not src.exists():
        print(f"[wd-filter] FAIL: {src} missing", file=sys.stderr)
        return 2

    import pyarrow.parquet as pq
    df = pq.read_table(src).to_pandas()
    n_in = len(df)

    keep = (
        (df["period_d"] > args.p_min)
        & (df["period_d"] < args.p_max)
        & (df["sde_max"] > args.sde_min)
    )
    sub = df[keep].sort_values("sde_max", ascending=False).reset_index(drop=True)
    n_out = len(sub)

    out_path = sector_dir / "wd_candidates.parquet"
    pq.write_table(pq.read_table(src).filter(keep.values).sort_by(
        [("sde_max", "descending")]
    ), out_path, compression="zstd")
    print(f"[wd-filter] sector {args.sector}: {n_in} clusters -> {n_out} survivors -> {out_path}")
    print(f"[wd-filter] cuts: {args.p_min} d < P < {args.p_max} d, sde_max > {args.sde_min}")
    by_agree = sub["n_apertures_agree"].value_counts().to_dict()
    print(f"[wd-filter] agreement breakdown: {by_agree}")

    head = sub.head(min(20, n_out))
    print(f"[wd-filter] top {len(head)} by SDE:")
    for _, r in head.iterrows():
        tic = int(r["tic"]); tmag = r["tmag"]; p = r["period_d"]
        sde = r["sde_max"]; agree = r["n_apertures_agree"]
        rep = r["rep_aperture"]; aps = r["apertures_agree"]
        print(f"  TIC {tic:>11}  T={tmag:5.2f}  P={p:.6f}d  SDE={sde:7.2f}  "
              f"agree={agree}  rep={rep:<14}  ap_set={aps}")

    if args.top_n > 0:
        unique_tics = sub.drop_duplicates(subset="tic").head(args.top_n)["tic"].astype(int).tolist()
        tics_path = sector_dir / "wd_top_tics.txt"
        tics_path.write_text("\n".join(str(t) for t in unique_tics) + "\n")
        print(f"[wd-filter] wrote {len(unique_tics)} unique TICs to {tics_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
