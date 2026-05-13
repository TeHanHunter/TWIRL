#!/usr/bin/env python3
"""Apply the TWIRL heuristic vetter to a Stage 2 sector-wide BLS run.

Reads ``consolidated.parquet`` (all BLS peaks for all WDs in the sector),
takes the best peak per TIC, computes vetting features, and writes a
ranked ``vetted_per_tic.parquet`` sibling. Reports the rank of WD 1856 (TIC
267574918) before and after vetting cuts to gate the Week-1 sprint
milestone.

Cuts at this stage are FP rejectors, not score modifiers — survivors are
re-ranked by SDE. No ML-style combined score (deferred to Year 2).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

WD_1856_TIC = 267574918


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sector", type=int, default=56)
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=REPO_ROOT / "data_local/stage2/bls_first_pass",
    )
    ap.add_argument(
        "--m-wd",
        type=float,
        default=0.4,
        help="Canonical WD mass (M_sun) for the duration envelope. Default "
             "is the low-mass tail (~1%% of WDs below this) so the envelope "
             "is conservatively generous. The cut barely changes from M=0.4 "
             "to M=1.0; per-WD precision is unnecessary at this stage.",
    )
    ap.add_argument(
        "--r-wd",
        type=float,
        default=0.018,
        help="Canonical WD radius (R_sun) for the duration envelope "
             "(paired with M_WD via the mass-radius relation).",
    )
    ap.add_argument(
        "--r-comp-max",
        type=float,
        default=2.0,
        help="Max companion radius in R_jup (above this is brown-dwarf+; "
             "default 2.0 covers the largest inflated super-Jupiters).",
    )
    ap.add_argument(
        "--grid-max-duration-min",
        type=float,
        default=20.0,
        help="BLS duration-grid ceiling for the run being vetted "
             "(see bls_default.yaml durations_min).",
    )
    ap.add_argument(
        "--p-alias-min-d",
        type=float,
        default=0.10,
        help="Reject candidates with P < this value (default 0.10 d ~ 2.4 h). "
             "Cuts the BLS p_min aliasing wall (~26% of WDs pile up at "
             "the 0.0833 d boundary) and orbits below the iron-core Roche "
             "limit for any plausible companion. Long-term fix is in the "
             "BLS p_min setting.",
    )
    ap.add_argument(
        "--max-period-cluster",
        type=int,
        default=50,
        help="Reject candidates whose fitted period sits in a cluster of "
             ">N TICs within 0.1%% (default 50). Catches instrumental alias "
             "walls beyond p_min that produce sharp pile-ups (e.g. the "
             "P=0.1065 d cluster in S56).",
    )
    args = ap.parse_args()

    import numpy as np
    import pandas as pd
    from twirl.vetting.heuristic import add_vetting_features

    sector_dir = args.out_dir / f"sector_{args.sector:04d}"
    src = sector_dir / "consolidated.parquet"
    if not src.exists():
        print(f"[vet] FAIL: {src} missing", file=sys.stderr)
        return 2

    df = pd.read_parquet(src)
    print(f"[vet] sector {args.sector}: read {len(df):,} BLS peaks "
          f"from {df['tic'].nunique():,} TICs")

    best = df.loc[df.groupby("tic")["sde_max"].idxmax()].copy()
    best = best.sort_values("sde_max", ascending=False).reset_index(drop=True)
    best["blind_rank"] = np.arange(1, len(best) + 1)
    print(f"[vet] per-TIC best peak: {len(best):,} rows")

    vetted = add_vetting_features(
        best,
        grid_max_duration_min=args.grid_max_duration_min,
        p_alias_min_d=args.p_alias_min_d,
        m_wd_msun=args.m_wd,
        r_wd_rsun=args.r_wd,
        r_comp_max_rjup=args.r_comp_max,
    )

    n_total = len(vetted)
    vetted["p_cluster_pass"] = (
        vetted["period_cluster_count"].values <= args.max_period_cluster
    )
    n_p = vetted["p_alias_pass"].sum()
    n_d = vetted["dur_envelope_pass"].sum()
    n_c = vetted["p_cluster_pass"].sum()
    print(
        f"[vet] cuts (independent):\n"
        f"  p_alias_pass (P > {args.p_alias_min_d:.2f} d):                "
        f"{n_p:>6,} / {n_total:,}  ({100 * n_p / n_total:.1f}%)\n"
        f"  dur_envelope_pass (R_comp <= {args.r_comp_max:.1f} R_jup):     "
        f"{n_d:>6,} / {n_total:,}  ({100 * n_d / n_total:.1f}%)\n"
        f"  p_cluster_pass (cluster <= {args.max_period_cluster} TICs):         "
        f"{n_c:>6,} / {n_total:,}  ({100 * n_c / n_total:.1f}%)"
    )
    keep = (
        vetted["p_alias_pass"]
        & vetted["dur_envelope_pass"]
        & vetted["p_cluster_pass"]
    ).values
    survivors = vetted[keep].copy().reset_index(drop=True)
    survivors["vetted_rank"] = np.arange(1, len(survivors) + 1)
    n_out = len(survivors)
    print(f"  combined (AND):                                  "
          f"{n_out:>6,} / {n_total:,}  ({100 * n_out / n_total:.1f}%)")

    # Class-based labeling — vet_class column now produced by
    # add_vetting_features. Non-rejecting: every survivor keeps a label;
    # consumers filter by class. Classes:
    #   planet_candidate        — passes Roche AND below grid ceiling
    #   sub_roche_pceb_suspect  — sub-Roche period (P < 0.216 d at iron density);
    #                             likely WD+M-dwarf / WD+WD binary / pulsator
    #   pceb_grid_ceiling       — dur pegged at BLS ceiling; long-duration PCEB
    #   alias_artifact          — period in BLS alias band or known cluster
    #   duration_violator       — observed dur exceeds chord envelope
    class_counts = survivors["vet_class"].value_counts()
    print(f"\n[vet] vet_class distribution:")
    for cls in ["planet_candidate", "sub_roche_pceb_suspect",
                "pceb_grid_ceiling", "alias_artifact", "duration_violator"]:
        n_c = int(class_counts.get(cls, 0))
        print(f"  {cls:<28} {n_c:>6,}")

    # Per-class ranking (by sde_max, descending) for downstream "top-N within class" reporting.
    survivors["class_rank"] = (
        survivors.groupby("vet_class")["sde_max"]
        .rank(method="dense", ascending=False)
        .astype(int)
    )

    wd_blind = vetted[vetted["tic"] == WD_1856_TIC]
    if not wd_blind.empty:
        r = wd_blind.iloc[0]
        print(f"\n[wd1856] BEFORE vetting:")
        print(f"  blind_rank      : {int(r['blind_rank']):,} / {len(vetted):,}")
        print(f"  period_d        : {r['period_d']:.6f}  (truth 1.40793909)")
        print(f"  duration_min    : {r['duration_min']}")
        print(f"  sde_max         : {r['sde_max']:.2f}")
        print(f"  dur_envelope_min: {r['dur_envelope_min']:.2f} min  "
              f"(b=0 chord for {args.r_comp_max:.1f} R_jup at this P)")
        print(f"  dur_envelope_pass: {bool(r['dur_envelope_pass'])}")
        print(f"  dur_grid_ceiling_hit: {bool(r['dur_grid_ceiling_hit'])}")
        print(f"  p_alias_pass     : {bool(r['p_alias_pass'])}")
        print(f"  period_cluster   : {int(r['period_cluster_count'])} TICs in same period band")

    wd_vetted = survivors[survivors["tic"] == WD_1856_TIC]
    print(f"\n[wd1856] AFTER vetting:")
    if not wd_vetted.empty:
        r = wd_vetted.iloc[0]
        print(f"  vetted_rank (all survivors): {int(r['vetted_rank']):>5,} / {len(survivors):,}")
        print(f"  vet_class                  : {r['vet_class']}")
        print(f"  class_rank                 : {int(r['class_rank'])}")
        print(f"  roche_pass                 : {bool(r['roche_pass'])} (P_Roche={r['roche_period_d']:.3f} d)")
    elif not wd_blind.empty:
        print("  REJECTED by hard vetter cut")
        return 1

    planet = survivors[survivors["vet_class"] == "planet_candidate"]
    print(f"\n[vet] top 20 planet_candidate class (sorted by SDE):")
    for _, r in planet.head(20).iterrows():
        flag = " *WD1856*" if int(r["tic"]) == WD_1856_TIC else ""
        print(
            f"  c#{int(r['class_rank']):>3} (blind #{int(r['blind_rank']):>5}) "
            f"TIC {int(r['tic']):>11}  T={r['tmag']:5.2f}  "
            f"P={r['period_d']:8.4f}d  dur={r['duration_min']:5.1f}min  "
            f"SDE={r['sde_max']:7.2f}  env={r['dur_envelope_min']:5.1f}min  "
            f"agree={int(r['n_apertures_agree'])}{flag}"
        )

    # Single classified output; backward-compat aliases for downstream
    # code that still expects vetted_planet_candidates.parquet.
    survivors.to_parquet(sector_dir / "vetted_per_tic.parquet", compression="zstd")
    planet.to_parquet(sector_dir / "vetted_planet_candidates.parquet", compression="zstd")
    print(f"\n[vet] wrote:")
    print(f"  {sector_dir / 'vetted_per_tic.parquet'}         ({len(survivors):,} rows, all classes)")
    print(f"  {sector_dir / 'vetted_planet_candidates.parquet'}   ({len(planet):,} planet_candidate rows — alias for downstream)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
