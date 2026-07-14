#!/usr/bin/env python3
"""Finalize fresh A2v1 BLS and Teacher-v1 recovery tables and figures."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import pandas as pd

from plot_s56_duration_aware_recovery import plot_publication_period_radius_recovery_map
from twirl.vetting.injection_teacher_recovery import (
    aggregate_teacher_injection_recovery,
)


SMALL_APERTURE = "DET_FLUX_ADP_SML"


def _bool(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values.fillna(False)
    return (
        values.fillna("").astype(str).str.lower().isin({"1", "true", "t", "yes", "y"})
    )


def _read_many(paths: list[Path]) -> pd.DataFrame:
    frames = [
        pd.read_parquet(path)
        if path.suffix.lower() == ".parquet"
        else pd.read_csv(path, low_memory=False)
        for path in sorted(paths)
    ]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _topk_table(manifest: pd.DataFrame, peaks: pd.DataFrame) -> pd.DataFrame:
    small = peaks.loc[peaks["aperture"].astype(str).eq(SMALL_APERTURE)].copy()
    small["peak_rank"] = pd.to_numeric(small["peak_rank"], errors="coerce")
    small["truth_match"] = _bool(small["is_injected_signal_peak"])
    out = manifest[["injection_id"]].copy()
    for k in (1, 3, 5):
        matched = (
            small.loc[small["peak_rank"].between(1, k)]
            .groupby("injection_id")["truth_match"]
            .any()
        )
        out[f"bls_top{k}_recovered"] = (
            out["injection_id"].map(matched).fillna(False).astype(bool)
        )
    return out


def _tmag_bin(values: pd.Series) -> pd.Series:
    return pd.cut(
        pd.to_numeric(values, errors="coerce"),
        [-np.inf, 17.0, 18.0, 19.0, np.inf],
        right=False,
        labels=["Tmag < 17", "17 <= Tmag < 18", "18 <= Tmag < 19", "Tmag >= 19"],
    )


def _raw_grid(outcomes: pd.DataFrame, outcome_columns: tuple[str, ...]) -> pd.DataFrame:
    work = outcomes.copy()
    work["tmag_bin"] = _tmag_bin(work["tmag"])
    aggregations: dict[str, tuple[str, str]] = {
        "n_injections": ("injection_id", "size")
    }
    for column in outcome_columns:
        aggregations[f"n_{column}"] = (column, "sum")
        aggregations[f"fraction_{column}"] = (column, "mean")
    return (
        work.groupby(
            ["tmag_bin", "grid_period_bin", "grid_radius_bin"],
            observed=False,
            dropna=False,
        )
        .agg(**aggregations)
        .reset_index()
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifests", type=Path, nargs="+", required=True)
    parser.add_argument("--peak-tables", type=Path, nargs="+", required=True)
    parser.add_argument("--teacher-scores", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--expect-injections", type=int, default=20_000)
    parser.add_argument("--preserve-threshold", type=float, default=0.5)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    manifest = _read_many(args.manifests)
    peaks = _read_many(args.peak_tables)
    scores = pd.read_parquet(args.teacher_scores)
    for name, frame in (("manifest", manifest), ("peaks", peaks), ("scores", scores)):
        if "injection_id" not in frame:
            raise KeyError(f"{name} has no injection_id column")
    if (
        len(manifest) != args.expect_injections
        or manifest["injection_id"].nunique() != args.expect_injections
    ):
        raise ValueError(
            "merged injection manifests do not contain exactly the expected unique rows"
        )
    if manifest["tic"].nunique() != args.expect_injections:
        raise ValueError("fresh injection manifests reuse host TICs")
    sectors = (
        pd.to_numeric(manifest["sector"], errors="coerce").dropna().astype(int).unique()
    )
    if len(sectors) != 1:
        raise ValueError("fresh injection manifests must contain exactly one sector")
    if scores["review_id"].duplicated().any():
        raise ValueError("Teacher score table has duplicate review_id values")

    peaks.to_parquet(
        args.out_dir / "adp_bls_top10_peaks.parquet", compression="zstd", index=False
    )
    topk = _topk_table(manifest, peaks)
    base = manifest.rename(
        columns={
            "period_d": "truth_period_d",
            "radius_reearth": "truth_radius_reearth",
            "duration_min": "truth_duration_min",
            "n_good_in_transit": "truth_n_good_in_transit",
            "tessmag": "tmag",
        }
    ).copy()
    base["plot_radius_rearth"] = pd.to_numeric(
        base["truth_radius_reearth"], errors="coerce"
    )
    base = base.merge(topk, on="injection_id", how="left", validate="one_to_one")
    outcomes, teacher_summary = aggregate_teacher_injection_recovery(
        base,
        scores,
        preserve_threshold=args.preserve_threshold,
        max_peak_rank=5,
    )
    if not outcomes["bls_adp_only_recovered"].equals(outcomes["bls_top5_recovered"]):
        mismatch = int(
            (outcomes["bls_adp_only_recovered"] != outcomes["bls_top5_recovered"]).sum()
        )
        raise RuntimeError(
            f"BLS and Teacher candidate sets disagree for {mismatch} injections"
        )
    outcomes.to_parquet(
        args.out_dir / "injection_recovery_outcomes.parquet",
        compression="zstd",
        index=False,
    )
    topk.to_csv(args.out_dir / "bls_recovery_at_1_3_5.csv", index=False)

    strict_planet = scores.loc[
        _bool(scores["is_injected_signal_peak"])
        & scores["predicted_morphology"].astype(str).eq("planet_like")
    ].copy()
    strict_planet.to_parquet(
        args.out_dir / "teacher_v1_strict_planet_truth_matches.parquet",
        compression="zstd",
        index=False,
    )
    raw_grid = _raw_grid(
        outcomes,
        (
            "bls_top1_recovered",
            "bls_top3_recovered",
            "bls_top5_recovered",
            "teacher_v1_preserve_recovered",
            "teacher_v1_planet_recovered",
        ),
    )
    raw_grid.to_csv(args.out_dir / "period_radius_recovery_unsmoothed.csv", index=False)
    conditional = (
        outcomes.assign(tmag_bin=_tmag_bin(outcomes["tmag"]))
        .groupby("tmag_bin", observed=False)
        .agg(
            n_injections=("injection_id", "size"),
            n_bls_top5=("bls_top5_recovered", "sum"),
            n_teacher_preserve=("teacher_v1_preserve_recovered", "sum"),
            bls_top5_fraction=("bls_top5_recovered", "mean"),
            teacher_end_to_end_fraction=("teacher_v1_preserve_recovered", "mean"),
        )
        .reset_index()
    )
    conditional["teacher_retention_given_bls"] = conditional[
        "n_teacher_preserve"
    ] / conditional["n_bls_top5"].replace(0, np.nan)
    conditional.to_csv(args.out_dir / "teacher_retention_by_tmag.csv", index=False)

    bls_figures = plot_publication_period_radius_recovery_map(
        outcomes,
        args.out_dir / "bls_top5_recovery_map",
        recovered_col="bls_top5_recovered",
        colorbar_label="Kernel-smoothed ADP BLS top-5 recovery fraction",
    )
    teacher_figures = plot_publication_period_radius_recovery_map(
        outcomes,
        args.out_dir / "teacher_v1_preserve_recovery_map",
        recovered_col="teacher_v1_preserve_recovered",
        colorbar_label="Kernel-smoothed BLS + Teacher-v1 recovery fraction",
    )
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sector": int(sectors[0]),
        "n_injections": int(len(outcomes)),
        "n_unique_tics": int(outcomes["tic"].nunique()),
        "n_peak_rows": int(len(peaks)),
        "n_teacher_candidate_rows": int(len(scores)),
        "bls_recovery": {
            f"top_{k}": {
                "n": int(outcomes[f"bls_top{k}_recovered"].sum()),
                "fraction": float(outcomes[f"bls_top{k}_recovered"].mean()),
            }
            for k in (1, 3, 5)
        },
        "teacher": teacher_summary,
        "n_strict_planet_truth_match_rows": int(len(strict_planet)),
        "figures": {"bls": bls_figures, "teacher": teacher_figures},
        "candidate_set_identity_verified": True,
    }
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
