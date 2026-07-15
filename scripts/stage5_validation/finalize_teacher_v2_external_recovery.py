#!/usr/bin/env python3
"""Finalize locked/external BLS and Teacher-v2 recovery tables and figures."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import pandas as pd

from plot_s56_duration_aware_recovery import plot_publication_period_radius_recovery_map
from twirl.injections.a2v1_recovery import normalize_fresh_injection_manifest_truth
from twirl.plotting.style import apply_twirl_style
from twirl.vetting.teacher_v2_recovery import (
    aggregate_compact_recovery,
    bls_topk_recovery_table,
    compare_recovery_models,
    period_radius_tmag_support,
    select_role_scoped_injection_truth,
)


def _read(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path) if path.suffix.lower() == ".parquet" else pd.read_csv(path, low_memory=False)


def _bool(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False).astype(bool)
    return values.fillna("").astype(str).str.lower().isin({"1", "1.0", "true", "t", "yes", "y"})


def _tmag_bin(values: pd.Series) -> pd.Series:
    return pd.cut(
        pd.to_numeric(values, errors="coerce"),
        [-np.inf, 17.0, 18.0, 19.0, np.inf],
        right=False,
        labels=["Tmag < 17", "17 <= Tmag < 18", "18 <= Tmag < 19", "Tmag >= 19"],
    )


def _unsmoothed_grid(outcomes: pd.DataFrame, columns: tuple[str, ...]) -> pd.DataFrame:
    work = outcomes.copy()
    work["tmag_bin"] = _tmag_bin(work["tmag"])
    aggregations: dict[str, tuple[str, str]] = {"n_injections": ("injection_id", "size")}
    for column in columns:
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


def _plot_difference_map(
    *,
    v2_grid_path: Path,
    v1_grid_path: Path,
    out_dir: Path,
) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    v2 = pd.read_csv(v2_grid_path)
    v1 = pd.read_csv(v1_grid_path)
    keys = ["tmag_bin", "period_d", "radius_rearth"]
    merged = v2.merge(v1, on=keys, suffixes=("_v2", "_v1"), validate="one_to_one")
    merged["recovery_difference_v2_minus_v1"] = (
        merged["kernel_recovery_fraction_v2"] - merged["kernel_recovery_fraction_v1"]
    )
    merged["masked_for_low_support"] = _bool(merged["masked_for_low_support_v2"]) | _bool(
        merged["masked_for_low_support_v1"]
    )
    apply_twirl_style("full_page")
    labels = ["Tmag < 17", "17 <= Tmag < 18", "18 <= Tmag < 19", "Tmag >= 19"]
    finite = merged.loc[
        ~merged["masked_for_low_support"], "recovery_difference_v2_minus_v1"
    ].to_numpy(dtype=float)
    limit = max(0.1, min(0.5, float(np.nanpercentile(np.abs(finite), 98)))) if len(finite) else 0.1
    fig, axes = plt.subplots(2, 2, figsize=(9.9, 8.1), sharex=True, sharey=True)
    mesh = None
    for ax, label in zip(axes.ravel(), labels):
        part = merged.loc[merged["tmag_bin"].astype(str).eq(label)].copy()
        periods = np.sort(part["period_d"].unique())
        radii = np.sort(part["radius_rearth"].unique())
        pivot = part.pivot(index="period_d", columns="radius_rearth", values="recovery_difference_v2_minus_v1").reindex(
            index=periods, columns=radii
        )
        mask = part.pivot(index="period_d", columns="radius_rearth", values="masked_for_low_support").reindex(
            index=periods, columns=radii
        )
        values = pivot.to_numpy(dtype=float, copy=True)
        values[mask.to_numpy(dtype=bool)] = np.nan
        mesh = ax.pcolormesh(
            periods,
            radii,
            values.T,
            shading="nearest",
            cmap="RdBu_r",
            vmin=-limit,
            vmax=limit,
            rasterized=True,
        )
        ax.axhline(1.0, color="0.25", linewidth=0.5, alpha=0.35)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(0.12, 13.0)
        ax.set_ylim(0.18, 18.0)
        ax.set_title(label, fontsize=10.8)
        ax.grid(False)
        ax.set_xlabel(r"Injected orbital period, $P$ [days]")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"Injected companion radius, $R_p$ [$R_\oplus$]")
    if mesh is not None:
        cax = fig.add_axes([0.885, 0.10, 0.022, 0.82])
        colorbar = fig.colorbar(mesh, cax=cax)
        colorbar.set_label("Teacher v2 - Teacher v1 recovery fraction")
        colorbar.ax.tick_params(direction="in", which="both", labelsize=8.5)
    fig.subplots_adjust(left=0.09, right=0.855, top=0.92, bottom=0.10, wspace=0.08, hspace=0.20)
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / "teacher_v2_minus_v1_recovery.png"
    pdf = out_dir / "teacher_v2_minus_v1_recovery.pdf"
    csv = out_dir / "teacher_v2_minus_v1_recovery_grid.csv"
    merged.to_csv(csv, index=False)
    fig.savefig(png, dpi=260, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return {"png": str(png), "pdf": str(pdf), "grid_csv": str(csv)}


def _verify_candidate_identity(candidates: pd.DataFrame, scores: pd.DataFrame, *, name: str) -> None:
    if candidates["review_id"].duplicated().any() or scores["review_id"].duplicated().any():
        raise ValueError(f"{name} candidates or scores contain duplicate review IDs")
    if set(candidates["review_id"].astype(str)) != set(scores["review_id"].astype(str)):
        raise ValueError(f"{name} scores do not cover the identical stored candidate set")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        required=True,
        help="Derived role manifest defining the locked/external injection IDs.",
    )
    parser.add_argument(
        "--truth-manifests",
        type=Path,
        nargs="+",
        required=True,
        help="Immutable post-injection shard manifests containing physical truth.",
    )
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--teacher-v2-scores", type=Path, required=True)
    parser.add_argument("--frozen-selection", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--expect-min-injections", type=int, default=10_000)
    parser.add_argument("--teacher-v1-scores", type=Path)
    parser.add_argument("--teacher-v1-threshold", type=float)
    parser.add_argument("--metadata-scores", type=Path)
    parser.add_argument("--metadata-threshold", type=float)
    args = parser.parse_args()

    frozen = json.loads(args.frozen_selection.read_text())
    if not frozen.get("architecture_frozen") or not frozen.get("threshold_frozen"):
        raise RuntimeError("Teacher-v2 architecture and threshold must be frozen before external evaluation")
    threshold = float(frozen["frozen_compact_threshold"])
    role_manifest = _read(args.manifest)
    missing_truth_manifests = [path for path in args.truth_manifests if not path.exists()]
    if missing_truth_manifests:
        raise FileNotFoundError(
            f"missing immutable injection truth manifests: {missing_truth_manifests}"
        )
    truth_manifest = pd.concat(
        [_read(path) for path in args.truth_manifests],
        ignore_index=True,
        sort=False,
    )
    manifest = select_role_scoped_injection_truth(role_manifest, truth_manifest)
    candidates = _read(args.candidates)
    allowed_ids = frozenset(manifest["injection_id"].astype(str))
    candidates = candidates.loc[candidates["injection_id"].astype(str).isin(allowed_ids)].copy()
    v2_scores = _read(args.teacher_v2_scores)
    _verify_candidate_identity(candidates, v2_scores, name="Teacher-v2")
    if len(manifest) < int(args.expect_min_injections):
        raise ValueError("external injection support is below the accepted minimum")
    if manifest["injection_id"].nunique() != len(manifest) or manifest["tic"].nunique() != len(manifest):
        raise ValueError("external recovery manifest must use unique injection IDs and host TICs")
    if manifest["grid_cell_id"].nunique() != 2500:
        raise ValueError("external recovery manifest does not represent every period-radius cell")
    sectors = pd.to_numeric(manifest["sector"], errors="coerce").dropna().astype(int).unique()
    if len(sectors) != 1:
        raise ValueError("external recovery manifest must contain exactly one sector")

    base = normalize_fresh_injection_manifest_truth(manifest)
    support = period_radius_tmag_support(base)
    topk = bls_topk_recovery_table(base, candidates)
    v2_outcomes, v2_summary = aggregate_compact_recovery(
        base,
        v2_scores,
        threshold=threshold,
        outcome_prefix="teacher_v2",
    )
    truth = candidates.copy()
    truth["_truth_match"] = _bool(truth["is_injected_signal_peak"])
    rank = pd.to_numeric(truth["rep_peak_rank"], errors="coerce")
    for k in (1, 3, 5):
        recovered_ids = set(
            truth.loc[truth["_truth_match"] & rank.le(k), "injection_id"].astype(str)
        )
        v2_outcomes[f"bls_top{k}_recovered"] = v2_outcomes["injection_id"].astype(str).isin(
            recovered_ids
        )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    v2_outcomes.to_parquet(
        args.out_dir / "injection_recovery_outcomes.parquet", compression="zstd", index=False
    )
    topk.to_csv(args.out_dir / "bls_recovery_at_1_3_5.csv", index=False)
    support.to_csv(args.out_dir / "period_radius_tmag_support.csv", index=False)
    _unsmoothed_grid(
        v2_outcomes,
        ("bls_top1_recovered", "bls_top3_recovered", "bls_top5_recovered", "teacher_v2_compact_recovered"),
    ).to_csv(args.out_dir / "period_radius_recovery_unsmoothed.csv", index=False)
    conditional = (
        v2_outcomes.assign(tmag_bin=_tmag_bin(v2_outcomes["tmag"]))
        .groupby("tmag_bin", observed=False)
        .agg(
            n_injections=("injection_id", "size"),
            n_bls_top5=("bls_top5_recovered", "sum"),
            n_teacher_v2=("teacher_v2_compact_recovered", "sum"),
            bls_top5_fraction=("bls_top5_recovered", "mean"),
            teacher_v2_end_to_end_fraction=("teacher_v2_compact_recovered", "mean"),
        )
        .reset_index()
    )
    conditional["teacher_v2_retention_given_bls"] = conditional["n_teacher_v2"] / conditional[
        "n_bls_top5"
    ].replace(0, np.nan)
    conditional.to_csv(args.out_dir / "teacher_retention_by_tmag.csv", index=False)
    strict_planet = v2_scores.loc[
        _bool(v2_scores["is_injected_signal_peak"])
        & v2_scores["predicted_morphology"].astype(str).eq("planet_like")
    ].copy()
    strict_planet.to_parquet(
        args.out_dir / "teacher_v2_strict_planet_truth_matches.parquet",
        compression="zstd",
        index=False,
    )

    bls_figures = plot_publication_period_radius_recovery_map(
        v2_outcomes,
        args.out_dir / "bls_top5_recovery_map",
        recovered_col="bls_top5_recovered",
        colorbar_label="Kernel-smoothed ADP BLS top-5 recovery fraction",
    )
    v2_figures = plot_publication_period_radius_recovery_map(
        v2_outcomes,
        args.out_dir / "teacher_v2_compact_recovery_map",
        recovered_col="teacher_v2_compact_recovered",
        colorbar_label="Kernel-smoothed BLS + Teacher-v2 recovery fraction",
    )

    model_outcomes = {"teacher_v2": v2_outcomes}
    outcome_columns = {"teacher_v2": "teacher_v2_compact_recovered"}
    baseline_figures: dict[str, object] = {}
    if args.teacher_v1_scores:
        if args.teacher_v1_threshold is None:
            raise ValueError("--teacher-v1-threshold is required with Teacher-v1 scores")
        v1_scores = _read(args.teacher_v1_scores)
        _verify_candidate_identity(candidates, v1_scores, name="Teacher-v1")
        v1_outcomes, _ = aggregate_compact_recovery(
            base,
            v1_scores,
            threshold=args.teacher_v1_threshold,
            score_column="p_preserve",
            outcome_prefix="teacher_v1",
        )
        model_outcomes["teacher_v1"] = v1_outcomes
        outcome_columns["teacher_v1"] = "teacher_v1_compact_recovered"
        baseline_figures["teacher_v1"] = plot_publication_period_radius_recovery_map(
            v1_outcomes,
            args.out_dir / "teacher_v1_workload_matched_recovery_map",
            recovered_col="teacher_v1_compact_recovered",
            colorbar_label="Kernel-smoothed BLS + Teacher-v1 recovery fraction",
        )
        baseline_figures["v2_minus_v1"] = _plot_difference_map(
            v2_grid_path=Path(v2_figures["period_radius_empirical_publication_csv"]),
            v1_grid_path=Path(
                baseline_figures["teacher_v1"]["period_radius_empirical_publication_csv"]
            ),
            out_dir=args.out_dir / "teacher_v2_minus_v1_difference_map",
        )
    if args.metadata_scores:
        if args.metadata_threshold is None:
            raise ValueError("--metadata-threshold is required with metadata scores")
        metadata_scores = _read(args.metadata_scores)
        _verify_candidate_identity(candidates, metadata_scores, name="metadata-only")
        metadata_outcomes, _ = aggregate_compact_recovery(
            base,
            metadata_scores,
            threshold=args.metadata_threshold,
            outcome_prefix="metadata_only",
        )
        model_outcomes["metadata_only"] = metadata_outcomes
        outcome_columns["metadata_only"] = "metadata_only_compact_recovered"
    comparison = compare_recovery_models(model_outcomes, outcome_columns=outcome_columns)
    comparison.to_csv(args.out_dir / "workload_matched_model_comparison.csv", index=False)
    comparison_by_name = comparison.set_index("model")
    v2_conditional = float(comparison_by_name.loc["teacher_v2", "conditional_on_bls"])
    acceptance = {
        "s57_bls_conditional_recovery_at_least_0p80": bool(v2_conditional >= 0.80),
        "improvement_over_teacher_v1_at_least_0p15": (
            bool(
                v2_conditional
                - float(comparison_by_name.loc["teacher_v1", "conditional_on_bls"])
                >= 0.15
            )
            if "teacher_v1" in comparison_by_name.index
            else False
        ),
        "improvement_over_metadata_at_least_0p10": (
            bool(
                v2_conditional
                - float(comparison_by_name.loc["metadata_only", "conditional_on_bls"])
                >= 0.10
            )
            if "metadata_only" in comparison_by_name.index
            else False
        ),
    }
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sector": int(sectors[0]),
        "n_injections": int(len(base)),
        "n_unique_tics": int(base["tic"].nunique()),
        "n_candidates": int(len(candidates)),
        "frozen_selection": frozen,
        "teacher_v2": v2_summary,
        "tmag_support": support.to_dict("records"),
        "n_strict_planet_truth_match_rows": int(len(strict_planet)),
        "model_comparison": comparison.to_dict("records"),
        "acceptance": acceptance,
        "acceptance_passed": all(acceptance.values()),
        "figures": {"bls": bls_figures, "teacher_v2": v2_figures, **baseline_figures},
        "candidate_set_identity_verified": True,
        "role_manifest": str(args.manifest),
        "truth_manifests": [str(path) for path in args.truth_manifests],
    }
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
