#!/usr/bin/env python3
"""Build a source-blinded S56 teacher queue with recovery-boundary injections."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
SRC_ROOT = REPO_ROOT / "src"
for path in (str(SCRIPT_DIR), str(SRC_ROOT)):
    if path not in sys.path:
        sys.path.insert(0, path)

from build_s56_mixed_teacher_queue import (  # noqa: E402
    DEFAULT_REAL_CANDIDATES,
    LABEL_HEADER,
    _blind_for_browser,
    _finalize_queue,
    _sample_balanced_cells,
    _value_counts,
    select_real_pool,
)


DEFAULT_INJECTION_RECOVERY_CSV = (
    REPO_ROOT
    / "reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo"
    / "small_pair_200k/injection_bls_recoveries.csv"
)
DEFAULT_OUT_DIR = REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue"
DEFAULT_TWIRL_BRANCH = "current_adp"

PERIOD_EDGES = [0.08, 0.15, 0.25, 0.40, 0.70, 1.20, 2.00, 4.00, 8.00, 13.10]
RADIUS_EDGES = [0.18, 0.35, 0.60, 1.00, 1.70, 3.00, 5.00, 8.50, 13.00, 18.10]
TMAG_EDGES = [0.0, 17.0, 18.0, 19.0, 21.0]


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=_json_default) + "\n")


def _coerce_numeric(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = frame.copy()
    for column in columns:
        if column in out:
            out[column] = pd.to_numeric(out[column], errors="coerce")
    return out


def _safe_sheet_name(review_id: Any, branch_name: str = DEFAULT_TWIRL_BRANCH) -> str:
    safe = str(review_id).replace("/", "_").replace(":", "_")
    return f"{safe}_twirl_twoap_{branch_name}.png"


def _add_recovery_cells(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    work = frame.copy()
    work["strict_top1_recovered"] = work["recovery_status"].fillna("").astype(str).eq("bls_recovered")
    work["any_exact_or_harmonic_recovered"] = (
        work["strict_top1_recovered"]
        | work["topn_recovery_status"]
        .fillna("")
        .astype(str)
        .isin({"bls_top1_recovered", "bls_topn_recovered", "bls_topn_harmonic_match"})
    )
    work["boundary_period_bin"] = pd.cut(work["truth_period_d"], PERIOD_EDGES, include_lowest=True)
    work["boundary_radius_bin"] = pd.cut(work["truth_radius_rearth"], RADIUS_EDGES, include_lowest=True)
    work["boundary_tmag_bin"] = pd.cut(work["tmag"], TMAG_EDGES, include_lowest=True)
    cell_cols = ["boundary_period_bin", "boundary_radius_bin", "boundary_tmag_bin"]
    cells = (
        work.groupby(cell_cols, observed=True)
        .agg(
            recovery50_cell_n=("injection_id", "size"),
            recovery50_cell_strict_n=("strict_top1_recovered", "sum"),
            recovery50_cell_any_n=("any_exact_or_harmonic_recovered", "sum"),
        )
        .reset_index()
    )
    cells["recovery50_cell_strict_frac"] = (
        cells["recovery50_cell_strict_n"] / cells["recovery50_cell_n"]
    )
    cells["recovery50_cell_any_frac"] = cells["recovery50_cell_any_n"] / cells["recovery50_cell_n"]
    cells["recovery50_cell_distance_to_0p5"] = (cells["recovery50_cell_any_frac"] - 0.5).abs()

    for column in cell_cols:
        work[f"{column}_key"] = work[column].astype(str)
        cells[f"{column}_key"] = cells[column].astype(str)
    joined = work.merge(
        cells[
            [
                *(f"{column}_key" for column in cell_cols),
                "recovery50_cell_n",
                "recovery50_cell_strict_frac",
                "recovery50_cell_any_frac",
                "recovery50_cell_distance_to_0p5",
            ]
        ],
        on=[f"{column}_key" for column in cell_cols],
        how="left",
    )
    return joined, cells


def select_recovery50_injections(
    path: Path,
    *,
    n_injected: int,
    random_state: int,
    min_cell_recovery: float,
    max_cell_recovery: float,
    min_good_in_transit: int,
    require_both_apertures: bool,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    df = pd.read_csv(path)
    df = _coerce_numeric(
        df,
        [
            "period_d",
            "t0_bjd",
            "duration_min",
            "depth",
            "depth_snr",
            "sde_max",
            "tmag",
            "truth_period_d",
            "truth_t0_bjd",
            "truth_duration_min",
            "truth_depth",
            "truth_sampled_model_depth",
            "truth_radius_rearth",
            "truth_impact_b",
            "truth_inclination_deg",
            "truth_n_good_in_transit",
            "sde_DET_FLUX_ADP_SML",
            "sde_DET_FLUX_SML",
            "period_rel_err_DET_FLUX_ADP_SML",
            "period_rel_err_DET_FLUX_SML",
        ],
    )
    df, cells = _add_recovery_cells(df)
    df["both_apertures_top1_recovered"] = (
        df.get("recovery_status_DET_FLUX_ADP_SML", pd.Series("", index=df.index))
        .fillna("")
        .astype(str)
        .eq("bls_recovered")
        & df.get("recovery_status_DET_FLUX_SML", pd.Series("", index=df.index))
        .fillna("")
        .astype(str)
        .eq("bls_recovered")
    )
    df["min_twoap_sde"] = df[["sde_DET_FLUX_ADP_SML", "sde_DET_FLUX_SML"]].min(axis=1)
    base = (
        df["strict_top1_recovered"]
        & df["period_d"].notna()
        & df["t0_bjd"].notna()
        & df["duration_min"].notna()
        & df["truth_period_d"].notna()
        & df["truth_radius_rearth"].notna()
        & df["truth_sampled_model_depth"].notna()
        & df["truth_impact_b"].notna()
        & df["truth_inclination_deg"].notna()
        & df["truth_n_good_in_transit"].fillna(0).ge(min_good_in_transit)
        & df["h5_group"].fillna("").astype(str).ne("")
        & df["source_h5"].fillna("").astype(str).ne("")
        & df["recovery50_cell_any_frac"].between(min_cell_recovery, max_cell_recovery)
    )
    if require_both_apertures:
        base &= df["both_apertures_top1_recovered"]
    candidates = df.loc[base].copy()
    if len(candidates) < n_injected:
        raise ValueError(
            f"only {len(candidates):,} recovery-boundary injected candidates; need {n_injected:,}"
        )

    candidates["selection_bucket"] = "inj_recovery50_strict_top1"
    candidates["selection_note"] = (
        "strict_top1_recovered; coarse period-radius-tmag cell any/topN/harmonic recovery "
        f"in [{min_cell_recovery:.2f},{max_cell_recovery:.2f}]"
    )
    selected = _sample_balanced_cells(
        candidates,
        n=n_injected,
        random_state=random_state,
        cell_columns=("boundary_period_bin_key", "boundary_radius_bin_key", "boundary_tmag_bin_key"),
    )
    selected = selected.sort_values(
        ["recovery50_cell_distance_to_0p5", "min_twoap_sde"],
        ascending=[True, False],
        kind="stable",
    ).head(n_injected)
    selected = selected.sample(frac=1.0, random_state=random_state + 1000).reset_index(drop=True)

    if "truth_source_kind" not in selected:
        selected["truth_source_kind"] = "injection_recovery"
    if "truth_source_bucket" not in selected:
        selected["truth_source_bucket"] = selected.get("source_bucket", "")
    selected["source_kind"] = "injection_recovery"
    selected["truth_source_kind"] = "injection_recovery"
    selected["source_bucket"] = selected["selection_bucket"]
    selected["review_id"] = selected.get("review_id", selected["injection_id"].map(lambda x: f"inj:{x}"))
    selected["twirl_vet_sheet_name"] = selected["review_id"].map(_safe_sheet_name)
    selected["twirl_vet_sheet_pdf_name"] = selected["twirl_vet_sheet_name"].str.replace(".png", ".pdf", regex=False)

    selected_cell_counts = (
        selected.groupby(
            ["boundary_period_bin_key", "boundary_radius_bin_key", "boundary_tmag_bin_key"],
            dropna=False,
        )
        .size()
        .reset_index(name="selected_cell_n")
    )
    available_cell_counts = (
        candidates.groupby(
            ["boundary_period_bin_key", "boundary_radius_bin_key", "boundary_tmag_bin_key"],
            dropna=False,
        )
        .size()
        .reset_index(name="available_cell_n")
    )
    weights = selected_cell_counts.merge(
        available_cell_counts,
        on=["boundary_period_bin_key", "boundary_radius_bin_key", "boundary_tmag_bin_key"],
        how="left",
    )
    weights["selection_weight"] = weights["available_cell_n"] / weights["selected_cell_n"].clip(lower=1)
    selected = selected.merge(
        weights[
            [
                "boundary_period_bin_key",
                "boundary_radius_bin_key",
                "boundary_tmag_bin_key",
                "selection_weight",
            ]
        ],
        on=["boundary_period_bin_key", "boundary_radius_bin_key", "boundary_tmag_bin_key"],
        how="left",
    )

    summary = {
        "input_rows": int(len(df)),
        "coarse_recovery_cells": int(len(cells)),
        "candidate_rows_after_filters": int(len(candidates)),
        "selected_rows": int(len(selected)),
        "target_cell_any_recovery_range": [float(min_cell_recovery), float(max_cell_recovery)],
        "require_both_apertures": bool(require_both_apertures),
        "min_good_in_transit": int(min_good_in_transit),
        "selected_cell_any_recovery_quantiles": {
            str(q): float(selected["recovery50_cell_any_frac"].quantile(q))
            for q in (0.0, 0.1, 0.5, 0.9, 1.0)
        },
        "selected_strict_cell_recovery_quantiles": {
            str(q): float(selected["recovery50_cell_strict_frac"].quantile(q))
            for q in (0.0, 0.1, 0.5, 0.9, 1.0)
        },
        "selected_truth_period_d_quantiles": {
            str(q): float(selected["truth_period_d"].quantile(q))
            for q in (0.0, 0.1, 0.5, 0.9, 1.0)
        },
        "selected_truth_radius_rearth_quantiles": {
            str(q): float(selected["truth_radius_rearth"].quantile(q))
            for q in (0.0, 0.1, 0.5, 0.9, 1.0)
        },
        "selected_depth_pct_quantiles": {
            str(q): float((100.0 * selected["truth_sampled_model_depth"]).quantile(q))
            for q in (0.0, 0.1, 0.5, 0.9, 1.0)
        },
        "selected_tmag_quantiles": {
            str(q): float(selected["tmag"].quantile(q))
            for q in (0.0, 0.1, 0.5, 0.9, 1.0)
        },
        "topn_recovery_counts": _value_counts(selected, "topn_recovery_status"),
        "recovery_status_counts": _value_counts(selected, "recovery_status"),
    }
    return selected, summary


def downsample_real_teacher_rows(
    real_pool: pd.DataFrame,
    *,
    n_real: int,
    random_state: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Downsample the broad real pool into a compact teacher-review mix."""

    quotas = [
        ("real_wd1856_benchmark", 1),
        ("real_high_sde_planet_like", 260),
        ("real_eb_pceb_like", 160),
        ("real_variability_broad_duration", 95),
        ("real_aperture_disagreement", 70),
        ("real_cadence_alias_systematic", 45),
        ("real_mid_sde_control", 45),
        ("real_low_sde_control", 24),
    ]
    available_counts = _value_counts(real_pool, "selection_bucket")
    pieces: list[pd.DataFrame] = []
    used_review_ids: set[str] = set()
    for idx, (bucket, quota) in enumerate(quotas):
        if sum(len(piece) for piece in pieces) >= n_real:
            break
        sub = real_pool.loc[real_pool["selection_bucket"].fillna("").astype(str).eq(bucket)].copy()
        sub = sub.loc[~sub["review_id"].fillna("").astype(str).isin(used_review_ids)]
        take = min(quota, n_real - sum(len(piece) for piece in pieces), len(sub))
        if take <= 0:
            continue
        if bucket in {"real_wd1856_benchmark", "real_high_sde_planet_like", "real_eb_pceb_like"}:
            sort_cols = [col for col in ("class_rank", "sde_max") if col in sub]
            if sort_cols:
                ascending = [True if col == "class_rank" else False for col in sort_cols]
                selected = sub.sort_values(sort_cols, ascending=ascending, kind="stable").head(take)
            else:
                selected = sub.head(take)
        else:
            selected = sub.sample(n=take, random_state=random_state + idx)
        pieces.append(selected)
        used_review_ids.update(selected["review_id"].fillna("").astype(str).tolist())

    selected = pd.concat(pieces, ignore_index=False) if pieces else real_pool.head(0)
    if len(selected) < n_real:
        remaining = real_pool.loc[
            ~real_pool["review_id"].fillna("").astype(str).isin(used_review_ids)
        ].copy()
        fill_n = min(n_real - len(selected), len(remaining))
        if fill_n > 0:
            selected = pd.concat(
                [selected, remaining.sample(n=fill_n, random_state=random_state + 10_000)],
                ignore_index=False,
            )
    if len(selected) != n_real:
        raise ValueError(f"selected {len(selected):,} real rows; expected {n_real:,}")
    selected = selected.sample(frac=1.0, random_state=random_state + 20_000).reset_index(drop=True)
    summary = {
        "broad_pool_rows": int(len(real_pool)),
        "selected_rows": int(len(selected)),
        "broad_pool_bucket_counts": available_counts,
        "selected_bucket_counts": _value_counts(selected, "selection_bucket"),
        "selection_rule": "fixed teacher-review quotas from broad stratified real pool",
    }
    return selected, summary


def _verify_queue(queue: pd.DataFrame, *, n_real: int, n_injected: int) -> dict[str, Any]:
    failures: list[str] = []
    source = queue["source_kind"].fillna("").astype(str)
    if len(queue) != n_real + n_injected:
        failures.append(f"queue rows {len(queue)} != expected {n_real + n_injected}")
    if int(source.eq("real_candidate").sum()) != n_real:
        failures.append("real row count mismatch")
    if int(source.eq("injection_recovery").sum()) != n_injected:
        failures.append("injected row count mismatch")
    if queue["review_id"].fillna("").astype(str).duplicated().any():
        failures.append("duplicate review_id")
    for column in ("period_d", "t0_bjd", "duration_min"):
        bad = int(pd.to_numeric(queue[column], errors="coerce").isna().sum())
        if bad:
            failures.append(f"{bad} non-finite {column} rows")
    visible_source = queue.get("source_bucket", pd.Series("", index=queue.index)).fillna("").astype(str)
    if visible_source.ne("review_candidate").any():
        failures.append("visible source_bucket is not fully blinded")
    visible_class = queue.get("vet_class", pd.Series("", index=queue.index)).fillna("").astype(str)
    if visible_class.ne("review_candidate").any():
        failures.append("visible vet_class is not fully blinded")
    return {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "passed": not failures,
        "failures": failures,
        "n_rows": int(len(queue)),
        "source_kind_counts": _value_counts(queue, "source_kind"),
        "truth_source_kind_counts": _value_counts(queue, "truth_source_kind"),
        "selection_bucket_counts": _value_counts(queue, "selection_bucket"),
    }


def build_queue(
    *,
    real_candidates: Path,
    injection_recovery_csv: Path,
    out_dir: Path,
    n_real: int,
    n_real_source_pool: int,
    n_injected: int,
    random_state: int,
    min_cell_recovery: float,
    max_cell_recovery: float,
    min_good_in_transit: int,
    require_both_apertures: bool,
    cadence_alias_tolerance: float,
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    real_pool, real_pool_summary = select_real_pool(
        real_candidates,
        n_real=n_real_source_pool,
        random_state=random_state,
        alias_tolerance=cadence_alias_tolerance,
    )
    real, real_summary = downsample_real_teacher_rows(
        real_pool,
        n_real=n_real,
        random_state=random_state + 500,
    )
    real_summary["source_pool_summary"] = real_pool_summary
    injected, injected_summary = select_recovery50_injections(
        injection_recovery_csv,
        n_injected=n_injected,
        random_state=random_state + 1000,
        min_cell_recovery=min_cell_recovery,
        max_cell_recovery=max_cell_recovery,
        min_good_in_transit=min_good_in_transit,
        require_both_apertures=require_both_apertures,
    )
    pool = _finalize_queue(real, injected)
    pool["teacher_pool_version"] = "s56_recovery50_teacher_1k_v1"
    pool["pool_random_state"] = int(random_state)
    pool["truth_vet_class"] = pool.get("vet_class", "")
    pool["twirl_vet_sheet_name"] = pool["review_id"].map(_safe_sheet_name)
    pool["twirl_vet_sheet_pdf_name"] = pool["twirl_vet_sheet_name"].str.replace(".png", ".pdf", regex=False)

    unblinded = pool.sample(frac=1.0, random_state=random_state + 2000).reset_index(drop=True)
    review = _blind_for_browser(unblinded)
    review["review_order_seed"] = int(random_state)
    review_csv = out_dir / "review_queue_1k.csv"
    review.to_csv(review_csv, index=False)
    unblinded.to_csv(out_dir / "review_queue_1k_unblinded.csv", index=False)
    review.loc[review["source_kind"].fillna("").astype(str).eq("real_candidate")].to_csv(
        out_dir / "review_queue_1k_real.csv",
        index=False,
    )
    review.loc[review["source_kind"].fillna("").astype(str).eq("injection_recovery")].to_csv(
        out_dir / "review_queue_1k_injected.csv",
        index=False,
    )
    label_path = out_dir / "human_labels_vetted.csv"
    if not label_path.exists():
        pd.DataFrame(columns=LABEL_HEADER).to_csv(label_path, index=False)

    verification = _verify_queue(review, n_real=n_real, n_injected=n_injected)
    _write_json(out_dir / "verification.json", verification)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "s56_recovery50_teacher_queue",
        "random_state": int(random_state),
        "n_review": int(len(review)),
        "n_real": int(n_real),
        "n_real_source_pool": int(n_real_source_pool),
        "n_injected": int(n_injected),
        "injected_fraction": float(n_injected / max(n_real + n_injected, 1)),
        "real_candidates": str(real_candidates),
        "injection_recovery_csv": str(injection_recovery_csv),
        "source_blinding": {
            "browser_visible_source_bucket": "review_candidate",
            "browser_visible_vet_class": "review_candidate",
            "hidden_truth_columns_retained": True,
        },
        "injected_selection": injected_summary,
        "real_selection": real_summary,
        "outputs": {
            "review_queue_1k_csv": str(review_csv),
            "review_queue_1k_real_csv": str(out_dir / "review_queue_1k_real.csv"),
            "review_queue_1k_injected_csv": str(out_dir / "review_queue_1k_injected.csv"),
            "human_labels_vetted_csv": str(label_path),
            "twirl_vet_sheets": str(out_dir / "twirl_vet_sheets"),
        },
        "verification_passed": bool(verification["passed"]),
        "verification_failures": verification["failures"],
    }
    _write_json(out_dir / "summary.json", summary)
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--real-candidates", type=Path, default=DEFAULT_REAL_CANDIDATES)
    parser.add_argument("--injection-recovery-csv", type=Path, default=DEFAULT_INJECTION_RECOVERY_CSV)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--n-real", type=int, default=700)
    parser.add_argument("--n-real-source-pool", type=int, default=9000)
    parser.add_argument("--n-injected", type=int, default=300)
    parser.add_argument("--random-state", type=int, default=561050)
    parser.add_argument("--min-cell-recovery", type=float, default=0.35)
    parser.add_argument("--max-cell-recovery", type=float, default=0.65)
    parser.add_argument("--min-good-in-transit", type=int, default=5)
    parser.add_argument("--allow-single-aperture-recovery", action="store_true")
    parser.add_argument("--cadence-alias-tolerance", type=float, default=0.02)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    summary = build_queue(
        real_candidates=args.real_candidates,
        injection_recovery_csv=args.injection_recovery_csv,
        out_dir=args.out_dir,
        n_real=int(args.n_real),
        n_real_source_pool=int(args.n_real_source_pool),
        n_injected=int(args.n_injected),
        random_state=int(args.random_state),
        min_cell_recovery=float(args.min_cell_recovery),
        max_cell_recovery=float(args.max_cell_recovery),
        min_good_in_transit=int(args.min_good_in_transit),
        require_both_apertures=not bool(args.allow_single_aperture_recovery),
        cadence_alias_tolerance=float(args.cadence_alias_tolerance),
    )
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return 0 if summary["verification_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
