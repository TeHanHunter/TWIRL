#!/usr/bin/env python3
"""Compare two injected BLS peak tables branch-by-branch.

This is intended for the next search-method iteration after the standard S56
peak table finishes. The standard branch establishes the baseline recall; an
alternate branch, such as a short-period capped search, should be credited only
for injections it actually rescues without hiding degradations.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from summarize_injection_peak_gate import (  # noqa: E402
    DEFAULT_PERIOD_BINS,
    DEFAULT_RADIUS_BINS,
    DEFAULT_TMAG_BINS,
    _parse_float_tuple,
    _read_table,
    add_recovery_bins,
    build_injection_level_table,
)


DEFAULT_BASELINE_TABLE = (
    Path(__file__).resolve().parents[2]
    / "reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/"
    / "peak_training/s56_20k_injection_bls_peaks_chunked.csv"
)
DEFAULT_OUT_DIR = (
    Path(__file__).resolve().parents[2]
    / "reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/"
    / "peak_training_branch_compare"
)


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    text = series.astype("string").str.strip().str.lower()
    return text.isin({"true", "1", "yes", "y"}).fillna(False).astype(bool)


def _filter_branch(df: pd.DataFrame, branch: str | None) -> pd.DataFrame:
    if branch is None or branch == "":
        return df.copy()
    if "search_branch" not in df.columns:
        raise KeyError(f"cannot filter branch={branch!r}; peak table lacks search_branch")
    return df.loc[df["search_branch"].fillna("").astype(str).eq(branch)].copy()


def _collapse_branch(
    peaks: pd.DataFrame,
    *,
    top_k: int,
    branch_label: str,
    branch_filter: str | None = None,
) -> pd.DataFrame:
    filtered = _filter_branch(peaks, branch_filter)
    if filtered.empty:
        raise ValueError(f"no peak rows remain for branch {branch_label!r}")
    collapsed = build_injection_level_table(filtered, top_k=top_k)
    keep_cols = [
        col
        for col in (
            "injection_id",
            "n_peak_rows",
            "n_candidate_peak_rows",
            "top1_match",
            f"top{top_k}_match",
            f"top{top_k}_exact_match",
            f"top{top_k}_harmonic_match",
            "best_signal_peak_rank",
            "ranking_loss",
            f"not_in_top{top_k}",
        )
        if col in collapsed.columns
    ]
    out = collapsed.loc[:, keep_cols].copy()
    rename = {col: f"{branch_label}_{col}" for col in keep_cols if col != "injection_id"}
    return out.rename(columns=rename)


def _truth_columns(peaks: pd.DataFrame) -> pd.DataFrame:
    truth_cols = [
        col
        for col in (
            "injection_id",
            "tic",
            "sector",
            "cam",
            "ccd",
            "tmag",
            "signal_family",
            "truth_period_d",
            "truth_t0_bjd",
            "truth_duration_min",
            "truth_depth",
            "truth_model_depth",
            "truth_sampled_model_depth",
            "truth_radius_rearth",
            "truth_radius_rwd",
            "truth_impact_b",
            "truth_n_in_transit",
            "truth_n_good_in_transit",
            "truth_grid_cell_id",
            "truth_grid_period_bin",
            "truth_grid_radius_bin",
        )
        if col in peaks.columns
    ]
    if not truth_cols:
        return pd.DataFrame(columns=["injection_id"])
    return peaks.groupby("injection_id", dropna=False)[truth_cols].first().reset_index(drop=True)


def compare_branches(
    baseline_peaks: pd.DataFrame,
    branch_peaks: pd.DataFrame,
    *,
    top_k: int,
    baseline_label: str = "baseline",
    branch_label: str = "branch",
    baseline_filter: str | None = None,
    branch_filter: str | None = None,
    period_bins: tuple[float, ...] = DEFAULT_PERIOD_BINS,
    radius_bins: tuple[float, ...] = DEFAULT_RADIUS_BINS,
    tmag_bins: tuple[float, ...] = DEFAULT_TMAG_BINS,
    comparison_scope: str = "intersection",
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    baseline = _collapse_branch(
        baseline_peaks,
        top_k=top_k,
        branch_label=baseline_label,
        branch_filter=baseline_filter,
    )
    branch = _collapse_branch(
        branch_peaks,
        top_k=top_k,
        branch_label=branch_label,
        branch_filter=branch_filter,
    )
    baseline_ids = set(baseline["injection_id"].dropna().astype(str))
    branch_ids = set(branch["injection_id"].dropna().astype(str))
    overlap_ids = baseline_ids & branch_ids
    truth = _truth_columns(pd.concat([baseline_peaks, branch_peaks], ignore_index=True, sort=False))
    merged = baseline.merge(branch, on="injection_id", how="outer", validate="one_to_one")
    merged[f"{baseline_label}_evaluated"] = merged["injection_id"].astype(str).isin(baseline_ids)
    merged[f"{branch_label}_evaluated"] = merged["injection_id"].astype(str).isin(branch_ids)
    if comparison_scope not in {"intersection", "union"}:
        raise ValueError("comparison_scope must be 'intersection' or 'union'")
    if comparison_scope == "intersection":
        merged = merged.loc[merged["injection_id"].astype(str).isin(overlap_ids)].copy()
    if not truth.empty:
        merged = merged.merge(truth, on="injection_id", how="left")

    base_top1_col = f"{baseline_label}_top1_match"
    base_topk_col = f"{baseline_label}_top{top_k}_match"
    branch_top1_col = f"{branch_label}_top1_match"
    branch_topk_col = f"{branch_label}_top{top_k}_match"
    for col in (base_top1_col, base_topk_col, branch_top1_col, branch_topk_col):
        if col not in merged:
            merged[col] = False
        merged[col] = _as_bool(merged[col])

    base_top1 = merged[base_top1_col]
    base_topk = merged[base_topk_col]
    alt_top1 = merged[branch_top1_col]
    alt_topk = merged[branch_topk_col]
    both_evaluated = (
        merged[f"{baseline_label}_evaluated"].fillna(False).astype(bool)
        & merged[f"{branch_label}_evaluated"].fillna(False).astype(bool)
    )
    merged["union_top1_match"] = base_top1 | alt_top1
    merged[f"union_top{top_k}_match"] = base_topk | alt_topk
    merged["branch_rescued_top1"] = both_evaluated & ~base_top1 & alt_top1
    merged[f"branch_rescued_top{top_k}"] = both_evaluated & ~base_topk & alt_topk
    merged["branch_reranked_to_top1"] = both_evaluated & base_topk & ~base_top1 & alt_top1
    merged["branch_degraded_top1"] = both_evaluated & base_top1 & ~alt_top1
    merged[f"branch_degraded_top{top_k}"] = both_evaluated & base_topk & ~alt_topk

    status = np.full(len(merged), "both_missed", dtype=object)
    baseline_only = merged[f"{baseline_label}_evaluated"] & ~merged[f"{branch_label}_evaluated"]
    branch_only = ~merged[f"{baseline_label}_evaluated"] & merged[f"{branch_label}_evaluated"]
    status[baseline_only.to_numpy()] = "baseline_only_not_evaluated_by_branch"
    status[branch_only.to_numpy()] = "branch_only_not_evaluated_by_baseline"
    status[base_top1 & alt_top1] = "both_top1"
    status[base_topk & alt_topk & ~(base_top1 & alt_top1)] = "both_topk_not_both_top1"
    status[merged[f"branch_rescued_top{top_k}"].to_numpy()] = f"branch_rescued_top{top_k}"
    status[merged["branch_rescued_top1"].to_numpy()] = "branch_rescued_top1"
    status[merged["branch_reranked_to_top1"].to_numpy()] = "branch_reranked_to_top1"
    status[merged[f"branch_degraded_top{top_k}"].to_numpy()] = f"branch_degraded_top{top_k}"
    status[merged["branch_degraded_top1"].to_numpy()] = "branch_degraded_top1"
    merged["branch_comparison_status"] = status

    binned = add_recovery_bins(
        merged,
        period_bins=period_bins,
        radius_bins=radius_bins,
        tmag_bins=tmag_bins,
    )
    cell_rows = []
    for keys, group in binned.groupby(["tmag_bin", "period_bin", "radius_bin"], dropna=False):
        n = int(len(group))
        cell_rows.append(
            {
                "tmag_bin": str(keys[0]),
                "period_bin": str(keys[1]),
                "radius_bin": str(keys[2]),
                "n": n,
                f"{baseline_label}_top1_n": int(group[base_top1_col].sum()),
                f"{baseline_label}_top{top_k}_n": int(group[base_topk_col].sum()),
                f"{branch_label}_top1_n": int(group[branch_top1_col].sum()),
                f"{branch_label}_top{top_k}_n": int(group[branch_topk_col].sum()),
                "union_top1_n": int(group["union_top1_match"].sum()),
                f"union_top{top_k}_n": int(group[f"union_top{top_k}_match"].sum()),
                "branch_rescued_top1_n": int(group["branch_rescued_top1"].sum()),
                f"branch_rescued_top{top_k}_n": int(group[f"branch_rescued_top{top_k}"].sum()),
                "branch_degraded_top1_n": int(group["branch_degraded_top1"].sum()),
                f"branch_degraded_top{top_k}_n": int(group[f"branch_degraded_top{top_k}"].sum()),
            }
        )
    cells = (
        pd.DataFrame(cell_rows)
        .sort_values([f"branch_rescued_top{top_k}_n", "n"], ascending=[False, False])
        .reset_index(drop=True)
        if cell_rows
        else pd.DataFrame()
    )
    summary = summarize_comparison(
        merged,
        top_k=top_k,
        baseline_label=baseline_label,
        branch_label=branch_label,
        comparison_scope=comparison_scope,
        n_baseline_injections=len(baseline_ids),
        n_branch_injections=len(branch_ids),
        n_overlap_injections=len(overlap_ids),
    )
    return binned, cells, summary


def _n(series: pd.Series) -> int:
    return int(_as_bool(series).sum())


def summarize_comparison(
    comparison: pd.DataFrame,
    *,
    top_k: int,
    baseline_label: str,
    branch_label: str,
    comparison_scope: str = "intersection",
    n_baseline_injections: int | None = None,
    n_branch_injections: int | None = None,
    n_overlap_injections: int | None = None,
) -> dict[str, Any]:
    n = int(len(comparison))
    base_top1 = f"{baseline_label}_top1_match"
    base_topk = f"{baseline_label}_top{top_k}_match"
    alt_top1 = f"{branch_label}_top1_match"
    alt_topk = f"{branch_label}_top{top_k}_match"
    keys = {
        f"{baseline_label}_top1": base_top1,
        f"{baseline_label}_top{top_k}": base_topk,
        f"{branch_label}_top1": alt_top1,
        f"{branch_label}_top{top_k}": alt_topk,
        "union_top1": "union_top1_match",
        f"union_top{top_k}": f"union_top{top_k}_match",
        "branch_rescued_top1": "branch_rescued_top1",
        f"branch_rescued_top{top_k}": f"branch_rescued_top{top_k}",
        "branch_reranked_to_top1": "branch_reranked_to_top1",
        "branch_degraded_top1": "branch_degraded_top1",
        f"branch_degraded_top{top_k}": f"branch_degraded_top{top_k}",
    }
    summary: dict[str, Any] = {
        "n_injections": n,
        "top_k": int(top_k),
        "baseline_label": baseline_label,
        "branch_label": branch_label,
        "comparison_scope": comparison_scope,
        "n_baseline_injections": int(n_baseline_injections if n_baseline_injections is not None else n),
        "n_branch_injections": int(n_branch_injections if n_branch_injections is not None else n),
        "n_overlap_injections": int(n_overlap_injections if n_overlap_injections is not None else n),
        "status_counts": {
            str(key): int(value)
            for key, value in comparison["branch_comparison_status"].value_counts().sort_index().items()
        },
    }
    for name, col in keys.items():
        value = _n(comparison[col]) if col in comparison else 0
        summary[name] = {
            "n": value,
            "fraction": float(value / n) if n else float("nan"),
        }
    return summary


def write_markdown_report(summary: dict[str, Any], *, out_path: Path) -> None:
    top_k = int(summary["top_k"])
    baseline = summary["baseline_label"]
    branch = summary["branch_label"]
    lines = [
        "# BLS Branch Comparison",
        "",
        f"- injections compared: `{summary['n_injections']:,}`",
        f"- comparison scope: `{summary['comparison_scope']}`",
        f"- baseline injections: `{summary['n_baseline_injections']:,}`",
        f"- branch injections: `{summary['n_branch_injections']:,}`",
        f"- overlap injections: `{summary['n_overlap_injections']:,}`",
        f"- `{baseline}` top-1: `{summary[f'{baseline}_top1']['n']:,}`",
        f"- `{baseline}` top-{top_k}: `{summary[f'{baseline}_top{top_k}']['n']:,}`",
        f"- `{branch}` top-1: `{summary[f'{branch}_top1']['n']:,}`",
        f"- `{branch}` top-{top_k}: `{summary[f'{branch}_top{top_k}']['n']:,}`",
        f"- union top-1: `{summary['union_top1']['n']:,}`",
        f"- union top-{top_k}: `{summary[f'union_top{top_k}']['n']:,}`",
        f"- branch rescued top-1: `{summary['branch_rescued_top1']['n']:,}`",
        f"- branch rescued top-{top_k}: `{summary[f'branch_rescued_top{top_k}']['n']:,}`",
        f"- branch degraded top-1: `{summary['branch_degraded_top1']['n']:,}`",
        f"- branch degraded top-{top_k}: `{summary[f'branch_degraded_top{top_k}']['n']:,}`",
        "",
        "Status counts:",
    ]
    for key, value in summary["status_counts"].items():
        lines.append(f"- `{key}`: `{value:,}`")
    lines.append("")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-peak-table", type=Path, default=DEFAULT_BASELINE_TABLE)
    parser.add_argument("--branch-peak-table", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument("--baseline-label", default="standard")
    parser.add_argument("--branch-label", default="short_pmax2")
    parser.add_argument("--baseline-search-branch", default="")
    parser.add_argument("--branch-search-branch", default="")
    parser.add_argument("--comparison-scope", choices=("intersection", "union"), default="intersection")
    parser.add_argument("--period-bins", default=",".join(str(v) for v in DEFAULT_PERIOD_BINS))
    parser.add_argument("--radius-bins", default=",".join(str(v) for v in DEFAULT_RADIUS_BINS))
    parser.add_argument("--tmag-bins", default=",".join(str(v) for v in DEFAULT_TMAG_BINS))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    baseline = _read_table(args.baseline_peak_table)
    branch = _read_table(args.branch_peak_table)
    comparison, cells, summary = compare_branches(
        baseline,
        branch,
        top_k=args.top_k,
        baseline_label=args.baseline_label,
        branch_label=args.branch_label,
        baseline_filter=args.baseline_search_branch or None,
        branch_filter=args.branch_search_branch or None,
        comparison_scope=args.comparison_scope,
        period_bins=_parse_float_tuple(args.period_bins),
        radius_bins=_parse_float_tuple(args.radius_bins),
        tmag_bins=_parse_float_tuple(args.tmag_bins),
    )
    summary.update(
        {
            "baseline_peak_table": str(args.baseline_peak_table),
            "branch_peak_table": str(args.branch_peak_table),
            "out_dir": str(args.out_dir),
            "baseline_search_branch": args.baseline_search_branch,
            "branch_search_branch": args.branch_search_branch,
            "comparison_scope": args.comparison_scope,
            "period_bins": list(_parse_float_tuple(args.period_bins)),
            "radius_bins": list(_parse_float_tuple(args.radius_bins)),
            "tmag_bins": list(_parse_float_tuple(args.tmag_bins)),
        }
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(args.out_dir / "branch_comparison_by_injection.csv", index=False)
    cells.to_csv(args.out_dir / "branch_comparison_by_cell.csv", index=False)
    rescued = comparison.loc[
        comparison[f"branch_rescued_top{args.top_k}"].fillna(False)
        | comparison["branch_rescued_top1"].fillna(False)
        | comparison["branch_reranked_to_top1"].fillna(False)
    ].copy()
    rescued.to_csv(args.out_dir / "branch_rescued_injections.csv", index=False)
    degraded = comparison.loc[
        comparison[f"branch_degraded_top{args.top_k}"].fillna(False)
        | comparison["branch_degraded_top1"].fillna(False)
    ].copy()
    degraded.to_csv(args.out_dir / "branch_degraded_injections.csv", index=False)
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )
    write_markdown_report(summary, out_path=args.out_dir / "summary.md")

    print("[bls-branch-compare] complete")
    print(f"  injections: {summary['n_injections']:,}")
    print(f"  {args.baseline_label} top{args.top_k}: {summary[f'{args.baseline_label}_top{args.top_k}']['n']:,}")
    print(f"  {args.branch_label} top{args.top_k}: {summary[f'{args.branch_label}_top{args.top_k}']['n']:,}")
    print(f"  rescued top{args.top_k}: {summary[f'branch_rescued_top{args.top_k}']['n']:,}")
    print(f"  degraded top{args.top_k}: {summary[f'branch_degraded_top{args.top_k}']['n']:,}")
    print(f"  out: {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
