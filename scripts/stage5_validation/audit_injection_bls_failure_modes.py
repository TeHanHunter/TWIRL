#!/usr/bin/env python3
"""Audit injected BLS failure modes after peak-table generation.

The peak gate answers whether injected truth is present in the top-N peak list.
This companion audit asks which non-recoveries look physically observable from
the injection metadata and therefore deserve search tuning before human triage.
It intentionally uses only truth columns for audit stratification, not for any
real-candidate inference.
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


DEFAULT_PEAK_TABLE = (
    Path(__file__).resolve().parents[2]
    / "reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/"
    / "peak_training/s56_20k_injection_bls_peaks_chunked.csv"
)
DEFAULT_OUT_DIR = (
    Path(__file__).resolve().parents[2]
    / "reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/"
    / "peak_training_gate/failure_modes"
)


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _first_existing(df: pd.DataFrame, candidates: tuple[str, ...]) -> str | None:
    for col in candidates:
        if col in df.columns:
            return col
    return None


def _numeric(df: pd.DataFrame, col: str | None, default: float = np.nan) -> pd.Series:
    if col is None or col not in df.columns:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").astype(float)


def _safe_log10(values: pd.Series) -> pd.Series:
    arr = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    return pd.Series(np.log10(np.where(arr > 0, arr, np.nan)), index=values.index)


def add_observability_proxy(injections: pd.DataFrame) -> pd.DataFrame:
    """Attach a simple physical proxy for injection detectability.

    The proxy is not a calibrated SNR. It is the sampled/model depth multiplied
    by sqrt(number of good in-transit cadences), then used only for relative
    ordering within broad Tmag bins.
    """

    out = injections.copy()
    depth_col = _first_existing(
        out,
        ("truth_sampled_model_depth", "truth_model_depth", "truth_depth", "depth"),
    )
    n_good_col = _first_existing(out, ("truth_n_good_in_transit", "truth_n_in_transit"))
    depth = _numeric(out, depth_col).abs()
    n_good = _numeric(out, n_good_col)
    n_good = n_good.where(n_good > 0, np.nan)
    out["truth_depth_for_audit"] = depth
    out["truth_n_good_in_transit_for_audit"] = n_good
    out["observability_proxy"] = depth * np.sqrt(n_good)
    out["log_observability_proxy"] = _safe_log10(out["observability_proxy"])
    out["low_in_transit_cadences"] = n_good.fillna(0).lt(3)
    return out


def _reference_thresholds(
    injections: pd.DataFrame,
    *,
    top_k: int,
    quantile: float,
    min_reference_count: int,
) -> dict[str, float]:
    topk_col = f"top{top_k}_match"
    recovered = injections.loc[injections[topk_col].fillna(False)].copy()
    finite = recovered["log_observability_proxy"].replace([np.inf, -np.inf], np.nan).dropna()
    global_threshold = float(finite.quantile(quantile)) if len(finite) else float("inf")
    thresholds: dict[str, float] = {}
    if "tmag_bin" not in injections.columns:
        thresholds[""] = global_threshold
        return thresholds
    for key, group in recovered.groupby("tmag_bin", dropna=False):
        values = group["log_observability_proxy"].replace([np.inf, -np.inf], np.nan).dropna()
        if len(values) >= min_reference_count:
            thresholds[str(key)] = float(values.quantile(quantile))
        else:
            thresholds[str(key)] = global_threshold
    for key in injections["tmag_bin"].astype(str).fillna("").unique():
        thresholds.setdefault(str(key), global_threshold)
    return thresholds


def _best_available_peak(peaks: pd.DataFrame) -> pd.DataFrame:
    work = peaks.copy()
    work["_is_candidate_peak"] = (
        work["is_candidate_peak"].fillna(False).astype(bool) if "is_candidate_peak" in work else True
    )
    work = work.loc[work["_is_candidate_peak"]].copy()
    if work.empty:
        return pd.DataFrame(columns=["injection_id"])
    for col in ("sde", "depth_snr", "peak_rank", "period_d", "duration_min", "period_ratio"):
        if col in work.columns:
            work[f"_{col}"] = pd.to_numeric(work[col], errors="coerce")
        else:
            work[f"_{col}"] = np.nan
    work["_sort_sde"] = work["_sde"].fillna(-np.inf)
    work["_sort_depth_snr"] = work["_depth_snr"].fillna(-np.inf)
    work["_sort_peak_rank"] = work["_peak_rank"].fillna(np.inf)
    idx = (
        work.sort_values(
            ["injection_id", "_sort_sde", "_sort_depth_snr", "_sort_peak_rank"],
            ascending=[True, False, False, True],
            kind="stable",
        )
        .groupby("injection_id", sort=False)
        .head(1)
        .index
    )
    keep = []
    for col in (
        "injection_id",
        "aperture",
        "peak_rank",
        "period_d",
        "duration_min",
        "depth",
        "depth_snr",
        "sde",
        "match_kind",
        "period_ratio",
        "nearest_harmonic_factor",
        "harmonic_period_rel_err",
        "t0_phase_err_min",
        "n_cad_quality",
        "n_cad_kept",
        "n_cad_edge_trimmed",
        "n_cad_sigma_clipped",
        "dropout_frac",
        "quality_dropout_frac",
    ):
        if col in work.columns:
            keep.append(col)
    best = work.loc[idx, keep].copy()
    rename = {col: f"best_peak_{col}" for col in keep if col != "injection_id"}
    return best.rename(columns=rename).reset_index(drop=True)


def build_failure_mode_table(
    peaks: pd.DataFrame,
    *,
    top_k: int,
    period_bins: tuple[float, ...],
    radius_bins: tuple[float, ...],
    tmag_bins: tuple[float, ...],
    observability_quantile: float,
    min_reference_count: int,
) -> pd.DataFrame:
    injections = build_injection_level_table(peaks, top_k=top_k)
    injections = add_recovery_bins(
        injections,
        period_bins=period_bins,
        radius_bins=radius_bins,
        tmag_bins=tmag_bins,
    )
    injections = add_observability_proxy(injections)
    thresholds = _reference_thresholds(
        injections,
        top_k=top_k,
        quantile=observability_quantile,
        min_reference_count=min_reference_count,
    )
    injections["observability_reference_log_threshold"] = injections["tmag_bin"].astype(str).map(thresholds)
    topk_col = f"top{top_k}_match"
    injections["high_observability_not_in_topk"] = (
        ~injections[topk_col].fillna(False)
        & ~injections["low_in_transit_cadences"].fillna(False)
        & injections["log_observability_proxy"].ge(injections["observability_reference_log_threshold"])
    )

    mode = np.full(len(injections), f"low_observability_not_in_top{top_k}", dtype=object)
    mode[injections["low_in_transit_cadences"].fillna(False).to_numpy()] = f"low_cadence_not_in_top{top_k}"
    mode[injections["high_observability_not_in_topk"].fillna(False).to_numpy()] = (
        f"high_observability_not_in_top{top_k}"
    )
    mode[injections["ranking_loss"].fillna(False).to_numpy()] = "ranker_fixable"
    mode[injections["top1_match"].fillna(False).to_numpy()] = "top1_recovered"
    injections["failure_mode"] = mode

    best = _best_available_peak(peaks)
    if not best.empty:
        injections = injections.merge(best, on="injection_id", how="left")
    return injections.sort_values(
        ["failure_mode", "log_observability_proxy"], ascending=[True, False], kind="stable"
    ).reset_index(drop=True)


def summarize_failure_modes(table: pd.DataFrame, *, top_k: int) -> dict[str, Any]:
    n = int(len(table))
    topk_col = f"top{top_k}_match"
    mode_counts = table["failure_mode"].value_counts().sort_index()
    summary: dict[str, Any] = {
        "n_injections": n,
        "top_k": int(top_k),
        "top1_n": int(table["top1_match"].fillna(False).sum()),
        f"top{top_k}_n": int(table[topk_col].fillna(False).sum()),
        "ranking_loss_n": int(table["ranking_loss"].fillna(False).sum()),
        f"not_in_top{top_k}_n": int((~table[topk_col].fillna(False)).sum()),
        "high_observability_not_in_topk_n": int(table["high_observability_not_in_topk"].fillna(False).sum()),
        "failure_mode_counts": {str(key): int(value) for key, value in mode_counts.items()},
    }
    if n:
        summary["top1_fraction"] = float(summary["top1_n"] / n)
        summary[f"top{top_k}_fraction"] = float(summary[f"top{top_k}_n"] / n)
        summary["high_observability_not_in_topk_fraction"] = float(
            summary["high_observability_not_in_topk_n"] / n
        )
    grouped_rows = []
    for keys, group in table.groupby(["tmag_bin", "period_bin", "radius_bin"], dropna=False):
        grouped_rows.append(
            {
                "tmag_bin": str(keys[0]),
                "period_bin": str(keys[1]),
                "radius_bin": str(keys[2]),
                "n": int(len(group)),
                "top1_n": int(group["top1_match"].fillna(False).sum()),
                f"top{top_k}_n": int(group[topk_col].fillna(False).sum()),
                "high_observability_not_in_topk_n": int(
                    group["high_observability_not_in_topk"].fillna(False).sum()
                ),
            }
        )
    summary["n_cells"] = len(grouped_rows)
    summary["n_cells_with_high_observability_misses"] = int(
        sum(row["high_observability_not_in_topk_n"] > 0 for row in grouped_rows)
    )
    return summary


def build_cell_failure_table(table: pd.DataFrame, *, top_k: int) -> pd.DataFrame:
    """Summarize failure modes in broad period/radius/Tmag cells."""

    topk_col = f"top{top_k}_match"
    rows = []
    for keys, group in table.groupby(["tmag_bin", "period_bin", "radius_bin"], dropna=False):
        n = int(len(group))
        top1_n = int(group["top1_match"].fillna(False).sum())
        topk_n = int(group[topk_col].fillna(False).sum())
        high_miss_n = int(group["high_observability_not_in_topk"].fillna(False).sum())
        rows.append(
            {
                "tmag_bin": str(keys[0]),
                "period_bin": str(keys[1]),
                "radius_bin": str(keys[2]),
                "n": n,
                "top1_n": top1_n,
                "top1_fraction": float(top1_n / n) if n else float("nan"),
                f"top{top_k}_n": topk_n,
                f"top{top_k}_fraction": float(topk_n / n) if n else float("nan"),
                "ranking_loss_n": int(group["ranking_loss"].fillna(False).sum()),
                "high_observability_not_in_topk_n": high_miss_n,
                "high_observability_not_in_topk_fraction": float(high_miss_n / n) if n else float("nan"),
                "median_log_observability_proxy": float(
                    group["log_observability_proxy"].replace([np.inf, -np.inf], np.nan).median()
                ),
            }
        )
    if not rows:
        return pd.DataFrame()
    return (
        pd.DataFrame(rows)
        .sort_values(
            ["high_observability_not_in_topk_n", "n", "tmag_bin", "period_bin", "radius_bin"],
            ascending=[False, False, True, True, True],
        )
        .reset_index(drop=True)
    )


def write_markdown_report(
    summary: dict[str, Any],
    hard_misses: pd.DataFrame,
    *,
    out_path: Path,
    top_k: int,
) -> None:
    lines = [
        "# Injected BLS Failure-Mode Audit",
        "",
        f"- injections audited: `{summary['n_injections']:,}`",
        f"- top-1 truth recovery: `{summary['top1_n']:,}/{summary['n_injections']:,}`",
        f"- top-{top_k} truth recovery: `{summary[f'top{top_k}_n']:,}/{summary['n_injections']:,}`",
        f"- ranker-fixable rows: `{summary['ranking_loss_n']:,}`",
        f"- high-observability not-in-top-{top_k} rows: "
        f"`{summary['high_observability_not_in_topk_n']:,}`",
        "",
        "Failure modes:",
    ]
    for key, value in summary["failure_mode_counts"].items():
        lines.append(f"- `{key}`: `{value:,}`")
    lines.extend(
        [
            "",
            "High-observability misses use a relative audit proxy: sampled/model depth "
            "times sqrt(good in-transit cadences), compared against recovered injections "
            "within broad Tmag bins. They are the rows to inspect first when tuning BLS.",
            "",
        ]
    )
    if not hard_misses.empty:
        lines.append("Top high-observability misses:")
        cols = [
            col
            for col in (
                "injection_id",
                "tic",
                "tmag",
                "truth_period_d",
                "truth_radius_rearth",
                "truth_duration_min",
                "truth_n_good_in_transit_for_audit",
                "observability_proxy",
                "best_peak_period_d",
                "best_peak_sde",
                "best_peak_period_ratio",
            )
            if col in hard_misses.columns
        ]
        for _, row in hard_misses.head(10).iterrows():
            bits = []
            for col in cols:
                value = row.get(col)
                if isinstance(value, float):
                    bits.append(f"{col}={value:.5g}")
                else:
                    bits.append(f"{col}={value}")
            lines.append("- " + ", ".join(bits))
        lines.append("")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peak-table", type=Path, default=DEFAULT_PEAK_TABLE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument("--period-bins", default=",".join(str(v) for v in DEFAULT_PERIOD_BINS))
    parser.add_argument("--radius-bins", default=",".join(str(v) for v in DEFAULT_RADIUS_BINS))
    parser.add_argument("--tmag-bins", default=",".join(str(v) for v in DEFAULT_TMAG_BINS))
    parser.add_argument("--observability-quantile", type=float, default=0.25)
    parser.add_argument("--min-reference-count", type=int, default=20)
    parser.add_argument("--max-hard-misses", type=int, default=200)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    peaks = _read_table(args.peak_table)
    table = build_failure_mode_table(
        peaks,
        top_k=args.top_k,
        period_bins=_parse_float_tuple(args.period_bins),
        radius_bins=_parse_float_tuple(args.radius_bins),
        tmag_bins=_parse_float_tuple(args.tmag_bins),
        observability_quantile=args.observability_quantile,
        min_reference_count=args.min_reference_count,
    )
    summary = summarize_failure_modes(table, top_k=args.top_k)
    summary.update(
        {
            "peak_table": str(args.peak_table),
            "out_dir": str(args.out_dir),
            "observability_quantile": float(args.observability_quantile),
            "min_reference_count": int(args.min_reference_count),
            "period_bins": list(_parse_float_tuple(args.period_bins)),
            "radius_bins": list(_parse_float_tuple(args.radius_bins)),
            "tmag_bins": list(_parse_float_tuple(args.tmag_bins)),
        }
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out_dir / "failure_modes_by_injection.csv", index=False)
    cell_table = build_cell_failure_table(table, top_k=args.top_k)
    cell_table.to_csv(args.out_dir / "failure_modes_by_cell.csv", index=False)
    counts = table["failure_mode"].value_counts().sort_index().rename_axis("failure_mode").reset_index(name="n")
    counts.to_csv(args.out_dir / "failure_mode_counts.csv", index=False)
    hard = (
        table.loc[table["high_observability_not_in_topk"].fillna(False)]
        .sort_values("log_observability_proxy", ascending=False, kind="stable")
        .head(args.max_hard_misses)
    )
    hard.to_csv(args.out_dir / "high_observability_misses.csv", index=False)
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )
    write_markdown_report(summary, hard, out_path=args.out_dir / "summary.md", top_k=args.top_k)

    print("[bls-failure-audit] complete")
    print(f"  injections: {summary['n_injections']:,}")
    print(f"  high-observability not-in-top{args.top_k}: {summary['high_observability_not_in_topk_n']:,}")
    print(f"  out: {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
