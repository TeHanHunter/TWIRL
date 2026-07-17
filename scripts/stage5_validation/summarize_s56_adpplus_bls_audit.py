#!/usr/bin/env python3
"""Summarize the S56 ADP+ two-aperture BLS audit into decision tables."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError


DEFAULT_AUDIT_DIR = Path("reports/stage5_validation/s56_adpplus_bls_audit_pdo")


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _num(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        return pd.Series(np.nan, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce").replace([np.inf, -np.inf], np.nan)


def _read_csv_if_present(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except EmptyDataError:
        return pd.DataFrame()


def _bool(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        return pd.Series(False, index=frame.index, dtype=bool)
    vals = frame[column]
    if vals.dtype == bool:
        return vals.fillna(False).astype(bool)
    text = vals.fillna("").astype(str).str.lower()
    return text.isin(["true", "1", "yes", "y"])


def _add_injection_bins(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    tmag = _num(out, "tmag")
    period = _num(out, "truth_period_d")
    radius = _num(out, "truth_radius_rearth")
    out["tmag_bin"] = pd.cut(
        tmag,
        [-np.inf, 17, 18, 19, np.inf],
        labels=["T<17", "17<=T<18", "18<=T<19", "T>=19"],
    ).astype(object).where(tmag.notna(), "missing")
    out["period_bin"] = pd.cut(
        period,
        [0.0, 0.12, 0.2, 0.35, 0.6, 1.0, 1.8, 3.2, 5.6, 10.0, np.inf],
        labels=["0.00-0.12", "0.12-0.20", "0.20-0.35", "0.35-0.60", "0.60-1.00", "1.00-1.80", "1.80-3.20", "3.20-5.60", "5.60-10.0", ">=10.0"],
    ).astype(object).where(period.notna(), "missing")
    out["radius_bin"] = pd.cut(
        radius,
        [0.0, 0.5, 1, 2, 4, 8, 12, np.inf],
        labels=["0.0-0.5", "0.5-1", "1-2", "2-4", "4-8", "8-12", ">=12"],
    ).astype(object).where(radius.notna(), "missing")
    return out


def _peak_rows(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty or "peak_rank" not in frame:
        return frame.head(0).copy()
    out = frame.copy()
    out["peak_rank"] = _num(out, "peak_rank").fillna(0).astype(int)
    return out[out["peak_rank"].gt(0)].copy()


def _summarize_injections(frame: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    peaks = _peak_rows(frame)
    columns = group_cols + [
        "n_injections",
        "top1_exact_n",
        "top1_exact_frac",
        "topn_exact_n",
        "topn_exact_frac",
        "topn_exact_or_harmonic_n",
        "topn_exact_or_harmonic_frac",
        "median_depth_retention",
        "p10_depth_retention",
        "median_sde_rank1",
        "median_trend_ptp",
    ]
    if peaks.empty:
        return pd.DataFrame(columns=columns)
    peaks["_exact"] = _bool(peaks, "exact_ephemeris_match")
    peaks["_harmonic"] = _bool(peaks, "harmonic_ephemeris_match")
    rows: list[dict[str, Any]] = []
    for keys, group in peaks.groupby(group_cols, dropna=False, sort=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        top1 = group[group["peak_rank"].eq(1)]
        by_object = group.groupby("injection_id", dropna=False)
        exact_any = by_object["_exact"].any()
        harmonic_any = by_object["_harmonic"].any()
        ret = _num(group, "depth_retention_frac")
        sde = _num(top1, "sde")
        trend = _num(group, "trend_ptp")
        row = dict(zip(group_cols, keys))
        row.update(
            {
                "n_injections": int(group["injection_id"].nunique()),
                "top1_exact_n": int(top1["_exact"].sum()),
                "top1_exact_frac": float(top1["_exact"].mean()) if len(top1) else np.nan,
                "topn_exact_n": int(exact_any.sum()),
                "topn_exact_frac": float(exact_any.mean()) if len(exact_any) else np.nan,
                "topn_exact_or_harmonic_n": int((exact_any | harmonic_any).sum()),
                "topn_exact_or_harmonic_frac": float((exact_any | harmonic_any).mean()) if len(exact_any) else np.nan,
                "median_depth_retention": float(np.nanmedian(ret)) if ret.notna().any() else np.nan,
                "p10_depth_retention": float(np.nanpercentile(ret.dropna(), 10)) if ret.notna().any() else np.nan,
                "median_sde_rank1": float(np.nanmedian(sde)) if sde.notna().any() else np.nan,
                "median_trend_ptp": float(np.nanmedian(trend)) if trend.notna().any() else np.nan,
            }
        )
        rows.append(row)
    return pd.DataFrame(rows, columns=columns)


def _summarize_real(frame: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    peaks = _peak_rows(frame)
    columns = group_cols + ["n_tics", "median_sde_rank1", "median_trend_ptp", "p90_trend_ptp"]
    if peaks.empty:
        return pd.DataFrame(columns=columns)
    rows: list[dict[str, Any]] = []
    for keys, group in peaks.groupby(group_cols, dropna=False, sort=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        top1 = group[group["peak_rank"].eq(1)]
        trend = _num(group, "trend_ptp")
        sde = _num(top1, "sde")
        row = dict(zip(group_cols, keys))
        row.update(
            {
                "n_tics": int(group["tic"].nunique()),
                "median_sde_rank1": float(np.nanmedian(sde)) if sde.notna().any() else np.nan,
                "median_trend_ptp": float(np.nanmedian(trend)) if trend.notna().any() else np.nan,
                "p90_trend_ptp": float(np.nanpercentile(trend.dropna(), 90)) if trend.notna().any() else np.nan,
            }
        )
        rows.append(row)
    return pd.DataFrame(rows, columns=columns)


def _score_branches(injection_summary: pd.DataFrame, real_summary: pd.DataFrame) -> pd.DataFrame:
    if injection_summary.empty:
        return pd.DataFrame()
    key = ["branch", "mask_mode", "aperture"]
    inj = injection_summary.copy()
    real = real_summary.copy() if not real_summary.empty else pd.DataFrame(columns=key)
    keep = key + [
        "n_injections",
        "top1_exact_frac",
        "topn_exact_or_harmonic_frac",
        "median_depth_retention",
        "p10_depth_retention",
        "median_sde_rank1",
    ]
    out = inj[keep].merge(
        real[key + ["median_trend_ptp", "p90_trend_ptp"]],
        on=key,
        how="left",
        suffixes=("", "_real"),
    )
    base = out[out["branch"].eq("current_adp")].copy()
    base = base.rename(
        columns={
            "top1_exact_frac": "base_top1_exact_frac",
            "topn_exact_or_harmonic_frac": "base_topn_exact_or_harmonic_frac",
            "median_depth_retention": "base_median_depth_retention",
            "median_trend_ptp": "base_median_trend_ptp",
        }
    )
    base = base[["aperture", "base_top1_exact_frac", "base_topn_exact_or_harmonic_frac", "base_median_depth_retention", "base_median_trend_ptp"]]
    out = out.merge(base, on="aperture", how="left")
    out["delta_top1_exact_frac_vs_current"] = out["top1_exact_frac"] - out["base_top1_exact_frac"]
    out["delta_topn_exact_or_harmonic_frac_vs_current"] = out["topn_exact_or_harmonic_frac"] - out["base_topn_exact_or_harmonic_frac"]
    out["delta_median_depth_retention_vs_current"] = out["median_depth_retention"] - out["base_median_depth_retention"]
    out["trend_ptp_ratio_vs_current"] = out["median_trend_ptp"] / out["base_median_trend_ptp"]
    out["acceptance_note"] = np.where(
        out["branch"].eq("current_adp"),
        "baseline",
        np.where(
            (out["delta_top1_exact_frac_vs_current"].fillna(-1) >= 0)
            & (out["delta_median_depth_retention_vs_current"].fillna(-1) > -0.05)
            & (out["trend_ptp_ratio_vs_current"].fillna(np.inf) <= 1.0),
            "candidate_improvement",
            "review_tradeoff",
        ),
    )
    return out.sort_values(["aperture", "acceptance_note", "branch", "mask_mode"])


def summarize_audit(audit_dir: Path, out_dir: Path) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    inj_path = audit_dir / "injection_branch_peak_table.csv"
    real_path = audit_dir / "real_branch_peak_table.csv"
    inj = _read_csv_if_present(inj_path)
    real = _read_csv_if_present(real_path)
    inj = _add_injection_bins(inj) if not inj.empty else inj

    inj_overall = _summarize_injections(inj, ["branch", "mask_mode", "aperture"])
    inj_tmag = _summarize_injections(inj, ["branch", "mask_mode", "aperture", "tmag_bin"])
    inj_period = _summarize_injections(inj, ["branch", "mask_mode", "aperture", "period_bin"])
    inj_radius = _summarize_injections(inj, ["branch", "mask_mode", "aperture", "radius_bin"])
    inj_period_radius = _summarize_injections(inj, ["branch", "mask_mode", "aperture", "period_bin", "radius_bin"])
    real_overall = _summarize_real(real, ["branch", "mask_mode", "aperture"])
    real_bucket = _summarize_real(real, ["branch", "mask_mode", "aperture", "real_audit_bucket"])
    score = _score_branches(inj_overall, real_overall)

    outputs = {
        "injection_overall": out_dir / "injection_recovery_overall.csv",
        "injection_by_tmag": out_dir / "injection_recovery_by_tmag.csv",
        "injection_by_period": out_dir / "injection_recovery_by_period.csv",
        "injection_by_radius": out_dir / "injection_recovery_by_radius.csv",
        "injection_by_period_radius": out_dir / "injection_recovery_by_period_radius.csv",
        "real_overall": out_dir / "real_trend_overall.csv",
        "real_by_bucket": out_dir / "real_trend_by_bucket.csv",
        "branch_score": out_dir / "branch_acceptance_score.csv",
    }
    inj_overall.to_csv(outputs["injection_overall"], index=False)
    inj_tmag.to_csv(outputs["injection_by_tmag"], index=False)
    inj_period.to_csv(outputs["injection_by_period"], index=False)
    inj_radius.to_csv(outputs["injection_by_radius"], index=False)
    inj_period_radius.to_csv(outputs["injection_by_period_radius"], index=False)
    real_overall.to_csv(outputs["real_overall"], index=False)
    real_bucket.to_csv(outputs["real_by_bucket"], index=False)
    score.to_csv(outputs["branch_score"], index=False)

    status_counts = {
        "injection": inj.get("status", pd.Series(dtype=object)).fillna("").astype(str).value_counts().to_dict(),
        "real": real.get("status", pd.Series(dtype=object)).fillna("").astype(str).value_counts().to_dict(),
    }
    summary = {
        "audit_dir": str(audit_dir),
        "out_dir": str(out_dir),
        "n_injection_rows": int(len(inj)),
        "n_real_rows": int(len(real)),
        "status_counts": status_counts,
        "outputs": {key: str(path) for key, path in outputs.items()},
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-dir", type=Path, default=DEFAULT_AUDIT_DIR)
    parser.add_argument("--out-dir", type=Path, default=None)
    args = parser.parse_args(argv)
    out_dir = args.out_dir or (args.audit_dir / "summary_tables")
    summary = summarize_audit(args.audit_dir, out_dir)
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
