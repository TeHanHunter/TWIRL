#!/usr/bin/env python3
"""Summarize aperture-dependent signal survival for pre-detrend injections."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


TMAG_BINS = (-np.inf, 16.0, 17.0, 18.0, 19.0, 20.0, np.inf)
TMAG_LABELS = ("<16", "16-17", "17-18", "18-19", "19-20", ">20")
PERIOD_BINS = (0.0, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, np.inf)
PERIOD_LABELS = ("<0.25", "0.25-0.5", "0.5-1", "1-2", "2-5", "5-10", ">10")
DEPTH_BINS = (0.0, 0.01, 0.03, 0.1, 0.3, 1.0, np.inf)
DEPTH_LABELS = ("<1%", "1-3%", "3-10%", "10-30%", "30-100%", ">100%")

DEFAULT_GROUPS = {
    "adp": ("DET_FLUX_ADP_SML", "DET_FLUX_ADP", "DET_FLUX_ADP_LAG"),
    "canonical": ("DET_FLUX_SML", "DET_FLUX", "DET_FLUX_LAG"),
    "detrended_all": (
        "DET_FLUX_ADP_SML",
        "DET_FLUX_ADP",
        "DET_FLUX_ADP_LAG",
        "DET_FLUX_SML",
        "DET_FLUX",
        "DET_FLUX_LAG",
    ),
    "raw": ("RAW_FLUX_Small", "RAW_FLUX_Primary", "RAW_FLUX_Large"),
}


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def _read_survival(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"injection_id", "aperture", "delta_signal_multi_snr_mad"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"missing required columns: {missing}")
    for col in (
        "delta_signal_multi_snr_mad",
        "delta_signal_snr_mad",
        "depth_retention_frac",
        "depth_retention_vs_model_frac",
        "truth_period_d",
        "truth_model_depth",
        "truth_depth",
        "truth_radius_rearth",
        "tmag",
        "n_good_ap_in_transit",
    ):
        if col in df:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _best_rows(df: pd.DataFrame, apertures: tuple[str, ...]) -> pd.DataFrame:
    sub = df[df["aperture"].isin(apertures)].copy()
    if sub.empty:
        return sub
    sub = sub.sort_values(["injection_id", "delta_signal_multi_snr_mad"], na_position="first")
    return sub.groupby("injection_id", sort=False).tail(1).copy()


def _aperture_summary(df: pd.DataFrame, recoverable_snr: float) -> pd.DataFrame:
    grouped = (
        df.groupby("aperture", observed=False)
        .agg(
            n=("injection_id", "nunique"),
            n_snr_ge_threshold=("delta_signal_multi_snr_mad", lambda x: int((x >= recoverable_snr).sum())),
            median_snr=("delta_signal_multi_snr_mad", "median"),
            p90_snr=("delta_signal_multi_snr_mad", lambda x: float(np.nanpercentile(x, 90))),
            median_depth_retention=("depth_retention_frac", "median"),
            median_tmag=("tmag", "median"),
        )
        .reset_index()
    )
    grouped["frac_snr_ge_threshold"] = grouped["n_snr_ge_threshold"] / grouped["n"].replace(0, np.nan)
    return grouped.sort_values(["n_snr_ge_threshold", "median_snr"], ascending=[False, False])


def _group_summary(df: pd.DataFrame, recoverable_snr: float) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    rows: list[dict[str, Any]] = []
    best_tables: dict[str, pd.DataFrame] = {}
    for name, apertures in DEFAULT_GROUPS.items():
        best = _best_rows(df, apertures)
        best_tables[name] = best
        if best.empty:
            rows.append({"aperture_group": name, "n": 0, "n_snr_ge_threshold": 0})
            continue
        rows.append(
            {
                "aperture_group": name,
                "apertures": ",".join(apertures),
                "n": int(best["injection_id"].nunique()),
                "n_snr_ge_threshold": int((best["delta_signal_multi_snr_mad"] >= recoverable_snr).sum()),
                "frac_snr_ge_threshold": float((best["delta_signal_multi_snr_mad"] >= recoverable_snr).mean()),
                "median_snr": float(best["delta_signal_multi_snr_mad"].median()),
                "median_depth_retention": float(best["depth_retention_frac"].median())
                if "depth_retention_frac" in best
                else np.nan,
                "best_aperture_counts": dict(best["aperture"].value_counts().sort_values(ascending=False)),
            }
        )
    out = pd.DataFrame(rows).sort_values(["n_snr_ge_threshold", "median_snr"], ascending=[False, False])
    return out, best_tables


def _current_vs_best(df: pd.DataFrame, current_aperture: str, best: pd.DataFrame, recoverable_snr: float) -> pd.DataFrame:
    current = df[df["aperture"].eq(current_aperture)].copy()
    keep_cols = [
        col
        for col in (
            "injection_id",
            "tic",
            "tmag",
            "truth_period_d",
            "truth_model_depth",
            "truth_radius_rearth",
            "delta_signal_multi_snr_mad",
            "depth_retention_frac",
        )
        if col in current.columns
    ]
    current = current.loc[:, keep_cols].rename(
        columns={
            "delta_signal_multi_snr_mad": "current_snr",
            "depth_retention_frac": "current_depth_retention",
        }
    )
    best_cols = [
        col
        for col in (
            "injection_id",
            "aperture",
            "delta_signal_multi_snr_mad",
            "depth_retention_frac",
            "n_good_ap_in_transit",
        )
        if col in best.columns
    ]
    best = best.loc[:, best_cols].rename(
        columns={
            "aperture": "best_aperture",
            "delta_signal_multi_snr_mad": "best_snr",
            "depth_retention_frac": "best_depth_retention",
            "n_good_ap_in_transit": "best_n_good_in_transit",
        }
    )
    joined = current.merge(best, on="injection_id", how="outer")
    joined["current_snr_ge_threshold"] = joined["current_snr"] >= recoverable_snr
    joined["best_snr_ge_threshold"] = joined["best_snr"] >= recoverable_snr
    joined["newly_recoverable_by_aperture"] = ~joined["current_snr_ge_threshold"] & joined["best_snr_ge_threshold"]
    joined["lost_by_aperture_choice"] = joined["current_snr_ge_threshold"] & ~joined["best_snr_ge_threshold"]
    joined["snr_gain"] = joined["best_snr"] - joined["current_snr"]
    joined["snr_gain_factor"] = joined["best_snr"] / joined["current_snr"].replace(0, np.nan)
    if "tmag" in joined:
        joined["tmag_bin"] = pd.cut(joined["tmag"], TMAG_BINS, labels=TMAG_LABELS)
    if "truth_period_d" in joined:
        joined["period_bin"] = pd.cut(joined["truth_period_d"], PERIOD_BINS, labels=PERIOD_LABELS)
    if "truth_model_depth" in joined:
        joined["depth_bin"] = pd.cut(joined["truth_model_depth"], DEPTH_BINS, labels=DEPTH_LABELS)
    return joined.sort_values(["newly_recoverable_by_aperture", "snr_gain"], ascending=[False, False])


def _write_bin_tables(joined: pd.DataFrame, out_dir: Path) -> dict[str, dict[str, dict[str, int]]]:
    tables: dict[str, dict[str, dict[str, int]]] = {}
    for col in ("tmag_bin", "period_bin", "depth_bin", "best_aperture"):
        if col not in joined:
            continue
        table = pd.crosstab(joined[col], joined["best_snr_ge_threshold"], dropna=False)
        table = table.rename(columns={False: "not_recoverable_best", True: "recoverable_best"})
        table.to_csv(out_dir / f"best_recoverable_by_{col}.csv")
        tables[col] = {
            str(idx): {str(key): int(value) for key, value in row.items()}
            for idx, row in table.iterrows()
        }
    return tables


def _write_plots(aperture_summary: pd.DataFrame, joined: pd.DataFrame, out_dir: Path) -> dict[str, str]:
    paths: dict[str, str] = {}
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        return {"plot_error": f"{type(exc).__name__}: {exc}"}

    fig, ax = plt.subplots(figsize=(8, 4))
    plot_df = aperture_summary.sort_values("n_snr_ge_threshold", ascending=True)
    ax.barh(plot_df["aperture"], plot_df["n_snr_ge_threshold"], color="#2c7fb8")
    ax.set_xlabel("Rows with empirical SNR >= threshold")
    ax.set_ylabel("Aperture/product")
    ax.figure.tight_layout()
    path = out_dir / "recoverable_count_by_aperture.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths["recoverable_count_by_aperture"] = str(path)

    if "tmag_bin" in joined:
        table = pd.crosstab(joined["tmag_bin"], joined["best_snr_ge_threshold"], dropna=False)
        table = table.rename(columns={False: "not recoverable", True: "recoverable"})
        frac = table.div(table.sum(axis=1), axis=0).fillna(0.0)
        fig, ax = plt.subplots(figsize=(7, 3.5))
        bottom = np.zeros(len(frac))
        x = np.arange(len(frac))
        for col, color in (("not recoverable", "#bdbdbd"), ("recoverable", "#41ab5d")):
            values = frac[col].to_numpy() if col in frac else np.zeros(len(frac))
            ax.bar(x, values, bottom=bottom, label=col, color=color, width=0.85)
            bottom += values
        ax.set_xticks(x, frac.index.astype(str))
        ax.set_ylim(0, 1)
        ax.set_ylabel("Fraction")
        ax.set_xlabel("Tmag bin")
        ax.legend(fontsize=8)
        ax.figure.tight_layout()
        path = out_dir / "best_recoverable_fraction_by_tmag.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths["best_recoverable_fraction_by_tmag"] = str(path)
    return paths


def summarize(
    survival_csv: Path,
    out_dir: Path,
    current_aperture: str = "DET_FLUX_ADP",
    recoverable_snr: float = 7.0,
    best_group: str = "detrended_all",
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    df = _read_survival(survival_csv)
    aperture_summary = _aperture_summary(df, recoverable_snr)
    aperture_summary.to_csv(out_dir / "aperture_survival_summary.csv", index=False)
    group_summary, best_tables = _group_summary(df, recoverable_snr)
    group_summary.to_csv(out_dir / "aperture_group_survival_summary.csv", index=False)

    if best_group not in best_tables:
        raise ValueError(f"unknown best group {best_group!r}; choices are {sorted(best_tables)}")
    joined = _current_vs_best(df, current_aperture, best_tables[best_group], recoverable_snr)
    joined.to_csv(out_dir / "best_aperture_per_injection.csv", index=False)
    newly = joined[joined["newly_recoverable_by_aperture"]].copy()
    newly.to_csv(out_dir / "newly_recoverable_by_best_aperture.csv", index=False)
    bin_tables = _write_bin_tables(joined, out_dir)

    current_rows = df[df["aperture"].eq(current_aperture)]
    best_rows = best_tables[best_group]
    best_single = aperture_summary.iloc[0].to_dict() if not aperture_summary.empty else {}
    summary: dict[str, Any] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "survival_csv": str(survival_csv),
        "current_aperture": current_aperture,
        "recoverable_snr": recoverable_snr,
        "best_group": best_group,
        "n_injections": int(df["injection_id"].nunique()),
        "current_n_snr_ge_threshold": int((current_rows["delta_signal_multi_snr_mad"] >= recoverable_snr).sum()),
        "current_median_snr": float(current_rows["delta_signal_multi_snr_mad"].median()),
        "best_group_n_snr_ge_threshold": int((best_rows["delta_signal_multi_snr_mad"] >= recoverable_snr).sum()),
        "best_group_median_snr": float(best_rows["delta_signal_multi_snr_mad"].median()),
        "newly_recoverable_by_best_group": int(joined["newly_recoverable_by_aperture"].sum()),
        "best_single_aperture": {str(key): _json_default(value) for key, value in best_single.items()},
        "best_aperture_counts": {
            str(key): int(value)
            for key, value in best_rows["aperture"].value_counts().sort_values(ascending=False).items()
        },
        "bin_tables": bin_tables,
    }
    summary["plot_paths"] = _write_plots(aperture_summary, joined, out_dir)
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--survival-csv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--current-aperture", default="DET_FLUX_ADP")
    parser.add_argument("--recoverable-snr", type=float, default=7.0)
    parser.add_argument("--best-group", choices=sorted(DEFAULT_GROUPS), default="detrended_all")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    summarize(
        args.survival_csv,
        args.out_dir,
        current_aperture=args.current_aperture,
        recoverable_snr=args.recoverable_snr,
        best_group=args.best_group,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
