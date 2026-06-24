#!/usr/bin/env python3
"""Diagnose where injected signals are lost by the BLS recovery step.

This joins the review queue's truth/recovery metadata to the pre-detrend
signal-survival diagnostics. The key diagnostic is empirical signal SNR after
the light-curve product and detrending, not just the injected model depth.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


MODE_ORDER = (
    "bls_top1_recovered",
    "bls_topn_recovered",
    "bls_topn_harmonic_match",
    "bls_peak_mismatch",
)
TOP1_MODES = {"bls_top1_recovered", "bls_recovered"}
EXACT_MODES = TOP1_MODES | {"bls_topn_recovered"}
ANY_MATCH_MODES = EXACT_MODES | {"bls_topn_harmonic_match"}

SNR_BINS = (-np.inf, 1.0, 3.0, 5.0, 7.0, 10.0, 20.0, np.inf)
SNR_LABELS = ("<1", "1-3", "3-5", "5-7", "7-10", "10-20", ">20")
TMAG_BINS = (-np.inf, 16.0, 17.0, 18.0, 19.0, 20.0, np.inf)
TMAG_LABELS = ("<16", "16-17", "17-18", "18-19", "19-20", ">20")
PERIOD_BINS = (0.0, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, np.inf)
PERIOD_LABELS = ("<0.25", "0.25-0.5", "0.5-1", "1-2", "2-5", "5-10", ">10")
DEPTH_BINS = (0.0, 0.01, 0.03, 0.1, 0.3, 1.0, np.inf)
DEPTH_LABELS = ("<1%", "1-3%", "3-10%", "10-30%", "30-100%", ">100%")
RADIUS_BINS = (0.0, 1.0, 2.0, 4.0, 8.0, 11.2, 16.0, 25.0, np.inf)
RADIUS_LABELS = ("<1", "1-2", "2-4", "4-8", "8-11.2", "11.2-16", "16-25", ">25")
DEFAULT_RECOVERABLE_SNR = 7.0


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def _mode_series(df: pd.DataFrame) -> pd.Series:
    if "topn_recovery_status" in df:
        topn = df["topn_recovery_status"].fillna("").astype(str)
        strict = df.get("recovery_status", pd.Series("", index=df.index)).fillna("").astype(str)
        return topn.where(topn.ne(""), strict)
    return df.get("recovery_status", pd.Series("", index=df.index)).fillna("").astype(str)


def _ordered_counts(series: pd.Series) -> dict[str, int]:
    counts = series.value_counts(dropna=False).to_dict()
    ordered: dict[str, int] = {}
    for key in MODE_ORDER:
        if key in counts:
            ordered[key] = int(counts.pop(key))
    for key in sorted(counts):
        ordered[str(key)] = int(counts[key])
    return ordered


def _read_survival(survival_csv: Path, aperture: str) -> pd.DataFrame:
    survival = pd.read_csv(survival_csv)
    if "aperture" in survival:
        survival = survival[survival["aperture"].astype(str).eq(aperture)].copy()
    keep = [
        col
        for col in (
            "injection_id",
            "aperture",
            "delta_signal_multi_snr_mad",
            "delta_signal_snr_mad",
            "depth_retention_frac",
            "depth_retention_vs_model_frac",
            "measured_delta_depth",
            "n_good_ap_in_transit",
            "n_good_ap_oot",
        )
        if col in survival.columns
    ]
    return survival.loc[:, keep].copy()


def _assign_failure_class(df: pd.DataFrame, recoverable_snr: float) -> pd.Series:
    snr = pd.to_numeric(df.get("delta_signal_multi_snr_mad"), errors="coerce")
    mode = df["recovery_mode"].fillna("").astype(str)
    out = pd.Series("unknown_or_missing_snr", index=df.index, dtype=object)
    out = out.where(~df["top1_recovered"], "matched_top1")
    out = out.where(~(df["exact_recovered"] & ~df["top1_recovered"]), "matched_exact_topn")
    out = out.where(~(mode.eq("bls_topn_harmonic_match")), "matched_harmonic_topn")
    unmatched = df["unmatched"].astype(bool)
    out = out.where(~(unmatched & snr.notna() & snr.lt(recoverable_snr)), "low_snr_unmatched")
    out = out.where(~(unmatched & snr.notna() & snr.ge(recoverable_snr)), "snr_qualified_bls_ranking_loss")
    return out


def load_joined_diagnostics(
    queue_csv: Path,
    survival_csv: Path,
    aperture: str,
    recoverable_snr: float = DEFAULT_RECOVERABLE_SNR,
) -> pd.DataFrame:
    queue = pd.read_csv(queue_csv)
    queue["recovery_mode"] = _mode_series(queue)
    survival = _read_survival(survival_csv, aperture)
    joined = queue.merge(survival, on="injection_id", how="left", suffixes=("", "_survival"))

    joined["top1_recovered"] = joined["recovery_mode"].isin(TOP1_MODES)
    joined["exact_recovered"] = joined["recovery_mode"].isin(EXACT_MODES)
    joined["any_period_match"] = joined["recovery_mode"].isin(ANY_MATCH_MODES)
    joined["unmatched"] = ~joined["any_period_match"]

    for col in (
        "delta_signal_multi_snr_mad",
        "delta_signal_snr_mad",
        "depth_retention_frac",
        "truth_period_d",
        "truth_model_depth",
        "truth_depth",
        "truth_radius_rearth",
        "tmag",
    ):
        if col in joined:
            joined[col] = pd.to_numeric(joined[col], errors="coerce")

    if "delta_signal_multi_snr_mad" in joined:
        joined["snr_bin"] = pd.cut(joined["delta_signal_multi_snr_mad"], SNR_BINS, labels=SNR_LABELS)
    if "tmag" in joined:
        joined["tmag_bin"] = pd.cut(joined["tmag"], TMAG_BINS, labels=TMAG_LABELS)
    if "truth_period_d" in joined:
        joined["period_bin"] = pd.cut(joined["truth_period_d"], PERIOD_BINS, labels=PERIOD_LABELS)
    depth_col = "truth_model_depth" if "truth_model_depth" in joined else "truth_depth"
    if depth_col in joined:
        joined["depth_bin"] = pd.cut(joined[depth_col], DEPTH_BINS, labels=DEPTH_LABELS)
    if "truth_radius_rearth" in joined:
        joined["radius_rearth_bin"] = pd.cut(joined["truth_radius_rearth"], RADIUS_BINS, labels=RADIUS_LABELS)
    joined["failure_class"] = _assign_failure_class(joined, recoverable_snr)
    return joined


def summarize_by_bin(df: pd.DataFrame, bin_col: str) -> pd.DataFrame:
    if bin_col not in df:
        return pd.DataFrame()
    grouped = (
        df.groupby(bin_col, observed=False)
        .agg(
            n=("recovery_mode", "size"),
            top1_recovered=("top1_recovered", "sum"),
            exact_recovered=("exact_recovered", "sum"),
            any_period_match=("any_period_match", "sum"),
            unmatched=("unmatched", "sum"),
            median_empirical_snr=("delta_signal_multi_snr_mad", "median"),
            median_tmag=("tmag", "median"),
            median_period_d=("truth_period_d", "median"),
            median_depth=("truth_model_depth", "median"),
            median_radius_rearth=("truth_radius_rearth", "median"),
            median_depth_retention=("depth_retention_frac", "median"),
        )
        .reset_index()
    )
    denom = grouped["n"].replace(0, np.nan)
    for col in ("top1_recovered", "exact_recovered", "any_period_match", "unmatched"):
        grouped[f"{col}_frac"] = grouped[col] / denom
    return grouped


def _mode_stats(df: pd.DataFrame) -> pd.DataFrame:
    stats = (
        df.groupby("recovery_mode", observed=False)
        .agg(
            n=("recovery_mode", "size"),
            median_empirical_snr=("delta_signal_multi_snr_mad", "median"),
            median_tmag=("tmag", "median"),
            median_period_d=("truth_period_d", "median"),
            median_depth=("truth_model_depth", "median"),
            median_radius_rearth=("truth_radius_rearth", "median"),
            median_depth_retention=("depth_retention_frac", "median"),
            median_good_in_transit=("n_good_ap_in_transit", "median"),
        )
        .reset_index()
    )
    stats["_order"] = stats["recovery_mode"].map({mode: idx for idx, mode in enumerate(MODE_ORDER)}).fillna(99)
    return stats.sort_values(["_order", "recovery_mode"]).drop(columns="_order")


def _failure_class_stats(df: pd.DataFrame) -> pd.DataFrame:
    order = {
        "matched_top1": 0,
        "matched_exact_topn": 1,
        "matched_harmonic_topn": 2,
        "low_snr_unmatched": 3,
        "snr_qualified_bls_ranking_loss": 4,
        "unknown_or_missing_snr": 5,
    }
    stats = (
        df.groupby("failure_class", observed=False)
        .agg(
            n=("failure_class", "size"),
            median_empirical_snr=("delta_signal_multi_snr_mad", "median"),
            median_tmag=("tmag", "median"),
            median_period_d=("truth_period_d", "median"),
            median_depth=("truth_model_depth", "median"),
            median_radius_rearth=("truth_radius_rearth", "median"),
            median_depth_retention=("depth_retention_frac", "median"),
            median_good_in_transit=("n_good_ap_in_transit", "median"),
        )
        .reset_index()
    )
    stats["_order"] = stats["failure_class"].map(order).fillna(99)
    return stats.sort_values(["_order", "failure_class"]).drop(columns="_order")


def _write_failure_crosstab(df: pd.DataFrame, index_col: str, out_path: Path) -> pd.DataFrame:
    if index_col not in df:
        table = pd.DataFrame()
    else:
        table = pd.crosstab(df[index_col], df["failure_class"], dropna=False)
    table.to_csv(out_path)
    return table


def _write_plots(df: pd.DataFrame, out_dir: Path) -> dict[str, str]:
    paths: dict[str, str] = {}
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        return {"plot_error": f"{type(exc).__name__}: {exc}"}

    if {"tmag", "delta_signal_multi_snr_mad", "recovery_mode"}.issubset(df.columns):
        colors = {
            "bls_top1_recovered": "#2c7fb8",
            "bls_topn_recovered": "#41ab5d",
            "bls_topn_harmonic_match": "#fdae61",
            "bls_peak_mismatch": "#969696",
            "bls_recovered": "#2c7fb8",
        }
        fig, ax = plt.subplots(figsize=(7.0, 4.5))
        for mode in list(MODE_ORDER) + [m for m in sorted(df["recovery_mode"].dropna().unique()) if m not in MODE_ORDER]:
            sub = df[df["recovery_mode"].eq(mode)]
            if sub.empty:
                continue
            y = sub["delta_signal_multi_snr_mad"].clip(lower=0.03)
            ax.scatter(sub["tmag"], y, s=18, alpha=0.75, label=mode, color=colors.get(mode, "#6a51a3"), linewidths=0)
        ax.axhline(7.0, color="black", linewidth=1.0, linestyle="--", alpha=0.7)
        ax.set_yscale("log")
        ax.set_xlabel("TESS magnitude")
        ax.set_ylabel("Empirical injected-signal SNR after detrending")
        ax.legend(fontsize=7, loc="upper right")
        ax.figure.tight_layout()
        path = out_dir / "bls_recovery_mode_snr_tmag.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths["snr_tmag_scatter"] = str(path)

    fraction_tables = []
    for bin_col, title in (
        ("snr_bin", "Empirical SNR bin"),
        ("tmag_bin", "Tmag bin"),
        ("period_bin", "Period bin [d]"),
        ("radius_rearth_bin", "Radius bin [Earth radii]"),
    ):
        table = summarize_by_bin(df, bin_col)
        if table.empty:
            continue
        fraction_tables.append((table, bin_col, title))
    if fraction_tables:
        fig, axes = plt.subplots(len(fraction_tables), 1, figsize=(7.0, 2.3 * len(fraction_tables)))
        if len(fraction_tables) == 1:
            axes = [axes]
        for ax, (table, bin_col, title) in zip(axes, fraction_tables, strict=True):
            labels = table[bin_col].astype(str)
            x = np.arange(len(labels))
            ax.plot(x, table["top1_recovered_frac"], marker="o", label="top-1 exact")
            ax.plot(x, table["exact_recovered_frac"], marker="o", label="top-N exact")
            ax.plot(x, table["any_period_match_frac"], marker="o", label="top-N incl. harmonic")
            ax.set_xticks(x, labels, rotation=0)
            ax.set_ylim(-0.03, 1.03)
            ax.set_ylabel("Fraction")
            ax.set_title(title, fontsize=10)
            ax.grid(axis="y", alpha=0.25)
        axes[0].legend(fontsize=8, ncol=3, loc="upper left", bbox_to_anchor=(0.0, 1.34))
        axes[-1].set_xlabel("Bin")
        fig.tight_layout()
        path = out_dir / "bls_recovery_fraction_by_bins.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths["fraction_by_bins"] = str(path)

    if {"snr_bin", "tmag_bin", "failure_class"}.issubset(df.columns):
        class_order = [
            "matched_top1",
            "matched_exact_topn",
            "matched_harmonic_topn",
            "low_snr_unmatched",
            "snr_qualified_bls_ranking_loss",
            "unknown_or_missing_snr",
        ]
        colors = {
            "matched_top1": "#2c7fb8",
            "matched_exact_topn": "#41ab5d",
            "matched_harmonic_topn": "#fdae61",
            "low_snr_unmatched": "#bdbdbd",
            "snr_qualified_bls_ranking_loss": "#d95f0e",
            "unknown_or_missing_snr": "#756bb1",
        }
        fig, axes = plt.subplots(2, 1, figsize=(8.0, 5.8))
        for ax, bin_col, title in (
            (axes[0], "snr_bin", "Failure class by empirical SNR"),
            (axes[1], "tmag_bin", "Failure class by Tmag"),
        ):
            table = pd.crosstab(df[bin_col], df["failure_class"], dropna=False)
            plot_cols = [col for col in class_order if col in table.columns]
            plot_cols.extend([col for col in table.columns if col not in plot_cols])
            frac = table[plot_cols].div(table[plot_cols].sum(axis=1), axis=0).fillna(0.0)
            bottom = np.zeros(len(frac))
            x = np.arange(len(frac))
            for col in plot_cols:
                ax.bar(x, frac[col].to_numpy(), bottom=bottom, label=col, color=colors.get(col), width=0.86)
                bottom += frac[col].to_numpy()
            ax.set_xticks(x, frac.index.astype(str))
            ax.set_ylim(0, 1)
            ax.set_ylabel("Fraction")
            ax.set_title(title, fontsize=10)
            ax.grid(axis="y", alpha=0.25)
        axes[-1].set_xlabel("Bin")
        axes[0].legend(fontsize=7, ncol=3, loc="upper left", bbox_to_anchor=(0.0, 1.35))
        fig.tight_layout()
        path = out_dir / "bls_failure_class_by_snr_tmag.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths["failure_class_by_snr_tmag"] = str(path)
    return paths


def diagnose(
    queue_csv: Path,
    survival_csv: Path,
    aperture: str,
    out_dir: Path,
    recoverable_snr: float = DEFAULT_RECOVERABLE_SNR,
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    joined = load_joined_diagnostics(queue_csv, survival_csv, aperture, recoverable_snr=recoverable_snr)
    joined.to_csv(out_dir / "bls_recovery_failure_joined.csv", index=False)

    bin_tables: dict[str, dict[str, Any]] = {}
    for bin_col in ("snr_bin", "tmag_bin", "period_bin", "depth_bin", "radius_rearth_bin"):
        table = summarize_by_bin(joined, bin_col)
        if table.empty:
            continue
        table.to_csv(out_dir / f"recovery_by_{bin_col}.csv", index=False)
        bin_tables[bin_col] = {
            str(row[bin_col]): {
                key: _json_default(value)
                for key, value in row.drop(labels=[bin_col]).dropna().to_dict().items()
            }
            for _, row in table.iterrows()
        }

    mode_stats = _mode_stats(joined)
    mode_stats.to_csv(out_dir / "recovery_mode_failure_stats.csv", index=False)

    failure_stats = _failure_class_stats(joined)
    failure_stats.to_csv(out_dir / "failure_class_stats.csv", index=False)
    failure_counts = joined["failure_class"].value_counts(dropna=False).rename_axis("failure_class").reset_index(name="count")
    failure_counts.to_csv(out_dir / "failure_class_counts.csv", index=False)
    failure_by_snr = _write_failure_crosstab(joined, "snr_bin", out_dir / "failure_class_by_snr_bin.csv")
    failure_by_tmag = _write_failure_crosstab(joined, "tmag_bin", out_dir / "failure_class_by_tmag_bin.csv")
    case_files: dict[str, str] = {}
    case_sort_cols = [
        col
        for col in (
            "delta_signal_multi_snr_mad",
            "tmag",
            "truth_period_d",
            "truth_model_depth",
            "truth_radius_rearth",
            "depth_retention_frac",
        )
        if col in joined.columns
    ]
    for failure_class in ("snr_qualified_bls_ranking_loss", "low_snr_unmatched"):
        cases = joined[joined["failure_class"].eq(failure_class)].copy()
        if not cases.empty and case_sort_cols:
            ascending = [False if col == "delta_signal_multi_snr_mad" else True for col in case_sort_cols]
            cases = cases.sort_values(case_sort_cols, ascending=ascending)
        path = out_dir / f"{failure_class}_cases.csv"
        cases.to_csv(path, index=False)
        case_files[failure_class] = str(path)

    tmag_snr_table = pd.crosstab(
        joined["tmag_bin"],
        joined["snr_bin"],
        values=joined["any_period_match"].astype(int),
        aggfunc="mean",
        dropna=False,
    )
    tmag_snr_table.to_csv(out_dir / "any_match_fraction_tmag_by_snr.csv")

    count_table = pd.crosstab(joined["tmag_bin"], joined["snr_bin"], dropna=False)
    count_table.to_csv(out_dir / "counts_tmag_by_snr.csv")

    summary: dict[str, Any] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "queue_csv": str(queue_csv),
        "survival_csv": str(survival_csv),
        "aperture": aperture,
        "recoverable_snr": float(recoverable_snr),
        "n_rows": int(len(joined)),
        "n_missing_survival_snr": int(joined["delta_signal_multi_snr_mad"].isna().sum())
        if "delta_signal_multi_snr_mad" in joined
        else int(len(joined)),
        "recovery_mode_counts": _ordered_counts(joined["recovery_mode"]),
        "top1_recovered_fraction": float(joined["top1_recovered"].mean()),
        "exact_recovered_fraction": float(joined["exact_recovered"].mean()),
        "any_period_match_fraction": float(joined["any_period_match"].mean()),
        "median_empirical_snr": float(joined["delta_signal_multi_snr_mad"].median()),
        "median_empirical_snr_by_mode": {
            str(row["recovery_mode"]): _json_default(row["median_empirical_snr"])
            for _, row in mode_stats.iterrows()
        },
        "median_tmag_by_mode": {
            str(row["recovery_mode"]): _json_default(row["median_tmag"])
            for _, row in mode_stats.iterrows()
        },
        "failure_class_counts": {
            str(row["failure_class"]): int(row["count"])
            for _, row in failure_counts.iterrows()
        },
        "failure_class_stats": {
            str(row["failure_class"]): {
                key: _json_default(value)
                for key, value in row.drop(labels=["failure_class"]).dropna().to_dict().items()
            }
            for _, row in failure_stats.iterrows()
        },
        "failure_class_by_snr_bin": {
            str(idx): {str(col): int(value) for col, value in row.items()}
            for idx, row in failure_by_snr.iterrows()
        },
        "failure_class_by_tmag_bin": {
            str(idx): {str(col): int(value) for col, value in row.items()}
            for idx, row in failure_by_tmag.iterrows()
        },
        "case_files": case_files,
        "bin_tables": bin_tables,
    }
    summary["plot_paths"] = _write_plots(joined, out_dir)
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue-csv", type=Path, required=True)
    parser.add_argument("--survival-csv", type=Path, required=True)
    parser.add_argument("--aperture", default="DET_FLUX_ADP")
    parser.add_argument(
        "--recoverable-snr",
        type=float,
        default=DEFAULT_RECOVERABLE_SNR,
        help="Empirical post-detrend signal SNR above which unmatched rows are counted as search/ranking losses.",
    )
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    diagnose(args.queue_csv, args.survival_csv, args.aperture, args.out_dir, recoverable_snr=args.recoverable_snr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
