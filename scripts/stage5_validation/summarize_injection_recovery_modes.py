#!/usr/bin/env python3
"""Summarize injected-row BLS recovery modes from a review queue.

The review queue keeps ``recovery_status`` as strict top-peak recovery. Newer
queues also carry ``topn_recovery_status`` so we can split misses into exact
top-N recoveries, harmonic/alias matches, and unmatched cases. When a signal
survival CSV is supplied, the script joins empirical injected-minus-original
SNR and writes mode-by-SNR tables.
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
SNR_BINS = (-np.inf, 1.0, 3.0, 5.0, 7.0, 10.0, 20.0, np.inf)
SNR_LABELS = ("<1", "1-3", "3-5", "5-7", "7-10", "10-20", ">20")


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
        mode = df["topn_recovery_status"].fillna("").astype(str)
        mode = mode.where(mode.ne(""), df["recovery_status"].fillna("").astype(str))
        return mode
    return df["recovery_status"].fillna("").astype(str)


def _ordered_counts(series: pd.Series) -> dict[str, int]:
    counts = series.value_counts(dropna=False).to_dict()
    ordered: dict[str, int] = {}
    for key in MODE_ORDER:
        if key in counts:
            ordered[key] = int(counts.pop(key))
    for key in sorted(counts):
        ordered[str(key)] = int(counts[key])
    return ordered


def _join_survival(queue: pd.DataFrame, survival_csv: Path | None, aperture: str) -> pd.DataFrame:
    if survival_csv is None:
        return queue
    survival = pd.read_csv(survival_csv)
    if "aperture" in survival:
        survival = survival[survival["aperture"].astype(str).eq(aperture)].copy()
    keep = [
        col
        for col in (
            "injection_id",
            "delta_signal_multi_snr_mad",
            "delta_signal_snr_mad",
            "depth_retention_frac",
            "truth_model_depth",
            "truth_period_d",
            "tmag",
        )
        if col in survival.columns
    ]
    return queue.merge(survival[keep], on="injection_id", how="left", suffixes=("", "_survival"))


def summarize(queue_csv: Path, survival_csv: Path | None, aperture: str, out_dir: Path) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    queue = pd.read_csv(queue_csv)
    queue["recovery_mode"] = _mode_series(queue)
    joined = _join_survival(queue, survival_csv, aperture)

    summary: dict[str, Any] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "queue_csv": str(queue_csv),
        "survival_csv": str(survival_csv) if survival_csv else "",
        "aperture": aperture,
        "n_rows": int(len(queue)),
        "strict_recovery_status_counts": _ordered_counts(queue["recovery_status"].fillna("").astype(str)),
        "recovery_mode_counts": _ordered_counts(queue["recovery_mode"]),
    }
    for col in ("n_apertures_agree", "n_topn_apertures_agree", "n_harmonic_apertures_agree"):
        if col in queue:
            summary[f"{col}_counts"] = {
                str(key): int(value)
                for key, value in queue[col].value_counts(dropna=False).sort_index().items()
            }

    mode_counts = queue["recovery_mode"].value_counts(dropna=False).rename_axis("recovery_mode").reset_index(name="count")
    mode_counts.to_csv(out_dir / "recovery_mode_counts.csv", index=False)

    if "delta_signal_multi_snr_mad" in joined:
        joined["snr_bin"] = pd.cut(joined["delta_signal_multi_snr_mad"], SNR_BINS, labels=SNR_LABELS)
        mode_by_snr = pd.crosstab(joined["snr_bin"], joined["recovery_mode"], dropna=False)
        mode_by_snr.to_csv(out_dir / "recovery_mode_by_snr.csv")
        summary["snr_bin_counts"] = {
            str(idx): {str(col): int(value) for col, value in row.items()}
            for idx, row in mode_by_snr.iterrows()
        }
        group_stats = (
            joined.groupby("recovery_mode", observed=False)
            .agg(
                n=("recovery_mode", "size"),
                median_snr=("delta_signal_multi_snr_mad", "median"),
                median_tmag=("tmag", "median"),
                median_period_d=("truth_period_d", "median"),
                median_depth_retention=("depth_retention_frac", "median"),
            )
            .reset_index()
        )
        group_stats.to_csv(out_dir / "recovery_mode_stats.csv", index=False)

        try:
            import matplotlib.pyplot as plt

            plot_cols = [col for col in MODE_ORDER if col in mode_by_snr.columns]
            plot_cols.extend([col for col in mode_by_snr.columns if col not in plot_cols])
            frac = mode_by_snr[plot_cols].div(mode_by_snr[plot_cols].sum(axis=1), axis=0).fillna(0.0)
            ax = frac.plot(kind="bar", stacked=True, figsize=(8, 4), width=0.86)
            ax.set_xlabel("Empirical multi-transit SNR after detrending")
            ax.set_ylabel("Fraction of injected rows")
            ax.legend(title="Recovery mode", fontsize=8, title_fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1.0))
            ax.figure.tight_layout()
            ax.figure.savefig(out_dir / "recovery_mode_by_snr.png", dpi=160)
            plt.close(ax.figure)
        except Exception as exc:
            summary["plot_error"] = f"{type(exc).__name__}: {exc}"

    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue-csv", type=Path, required=True)
    parser.add_argument("--survival-csv", type=Path, default=None)
    parser.add_argument("--aperture", default="DET_FLUX_ADP")
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    summary = summarize(args.queue_csv, args.survival_csv, args.aperture, args.out_dir)
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
