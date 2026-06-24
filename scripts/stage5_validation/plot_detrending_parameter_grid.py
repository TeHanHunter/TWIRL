#!/usr/bin/env python3
"""Plot detrending parameter-grid performance from a sweep detail CSV."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.plotting.style import apply_twirl_style  # noqa: E402


DEFAULT_SWEEP_DIR = (
    REPO_ROOT
    / "reports/stage5_validation/s56_1k_predetrend_detrending_method_sweep_finegrid_pdo"
)


def summarize(detail: pd.DataFrame, *, snr_threshold: float) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for method, group in detail.groupby("method", dropna=False):
        snr = pd.to_numeric(group["delta_signal_multi_snr_mad"], errors="coerce")
        trend = pd.to_numeric(group["det_trend_ptp"], errors="coerce")
        retention = pd.to_numeric(group["depth_retention_frac"], errors="coerce")
        rows.append(
            {
                "method": method,
                "method_kind": str(group["method_kind"].dropna().iloc[0]) if group["method_kind"].notna().any() else "",
                "median_bkspace_d": float(np.nanmedian(pd.to_numeric(group.get("bkspace_d"), errors="coerce"))),
                "median_window_d": float(np.nanmedian(pd.to_numeric(group.get("window_d"), errors="coerce"))),
                "n": int(len(group)),
                "n_snr_ge_threshold": int((snr >= snr_threshold).sum()),
                "frac_snr_ge_threshold": float((snr >= snr_threshold).mean()),
                "median_multi_snr": float(np.nanmedian(snr)),
                "median_depth_retention": float(np.nanmedian(retention)),
                "median_det_trend_ptp": float(np.nanmedian(trend)),
                "median_trend_reduction_frac": float(
                    np.nanmedian(pd.to_numeric(group["trend_reduction_frac"], errors="coerce"))
                ),
            }
        )
    return pd.DataFrame(rows).sort_values(["method_kind", "median_window_d", "median_bkspace_d", "method"])


def plot_grid(summary: pd.DataFrame, out_dir: Path) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    apply_twirl_style("full_page")
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.25), sharey=True)
    configs = [
        ("median_filter", "median_window_d", "Rolling-median window [days]"),
        ("spline", "median_bkspace_d", "Quantile-spline knot spacing [days]"),
    ]
    for ax, (kind, xcol, xlabel) in zip(axes, configs):
        sub = summary[summary["method_kind"].eq(kind) & np.isfinite(summary[xcol])].copy()
        if sub.empty:
            ax.set_visible(False)
            continue
        sub = sub.sort_values(xcol)
        ax.plot(
            sub[xcol],
            sub["n_snr_ge_threshold"],
            marker="o",
            color="#2f6f9f",
            label="SNR-qualified rows",
        )
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Rows with empirical signal SNR >= 7")
        ax.grid(True, alpha=0.3)
        ax2 = ax.twinx()
        ax2.plot(
            sub[xcol],
            sub["median_det_trend_ptp"],
            marker="s",
            color="#b85c38",
            label="Residual trend",
        )
        ax2.set_ylabel("Residual trend 90-10% range")
        best = sub.sort_values(["n_snr_ge_threshold", "median_det_trend_ptp"], ascending=[False, True]).iloc[0]
        ax.axvline(best[xcol], color="0.25", linestyle=":", linewidth=1.0)
        ax.text(
            best[xcol],
            float(best["n_snr_ge_threshold"]),
            str(best["method"]),
            fontsize=7,
            ha="left",
            va="bottom",
            rotation=25,
        )
        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines + lines2, labels + labels2, loc="best", fontsize=7)
    fig.suptitle("S56 pre-detrend detrending fine-grid sweep")
    fig.tight_layout()
    png = out_dir / "detrending_finegrid_parameter_performance.png"
    pdf = out_dir / "detrending_finegrid_parameter_performance.pdf"
    fig.savefig(png, dpi=220)
    fig.savefig(pdf)
    plt.close(fig)
    return {"finegrid_parameter_png": str(png), "finegrid_parameter_pdf": str(pdf)}


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sweep-dir", type=Path, default=DEFAULT_SWEEP_DIR)
    parser.add_argument("--snr-threshold", type=float, default=7.0)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    detail_path = args.sweep_dir / "detrending_method_detail.csv"
    detail = pd.read_csv(detail_path)
    summary = summarize(detail, snr_threshold=args.snr_threshold)
    out_csv = args.sweep_dir / "finegrid_parameter_summary.csv"
    summary.to_csv(out_csv, index=False)
    paths = plot_grid(summary, args.sweep_dir)
    print({"summary_csv": str(out_csv), **paths})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
