#!/usr/bin/env python3
"""Visualize the LEO-Vetter (WD-tuned) run on the S56 top-50 candidates.

Reads ``leo_metrics_top50.parquet`` and produces a single full-page figure:
  * Top-left:    scatter of fitted P vs SDE_bls, colored by class (PC/FA/FP)
                 with WD 1856 marked
  * Top-right:   per-rank class bar (rank 1 - 50 as horizontal stripes)
  * Bottom row:  metric distributions (MES, SHP, sine_sig, transit_RpRs)
                 with WD 1856 marked as a vertical line
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.plotting.style import BASE_RC, PLOT_TEMPLATES  # noqa: E402

WD_1856_TIC = 267574918
DATA_ROOT = Path("/Users/tehan/PycharmProjects/TWIRL")
BENCHMARK_DIR = DATA_ROOT / "benchmark" / "leo_vetter_s56_top50"
SRC_PARQUET = BENCHMARK_DIR / "leo_metrics.parquet"
OUT_FIG = BENCHMARK_DIR / "summary.pdf"

# Class colors
CMAP = {
    "PC": "#1f77b4",   # blue
    "FA": "#888888",   # gray
    "FP": "#d62728",   # red
    "ERR": "#cccccc",  # light gray
}


def classify(row: pd.Series) -> str:
    if row.get("error") and str(row["error"]) != "None" and str(row["error"]) != "nan":
        return "ERR"
    if bool(row.get("leo_PC", False)):
        return "PC"
    if bool(row.get("leo_FA", False)):
        return "FA"
    if bool(row.get("leo_FP", False)):
        return "FP"
    return "ERR"


def main() -> int:
    df = pd.read_parquet(SRC_PARQUET)
    df["class"] = df.apply(classify, axis=1)

    plt.rcParams.update(BASE_RC)
    tpl = PLOT_TEMPLATES["full_page"]
    fig = plt.figure(figsize=(8.5, 9.5))
    gs = fig.add_gridspec(
        4, 2,
        height_ratios=[2.4, 0.9, 2.0, 2.0],
        hspace=0.55, wspace=0.30,
        left=0.085, right=0.97, top=0.965, bottom=0.06,
    )

    # ---- Top-left: P vs SDE_bls scatter ----
    ax = fig.add_subplot(gs[0, 0])
    for cls in ["FA", "FP", "PC"]:
        sub = df[df["class"] == cls]
        if not len(sub):
            continue
        ax.scatter(
            sub["P_bls"], sub["sde_bls"],
            s=24, c=CMAP[cls], alpha=0.85,
            edgecolor="white", linewidth=0.5,
            label=f"{cls} (n={len(sub)})",
            zorder=3,
        )
    wd = df[df["tic"] == WD_1856_TIC]
    if len(wd):
        r = wd.iloc[0]
        ax.scatter(
            r["P_bls"], r["sde_bls"],
            marker="*", s=240, c="#fbbf24",
            edgecolor="black", linewidth=0.8,
            zorder=5, label="WD 1856",
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("BLS fitted period (d)", fontsize=tpl["label_size"])
    ax.set_ylabel("BLS SDE", fontsize=tpl["label_size"])
    ax.set_title("S56 top-50 — LEO-Vetter (WD-tuned) class",
                 fontsize=tpl["title_size"])
    ax.grid(True, which="both", alpha=0.3, linewidth=tpl["grid_linewidth"])
    ax.tick_params(labelsize=tpl["tick_size"])
    ax.legend(loc="upper right", fontsize=tpl["legend_size"], framealpha=0.9)

    # ---- Top-right: per-rank class strip ----
    ax = fig.add_subplot(gs[0, 1])
    for i, r in df.iterrows():
        ax.barh(
            -i, 1, height=0.85,
            color=CMAP[r["class"]],
            edgecolor="white", linewidth=0.3,
        )
        if int(r["tic"]) == WD_1856_TIC:
            ax.text(
                1.05, -i, f"WD 1856 (rank {int(r['rank_in'])})",
                fontsize=tpl["annotation_size"], va="center",
                color="black", fontweight="bold",
            )
        elif r["class"] == "PC":
            ax.text(
                1.05, -i, f"TIC {int(r['tic'])} (rank {int(r['rank_in'])})",
                fontsize=tpl["annotation_size"], va="center",
                color=CMAP["PC"],
            )
    ax.set_xlim(0, 4.0)
    ax.set_ylim(-len(df) + 0.5, 0.5)
    ax.set_yticks([-i for i in range(0, len(df), 5)])
    ax.set_yticklabels([f"#{i+1}" for i in range(0, len(df), 5)],
                       fontsize=tpl["tick_size"])
    ax.set_xticks([])
    ax.set_title("Per-rank class (top → bottom)", fontsize=tpl["title_size"])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    legend_handles = [
        Line2D([0], [0], marker="s", color="w", markerfacecolor=CMAP[c],
               markersize=8, label=f"{c}")
        for c in ["PC", "FA", "FP", "ERR"]
    ]
    ax.legend(handles=legend_handles, loc="lower right",
              fontsize=tpl["legend_size"], framealpha=0.9)

    # ---- Middle row: summary text ----
    ax = fig.add_subplot(gs[1, :])
    ax.axis("off")
    n_pc = int((df["class"] == "PC").sum())
    n_fa = int((df["class"] == "FA").sum())
    n_fp = int((df["class"] == "FP").sum())
    n_err = int((df["class"] == "ERR").sum())
    pcs = df[df["class"] == "PC"][["tic", "rank_in", "P_bls", "dur_bls_min", "sde_bls", "tmag"]]
    text = (
        f"$\\bf{{Class\\ tally}}$  (n=50):  PC={n_pc}    FA={n_fa}    FP={n_fp}    ERR={n_err}\n\n"
        f"$\\bf{{Planet\\ candidates}}$:\n"
    )
    for _, r in pcs.iterrows():
        flag = "  ← WD 1856" if int(r["tic"]) == WD_1856_TIC else ""
        text += (
            f"   rank {int(r['rank_in'])}: TIC {int(r['tic'])}  "
            f"T={r['tmag']:5.2f}  P={r['P_bls']:7.4f} d  dur={r['dur_bls_min']:4.1f} min  "
            f"SDE={r['sde_bls']:7.2f}{flag}\n"
        )
    ax.text(0.01, 0.98, text, fontsize=tpl["annotation_size"], family="monospace",
            verticalalignment="top", transform=ax.transAxes)

    # ---- Bottom: metric histograms ----
    metric_specs = [
        ("MES", "LEO MES (multi-event)", None, 6.2),
        ("SHP", "Shape parameter", (0, 1), 0.6),
        ("sine_sig", "Sinusoidal significance", None, 15),
        ("transit_RpRs", "Transit fit R_p / R_*", None, None),
    ]
    for i, (col, label, xlim, threshold) in enumerate(metric_specs):
        ax = fig.add_subplot(gs[2 + i // 2, i % 2])
        # Skip NaN/missing
        m = df[col].astype(float)
        m_clean = m[np.isfinite(m)]
        # Split by class for stacked hist
        for cls in ["FA", "FP", "PC"]:
            sub = df[df["class"] == cls][col].astype(float)
            sub = sub[np.isfinite(sub)]
            if len(sub):
                ax.hist(
                    sub, bins=20, alpha=0.7,
                    color=CMAP[cls], label=f"{cls}",
                    edgecolor="white", linewidth=0.4,
                )
        # WD 1856 marker
        if len(wd) and col in wd.columns:
            v = float(wd.iloc[0][col])
            if np.isfinite(v):
                ax.axvline(v, color="#fbbf24", linewidth=1.5,
                           linestyle="--", zorder=5,
                           label=f"WD 1856 = {v:.2f}")
        if threshold is not None:
            ax.axvline(threshold, color="black", linewidth=0.7,
                       linestyle=":", alpha=0.6,
                       label=f"FGKM thresh = {threshold}")
        if xlim:
            ax.set_xlim(xlim)
        ax.set_xlabel(label, fontsize=tpl["label_size"])
        ax.set_ylabel("Count", fontsize=tpl["label_size"])
        ax.tick_params(labelsize=tpl["tick_size"])
        ax.grid(True, alpha=0.3, linewidth=tpl["grid_linewidth"])
        ax.legend(fontsize=tpl["legend_size"], framealpha=0.9)

    fig.savefig(OUT_FIG, bbox_inches="tight")
    fig.savefig(OUT_FIG.with_suffix(".png"), bbox_inches="tight", dpi=160)
    print(f"wrote {OUT_FIG}")
    print(f"wrote {OUT_FIG.with_suffix('.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
