#!/usr/bin/env python3
"""Smoke-test WD-motivated LEO pass-through threshold presets.

This does not replace LEO-Vetter's threshold implementation. It audits the
existing LEO metric table and asks which BLS-recovered injected rows would be
rescued by physically motivated WD review rules.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
import json
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


DEFAULT_QUEUE_CSV = (
    REPO_ROOT
    / "reports/stage5_validation/s56_10k_predetrend_small_pair_200k_review_queue_pdo/"
    / "review_queue.csv"
)
DEFAULT_METRICS_CSV = (
    REPO_ROOT
    / "reports/stage5_validation/s56_10k_predetrend_small_pair_200k_review_queue_pdo/"
    / "leo_metrics.csv"
)
DEFAULT_OUT_DIR = (
    REPO_ROOT
    / "reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/"
    / "leo_wd_tuning_smoke"
)


@dataclass(frozen=True)
class Preset:
    name: str
    description: str
    new_mes_min: float
    shp_max: float
    primary_margin_min: float
    n_in_min: int
    new_n_transit_min: int = 2
    sig_pri_min: float = 5.0
    fred_max: float = 3.0
    dmm_min: float = 0.5
    dmm_max: float = 2.5
    rp_max_rearth: float = 22.4


PRESETS = (
    Preset(
        name="wd_review_high_purity",
        description=(
            "Under-resolved WD review pass: keep two or more surviving transits, "
            "allow modest MES, but require a compact shape statistic and a primary "
            "event at least 4 sigma stronger than secondary/tertiary/positive events."
        ),
        new_mes_min=5.0,
        shp_max=0.30,
        primary_margin_min=4.0,
        n_in_min=6,
    ),
    Preset(
        name="wd_review_balanced",
        description=(
            "Slightly looser shape pass for 200 s WD transits whose trapezoid/planet "
            "fits are cadence-limited; still requires primary dominance."
        ),
        new_mes_min=5.0,
        shp_max=0.40,
        primary_margin_min=4.0,
        n_in_min=6,
    ),
    Preset(
        name="wd_review_aggressive",
        description=(
            "Aggressive human-review catchment for exploratory triage; useful for "
            "estimating the purity cost of relaxing both the inherited stock LEO "
            "shape cut and the primary-dominance guard."
        ),
        new_mes_min=5.0,
        shp_max=0.45,
        primary_margin_min=3.0,
        n_in_min=6,
    ),
)


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _bool_series(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    text = series.astype(str).str.strip().str.lower()
    return text.isin({"true", "1", "yes", "y"})


def classify_recovery(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["strict_top1_recovered"] = out["recovery_status"].astype(str).eq("bls_recovered")
    out["topn_exact_recovered"] = out["topn_recovery_status"].astype(str).eq("bls_topn_recovered")
    out["topn_harmonic_recovered"] = out["topn_recovery_status"].astype(str).eq("bls_topn_harmonic_match")
    for col in [c for c in out.columns if c.startswith("topn_exact_recovered_")]:
        out["topn_exact_recovered"] |= _bool_series(out[col])
    for col in [c for c in out.columns if c.startswith("topn_harmonic_match_")]:
        out["topn_harmonic_recovered"] |= _bool_series(out[col])
    out["bls_recovered"] = (
        out["strict_top1_recovered"] | out["topn_exact_recovered"] | out["topn_harmonic_recovered"]
    )
    return out


def load_metric_frame(queue_csv: Path, metrics_csv: Path) -> pd.DataFrame:
    queue = classify_recovery(pd.read_csv(queue_csv))
    if "source_kind" in queue:
        queue = queue[queue["source_kind"].astype(str).eq("injection_recovery")].copy()
    metrics = pd.read_csv(metrics_csv)
    if "source_kind" in metrics:
        metrics = metrics[metrics["source_kind"].astype(str).eq("injection_recovery")].copy()
    merged = queue.merge(metrics.drop_duplicates("review_id"), on="review_id", how="inner", suffixes=("", "_metric"))
    leo_class = merged.get("leo_class_metric", merged.get("leo_class", pd.Series("", index=merged.index)))
    merged["leo_class_effective"] = leo_class.fillna(merged.get("leo_class", "")).astype(str)
    merged["current_leo_pc_or_fp"] = merged["leo_class_effective"].isin(["PC", "FP"])
    numeric_cols = [
        "MES",
        "new_MES",
        "SHP",
        "N_transit",
        "new_N_transit",
        "n_in",
        "sig_pri",
        "sig_sec",
        "sig_ter",
        "sig_pos",
        "DMM",
        "Fred",
        "Rp",
        "sde_max",
        "depth_snr",
        "truth_period_d",
        "truth_duration_min",
        "truth_radius_rearth",
        "tmag",
    ]
    for col in numeric_cols:
        if col in merged:
            merged[col] = pd.to_numeric(merged[col], errors="coerce")
    alt_cols = [col for col in ("sig_sec", "sig_ter", "sig_pos") if col in merged]
    if alt_cols:
        merged["max_alternate_sig"] = np.nanmax(
            np.vstack([merged[col].fillna(-np.inf).to_numpy(dtype=float) for col in alt_cols]),
            axis=0,
        )
    else:
        merged["max_alternate_sig"] = np.nan
    merged["primary_margin"] = merged["sig_pri"] - merged["max_alternate_sig"]
    return merged


def preset_mask(df: pd.DataFrame, preset: Preset) -> pd.Series:
    size_ok = df["Rp"].isna() | df["Rp"].le(preset.rp_max_rearth)
    dmm_ok = df["DMM"].isna() | df["DMM"].between(preset.dmm_min, preset.dmm_max)
    fred_ok = df["Fred"].isna() | df["Fred"].lt(preset.fred_max)
    return (
        df["new_MES"].ge(preset.new_mes_min)
        & df["SHP"].le(preset.shp_max)
        & df["new_N_transit"].ge(preset.new_n_transit_min)
        & df["n_in"].ge(preset.n_in_min)
        & df["sig_pri"].ge(preset.sig_pri_min)
        & df["primary_margin"].ge(preset.primary_margin_min)
        & dmm_ok
        & fred_ok
        & size_ok
    ).fillna(False)


def score_prediction(df: pd.DataFrame, pred: pd.Series, *, label: str) -> dict[str, Any]:
    y = df["bls_recovered"].fillna(False).astype(bool)
    pred = pred.fillna(False).astype(bool)
    current = df["current_leo_pc_or_fp"].fillna(False).astype(bool)
    tp = int((pred & y).sum())
    fp = int((pred & ~y).sum())
    fn = int((~pred & y).sum())
    tn = int((~pred & ~y).sum())
    return {
        "preset": label,
        "n": int(len(df)),
        "predicted_pass_n": int(pred.sum()),
        "tp_vs_bls": tp,
        "fp_vs_bls": fp,
        "fn_vs_bls": fn,
        "tn_vs_bls": tn,
        "precision_vs_bls": tp / (tp + fp) if tp + fp else float("nan"),
        "recall_vs_bls": tp / (tp + fn) if tp + fn else float("nan"),
        "specificity_vs_bls": tn / (tn + fp) if tn + fp else float("nan"),
        "added_tp_vs_current": int((pred & ~current & y).sum()),
        "added_fp_vs_current": int((pred & ~current & ~y).sum()),
        "added_review_n": int((pred & ~current).sum()),
    }


def evaluate_presets(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = [score_prediction(df, df["current_leo_pc_or_fp"], label="current_leo_pc_or_fp")]
    relabel_pieces: list[pd.DataFrame] = []
    current = df["current_leo_pc_or_fp"].fillna(False).astype(bool)
    for preset in PRESETS:
        candidate = preset_mask(df, preset)
        pred = current | candidate
        rows.append(score_prediction(df, pred, label=preset.name) | asdict(preset))
        added = df[pred & ~current].copy()
        if not added.empty:
            added["preset"] = preset.name
            relabel_pieces.append(
                added[
                    [
                        "preset",
                        "review_id",
                        "tic",
                        "bls_recovered",
                        "leo_class_effective",
                        "new_MES",
                        "MES",
                        "SHP",
                        "new_N_transit",
                        "n_in",
                        "sig_pri",
                        "max_alternate_sig",
                        "primary_margin",
                        "truth_period_d",
                        "truth_duration_min",
                        "truth_radius_rearth",
                        "tmag",
                    ]
                ]
            )
    relabels = pd.concat(relabel_pieces, ignore_index=True) if relabel_pieces else pd.DataFrame()
    return pd.DataFrame(rows), relabels


def plot_scores(scores: pd.DataFrame, out_dir: Path) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    apply_twirl_style("full_page")
    labels = scores["preset"].tolist()
    x = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    ax.bar(x - width / 2, scores["precision_vs_bls"], width=width, color="#2563eb", label="precision")
    ax.bar(x + width / 2, scores["recall_vs_bls"], width=width, color="#d97706", label="recall")
    for idx, row in scores.iterrows():
        ax.text(
            idx,
            min(1.04, float(row["recall_vs_bls"]) + 0.035),
            f"+{int(row['added_tp_vs_current'])} TP\n+{int(row['added_fp_vs_current'])} FP",
            ha="center",
            va="bottom",
            fontsize=7,
        )
    ax.set_ylim(0, 1.12)
    ax.set_ylabel("Metric relative to BLS exact/top-N/harmonic recovery")
    ax.set_xticks(x, labels=labels, rotation=18, ha="right")
    ax.legend(loc="upper right", frameon=True)
    ax.grid(True, axis="y", color="0.90", linewidth=0.7)
    fig.subplots_adjust(bottom=0.27, left=0.10, right=0.98, top=0.95)
    fig.text(
        0.5,
        0.03,
        "These are smoke-test pass-through presets for human review, not automatic planet-candidate labels.",
        ha="center",
        va="bottom",
        fontsize=8,
    )
    png = out_dir / "leo_wd_threshold_smoke_precision_recall.png"
    pdf = out_dir / "leo_wd_threshold_smoke_precision_recall.pdf"
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return {"precision_recall_png": str(png), "precision_recall_pdf": str(pdf)}


def markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    cols = [
        "preset",
        "predicted_pass_n",
        "tp_vs_bls",
        "fp_vs_bls",
        "precision_vs_bls",
        "recall_vs_bls",
        "added_tp_vs_current",
        "added_fp_vs_current",
    ]
    table = df[cols].copy()
    for col in ("precision_vs_bls", "recall_vs_bls"):
        table[col] = table[col].map(lambda value: f"{value:.1%}" if np.isfinite(value) else "")
    headers = table.columns.tolist()
    rows = table.astype(str).values.tolist()
    widths = [max(len(headers[i]), *(len(row[i]) for row in rows)) for i in range(len(headers))]
    lines = [
        "| " + " | ".join(headers[i].ljust(widths[i]) for i in range(len(headers))) + " |",
        "| " + " | ".join("-" * widths[i] for i in range(len(headers))) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row[i].ljust(widths[i]) for i in range(len(row))) + " |")
    return "\n".join(lines)


def write_summary(
    out_dir: Path,
    *,
    queue_csv: Path,
    metrics_csv: Path,
    scores: pd.DataFrame,
    paths: dict[str, str],
) -> None:
    lines = [
        "# LEO WD Threshold Smoke Test",
        "",
        f"Queue: `{queue_csv}`",
        f"Metrics: `{metrics_csv}`",
        "",
        "## Result",
        "",
        markdown_table(scores),
        "",
        "## Interpretation",
        "",
        "The current WD-tuned LEO PC/FP label is very high precision but low recall relative to the BLS exact/top-N/harmonic recovery label. A physically motivated WD review pass-through can recover additional BLS-tracked injections, but the purity cost rises quickly if the inherited shape tests are loosened too far.",
        "",
        "The most defensible near-term tune is `wd_review_high_purity`: keep `new_N_transit >= 2`, keep the large-object ceiling, require the primary modshift event to dominate secondary/tertiary/positive events, and only relax the under-resolved shape/MES behavior enough to route cases to human review. Do not promote these directly to PC without human labels.",
        "",
        "Physical rationale: WD transits last only minutes at 200 s cadence, so trapezoid/planet-shape diagnostics and pruned MES can be unstable when only a few cadences sample each event. Secondary/tertiary modshift peaks, however, remain physically meaningful for rejecting wrong ephemerides, so the tune should use primary dominance rather than dropping uniqueness checks.",
        "",
        "## Artifacts",
        "",
    ]
    lines.extend(f"- `{key}`: `{path}`" for key, path in paths.items())
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue-csv", type=Path, default=DEFAULT_QUEUE_CSV)
    parser.add_argument("--metrics-csv", type=Path, default=DEFAULT_METRICS_CSV)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    df = load_metric_frame(args.queue_csv, args.metrics_csv)
    scores, relabels = evaluate_presets(df)
    score_csv = args.out_dir / "leo_wd_threshold_smoke_scores.csv"
    relabel_csv = args.out_dir / "leo_wd_threshold_smoke_added_review_rows.csv"
    scores.to_csv(score_csv, index=False)
    relabels.to_csv(relabel_csv, index=False)
    paths = {
        "scores_csv": str(score_csv),
        "added_review_rows_csv": str(relabel_csv),
    }
    paths.update(plot_scores(scores, args.out_dir))
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "queue_csv": str(args.queue_csv),
        "metrics_csv": str(args.metrics_csv),
        "n_injection_rows": int(len(df)),
        "presets": [asdict(preset) for preset in PRESETS],
        "scores": scores.to_dict(orient="records"),
        "paths": paths,
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    write_summary(args.out_dir, queue_csv=args.queue_csv, metrics_csv=args.metrics_csv, scores=scores, paths=paths)
    print(json.dumps({"out_dir": str(args.out_dir), "scores": scores.to_dict(orient="records")}, indent=2, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
