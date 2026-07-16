#!/usr/bin/env python3
"""Render S56 Teacher-v2 recovery as a function of real-TIC review load."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from twirl.plotting.style import apply_twirl_style, get_ordered_palette


ROOT = Path(__file__).resolve().parent
INPUT_ROOT = ROOT / "workload_inputs"
OUT_ROOT = ROOT / "workload_performance"


def _as_bool(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False).astype(bool)
    return values.fillna("").astype(str).str.lower().isin({"1", "true", "t", "yes", "y"})


def _workload_curve(
    *,
    real_path: Path,
    holdout_path: Path,
    score_column: str,
    fractions: np.ndarray,
) -> tuple[pd.DataFrame, dict[str, float | int]]:
    real = pd.read_parquet(real_path, columns=["tic", score_column])
    real[score_column] = pd.to_numeric(real[score_column], errors="coerce")
    tic_scores = (
        real.dropna(subset=["tic", score_column])
        .groupby("tic", sort=False)[score_column]
        .max()
        .sort_values(ascending=False)
    )

    holdout = pd.read_parquet(
        holdout_path,
        columns=["injection_id", "is_injected_signal_peak", score_column],
    )
    holdout[score_column] = pd.to_numeric(holdout[score_column], errors="coerce")
    truth_scores = (
        holdout.loc[_as_bool(holdout["is_injected_signal_peak"])]
        .dropna(subset=["injection_id", score_column])
        .groupby("injection_id", sort=False)[score_column]
        .max()
    )
    n_injections = int(holdout["injection_id"].nunique())
    n_bls_recovered = int(len(truth_scores))
    if not len(tic_scores) or not n_bls_recovered:
        raise ValueError("workload curve requires nonempty real scores and BLS-recovered injections")

    rows: list[dict[str, float | int]] = []
    for requested_fraction in fractions:
        max_pass = max(1, int(np.floor(float(requested_fraction) * len(tic_scores))))
        threshold = float(tic_scores.iloc[max_pass - 1])
        n_pass = int(tic_scores.ge(threshold).sum())
        n_recovered = int(truth_scores.ge(threshold).sum())
        rows.append(
            {
                "requested_workload_fraction": float(requested_fraction),
                "actual_workload_fraction": n_pass / len(tic_scores),
                "threshold": threshold,
                "n_real_tics": int(len(tic_scores)),
                "n_real_tics_passed": n_pass,
                "n_injections": n_injections,
                "n_bls_recovered": n_bls_recovered,
                "n_model_recovered": n_recovered,
                "conditional_retention": n_recovered / n_bls_recovered,
                "end_to_end_recovery": n_recovered / n_injections,
            }
        )

    curve = pd.DataFrame(rows)
    at_five = curve.iloc[(curve["requested_workload_fraction"] - 0.05).abs().argmin()]
    summary = {
        "threshold_at_5pct": float(at_five["threshold"]),
        "actual_workload_fraction_at_5pct": float(at_five["actual_workload_fraction"]),
        "n_real_tics": int(at_five["n_real_tics"]),
        "n_real_tics_passed_at_5pct": int(at_five["n_real_tics_passed"]),
        "n_injections": int(at_five["n_injections"]),
        "n_bls_recovered": int(at_five["n_bls_recovered"]),
        "n_model_recovered_at_5pct": int(at_five["n_model_recovered"]),
        "conditional_retention_at_5pct": float(at_five["conditional_retention"]),
        "end_to_end_recovery_at_5pct": float(at_five["end_to_end_recovery"]),
    }
    return curve, summary


def main() -> int:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fractions = np.unique(np.r_[np.linspace(0.005, 0.30, 60), 0.05])
    models = {
        "Teacher v2": {
            "real": INPUT_ROOT / "shape_plus_raw_chronology.parquet",
            "holdout": INPUT_ROOT / "s56_holdout_teacher_v2.parquet",
            "score": "p_compact_transit",
        },
        "Metadata only": {
            "real": INPUT_ROOT / "metadata_only.parquet",
            "holdout": INPUT_ROOT / "s56_holdout_metadata_only.parquet",
            "score": "p_compact_transit",
        },
        "Teacher v1": {
            "real": INPUT_ROOT / "s56_real_teacher_v1.parquet",
            "holdout": INPUT_ROOT / "s56_holdout_teacher_v1.parquet",
            "score": "p_preserve",
        },
    }

    curves: list[pd.DataFrame] = []
    summaries: dict[str, dict[str, float | int]] = {}
    for name, spec in models.items():
        curve, summary = _workload_curve(
            real_path=spec["real"],
            holdout_path=spec["holdout"],
            score_column=spec["score"],
            fractions=fractions,
        )
        curve.insert(0, "model", name)
        curves.append(curve)
        summaries[name] = summary
    all_curves = pd.concat(curves, ignore_index=True)

    expected = pd.read_csv(ROOT / "workload_matched_model_comparison.csv").set_index("model")
    expected_names = {
        "Teacher v2": "teacher_v2",
        "Metadata only": "metadata_only",
        "Teacher v1": "teacher_v1",
    }
    for name, expected_name in expected_names.items():
        measured = summaries[name]
        reference = expected.loc[expected_name]
        if int(measured["n_model_recovered_at_5pct"]) != int(reference["n_model_recovered"]):
            raise RuntimeError(f"{name} workload reconstruction does not match the frozen summary")

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    all_curves.to_csv(OUT_ROOT / "recovery_vs_review_workload.csv", index=False)
    (OUT_ROOT / "summary.json").write_text(json.dumps(summaries, indent=2, sort_keys=True) + "\n")

    apply_twirl_style("full_page")
    colors = dict(zip(models, get_ordered_palette(len(models))))
    fig, (ax_curve, ax_bar) = plt.subplots(1, 2, figsize=(9.4, 3.8))

    for name in models:
        part = all_curves.loc[all_curves["model"].eq(name)]
        ax_curve.plot(
            100.0 * part["actual_workload_fraction"],
            100.0 * part["conditional_retention"],
            color=colors[name],
            linewidth=2.0,
            label=name,
        )
        point = summaries[name]
        ax_curve.scatter(
            100.0 * float(point["actual_workload_fraction_at_5pct"]),
            100.0 * float(point["conditional_retention_at_5pct"]),
            s=34,
            color=colors[name],
            edgecolor="black",
            linewidth=0.6,
            zorder=4,
        )
    ax_curve.axvline(5.0, color="0.3", linestyle="--", linewidth=1.0)
    ax_curve.text(5.35, 97.0, "Frozen 5% workload", color="0.25", fontsize=8.2, va="top")
    ax_curve.set_xlim(0.0, 30.0)
    ax_curve.set_ylim(0.0, 100.0)
    ax_curve.set_xlabel("Real TICs passed to review [%]")
    ax_curve.set_ylabel("Retention of BLS-recovered injections [%]")
    ax_curve.legend(frameon=True, edgecolor="black", loc="lower right")
    ax_curve.text(0.02, 0.96, "a", transform=ax_curve.transAxes, ha="left", va="top", fontweight="bold")

    n_injections = int(summaries["Teacher v2"]["n_injections"])
    n_bls = int(summaries["Teacher v2"]["n_bls_recovered"])
    bar_labels = ["BLS top 5", "Teacher v2", "Metadata only", "Teacher v1"]
    bar_counts = [
        n_bls,
        int(summaries["Teacher v2"]["n_model_recovered_at_5pct"]),
        int(summaries["Metadata only"]["n_model_recovered_at_5pct"]),
        int(summaries["Teacher v1"]["n_model_recovered_at_5pct"]),
    ]
    bar_colors = ["0.55", colors["Teacher v2"], colors["Metadata only"], colors["Teacher v1"]]
    positions = np.arange(len(bar_labels))
    bars = ax_bar.barh(positions, 100.0 * np.asarray(bar_counts) / n_injections, color=bar_colors)
    ax_bar.set_yticks(positions, labels=bar_labels)
    ax_bar.invert_yaxis()
    ax_bar.set_xlim(0.0, 22.0)
    ax_bar.set_xlabel("Recovery of all locked injections [%]")
    ax_bar.text(0.02, 0.96, "b", transform=ax_bar.transAxes, ha="left", va="top", fontweight="bold")
    for index, (bar, count) in enumerate(zip(bars, bar_counts)):
        if index == 0:
            detail = f"{count}/{n_injections}"
        else:
            detail = f"{count}/{n_injections}  ({100.0 * count / n_bls:.1f}% of BLS)"
        ax_bar.text(
            bar.get_width() + 0.25,
            bar.get_y() + bar.get_height() / 2.0,
            detail,
            va="center",
            ha="left",
            fontsize=8.1,
        )

    fig.subplots_adjust(left=0.09, right=0.985, bottom=0.17, top=0.97, wspace=0.36)
    fig.savefig(OUT_ROOT / "recovery_vs_review_workload.png", dpi=260, bbox_inches="tight")
    fig.savefig(OUT_ROOT / "recovery_vs_review_workload.pdf", bbox_inches="tight")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
