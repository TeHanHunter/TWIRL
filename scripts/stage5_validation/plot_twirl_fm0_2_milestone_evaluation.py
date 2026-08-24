#!/usr/bin/env python3
"""Plot one development-only TWIRL-FM0.2 milestone evaluation."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re
import sys

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.plotting.style import apply_twirl_style, get_ordered_palette  # noqa: E402


STEP_PATTERN = re.compile(
    r"\[fm0-train\] step=(?P<step>\d+)/\d+ loss=(?P<loss>[0-9.eE+-]+) "
    r"lr=(?P<lr>[0-9.eE+-]+) elapsed_s=(?P<elapsed>[0-9.eE+-]+)"
)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--representation-health", type=Path, required=True)
    parser.add_argument("--training-summary", type=Path, required=True)
    parser.add_argument("--training-log", type=Path, action="append", required=True)
    parser.add_argument("--output-stem", type=Path, required=True)
    parser.add_argument("--metrics-csv", type=Path, required=True)
    parser.add_argument("--training-trace-csv", type=Path, required=True)
    return parser


def _interval(payload: dict[str, object], key: str) -> tuple[float, float, float]:
    interval = payload[key]
    assert isinstance(interval, dict)
    return tuple(float(interval[name]) for name in ("mean", "lower", "upper"))


def _training_trace(paths: list[Path]) -> list[dict[str, float]]:
    by_step: dict[int, dict[str, float]] = {}
    for path in paths:
        for line in path.read_text(encoding="utf-8").splitlines():
            match = STEP_PATTERN.search(line)
            if match is None:
                continue
            step = int(match.group("step"))
            by_step[step] = {
                "step": float(step),
                "total_loss": float(match.group("loss")),
                "learning_rate": float(match.group("lr")),
                "elapsed_seconds_this_invocation": float(match.group("elapsed")),
            }
    if not by_step:
        raise ValueError("training logs contain no FM0 progress records")
    return [by_step[step] for step in sorted(by_step)]


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = _parser().parse_args()
    health = json.loads(args.representation_health.read_text(encoding="utf-8"))
    summary = json.loads(args.training_summary.read_text(encoding="utf-8"))
    if health.get("schema_version") != "twirl_fm0_representation_health_v2":
        raise ValueError("representation-health schema differs")
    if health.get("data_access", {}).get("sealed_test_access_count") != 0:
        raise ValueError("sealed-test access is nonzero")
    step = int(health["run"]["global_step"])
    if int(summary["global_step"]) != step:
        raise ValueError("training summary and evaluation milestone differ")

    trace = _training_trace(args.training_log)
    _write_csv(args.training_trace_csv, trace)
    trained = health["trained_encoder"]["representations"]
    random = health["random_encoder_control"]["representations"]
    trained_minus_random = health["source_paired_trained_minus_random"]
    masked = health["masked_reconstruction"]

    z_rank = float(trained["z_window"]["embedding_health"]["effective_rank"])
    z_constant = int(
        trained["z_window"]["embedding_health"]["zero_or_constant_dimensions"]
    )
    z_safe_mean, z_safe_lower, z_safe_upper = _interval(
        trained["z_window"]["safe_mask_pair_separation"],
        "paired_minus_unrelated_source_clustered_95_interval",
    )
    z_delta_mean, z_delta_lower, z_delta_upper = _interval(
        trained_minus_random["z_window"],
        "trained_minus_random_source_clustered_95_interval",
    )
    masked_huber = float(masked["masked_huber_mean"])
    zero_huber = float(masked["zero_prediction_masked_huber_mean"])
    gate_rows = [
        {
            "quantity": "z_window_effective_rank",
            "step_500_value": z_rank,
            "future_step_2000_rule": ">= 26",
            "would_meet_threshold_now": z_rank >= 26.0,
            "formal_gate_applied": False,
        },
        {
            "quantity": "z_window_zero_or_constant_dimensions",
            "step_500_value": z_constant,
            "future_step_2000_rule": "<= 0",
            "would_meet_threshold_now": z_constant <= 0,
            "formal_gate_applied": False,
        },
        {
            "quantity": "z_window_paired_minus_unrelated_ci_lower",
            "step_500_value": z_safe_lower,
            "future_step_2000_rule": "> 0",
            "would_meet_threshold_now": z_safe_lower > 0.0,
            "formal_gate_applied": False,
        },
        {
            "quantity": "z_window_trained_minus_random_ci_lower",
            "step_500_value": z_delta_lower,
            "future_step_2000_rule": "> 0",
            "would_meet_threshold_now": z_delta_lower > 0.0,
            "formal_gate_applied": False,
        },
        {
            "quantity": "development_masked_huber",
            "step_500_value": masked_huber,
            "future_step_2000_rule": "<= 0.004864579067",
            "would_meet_threshold_now": masked_huber <= 0.004864579067,
            "formal_gate_applied": False,
        },
        {
            "quantity": "sealed_test_access_count",
            "step_500_value": 0,
            "future_step_2000_rule": "== 0",
            "would_meet_threshold_now": True,
            "formal_gate_applied": False,
        },
    ]
    _write_csv(args.metrics_csv, gate_rows)

    template = apply_twirl_style("full_page")
    colors = get_ordered_palette(5)
    fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.6))
    ax_loss, ax_rank, ax_sep, ax_recon = axes.flat

    steps = np.asarray([row["step"] for row in trace])
    losses = np.asarray([row["total_loss"] for row in trace])
    ax_loss.plot(steps, losses, color=colors[2], lw=1.25, marker="o", ms=2.2)
    ax_loss.set(xlabel="Optimizer step", ylabel="Training objective")
    ax_loss.text(0.02, 0.94, "a", transform=ax_loss.transAxes, va="top", fontweight="bold")
    ax_loss.annotate(
        f"step {step}: {losses[-1]:.4f}",
        xy=(steps[-1], losses[-1]),
        xytext=(-72, 17),
        textcoords="offset points",
        arrowprops={"arrowstyle": "-", "color": "0.35", "lw": 0.7},
    )

    rank_labels = ["trained\n$h$", "trained\n$z$", "random\n$h$", "random\n$z$", "PCA\n(32-D)"]
    rank_values = [
        float(trained["h_window"]["embedding_health"]["effective_rank"]),
        z_rank,
        float(random["h_window"]["embedding_health"]["effective_rank"]),
        float(random["z_window"]["embedding_health"]["effective_rank"]),
        float(
            health["baselines"]["pca"]["development_evaluation"]
            ["embedding_health"]["effective_rank"]
        ),
    ]
    ax_rank.bar(np.arange(5), rank_values, color=[colors[3], colors[3], colors[0], colors[0], colors[1]])
    ax_rank.axhline(26.0, color="0.25", ls="--", lw=0.9, label="future step-2000 threshold")
    ax_rank.set(xticks=np.arange(5), xticklabels=rank_labels, ylabel="Effective rank", ylim=(0, 29))
    ax_rank.legend(loc="upper left")
    ax_rank.text(
        0.98, 0.94, "b", transform=ax_rank.transAxes,
        ha="right", va="top", fontweight="bold",
    )

    sep_labels = ["paired−unrelated\n$h$", "paired−unrelated\n$z$", "trained−random\n$h$", "trained−random\n$z$"]
    sep_payloads = [
        _interval(trained[name]["safe_mask_pair_separation"], "paired_minus_unrelated_source_clustered_95_interval")
        for name in ("h_window", "z_window")
    ] + [
        _interval(trained_minus_random[name], "trained_minus_random_source_clustered_95_interval")
        for name in ("h_window", "z_window")
    ]
    means = np.asarray([value[0] for value in sep_payloads])
    lower = np.asarray([value[1] for value in sep_payloads])
    upper = np.asarray([value[2] for value in sep_payloads])
    ax_sep.errorbar(
        np.arange(4), means, yerr=np.vstack((means - lower, upper - means)),
        fmt="o", ms=4, capsize=2.5, color=colors[4], ecolor=colors[4], lw=1.0,
    )
    ax_sep.axhline(0.0, color="0.25", ls="--", lw=0.9)
    ax_sep.set(xticks=np.arange(4), xticklabels=sep_labels, ylabel="Cosine-similarity difference")
    ax_sep.tick_params(axis="x", rotation=16)
    ax_sep.text(0.02, 0.94, "c", transform=ax_sep.transAxes, va="top", fontweight="bold")

    recon_values = np.asarray([masked_huber, zero_huber]) * 1e3
    ax_recon.bar([0, 1], recon_values, color=[colors[3], "0.72"])
    ax_recon.set(xticks=[0, 1], xticklabels=["trained model", "zero prediction"], ylabel=r"Development masked Huber ($\times 10^{-3}$)")
    ax_recon.set_ylim(0, max(recon_values) * 1.18)
    relative = 100.0 * (masked_huber / zero_huber - 1.0)
    ax_recon.text(0.5, max(recon_values) * 1.06, f"model is {relative:.2f}% worse", ha="center")
    ax_recon.text(0.02, 0.94, "d", transform=ax_recon.transAxes, va="top", fontweight="bold")

    fig.subplots_adjust(
        left=0.09, right=0.985, bottom=0.11, top=0.98,
        hspace=0.36, wspace=0.28,
    )
    args.output_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_stem.with_suffix(".png"), dpi=240)
    fig.savefig(args.output_stem.with_suffix(".pdf"))
    plt.close(fig)
    print(args.output_stem.with_suffix(".png"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
