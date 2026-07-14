#!/usr/bin/env python3
"""Summarize audit-only Teacher-v2 scores for Franklin's legacy labels."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--frozen-selection", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    scores = (
        pd.read_parquet(args.scores)
        if args.scores.suffix.lower() == ".parquet"
        else pd.read_csv(args.scores, low_memory=False)
    )
    frozen = json.loads(args.frozen_selection.read_text())
    threshold = float(frozen["frozen_compact_threshold"])
    required = {
        "review_id",
        "tic",
        "legacy_human_label",
        "a2v1_transfer_ok",
        "p_compact_transit",
        "p_preserve",
        "predicted_morphology",
        "teacher_v2_training_include",
    }
    missing = sorted(required - set(scores.columns))
    if missing:
        raise KeyError(f"Franklin score audit is missing columns: {missing}")
    include = scores["teacher_v2_training_include"]
    if include.dtype != bool:
        include = include.fillna("").astype(str).str.lower().isin(
            {"1", "true", "t", "yes", "y"}
        )
    if include.any():
        raise ValueError("Franklin audit rows were incorrectly activated for training")
    transfer = scores["a2v1_transfer_ok"]
    if transfer.dtype != bool:
        transfer = transfer.fillna("").astype(str).str.lower().isin(
            {"1", "true", "t", "yes", "y"}
        )
    work = scores.copy()
    work["a2v1_transfer_ok"] = transfer
    work["passes_frozen_compact_threshold"] = pd.to_numeric(
        work["p_compact_transit"], errors="raise"
    ).ge(threshold)
    rows: list[dict[str, object]] = []
    for (label, compatible), group in work.groupby(
        ["legacy_human_label", "a2v1_transfer_ok"], sort=True, dropna=False
    ):
        compact = pd.to_numeric(group["p_compact_transit"], errors="raise")
        preserve = pd.to_numeric(group["p_preserve"], errors="raise")
        rows.append(
            {
                "legacy_human_label": str(label),
                "a2v1_ephemeris_compatible": bool(compatible),
                "n": int(len(group)),
                "compact_pass_fraction": float(
                    group["passes_frozen_compact_threshold"].mean()
                ),
                "p_compact_q50": float(compact.quantile(0.50)),
                "p_compact_q90": float(compact.quantile(0.90)),
                "p_compact_q95": float(compact.quantile(0.95)),
                "p_compact_q99": float(compact.quantile(0.99)),
                "p_preserve_q50": float(preserve.quantile(0.50)),
                "p_preserve_q95": float(preserve.quantile(0.95)),
            }
        )
    distribution = pd.DataFrame(rows)
    background = work["legacy_human_label"].isin(
        {"instrumental_or_systematic", "uncertain", "no_visible_signal"}
    )
    compatible_background = background & transfer
    args.out_dir.mkdir(parents=True, exist_ok=True)
    distribution.to_csv(args.out_dir / "score_distribution_by_legacy_label.csv", index=False)
    work.sort_values(
        ["p_compact_transit", "p_preserve"], ascending=False, kind="stable"
    ).head(250).to_csv(args.out_dir / "highest_scoring_legacy_rows.csv", index=False)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "n_scored": int(len(work)),
        "n_unique_tics": int(work["tic"].nunique()),
        "n_legacy_background": int(background.sum()),
        "n_ephemeris_compatible_background": int(compatible_background.sum()),
        "frozen_compact_threshold": threshold,
        "legacy_background_pass_fraction_all_current_top1": (
            float(work.loc[background, "passes_frozen_compact_threshold"].mean())
            if background.any()
            else np.nan
        ),
        "legacy_background_pass_fraction_compatible_only": (
            float(
                work.loc[
                    compatible_background, "passes_frozen_compact_threshold"
                ].mean()
            )
            if compatible_background.any()
            else np.nan
        ),
        "training_include": False,
        "interpretation": (
            "Legacy labels are descriptive audit strata. Only the compatible subset "
            "refers to the same A2v1 event; neither subset is a Teacher-v2 target."
        ),
    }
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
