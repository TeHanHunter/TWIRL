#!/usr/bin/env python3
"""Compare old and retrained CNNs on the same adjudicated TIC-grouped holdout."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def _read(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)


def _metrics(truth: pd.Series, prediction: pd.Series) -> dict:
    truth = truth.fillna("").astype(str)
    prediction = prediction.fillna("").astype(str)
    labels = sorted(set(truth))
    recalls = []
    for label in labels:
        mask = truth.eq(label)
        recalls.append(float(prediction.loc[mask].eq(label).mean()))
    return {
        "n": int(len(truth)),
        "accuracy": float(prediction.eq(truth).mean()) if len(truth) else float("nan"),
        "balanced_accuracy": float(np.mean(recalls)) if recalls else float("nan"),
        "labels": labels,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--new-root", type=Path, required=True)
    parser.add_argument("--old-score-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for profile in ("cnn_shape_only", "cnn_shape_plus_bls"):
        new = _read(args.new_root / profile / "cnn_predictions.parquet")
        old = _read(args.old_score_root / profile / "cnn_scored_candidates.parquet")
        key = "source_uid" if "source_uid" in new and "source_uid" in old else "review_id"
        holdout = new.loc[new["cnn_training_split"].fillna("").astype(str).eq("test")].copy()
        if holdout[key].duplicated().any() or old[key].duplicated().any():
            raise ValueError(f"duplicate {key} values prevent checkpoint comparison")
        merged = holdout.merge(
            old.loc[:, [key, "cnn_label"]].rename(columns={"cnn_label": "old_cnn_label"}),
            on=key,
            how="left",
            validate="one_to_one",
        )
        old_classes = set(old.filter(regex=r"^cnn_p_").columns.str.replace("cnn_p_", "", regex=False))
        eligible_old = merged["main_teacher_target"].fillna("").astype(str).isin(old_classes)
        new_metrics = _metrics(merged["main_teacher_target"], merged["cnn_label"])
        old_metrics = _metrics(
            merged.loc[eligible_old, "main_teacher_target"],
            merged.loc[eligible_old, "old_cnn_label"],
        )
        records.append(
            {
                "profile": profile,
                "holdout_rows": int(len(merged)),
                "holdout_tics": int(merged["tic"].nunique()),
                "old_class_coverage_rows": int(eligible_old.sum()),
                "new_accuracy": new_metrics["accuracy"],
                "new_balanced_accuracy": new_metrics["balanced_accuracy"],
                "old_accuracy_on_supported_classes": old_metrics["accuracy"],
                "old_balanced_accuracy_on_supported_classes": old_metrics["balanced_accuracy"],
            }
        )
    comparison = pd.DataFrame(records)
    comparison.to_csv(args.out_dir / "checkpoint_comparison.csv", index=False)
    summary = {"profiles": records, "holdout_policy": "retrained TIC-grouped test split"}
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
