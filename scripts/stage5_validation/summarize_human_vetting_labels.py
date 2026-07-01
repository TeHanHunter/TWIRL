#!/usr/bin/env python3
"""Summarize human labels for a TWIRL review queue.

The vetting app writes labels separately from the review queue. This script
joins the latest human label per ``row_id`` back onto the queue and summarizes
label outcomes by hidden recovery/truth strata such as ``review_stratum``,
``recovery_mode``, empirical SNR, period, depth, and Tmag bins.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


SNR_BINS = (-np.inf, 1.0, 3.0, 5.0, 7.0, 10.0, 20.0, np.inf)
SNR_LABELS = ("<1", "1-3", "3-5", "5-7", "7-10", "10-20", ">20")
PERIOD_BINS = (0.0, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, np.inf)
PERIOD_LABELS = ("<0.25", "0.25-0.5", "0.5-1", "1-2", "2-5", "5-10", ">10")
TMAG_BINS = (-np.inf, 16.0, 17.0, 18.0, 19.0, 20.0, np.inf)
TMAG_LABELS = ("<16", "16-17", "17-18", "18-19", "19-20", ">20")
DEPTH_BINS = (0.0, 0.01, 0.03, 0.1, 0.3, 1.0, np.inf)
DEPTH_LABELS = ("<1%", "1-3%", "3-10%", "10-30%", "30-100%", ">100%")


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def _latest_labels(labels_path: Path) -> pd.DataFrame:
    columns = ["row_id", "candidate_key", "label", "label_source", "labeler", "notes", "updated_utc"]
    if not labels_path.exists() or labels_path.stat().st_size == 0:
        return pd.DataFrame(columns=columns)
    labels = pd.read_csv(labels_path)
    for col in columns:
        if col not in labels:
            labels[col] = ""
    if labels.empty:
        return labels.loc[:, columns]
    labels["row_id"] = pd.to_numeric(labels["row_id"], errors="coerce")
    labels = labels.dropna(subset=["row_id"]).copy()
    labels["row_id"] = labels["row_id"].astype(int)
    if "updated_utc" in labels:
        labels["_updated_sort"] = pd.to_datetime(labels["updated_utc"], errors="coerce")
        labels = labels.sort_values(["row_id", "_updated_sort"], na_position="first")
    return labels.drop_duplicates("row_id", keep="last").loc[:, columns]


def _with_row_ids(queue: pd.DataFrame) -> pd.DataFrame:
    out = queue.copy()
    if "row_id" not in out:
        out.insert(0, "row_id", np.arange(len(out), dtype=int))
    out["row_id"] = pd.to_numeric(out["row_id"], errors="coerce").astype(int)
    return out


def _merge_labels(queue_csv: Path, labels_csv: Path) -> pd.DataFrame:
    queue = _with_row_ids(pd.read_csv(queue_csv))
    for col in ("label", "label_source", "labeler", "notes", "updated_utc"):
        if col in queue:
            queue = queue.rename(columns={col: f"queue_{col}"})
    labels = _latest_labels(labels_csv)
    merged = queue.merge(labels, on="row_id", how="left", validate="one_to_one")
    for col in ("label", "label_source", "labeler", "notes", "updated_utc"):
        if col not in merged:
            merged[col] = ""
        merged[col] = merged[col].fillna("").astype(str)
    merged["is_labeled"] = merged["label"].ne("")
    return merged


def _write_crosstab(df: pd.DataFrame, index_col: str, out_path: Path) -> dict[str, dict[str, int]]:
    if index_col not in df:
        return {}
    labeled = df[df["is_labeled"]].copy()
    if labeled.empty:
        table = pd.DataFrame()
    else:
        index_values = labeled[index_col].astype(object).where(labeled[index_col].notna(), "").astype(str)
        label_values = labeled["label"].astype(object).where(labeled["label"].notna(), "").astype(str)
        table = pd.crosstab(index_values, label_values)
    table.to_csv(out_path)
    return {
        str(idx): {str(col): int(value) for col, value in row.items()}
        for idx, row in table.iterrows()
    }


def _add_bins(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "delta_signal_multi_snr_mad" in out:
        out["snr_bin"] = pd.cut(pd.to_numeric(out["delta_signal_multi_snr_mad"], errors="coerce"), SNR_BINS, labels=SNR_LABELS)
    if "truth_period_d" in out:
        out["period_bin"] = pd.cut(pd.to_numeric(out["truth_period_d"], errors="coerce"), PERIOD_BINS, labels=PERIOD_LABELS)
    if "tmag" in out:
        out["tmag_bin"] = pd.cut(pd.to_numeric(out["tmag"], errors="coerce"), TMAG_BINS, labels=TMAG_LABELS)
    depth_col = "truth_model_depth" if "truth_model_depth" in out else "truth_depth"
    if depth_col in out:
        out["depth_bin"] = pd.cut(pd.to_numeric(out[depth_col], errors="coerce"), DEPTH_BINS, labels=DEPTH_LABELS)
    return out


def summarize_labels(queue_csv: Path, labels_csv: Path, out_dir: Path) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    merged = _add_bins(_merge_labels(queue_csv, labels_csv))
    merged.to_csv(out_dir / "labeled_queue_joined.csv", index=False)

    labeled = merged[merged["is_labeled"]].copy()
    summary: dict[str, Any] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "queue_csv": str(queue_csv),
        "labels_csv": str(labels_csv),
        "n_queue_rows": int(len(merged)),
        "n_labeled": int(len(labeled)),
        "n_unlabeled": int(len(merged) - len(labeled)),
        "label_counts": labeled["label"].value_counts().sort_index().to_dict() if len(labeled) else {},
        "tables": {},
    }
    for index_col in (
        "review_stratum",
        "recovery_mode",
        "snr_bin",
        "period_bin",
        "tmag_bin",
        "depth_bin",
        "leo_class",
        "recovery_status",
        "topn_recovery_status",
    ):
        summary["tables"][index_col] = _write_crosstab(merged, index_col, out_dir / f"label_by_{index_col}.csv")

    if len(labeled):
        numeric_cols = [
            col
            for col in (
                "delta_signal_multi_snr_mad",
                "depth_retention_frac",
                "truth_period_d",
                "truth_model_depth",
                "truth_radius_rearth",
                "tmag",
            )
            if col in labeled
        ]
        if numeric_cols:
            stats = labeled.groupby("label", observed=False)[numeric_cols].median(numeric_only=True)
            stats.to_csv(out_dir / "label_median_truth_stats.csv")
            summary["label_median_truth_stats"] = {
                str(idx): {str(col): float(value) for col, value in row.dropna().items()}
                for idx, row in stats.iterrows()
            }

    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue-csv", type=Path, required=True)
    parser.add_argument("--labels-csv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    summarize_labels(args.queue_csv, args.labels_csv, args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
