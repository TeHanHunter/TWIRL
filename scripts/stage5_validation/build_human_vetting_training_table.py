#!/usr/bin/env python3
"""Build a training/audit table from browser human-vetting labels.

The browser app deliberately writes labels separately from the review queue.
This script joins the latest human label back to the row metadata without
modifying the label CSV, preserves hidden injection truth, and marks which rows
are usable as strong teacher labels versus audit-only examples.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd

SRC_ROOT = Path(__file__).resolve().parents[2] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.vetting.label_schema import (  # noqa: E402
    AUDIT_LABELS,
    BLS_TRUTH_MATCH_MODES,
    EXCLUDE_LABELS,
    STRONG_LABELS,
    teacher_target,
)
from twirl.vetting.label_io import (  # noqa: E402
    BASE_LABEL_COLUMNS,
    candidate_key as shared_candidate_key,
    latest_label_records,
    normalize_review_queue,
    validate_label_records,
)


DEFAULT_QUEUE = (
    Path(__file__).resolve().parents[2]
    / "reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_pdo/review_queue.csv"
)
DEFAULT_LABELS = DEFAULT_QUEUE.with_name("human_labels_vetted.csv")
DEFAULT_OUT_DIR = DEFAULT_QUEUE.with_name("human_training_table")

LABEL_COLUMNS = ("row_id", "candidate_key", "label", "label_source", "labeler", "notes", "updated_utc")


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".json", ".jsonl"}:
        return pd.read_json(path, lines=suffix == ".jsonl")
    raise ValueError(f"unsupported table format: {path}")


def _candidate_key(row: pd.Series) -> str:
    return shared_candidate_key(row)


def _with_row_ids(queue: pd.DataFrame) -> pd.DataFrame:
    return normalize_review_queue(queue)


def latest_labels(labels_csv: Path) -> pd.DataFrame:
    """Return the latest label per row_id, preserving the app schema."""

    return latest_label_records(Path(labels_csv), columns=BASE_LABEL_COLUMNS)


def join_queue_and_labels(queue_csv: Path, labels_csv: Path) -> pd.DataFrame:
    queue = _with_row_ids(_read_table(queue_csv))
    for col in ("label", "label_source", "labeler", "notes", "updated_utc"):
        if col in queue:
            queue = queue.rename(columns={col: f"queue_{col}"})
    raw_labels = latest_labels(labels_csv)
    validate_label_records(queue, raw_labels)
    labels = raw_labels.rename(
        columns={
            "candidate_key": "human_candidate_key",
            "label": "human_label",
            "label_source": "human_label_source",
            "labeler": "human_labeler",
            "notes": "human_notes",
            "updated_utc": "human_updated_utc",
        }
    )
    merged = queue.merge(labels, on="row_id", how="left", validate="one_to_one")
    for col in ("human_label", "human_label_source", "human_labeler", "human_notes", "human_updated_utc"):
        if col not in merged:
            merged[col] = ""
        merged[col] = merged[col].fillna("").astype(str)
    if "human_candidate_key" not in merged:
        merged["human_candidate_key"] = ""
    merged["human_candidate_key"] = merged["human_candidate_key"].fillna("").astype(str)
    merged["candidate_key_matches_label"] = (
        merged["human_candidate_key"].eq("")
        | merged["candidate_key"].fillna("").astype(str).eq(merged["human_candidate_key"])
    )
    return merged


def _truth_mode(df: pd.DataFrame) -> pd.Series:
    if "topn_recovery_status" in df:
        return df["topn_recovery_status"].fillna("").astype(str)
    if "recovery_status" in df:
        return df["recovery_status"].fillna("").astype(str)
    return pd.Series("", index=df.index)


def _teacher_target(label: str) -> str:
    return teacher_target(label)


def add_training_flags(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    label = out["human_label"].fillna("").astype(str)
    out["is_labeled"] = label.ne("")
    source = out["truth_source_kind"] if "truth_source_kind" in out else out.get("source_kind", pd.Series("", index=out.index))
    out["is_injected_row"] = source.fillna("").astype(str).eq("injection_recovery")
    mode = _truth_mode(out)
    out["bls_truth_match"] = mode.isin(BLS_TRUTH_MATCH_MODES)
    out["bls_truth_match_definition"] = "topn_recovery_status in " + ",".join(sorted(BLS_TRUTH_MATCH_MODES))
    out["leo_pc_or_fp"] = out.get("leo_class", "").fillna("").astype(str).isin(["PC", "FP"])
    out["human_signal_like"] = label.eq("planet_like")
    out["human_audit_label"] = label.isin(AUDIT_LABELS)
    out["human_exclude_label"] = label.isin(EXCLUDE_LABELS)
    out["human_strong_label"] = label.isin(STRONG_LABELS)
    out["teacher_target"] = label.map(_teacher_target).fillna("")
    out["teacher_include"] = out["human_strong_label"]
    out["audit_include"] = out["is_labeled"] & ~out["human_exclude_label"]
    return out


def attach_recovery_gate(df: pd.DataFrame, gate_csv: Path | None) -> pd.DataFrame:
    if gate_csv is None or not gate_csv.exists():
        return df
    gate = _read_table(gate_csv)
    if "injection_id" not in gate:
        raise KeyError(f"gate table missing injection_id: {gate_csv}")
    keep = [
        col
        for col in (
            "injection_id",
            "top1_match",
            "top20_match",
            "ranking_loss",
            "not_in_top20",
            "best_signal_peak_rank",
            "recovery_cell_n",
            "recovery_cell_top20_fraction",
            "recovery_gate",
        )
        if col in gate.columns
    ]
    gate = gate.loc[:, keep].drop_duplicates("injection_id")
    return df.merge(gate, on="injection_id", how="left", suffixes=("", "_gate"))


def assign_splits(
    df: pd.DataFrame,
    *,
    validation_fraction: float,
    test_fraction: float,
    random_state: int,
) -> pd.Series:
    split = pd.Series("unlabeled_or_audit", index=df.index, dtype=object)
    eligible = df["teacher_include"].fillna(False).astype(bool)
    if not eligible.any():
        return split
    rng = np.random.default_rng(random_state)
    keys = (
        df.loc[eligible, "teacher_target"].fillna("").astype(str)
        + "|"
        + df.loc[eligible, "source_kind"].fillna("").astype(str)
    )
    for _, idx_values in keys.groupby(keys).groups.items():
        idx = np.asarray(list(idx_values), dtype=int)
        rng.shuffle(idx)
        n = len(idx)
        n_test = int(round(test_fraction * n))
        n_val = int(round(validation_fraction * n))
        if n >= 3 and test_fraction > 0:
            n_test = max(1, n_test)
        if n - n_test >= 3 and validation_fraction > 0:
            n_val = max(1, n_val)
        n_test = min(n_test, max(0, n - 1))
        n_val = min(n_val, max(0, n - n_test - 1))
        split.loc[idx] = "train"
        if n_test:
            split.loc[idx[:n_test]] = "test"
        if n_val:
            split.loc[idx[n_test:n_test + n_val]] = "validation"
    return split


def build_training_table(
    *,
    queue_csv: Path,
    labels_csv: Path,
    out_dir: Path,
    recovery_gate_csv: Path | None = None,
    validation_fraction: float = 0.20,
    test_fraction: float = 0.20,
    random_state: int = 56,
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    table = join_queue_and_labels(queue_csv, labels_csv)
    table = attach_recovery_gate(table, recovery_gate_csv)
    table = add_training_flags(table)
    table["training_split"] = assign_splits(
        table,
        validation_fraction=validation_fraction,
        test_fraction=test_fraction,
        random_state=random_state,
    )

    labeled = table[table["is_labeled"]].copy()
    teacher = table[table["teacher_include"]].copy()
    audit = table[table["audit_include"]].copy()
    table.to_csv(out_dir / "human_vetting_training_table.csv", index=False)
    teacher.to_csv(out_dir / "teacher_labeled_rows.csv", index=False)
    audit.to_csv(out_dir / "audit_labeled_rows.csv", index=False)

    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "queue_csv": str(queue_csv),
        "labels_csv": str(labels_csv),
        "recovery_gate_csv": str(recovery_gate_csv) if recovery_gate_csv is not None else "",
        "n_rows": int(len(table)),
        "n_labeled": int(table["is_labeled"].sum()),
        "n_teacher_rows": int(table["teacher_include"].sum()),
        "n_audit_rows": int(table["audit_include"].sum()),
        "label_counts": {
            str(k): int(v) for k, v in labeled["human_label"].value_counts().sort_index().items()
        },
        "teacher_target_counts": {
            str(k): int(v) for k, v in teacher["teacher_target"].value_counts().sort_index().items()
        },
        "training_split_counts": {
            str(k): int(v) for k, v in table["training_split"].value_counts().sort_index().items()
        },
        "human_label_by_bls_truth_match": pd.crosstab(
            labeled["bls_truth_match"],
            labeled["human_label"],
        ).to_dict() if len(labeled) else {},
        "human_label_by_leo_class": pd.crosstab(
            labeled.get("leo_class", pd.Series("", index=labeled.index)).fillna("").astype(str),
            labeled["human_label"],
        ).to_dict() if len(labeled) else {},
        "outputs": {
            "training_table": str(out_dir / "human_vetting_training_table.csv"),
            "teacher_rows": str(out_dir / "teacher_labeled_rows.csv"),
            "audit_rows": str(out_dir / "audit_labeled_rows.csv"),
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue-csv", type=Path, default=DEFAULT_QUEUE)
    parser.add_argument("--labels-csv", type=Path, default=DEFAULT_LABELS)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--recovery-gate-csv", type=Path, default=None)
    parser.add_argument("--validation-fraction", type=float, default=0.20)
    parser.add_argument("--test-fraction", type=float, default=0.20)
    parser.add_argument("--random-state", type=int, default=56)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    build_training_table(
        queue_csv=args.queue_csv,
        labels_csv=args.labels_csv,
        out_dir=args.out_dir,
        recovery_gate_csv=args.recovery_gate_csv,
        validation_fraction=args.validation_fraction,
        test_fraction=args.test_fraction,
        random_state=args.random_state,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
