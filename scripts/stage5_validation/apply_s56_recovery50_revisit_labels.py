#!/usr/bin/env python3
"""Overlay revisit labels and harmonic notes onto the S56 recovery50 training table."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.vetting.recovery50_teacher import (  # noqa: E402
    LabelPolicy,
    add_deterministic_splits,
    add_harmonic_ephemeris_annotations,
    add_label_roles,
    json_default,
    latest_labels,
    read_table,
)
from twirl.vetting.label_io import normalize_review_queue, validate_label_records  # noqa: E402

DEFAULT_ROOT = REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue_2k"
DEFAULT_BASE_TABLE = DEFAULT_ROOT / "human_training_table/human_vetting_training_table.csv"
DEFAULT_REVISIT_DIR = DEFAULT_ROOT / "real_planet_revisit_wide"
DEFAULT_REVISIT_QUEUE = DEFAULT_REVISIT_DIR / "review_queue_real_planet_revisit.csv"
DEFAULT_REVISIT_LABELS = DEFAULT_REVISIT_DIR / "human_labels_revisit.csv"
DEFAULT_OUT_DIR = DEFAULT_ROOT / "human_training_table_revisit_harmonic"


def _as_bool_counts(series: pd.Series) -> dict[str, int]:
    return {str(k): int(v) for k, v in series.value_counts(dropna=False).sort_index().items()}


def load_revisit_labels(revisit_queue: Path, revisit_labels: Path) -> pd.DataFrame:
    queue = normalize_review_queue(read_table(revisit_queue))
    raw_labels = latest_labels(revisit_labels)
    validate_label_records(queue, raw_labels)
    labels = raw_labels.rename(
        columns={
            "label": "revisit_label",
            "label_source": "revisit_label_source",
            "labeler": "revisit_labeler",
            "notes": "revisit_notes",
            "updated_utc": "revisit_updated_utc",
            "candidate_key": "revisit_candidate_key",
        }
    )
    labels = labels.drop(columns=[c for c in ("tic", "sector") if c in labels], errors="ignore")
    merged = queue.merge(labels, on="row_id", how="left", validate="one_to_one")
    for col in (
        "revisit_label",
        "revisit_label_source",
        "revisit_labeler",
        "revisit_notes",
        "revisit_updated_utc",
        "revisit_candidate_key",
    ):
        if col not in merged:
            merged[col] = ""
        merged[col] = merged[col].fillna("").astype(str)
    return merged


def apply_revisit_overlay(
    *,
    base_table: Path,
    revisit_queue: Path,
    revisit_labels: Path,
    out_dir: Path,
    min_multiclass_count: int,
    validation_fraction: float,
    test_fraction: float,
    random_state: int,
) -> dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    table = read_table(base_table).copy()
    revisit = load_revisit_labels(revisit_queue, revisit_labels)
    if "review_id" not in table:
        raise ValueError("base table is missing review_id")
    if "revisit_source_review_id" not in revisit:
        raise ValueError("revisit queue is missing revisit_source_review_id")

    for col in ("human_label", "human_label_source", "human_labeler", "human_notes", "human_updated_utc"):
        if col in table and f"pre_revisit_{col}" not in table:
            table[f"pre_revisit_{col}"] = table[col]
    for col in (
        "revisit_row_id",
        "revisit_label",
        "revisit_label_source",
        "revisit_labeler",
        "revisit_notes",
        "revisit_updated_utc",
    ):
        if col not in table:
            table[col] = "" if col != "revisit_row_id" else np.nan

    applied = 0
    matched = 0
    labeled_revisit = revisit[revisit["revisit_label"].astype(str).ne("")].copy()
    for _, row in labeled_revisit.iterrows():
        source_review_id = str(row.get("revisit_source_review_id", "") or "")
        if not source_review_id:
            continue
        mask = table["review_id"].astype(str).eq(source_review_id)
        if not mask.any():
            continue
        matched += int(mask.sum())
        table.loc[mask, "revisit_row_id"] = int(row["row_id"])
        table.loc[mask, "revisit_label"] = row["revisit_label"]
        table.loc[mask, "revisit_label_source"] = row["revisit_label_source"]
        table.loc[mask, "revisit_labeler"] = row["revisit_labeler"]
        table.loc[mask, "revisit_notes"] = row["revisit_notes"]
        table.loc[mask, "revisit_updated_utc"] = row["revisit_updated_utc"]
        table.loc[mask, "human_label"] = row["revisit_label"]
        table.loc[mask, "human_label_source"] = "revisit:" + str(row["revisit_label_source"] or "")
        table.loc[mask, "human_labeler"] = row["revisit_labeler"]
        table.loc[mask, "human_notes"] = row["revisit_notes"]
        table.loc[mask, "human_updated_utc"] = row["revisit_updated_utc"]
        applied += 1

    table = add_harmonic_ephemeris_annotations(table)
    table = add_label_roles(table, LabelPolicy(min_multiclass_count=min_multiclass_count))
    table = add_deterministic_splits(
        table,
        validation_fraction=validation_fraction,
        test_fraction=test_fraction,
        random_state=random_state,
    )

    teacher = table[table["main_teacher_include"].fillna(False).astype(bool)].copy()
    audit = table[table["audit_include"].fillna(False).astype(bool)].copy()
    harmonic = table[table["harmonic_suspect"].fillna(False).astype(bool)].copy()
    table_path = out_dir / "human_vetting_training_table.csv"
    teacher_path = out_dir / "teacher_labeled_rows.csv"
    audit_path = out_dir / "audit_labeled_rows.csv"
    harmonic_path = out_dir / "harmonic_suspect_rows.csv"
    table.to_csv(table_path, index=False)
    teacher.to_csv(teacher_path, index=False)
    audit.to_csv(audit_path, index=False)
    harmonic.to_csv(harmonic_path, index=False)

    summary = {
        "base_table": str(base_table),
        "revisit_queue": str(revisit_queue),
        "revisit_labels": str(revisit_labels),
        "out_dir": str(out_dir),
        "n_rows": int(len(table)),
        "n_revisit_rows": int(len(revisit)),
        "n_revisit_labeled": int(len(labeled_revisit)),
        "n_revisit_applied": int(applied),
        "n_base_rows_matched": int(matched),
        "n_teacher_rows": int(table["main_teacher_include"].sum()),
        "n_harmonic_suspect": int(table["harmonic_suspect"].sum()),
        "label_counts": _as_bool_counts(table.loc[table["is_labeled"], "human_label"]),
        "revisit_label_counts": _as_bool_counts(labeled_revisit["revisit_label"]),
        "teacher_target_counts": _as_bool_counts(teacher["main_teacher_target"]),
        "ephemeris_status_counts": _as_bool_counts(table["ephemeris_status"]),
        "refold_factor_counts": _as_bool_counts(harmonic["refold_factor"].astype(str) if len(harmonic) else pd.Series(dtype=str)),
        "outputs": {
            "training_table": str(table_path),
            "teacher_rows": str(teacher_path),
            "audit_rows": str(audit_path),
            "harmonic_rows": str(harmonic_path),
            "summary": str(out_dir / "summary.json"),
        },
        "policy": (
            "revisit labels override prior human labels for matched rows; reviewer harmonic notes set model_* "
            "fold ephemerides but remain leakage-excluded audit columns"
        ),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=json_default) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True, default=json_default))
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-table", type=Path, default=DEFAULT_BASE_TABLE)
    parser.add_argument("--revisit-queue", type=Path, default=DEFAULT_REVISIT_QUEUE)
    parser.add_argument("--revisit-labels", type=Path, default=DEFAULT_REVISIT_LABELS)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--min-multiclass-count", type=int, default=40)
    parser.add_argument("--validation-fraction", type=float, default=0.20)
    parser.add_argument("--test-fraction", type=float, default=0.20)
    parser.add_argument("--random-state", type=int, default=56)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    apply_revisit_overlay(
        base_table=args.base_table,
        revisit_queue=args.revisit_queue,
        revisit_labels=args.revisit_labels,
        out_dir=args.out_dir,
        min_multiclass_count=args.min_multiclass_count,
        validation_fraction=args.validation_fraction,
        test_fraction=args.test_fraction,
        random_state=args.random_state,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
