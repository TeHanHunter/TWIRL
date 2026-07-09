#!/usr/bin/env python3
"""Build the S56 recovery50 human-training table with separate truth roles."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.vetting.recovery50_teacher import (  # noqa: E402
    LabelPolicy,
    add_display_ephemeris,
    add_deterministic_splits,
    add_label_roles,
    join_queue_labels,
    json_default,
)

DEFAULT_QUEUE = REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue/review_queue_1k.csv"
DEFAULT_LABELS = DEFAULT_QUEUE.with_name("human_labels_vetted.csv")
DEFAULT_OUT_DIR = DEFAULT_QUEUE.with_name("human_training_table")
DEFAULT_METRICS = (
    DEFAULT_QUEUE.with_name("twirl_vet_metrics_real_fullphase_binmatch.csv"),
    DEFAULT_QUEUE.with_name("twirl_vet_metrics_injected_fullphase_binmatch.csv"),
)


def build_training_table(
    *,
    queue_csv: Path,
    labels_csv: Path,
    out_dir: Path,
    min_multiclass_count: int,
    validation_fraction: float,
    test_fraction: float,
    random_state: int,
    metrics_tables: tuple[Path, ...] = (),
) -> dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    table = join_queue_labels(queue_csv, labels_csv)
    table = add_display_ephemeris(table, metrics_tables=metrics_tables)
    table = add_label_roles(table, LabelPolicy(min_multiclass_count=min_multiclass_count))
    table = add_deterministic_splits(
        table,
        validation_fraction=validation_fraction,
        test_fraction=test_fraction,
        random_state=random_state,
    )
    teacher = table[table["main_teacher_include"].fillna(False).astype(bool)].copy()
    audit = table[table["audit_include"].fillna(False).astype(bool)].copy()
    table.to_csv(out_dir / "human_vetting_training_table.csv", index=False)
    teacher.to_csv(out_dir / "teacher_labeled_rows.csv", index=False)
    audit.to_csv(out_dir / "audit_labeled_rows.csv", index=False)

    summary = {
        "queue_csv": str(queue_csv),
        "labels_csv": str(labels_csv),
        "metrics_tables": [str(path) for path in metrics_tables if path.exists()],
        "out_dir": str(out_dir),
        "n_rows": int(len(table)),
        "n_labeled": int(table["is_labeled"].sum()),
        "n_teacher_rows": int(table["main_teacher_include"].sum()),
        "n_audit_rows": int(table["audit_include"].sum()),
        "label_counts": {
            str(k): int(v)
            for k, v in table.loc[table["is_labeled"], "human_label"].value_counts().sort_index().items()
        },
        "teacher_target_counts": {
            str(k): int(v)
            for k, v in teacher["main_teacher_target"].value_counts().sort_index().items()
        },
        "teacher_source_kind_counts": {
            str(k): int(v) for k, v in teacher["source_kind"].value_counts().sort_index().items()
        },
        "training_split_counts": {
            str(k): int(v) for k, v in table["training_split"].value_counts().sort_index().items()
        },
        "display_ephemeris_source_counts": {
            str(k): int(v) for k, v in table["display_ephemeris_source"].value_counts().sort_index().items()
        },
        "n_display_ephemeris_anchor": int(table["display_ephemeris_used_anchor"].sum()),
        "outputs": {
            "training_table": str(out_dir / "human_vetting_training_table.csv"),
            "teacher_rows": str(out_dir / "teacher_labeled_rows.csv"),
            "audit_rows": str(out_dir / "audit_labeled_rows.csv"),
        },
        "label_policy": "main human-visible teacher target; injection truth auxiliary/evaluation only",
        "min_multiclass_count": int(min_multiclass_count),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=json_default) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True, default=json_default))
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue-csv", type=Path, default=DEFAULT_QUEUE)
    parser.add_argument("--labels-csv", type=Path, default=DEFAULT_LABELS)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--min-multiclass-count", type=int, default=40)
    parser.add_argument("--validation-fraction", type=float, default=0.20)
    parser.add_argument("--test-fraction", type=float, default=0.20)
    parser.add_argument("--random-state", type=int, default=56)
    parser.add_argument("--metrics-table", type=Path, action="append", default=list(DEFAULT_METRICS))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    build_training_table(
        queue_csv=args.queue_csv,
        labels_csv=args.labels_csv,
        out_dir=args.out_dir,
        min_multiclass_count=args.min_multiclass_count,
        validation_fraction=args.validation_fraction,
        test_fraction=args.test_fraction,
        random_state=args.random_state,
        metrics_tables=tuple(path for path in args.metrics_table if path.exists()),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
