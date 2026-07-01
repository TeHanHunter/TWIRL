#!/usr/bin/env python3
"""Audit whether human labels are ready for teacher-model training.

The joined human-vetting table can exist well before it is scientifically
usable for a model. This audit keeps three gates separate:

1. injection-visibility smoke: enough labels to ask whether humans see injected
   signals under the current BLS/LEO evidence surface;
2. binary teacher smoke: enough strong labels to run a small visible-signal vs
   not-visible smoke model;
3. object teacher training: enough real-candidate labels across the false
   positive taxonomy to train the astrophysical/object-class teacher.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
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
    BLS_TRUTH_MATCH_MODES,
    OBJECT_TEACHER_LABELS,
    VISIBILITY_NEGATIVE_LABELS,
    VISIBILITY_POSITIVE_LABELS,
)


DEFAULT_TRAINING_TABLE = (
    Path(__file__).resolve().parents[2]
    / "reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_pdo/"
    / "human_training_table_quicktest/human_vetting_training_table.csv"
)
DEFAULT_PRIORITY_TABLE = DEFAULT_TRAINING_TABLE.parent.with_name("human_label_priority_next") / "next_label_priority.csv"
DEFAULT_OUT_DIR = DEFAULT_TRAINING_TABLE.parent.with_name("human_training_readiness_quicktest")

DEFAULT_REQUIRED_OBJECT_CLASSES = OBJECT_TEACHER_LABELS

TMAG_BINS = (-np.inf, 17.0, 18.0, 19.0, np.inf)
TMAG_LABELS = ("<17", "17-18", "18-19", ">19")
PERIOD_BINS = (0.0, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, np.inf)
PERIOD_LABELS = ("<0.25", "0.25-0.5", "0.5-1", "1-2", "2-5", "5-10", ">10")
RADIUS_BINS = (0.0, 1.0, 2.0, 4.0, 8.0, 12.0, np.inf)
RADIUS_LABELS = ("<1", "1-2", "2-4", "4-8", "8-12", ">12")
COVERAGE_COLUMNS = ("topn_recovery_status", "leo_class", "tmag_bin", "period_bin", "radius_bin")


@dataclass(frozen=True)
class ReadinessConfig:
    min_smoke_labeled: int = 100
    min_smoke_positive: int = 20
    min_smoke_negative: int = 20
    min_smoke_bls_truth_match: int = 10
    min_smoke_bls_mismatch: int = 25
    min_binary_teacher_total: int = 150
    min_binary_teacher_per_class: int = 40
    min_teacher_total: int = 300
    min_teacher_per_class: int = 40
    min_train_per_class: int = 25
    min_validation_per_class: int = 5
    min_test_per_class: int = 5
    min_real_teacher: int = 100
    target_per_coverage_cell: int = 25
    required_object_classes: tuple[str, ...] = DEFAULT_REQUIRED_OBJECT_CLASSES


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


def _as_bool(series: pd.Series | bool, index: pd.Index) -> pd.Series:
    if isinstance(series, bool):
        return pd.Series(series, index=index)
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    text = series.fillna("").astype(str).str.strip().str.lower()
    return text.isin({"true", "1", "yes", "y"})


def _safe_series(df: pd.DataFrame, col: str, default: str = "") -> pd.Series:
    if col in df:
        return df[col]
    return pd.Series(default, index=df.index)


def _int_counts(series: pd.Series) -> dict[str, int]:
    return {str(k): int(v) for k, v in series.fillna("").astype(str).value_counts().sort_index().items()}


def _parse_csv_tuple(raw: str) -> tuple[str, ...]:
    return tuple(part.strip() for part in str(raw).split(",") if part.strip())


def _add_bins(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "tmag_bin" not in out and "tmag" in out:
        out["tmag_bin"] = pd.cut(pd.to_numeric(out["tmag"], errors="coerce"), TMAG_BINS, labels=TMAG_LABELS).astype(object)
    if "period_bin" not in out:
        period_col = "truth_period_d" if "truth_period_d" in out else "period_d"
        if period_col in out:
            out["period_bin"] = pd.cut(
                pd.to_numeric(out[period_col], errors="coerce"),
                PERIOD_BINS,
                labels=PERIOD_LABELS,
            ).astype(object)
    if "radius_bin" not in out and "truth_radius_rearth" in out:
        out["radius_bin"] = pd.cut(
            pd.to_numeric(out["truth_radius_rearth"], errors="coerce"),
            RADIUS_BINS,
            labels=RADIUS_LABELS,
        ).astype(object)
    return out


def _append_requirement(
    rows: list[dict[str, Any]],
    *,
    gate: str,
    dimension: str,
    value: str,
    current_n: int,
    required_n: int,
    note: str = "",
) -> bool:
    deficit = max(0, int(required_n) - int(current_n))
    passed = deficit == 0
    rows.append(
        {
            "gate": gate,
            "dimension": dimension,
            "value": value,
            "current_n": int(current_n),
            "required_n": int(required_n),
            "deficit": int(deficit),
            "status": "pass" if passed else "fail",
            "note": note,
        }
    )
    return passed


def _split_counts(teacher: pd.DataFrame) -> pd.DataFrame:
    if teacher.empty:
        return pd.DataFrame()
    return pd.crosstab(
        teacher["teacher_target"].fillna("").astype(str),
        teacher["training_split"].fillna("").astype(str),
    )


def _gate_status(deficits: pd.DataFrame, gate: str) -> str:
    subset = deficits[deficits["gate"].eq(gate)]
    if subset.empty:
        return "not_applicable"
    return "pass" if subset["deficit"].sum() == 0 else "fail"


def audit_training_readiness(
    table: pd.DataFrame,
    *,
    priority_table: pd.DataFrame | None = None,
    config: ReadinessConfig = ReadinessConfig(),
) -> tuple[dict[str, Any], pd.DataFrame, str]:
    """Return summary, requirement/coverage deficits, and a markdown report."""

    table = _add_bins(table)
    is_labeled = _as_bool(_safe_series(table, "is_labeled", False), table.index)
    teacher_include = _as_bool(_safe_series(table, "teacher_include", False), table.index)
    audit_include = _as_bool(_safe_series(table, "audit_include", False), table.index)

    source_kind = _safe_series(table, "source_kind").fillna("").astype(str)
    teacher_target = _safe_series(table, "teacher_target").fillna("").astype(str)
    human_label = _safe_series(table, "human_label").fillna("").astype(str)
    split = _safe_series(table, "training_split", "unlabeled_or_audit").fillna("").astype(str)
    topn = _safe_series(table, "topn_recovery_status").fillna("").astype(str)
    if not topn.ne("").any() and "recovery_status" in table:
        topn = table["recovery_status"].fillna("").astype(str)

    work = table.copy()
    work["source_kind"] = source_kind
    work["teacher_target"] = teacher_target
    work["human_label"] = human_label
    work["training_split"] = split
    work["_is_labeled"] = is_labeled
    work["_teacher_include"] = teacher_include
    work["_audit_include"] = audit_include
    work["_topn_recovery_status"] = topn
    if "bls_truth_match" in work:
        work["_bls_truth_match"] = _as_bool(work["bls_truth_match"], work.index)
    else:
        work["_bls_truth_match"] = topn.isin(BLS_TRUTH_MATCH_MODES)

    labeled = work[work["_is_labeled"]].copy()
    teacher = work[work["_teacher_include"]].copy()
    audit = work[work["_audit_include"]].copy()
    real_teacher = teacher[teacher["source_kind"].eq("real_candidate")]
    injected_teacher = teacher[teacher["source_kind"].eq("injection_recovery")]

    rows: list[dict[str, Any]] = []

    smoke_ok = []
    smoke_ok.append(
        _append_requirement(
            rows,
            gate="injection_visibility_smoke",
            dimension="labeled_rows",
            value="all",
            current_n=len(labeled),
            required_n=config.min_smoke_labeled,
            note="Enough rows to ask whether injected events are visually recognizable.",
        )
    )
    smoke_positive_n = int(teacher["teacher_target"].isin(VISIBILITY_POSITIVE_LABELS).sum())
    smoke_negative_n = int(teacher["teacher_target"].isin(VISIBILITY_NEGATIVE_LABELS).sum())
    smoke_ok.append(
        _append_requirement(
            rows,
            gate="injection_visibility_smoke",
            dimension="visibility_class",
            value="positive_planet_like",
            current_n=smoke_positive_n,
            required_n=config.min_smoke_positive,
        )
    )
    smoke_ok.append(
        _append_requirement(
            rows,
            gate="injection_visibility_smoke",
            dimension="visibility_class",
            value="negative_instrumental_or_systematic",
            current_n=smoke_negative_n,
            required_n=config.min_smoke_negative,
        )
    )
    smoke_ok.append(
        _append_requirement(
            rows,
            gate="injection_visibility_smoke",
            dimension="bls_truth_match",
            value="true",
            current_n=int(labeled["_bls_truth_match"].sum()),
            required_n=config.min_smoke_bls_truth_match,
        )
    )
    smoke_ok.append(
        _append_requirement(
            rows,
            gate="injection_visibility_smoke",
            dimension="bls_truth_match",
            value="false",
            current_n=int((~labeled["_bls_truth_match"]).sum()),
            required_n=config.min_smoke_bls_mismatch,
        )
    )

    binary_ok = []
    binary_ok.append(
        _append_requirement(
            rows,
            gate="binary_teacher_smoke",
            dimension="teacher_rows",
            value="all",
            current_n=len(teacher),
            required_n=config.min_binary_teacher_total,
        )
    )
    for label in (*VISIBILITY_POSITIVE_LABELS, *VISIBILITY_NEGATIVE_LABELS):
        binary_ok.append(
            _append_requirement(
                rows,
                gate="binary_teacher_smoke",
                dimension="teacher_target",
                value=label,
                current_n=int(teacher["teacher_target"].eq(label).sum()),
                required_n=config.min_binary_teacher_per_class,
            )
        )

    object_ok = []
    object_ok.append(
        _append_requirement(
            rows,
            gate="object_teacher_training",
            dimension="teacher_rows",
            value="all",
            current_n=len(teacher),
            required_n=config.min_teacher_total,
            note="Final teacher needs a broader, more stable label set than the 1k quick visual check.",
        )
    )
    object_ok.append(
        _append_requirement(
            rows,
            gate="object_teacher_training",
            dimension="real_candidate_teacher_rows",
            value="all",
            current_n=len(real_teacher),
            required_n=config.min_real_teacher,
            note="Real false positives cannot be replaced by synthetic injections.",
        )
    )
    for label in config.required_object_classes:
        object_ok.append(
            _append_requirement(
                rows,
                gate="object_teacher_training",
                dimension="teacher_target",
                value=label,
                current_n=int(teacher["teacher_target"].eq(label).sum()),
                required_n=config.min_teacher_per_class,
            )
        )
    splits = _split_counts(teacher)
    for label in config.required_object_classes:
        label_splits = splits.loc[label] if label in splits.index else pd.Series(dtype=int)
        for split_name, required in (
            ("train", config.min_train_per_class),
            ("validation", config.min_validation_per_class),
            ("test", config.min_test_per_class),
        ):
            object_ok.append(
                _append_requirement(
                    rows,
                    gate="object_teacher_training",
                    dimension=f"training_split:{split_name}",
                    value=label,
                    current_n=int(label_splits.get(split_name, 0)),
                    required_n=required,
                )
            )

    coverage_rows = []
    for col in COVERAGE_COLUMNS:
        if col not in work:
            continue
        labeled_counts = _int_counts(labeled[col])
        values = sorted({str(v) for v in work[col].fillna("").astype(str) if str(v)})
        for value in values:
            _append_requirement(
                coverage_rows,
                gate="label_coverage",
                dimension=col,
                value=value,
                current_n=int(labeled_counts.get(value, 0)),
                required_n=config.target_per_coverage_cell,
                note="Coverage target for selecting the next human-label batch.",
            )
    rows.extend(coverage_rows)
    deficits = pd.DataFrame(rows)

    gates = {
        "injection_visibility_smoke": {
            "status": _gate_status(deficits, "injection_visibility_smoke"),
            "passed": bool(all(smoke_ok)),
        },
        "binary_teacher_smoke": {
            "status": _gate_status(deficits, "binary_teacher_smoke"),
            "passed": bool(all(binary_ok)),
        },
        "object_teacher_training": {
            "status": _gate_status(deficits, "object_teacher_training"),
            "passed": bool(all(object_ok)),
        },
    }
    if gates["object_teacher_training"]["passed"]:
        recommendation = "ready_for_object_teacher_training"
    elif gates["binary_teacher_smoke"]["passed"]:
        recommendation = "ready_for_binary_teacher_smoke"
    elif gates["injection_visibility_smoke"]["passed"]:
        recommendation = "ready_for_injection_visibility_smoke_only"
    else:
        recommendation = "not_ready_label_more"

    priority_summary: dict[str, Any] = {}
    if priority_table is not None:
        priority_summary = {
            "n_priority_rows": int(len(priority_table)),
            "priority_by_topn_recovery_status": _int_counts(priority_table["topn_recovery_status"])
            if "topn_recovery_status" in priority_table
            else {},
            "priority_by_leo_class": _int_counts(priority_table["leo_class"]) if "leo_class" in priority_table else {},
        }

    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "recommendation": recommendation,
        "method_note": (
            "Injected truth is appropriate for BLS/ephemeris ranker training. "
            "Human labels from real candidates are still required for the final object-class teacher."
        ),
        "n_rows": int(len(work)),
        "n_labeled": int(len(labeled)),
        "n_teacher_rows": int(len(teacher)),
        "n_audit_rows": int(len(audit)),
        "n_real_teacher_rows": int(len(real_teacher)),
        "n_injected_teacher_rows": int(len(injected_teacher)),
        "label_counts": _int_counts(labeled["human_label"]) if len(labeled) else {},
        "teacher_target_counts": _int_counts(teacher["teacher_target"]) if len(teacher) else {},
        "source_kind_counts": _int_counts(work["source_kind"]),
        "teacher_source_kind_counts": _int_counts(teacher["source_kind"]) if len(teacher) else {},
        "training_split_counts": _int_counts(work["training_split"]),
        "teacher_split_by_target": {
            str(index): {str(col): int(value) for col, value in row.items()}
            for index, row in splits.to_dict(orient="index").items()
        },
        "labeled_by_bls_truth_match": {
            "true": int(labeled["_bls_truth_match"].sum()),
            "false": int((~labeled["_bls_truth_match"]).sum()),
        },
        "labeled_by_topn_recovery_status": _int_counts(labeled["_topn_recovery_status"]) if len(labeled) else {},
        "labeled_by_leo_class": _int_counts(labeled["leo_class"]) if "leo_class" in labeled else {},
        "gates": gates,
        "thresholds": config.__dict__,
        "priority_summary": priority_summary,
    }
    markdown = _render_markdown(summary, deficits)
    return summary, deficits, markdown


def _render_markdown(summary: dict[str, Any], deficits: pd.DataFrame) -> str:
    lines = [
        "# Human Training Readiness Audit",
        "",
        f"Recommendation: `{summary['recommendation']}`.",
        "",
        summary["method_note"],
        "",
        "## Counts",
        "",
        f"- labeled rows: `{summary['n_labeled']}` / `{summary['n_rows']}`",
        f"- strong teacher rows: `{summary['n_teacher_rows']}`",
        f"- audit rows including uncertain labels: `{summary['n_audit_rows']}`",
        f"- real-candidate teacher rows: `{summary['n_real_teacher_rows']}`",
        f"- injected teacher rows: `{summary['n_injected_teacher_rows']}`",
        "",
        "## Gates",
        "",
        "| Gate | Status |",
        "|---|---|",
    ]
    for gate, payload in summary["gates"].items():
        lines.append(f"| `{gate}` | `{payload['status']}` |")
    lines.extend(["", "## Teacher Targets", "", "| Label | Count |", "|---|---:|"])
    for label, count in summary["teacher_target_counts"].items():
        lines.append(f"| `{label}` | {count} |")

    failing = deficits[(deficits["status"].eq("fail")) & (~deficits["gate"].eq("label_coverage"))].copy()
    if len(failing):
        lines.extend(["", "## Blocking Deficits", "", "| Gate | Requirement | Current | Required | Deficit |", "|---|---|---:|---:|---:|"])
        for _, row in failing.iterrows():
            req = f"{row['dimension']}={row['value']}"
            lines.append(
                f"| `{row['gate']}` | `{req}` | {int(row['current_n'])} | {int(row['required_n'])} | {int(row['deficit'])} |"
            )

    coverage = deficits[(deficits["gate"].eq("label_coverage")) & (deficits["status"].eq("fail"))].copy()
    if len(coverage):
        coverage = coverage.sort_values(["deficit", "dimension", "value"], ascending=[False, True, True]).head(15)
        lines.extend(["", "## Largest Coverage Deficits", "", "| Dimension | Value | Current | Target | Deficit |", "|---|---|---:|---:|---:|"])
        for _, row in coverage.iterrows():
            lines.append(
                f"| `{row['dimension']}` | `{row['value']}` | {int(row['current_n'])} | {int(row['required_n'])} | {int(row['deficit'])} |"
            )
    lines.append("")
    return "\n".join(lines)


def build_readiness_audit(
    *,
    training_table: Path,
    out_dir: Path,
    priority_table: Path | None = None,
    config: ReadinessConfig = ReadinessConfig(),
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    table = _read_table(training_table)
    priority = _read_table(priority_table) if priority_table is not None and priority_table.exists() else None
    summary, deficits, markdown = audit_training_readiness(table, priority_table=priority, config=config)
    summary = {
        **summary,
        "training_table": str(training_table),
        "priority_table": str(priority_table) if priority_table is not None and priority_table.exists() else "",
        "out_dir": str(out_dir),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    deficits.to_csv(out_dir / "label_deficits.csv", index=False)
    (out_dir / "summary.md").write_text(markdown)
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training-table", type=Path, default=DEFAULT_TRAINING_TABLE)
    parser.add_argument("--priority-table", type=Path, default=DEFAULT_PRIORITY_TABLE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--min-smoke-labeled", type=int, default=ReadinessConfig.min_smoke_labeled)
    parser.add_argument("--min-smoke-positive", type=int, default=ReadinessConfig.min_smoke_positive)
    parser.add_argument("--min-smoke-negative", type=int, default=ReadinessConfig.min_smoke_negative)
    parser.add_argument("--min-smoke-bls-truth-match", type=int, default=ReadinessConfig.min_smoke_bls_truth_match)
    parser.add_argument("--min-smoke-bls-mismatch", type=int, default=ReadinessConfig.min_smoke_bls_mismatch)
    parser.add_argument("--min-binary-teacher-total", type=int, default=ReadinessConfig.min_binary_teacher_total)
    parser.add_argument("--min-binary-teacher-per-class", type=int, default=ReadinessConfig.min_binary_teacher_per_class)
    parser.add_argument("--min-teacher-total", type=int, default=ReadinessConfig.min_teacher_total)
    parser.add_argument("--min-teacher-per-class", type=int, default=ReadinessConfig.min_teacher_per_class)
    parser.add_argument("--min-train-per-class", type=int, default=ReadinessConfig.min_train_per_class)
    parser.add_argument("--min-validation-per-class", type=int, default=ReadinessConfig.min_validation_per_class)
    parser.add_argument("--min-test-per-class", type=int, default=ReadinessConfig.min_test_per_class)
    parser.add_argument("--min-real-teacher", type=int, default=ReadinessConfig.min_real_teacher)
    parser.add_argument("--target-per-coverage-cell", type=int, default=ReadinessConfig.target_per_coverage_cell)
    parser.add_argument("--required-object-classes", default=",".join(DEFAULT_REQUIRED_OBJECT_CLASSES))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    config = ReadinessConfig(
        min_smoke_labeled=args.min_smoke_labeled,
        min_smoke_positive=args.min_smoke_positive,
        min_smoke_negative=args.min_smoke_negative,
        min_smoke_bls_truth_match=args.min_smoke_bls_truth_match,
        min_smoke_bls_mismatch=args.min_smoke_bls_mismatch,
        min_binary_teacher_total=args.min_binary_teacher_total,
        min_binary_teacher_per_class=args.min_binary_teacher_per_class,
        min_teacher_total=args.min_teacher_total,
        min_teacher_per_class=args.min_teacher_per_class,
        min_train_per_class=args.min_train_per_class,
        min_validation_per_class=args.min_validation_per_class,
        min_test_per_class=args.min_test_per_class,
        min_real_teacher=args.min_real_teacher,
        target_per_coverage_cell=args.target_per_coverage_cell,
        required_object_classes=_parse_csv_tuple(args.required_object_classes),
    )
    priority_table = args.priority_table if args.priority_table.exists() else None
    build_readiness_audit(
        training_table=args.training_table,
        priority_table=priority_table,
        out_dir=args.out_dir,
        config=config,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
