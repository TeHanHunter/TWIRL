#!/usr/bin/env python3
"""Select high-value unlabeled rows for the next human-vetting pass.

This script does not modify the review queue or label CSV. It consumes the
joined human-vetting training table and writes a priority list of unlabeled
rows that improve coverage over BLS recovery mode, LEO class, magnitude, period,
and radius. Hidden truth columns are used only for planning the injected
training sample, not as browser-visible labels.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


DEFAULT_TRAINING_TABLE = (
    Path(__file__).resolve().parents[2]
    / "reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_pdo/"
    / "human_training_table_quicktest/human_vetting_training_table.csv"
)
DEFAULT_OUT_DIR = DEFAULT_TRAINING_TABLE.parent.with_name("human_label_priority_next")

TMAG_BINS = (-np.inf, 17.0, 18.0, 19.0, np.inf)
TMAG_LABELS = ("<17", "17-18", "18-19", ">19")
PERIOD_BINS = (0.0, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, np.inf)
PERIOD_LABELS = ("<0.25", "0.25-0.5", "0.5-1", "1-2", "2-5", "5-10", ">10")
RADIUS_BINS = (0.0, 1.0, 2.0, 4.0, 8.0, 12.0, np.inf)
RADIUS_LABELS = ("<1", "1-2", "2-4", "4-8", "8-12", ">12")
DEFAULT_COVERAGE_COLUMNS = (
    "topn_recovery_status",
    "leo_class",
    "tmag_bin",
    "period_bin",
    "radius_bin",
)
DEFAULT_KEEP_COLUMNS = (
    "priority_rank",
    "priority_score",
    "priority_reasons",
    "row_id",
    "review_id",
    "tic",
    "sector",
    "tmag",
    "period_d",
    "duration_min",
    "sde_max",
    "leo_class",
    "leo_report_name",
    "topn_recovery_status",
    "recovery_status",
    "truth_period_d",
    "truth_radius_rearth",
    "truth_sampled_model_depth",
    "truth_duration_min",
    "teacher_target",
    "human_label",
)


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


def add_priority_bins(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "tmag_bin" not in out and "tmag" in out:
        out["tmag_bin"] = pd.cut(
            pd.to_numeric(out["tmag"], errors="coerce"),
            TMAG_BINS,
            labels=TMAG_LABELS,
        ).astype(object)
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


def _coverage_deficits(
    table: pd.DataFrame,
    labeled: pd.DataFrame,
    *,
    columns: tuple[str, ...],
    target_per_cell: int,
) -> tuple[dict[str, dict[str, int]], pd.DataFrame]:
    deficits: dict[str, dict[str, int]] = {}
    rows = []
    for col in columns:
        if col not in table:
            continue
        counts = (
            labeled[col].fillna("").astype(str).value_counts().to_dict()
            if col in labeled
            else {}
        )
        values = sorted({str(v) for v in table[col].fillna("").astype(str).tolist() if str(v)})
        deficits[col] = {}
        for value in values:
            count = int(counts.get(value, 0))
            deficit = max(0, int(target_per_cell) - count)
            deficits[col][value] = deficit
            rows.append({"coverage_column": col, "value": value, "labeled_n": count, "deficit": deficit})
    return deficits, pd.DataFrame(rows)


def _score_row(
    row: pd.Series,
    deficits: dict[str, dict[str, int]],
    *,
    columns: tuple[str, ...],
) -> tuple[float, list[str]]:
    score = 0.0
    reasons: list[str] = []
    for col in columns:
        if col not in row:
            continue
        value = "" if pd.isna(row[col]) else str(row[col])
        deficit = int(deficits.get(col, {}).get(value, 0))
        if deficit > 0:
            score += float(deficit)
            reasons.append(f"{col}:{value}:{deficit}")

    leo = "" if "leo_class" not in row or pd.isna(row.get("leo_class")) else str(row.get("leo_class"))
    if leo in {"PC", "FP"}:
        score += 8.0
        reasons.append(f"rare_leo:{leo}")

    mode = (
        ""
        if "topn_recovery_status" not in row or pd.isna(row.get("topn_recovery_status"))
        else str(row.get("topn_recovery_status"))
    )
    if mode in {"bls_top1_recovered", "bls_topn_recovered", "bls_topn_harmonic_match"}:
        score += 5.0
        reasons.append(f"bls_truth_mode:{mode}")
    elif mode == "bls_peak_mismatch":
        score += 1.0
        reasons.append("bls_mismatch_control")

    gate = "" if "recovery_gate" not in row or pd.isna(row.get("recovery_gate")) else str(row.get("recovery_gate"))
    if gate == "below_empirical_threshold":
        score += 2.0
        reasons.append("below_empirical_threshold")
    elif gate == "above_empirical_threshold":
        score += 1.0
        reasons.append("above_empirical_threshold")

    return score, reasons


def select_next_rows(
    table: pd.DataFrame,
    *,
    n_rows: int,
    target_per_cell: int,
    coverage_columns: tuple[str, ...] = DEFAULT_COVERAGE_COLUMNS,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    table = add_priority_bins(table)
    is_labeled = table.get("is_labeled", False)
    if not isinstance(is_labeled, pd.Series):
        is_labeled = pd.Series(False, index=table.index)
    labeled = table[is_labeled.fillna(False).astype(bool)].copy()
    unlabeled = table[~is_labeled.fillna(False).astype(bool)].copy()
    deficits, coverage = _coverage_deficits(
        table,
        labeled,
        columns=coverage_columns,
        target_per_cell=target_per_cell,
    )
    scores = []
    reasons = []
    for _, row in unlabeled.iterrows():
        score, row_reasons = _score_row(row, deficits, columns=coverage_columns)
        scores.append(score)
        reasons.append(";".join(row_reasons))
    selected = unlabeled.copy()
    selected["priority_score"] = scores
    selected["priority_reasons"] = reasons
    secondary_cols = [col for col in ("sde_max", "tmag") if col in selected.columns]
    ascending = [False] + ([False] if "sde_max" in secondary_cols else []) + ([True] if "tmag" in secondary_cols else [])
    selected = selected.sort_values(
        ["priority_score", *secondary_cols],
        ascending=ascending,
        na_position="last",
        kind="stable",
    )
    selected = selected.head(int(n_rows)).copy().reset_index(drop=True)
    selected.insert(0, "priority_rank", np.arange(1, len(selected) + 1, dtype=int))
    summary = {
        "n_rows_total": int(len(table)),
        "n_labeled": int(len(labeled)),
        "n_unlabeled": int(len(unlabeled)),
        "n_selected": int(len(selected)),
        "target_per_cell": int(target_per_cell),
        "coverage_columns": list(coverage_columns),
        "selected_by_topn_recovery_status": selected["topn_recovery_status"].fillna("").astype(str).value_counts().to_dict()
        if "topn_recovery_status" in selected
        else {},
        "selected_by_leo_class": selected["leo_class"].fillna("").astype(str).value_counts().to_dict()
        if "leo_class" in selected
        else {},
        "selected_by_tmag_bin": selected["tmag_bin"].fillna("").astype(str).value_counts().to_dict()
        if "tmag_bin" in selected
        else {},
    }
    return selected, coverage, summary


def write_outputs(
    *,
    selected: pd.DataFrame,
    coverage: pd.DataFrame,
    summary: dict[str, Any],
    out_dir: Path,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    keep = [col for col in DEFAULT_KEEP_COLUMNS if col in selected.columns]
    extra = [col for col in ("tmag_bin", "period_bin", "radius_bin", "recovery_gate") if col in selected.columns]
    selected.loc[:, list(dict.fromkeys([*keep, *extra]))].to_csv(out_dir / "next_label_priority.csv", index=False)
    selected.to_csv(out_dir / "next_label_priority_full.csv", index=False)
    if "row_id" in selected:
        selected["row_id"].astype(int).to_csv(out_dir / "next_label_row_ids.txt", index=False, header=False)
    coverage.to_csv(out_dir / "label_coverage_deficits.csv", index=False)
    summary = {**summary, "created_utc": datetime.now(timezone.utc).isoformat()}
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")


def _parse_columns(raw: str) -> tuple[str, ...]:
    return tuple(part.strip() for part in str(raw).split(",") if part.strip())


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training-table", type=Path, default=DEFAULT_TRAINING_TABLE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--n-rows", type=int, default=200)
    parser.add_argument("--target-per-cell", type=int, default=25)
    parser.add_argument("--coverage-columns", default=",".join(DEFAULT_COVERAGE_COLUMNS))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    table = _read_table(args.training_table)
    selected, coverage, summary = select_next_rows(
        table,
        n_rows=args.n_rows,
        target_per_cell=args.target_per_cell,
        coverage_columns=_parse_columns(args.coverage_columns),
    )
    summary.update({"training_table": str(args.training_table), "out_dir": str(args.out_dir)})
    write_outputs(selected=selected, coverage=coverage, summary=summary, out_dir=args.out_dir)
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
