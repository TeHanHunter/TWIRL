#!/usr/bin/env python3
"""Verify a ranker-selected real-candidate review queue before human triage.

The injected-truth peak ranker produces candidate ephemerides for real S56
light curves. This verifier is the gate between that cheap queue construction
and expensive WD-tuned LEO PDF rendering.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


DEFAULT_QUEUE = (
    Path(__file__).resolve().parents[2]
    / "reports/stage5_validation/s56_ranker_selected_real_review_queue_pdo/review_queue.csv"
)
DEFAULT_SUMMARY = DEFAULT_QUEUE.with_name("summary.json")
KNOWN_APERTURES = {
    "DET_FLUX_SML",
    "DET_FLUX",
    "DET_FLUX_LAG",
    "DET_FLUX_ADP",
    "DET_FLUX_ADP_SML",
    "DET_FLUX_ADP_LAG",
}


def _read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix in {".json", ".jsonl"}:
        return pd.read_json(path, lines=suffix == ".jsonl")
    raise ValueError(f"unsupported table format: {path}")


def _as_text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str)


def _blank_or_missing(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df:
        return pd.Series(True, index=df.index)
    return _as_text(df[col]).eq("")


def _finite_positive(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df:
        return pd.Series(False, index=df.index)
    value = pd.to_numeric(df[col], errors="coerce")
    return np.isfinite(value) & value.gt(0)


def _finite_any(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df:
        return pd.Series(False, index=df.index)
    value = pd.to_numeric(df[col], errors="coerce")
    return np.isfinite(value)


def _count_bad(mask: pd.Series) -> int:
    return int((~mask.fillna(False)).sum())


def _failures_for_queue(
    queue: pd.DataFrame,
    *,
    min_rows: int,
    expect_real: int,
    expect_injected: int,
    expect_skip_leo: bool,
    require_reports: bool,
    reports_dir: Path | None,
    summary: dict[str, Any] | None,
) -> list[str]:
    failures: list[str] = []
    if len(queue) < min_rows:
        failures.append(f"queue has {len(queue)} rows; expected at least {min_rows}")

    required_columns = (
        "review_id",
        "source_kind",
        "source_bucket",
        "tic",
        "period_d",
        "t0_bjd",
        "duration_min",
        "rep_aperture",
        "sde_max",
        "recovery_status",
        "leo_class",
        "leo_report_name",
    )
    missing = [col for col in required_columns if col not in queue]
    if missing:
        failures.append("missing required columns: " + ", ".join(missing))
        return failures

    source_kind = _as_text(queue["source_kind"])
    n_real = int(source_kind.eq("real_candidate").sum())
    n_injected = int(source_kind.eq("injection_recovery").sum())
    if expect_real >= 0 and n_real != expect_real:
        failures.append(f"real_candidate rows = {n_real}; expected {expect_real}")
    if expect_injected >= 0 and n_injected != expect_injected:
        failures.append(f"injection_recovery rows = {n_injected}; expected {expect_injected}")
    if n_injected:
        failures.append("ranker-selected real queue should not contain injection rows")

    review_id = _as_text(queue["review_id"])
    if int(review_id.duplicated().sum()):
        failures.append(f"duplicate review_id rows: {int(review_id.duplicated().sum())}")
    bad_review_id = ~review_id.str.contains(r"^real:\d+:ranker:\d+$", regex=True)
    if bool(bad_review_id.any()):
        failures.append(f"{int(bad_review_id.sum())} review_id values are not real:<tic>:ranker:<rank>")

    recovery_status = _as_text(queue["recovery_status"])
    bad_recovery = ~recovery_status.eq("real_ranker_selected")
    if bool(bad_recovery.any()):
        failures.append(f"{int(bad_recovery.sum())} rows are not recovery_status=real_ranker_selected")

    source_bucket = _as_text(queue["source_bucket"])
    bad_source_bucket = ~source_bucket.eq("real_ranker_selected")
    if bool(bad_source_bucket.any()):
        failures.append(f"{int(bad_source_bucket.sum())} rows are not source_bucket=real_ranker_selected")

    for col in ("tic", "period_d", "t0_bjd", "sde_max"):
        bad = _count_bad(_finite_any(queue, col))
        if bad:
            failures.append(f"{bad} rows have non-finite {col}")
    bad_duration = _count_bad(_finite_positive(queue, "duration_min"))
    if bad_duration:
        failures.append(f"{bad_duration} rows have non-positive/non-finite duration_min")

    aperture = _as_text(queue["rep_aperture"])
    bad_aperture = ~aperture.isin(KNOWN_APERTURES)
    if bool(bad_aperture.any()):
        failures.append(f"{int(bad_aperture.sum())} rows have unknown rep_aperture")

    for col in ("injection_id", "truth_period_d", "truth_t0_bjd", "truth_radius_rearth"):
        bad_truth = ~_blank_or_missing(queue, col)
        if bool(bad_truth.any()):
            failures.append(f"{int(bad_truth.sum())} real rows have non-empty {col}")
    if "truth_source_kind" in queue:
        bad_truth_kind = ~_as_text(queue["truth_source_kind"]).eq("real_candidate")
        if bool(bad_truth_kind.any()):
            failures.append(f"{int(bad_truth_kind.sum())} rows are not truth_source_kind=real_candidate")

    if "ranker_selection_rank" in queue:
        rank = pd.to_numeric(queue["ranker_selection_rank"], errors="coerce")
        bad_rank = ~(np.isfinite(rank) & rank.ge(1))
        if bool(bad_rank.any()):
            failures.append(f"{int(bad_rank.sum())} rows have invalid ranker_selection_rank")
        if "tic" in queue:
            duplicate_units = queue.assign(_rank=rank.astype("Int64")).duplicated(["tic", "_rank"]).sum()
            if int(duplicate_units):
                failures.append(f"duplicate tic/ranker_selection_rank units: {int(duplicate_units)}")

    leo_names = _as_text(queue["leo_report_name"])
    n_named_reports = int(leo_names.ne("").sum())
    if expect_skip_leo and n_named_reports:
        failures.append(f"expected skip-LEO queue but {n_named_reports} rows have leo_report_name")
    if require_reports:
        if reports_dir is None:
            failures.append("--require-reports needs --reports-dir")
        elif n_named_reports != len(queue):
            failures.append(f"rows with leo_report_name = {n_named_reports}; expected {len(queue)}")
        elif reports_dir is not None:
            missing_reports = [name for name in leo_names if name and not (reports_dir / name).exists()]
            if missing_reports:
                failures.append(f"{len(missing_reports)} referenced LEO reports are missing; first={missing_reports[0]}")

    if summary is not None:
        if summary.get("real_selection") != "ranker_selected":
            failures.append(f"summary real_selection={summary.get('real_selection')!r}; expected ranker_selected")
        if int(summary.get("n_review_rows", -1)) != len(queue):
            failures.append(f"summary n_review_rows={summary.get('n_review_rows')}; queue has {len(queue)}")
        if expect_real >= 0 and int(summary.get("n_real", -1)) != expect_real:
            failures.append(f"summary n_real={summary.get('n_real')}; expected {expect_real}")
        if expect_injected >= 0 and int(summary.get("n_injections", -1)) != expect_injected:
            failures.append(f"summary n_injections={summary.get('n_injections')}; expected {expect_injected}")
        if expect_skip_leo and int(summary.get("n_leo_reports_attempted", -1)) != 0:
            failures.append(
                "expected skip-LEO summary but n_leo_reports_attempted="
                f"{summary.get('n_leo_reports_attempted')}"
            )
        if int(summary.get("n_leo_errors", 0)) > 0:
            failures.append(f"summary n_leo_errors={summary.get('n_leo_errors')}; expected 0")
        if int(summary.get("n_leo_plot_errors", 0)) > 0:
            failures.append(f"summary n_leo_plot_errors={summary.get('n_leo_plot_errors')}; expected 0")

    return failures


def verify_ranker_queue(
    *,
    queue_path: Path,
    summary_path: Path | None = None,
    reports_dir: Path | None = None,
    min_rows: int = 1,
    expect_real: int = -1,
    expect_injected: int = 0,
    expect_skip_leo: bool = False,
    require_reports: bool = False,
) -> dict[str, Any]:
    queue = _read_table(queue_path)
    summary = json.loads(summary_path.read_text()) if summary_path is not None and summary_path.exists() else None
    failures = _failures_for_queue(
        queue,
        min_rows=min_rows,
        expect_real=expect_real,
        expect_injected=expect_injected,
        expect_skip_leo=expect_skip_leo,
        require_reports=require_reports,
        reports_dir=reports_dir,
        summary=summary,
    )
    payload = {
        "queue": str(queue_path),
        "summary_json": str(summary_path) if summary_path is not None else "",
        "reports_dir": str(reports_dir) if reports_dir is not None else "",
        "n_rows": int(len(queue)),
        "source_kind_counts": _as_text(queue.get("source_kind", pd.Series("", index=queue.index))).value_counts().to_dict(),
        "recovery_status_counts": _as_text(queue.get("recovery_status", pd.Series("", index=queue.index))).value_counts().to_dict(),
        "leo_class_counts": _as_text(queue.get("leo_class", pd.Series("", index=queue.index))).value_counts().to_dict(),
        "n_unique_tic": int(queue["tic"].nunique()) if "tic" in queue else 0,
        "n_unique_review_id": int(queue["review_id"].nunique()) if "review_id" in queue else 0,
        "expect_skip_leo": bool(expect_skip_leo),
        "require_reports": bool(require_reports),
        "passed": not failures,
        "failures": failures,
    }
    if summary is not None:
        payload["run_summary"] = summary
    return payload


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue", type=Path, default=DEFAULT_QUEUE)
    parser.add_argument("--summary-json", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--reports-dir", type=Path, default=None)
    parser.add_argument("--out-json", type=Path, default=None)
    parser.add_argument("--min-rows", type=int, default=1)
    parser.add_argument("--expect-real", type=int, default=-1)
    parser.add_argument("--expect-injected", type=int, default=0)
    parser.add_argument("--expect-skip-leo", action="store_true")
    parser.add_argument("--require-reports", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    result = verify_ranker_queue(
        queue_path=args.queue,
        summary_path=args.summary_json,
        reports_dir=args.reports_dir,
        min_rows=args.min_rows,
        expect_real=args.expect_real,
        expect_injected=args.expect_injected,
        expect_skip_leo=args.expect_skip_leo,
        require_reports=args.require_reports,
    )
    if args.out_json is not None:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
