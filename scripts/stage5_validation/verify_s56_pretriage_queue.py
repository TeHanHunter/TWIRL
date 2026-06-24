#!/usr/bin/env python3
"""Verify the S56 pre-human-labeling review queue gates.

This is a lightweight guard for the checklist in ``doc/twirl_plan.md``:
human labels should start only after the queue is an end-to-end product with
multi-aperture injection recovery, WD-tuned LEO reports, and a passing WD 1856
benchmark row.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd


WD1856_TIC = 267574918
APERTURES = ("DET_FLUX_SML", "DET_FLUX", "DET_FLUX_LAG")


def _read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix == ".parquet":
        return pd.read_parquet(path)
    raise ValueError(f"unsupported table format: {path}")


def _failures_for_queue(
    queue: pd.DataFrame,
    *,
    reports_dir: Path,
    min_rows: int,
    expect_real: int,
    expect_injected: int,
    min_reports: int,
    apertures: tuple[str, ...],
    require_wd1856: bool,
) -> list[str]:
    failures: list[str] = []
    if len(queue) < min_rows:
        failures.append(f"queue has {len(queue)} rows; expected at least {min_rows}")

    source_kind = queue.get("source_kind", pd.Series("", index=queue.index)).astype(str)
    n_real = int(source_kind.eq("real_candidate").sum())
    n_inj = int(source_kind.eq("injection_recovery").sum())
    if expect_real >= 0 and n_real != expect_real:
        failures.append(f"real_candidate rows = {n_real}; expected {expect_real}")
    if expect_injected >= 0 and n_inj != expect_injected:
        failures.append(f"injection_recovery rows = {n_inj}; expected {expect_injected}")

    if require_wd1856 and "tic" not in queue.columns:
        failures.append("queue is missing tic column")
    elif require_wd1856:
        wd = queue[queue["tic"].astype(str).eq(str(WD1856_TIC))]
        if len(wd) != 1:
            failures.append(f"WD 1856 appears {len(wd)} times; expected exactly once")
        else:
            row = wd.iloc[0]
            if str(row.get("source_kind", "")) != "real_candidate":
                failures.append("WD 1856 row is not source_kind=real_candidate")
            if str(row.get("leo_class", "")) != "PC":
                failures.append(f"WD 1856 leo_class={row.get('leo_class', '')!r}; expected 'PC'")
            report_name = str(row.get("leo_report_name", ""))
            if not report_name:
                failures.append("WD 1856 is missing leo_report_name")
            elif not (reports_dir / report_name).exists():
                failures.append(f"WD 1856 report is missing on disk: {reports_dir / report_name}")

    missing_columns = []
    for aperture in apertures:
        for prefix in ("sde", "recovery_status"):
            col = f"{prefix}_{aperture}"
            if col not in queue.columns:
                missing_columns.append(col)
    if missing_columns:
        failures.append("queue is missing multi-aperture columns: " + ", ".join(missing_columns))

    injected = queue[source_kind.eq("injection_recovery")]
    if not injected.empty:
        if "rep_aperture" not in injected.columns:
            failures.append("injected rows are missing rep_aperture")
        else:
            bad = ~injected["rep_aperture"].astype(str).isin(apertures)
            if bool(bad.any()):
                failures.append(f"{int(bad.sum())} injected rows have unknown rep_aperture")
        for aperture in apertures:
            col = f"sde_{aperture}"
            if col in injected.columns and pd.to_numeric(injected[col], errors="coerce").notna().sum() == 0:
                failures.append(f"injected rows have no finite values in {col}")

    leo_names = queue.get("leo_report_name", pd.Series("", index=queue.index)).fillna("").astype(str)
    named_reports = leo_names[leo_names.ne("")]
    if len(named_reports) < min_reports:
        failures.append(f"rows with leo_report_name = {len(named_reports)}; expected at least {min_reports}")
    else:
        missing = [name for name in named_reports if not (reports_dir / name).exists()]
        if missing:
            failures.append(f"{len(missing)} referenced LEO reports are missing; first={missing[0]}")
    return failures


def _build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--queue", type=Path, required=True)
    ap.add_argument("--reports-dir", type=Path, required=True)
    ap.add_argument("--summary-json", type=Path, default=None)
    ap.add_argument("--min-rows", type=int, default=1000)
    ap.add_argument("--expect-real", type=int, default=100)
    ap.add_argument("--expect-injected", type=int, default=900)
    ap.add_argument("--min-reports", type=int, default=-1,
                    help="Minimum rows with LEO report names. Default: --min-rows.")
    ap.add_argument("--max-leo-errors", type=int, default=0)
    ap.add_argument("--apertures", nargs="+", default=list(APERTURES),
                    help="Apertures expected in per-aperture BLS columns.")
    ap.add_argument("--skip-wd1856-check", action="store_true",
                    help="Use for injection-only queues with no real WD 1856 benchmark row.")
    ap.add_argument("--out-json", type=Path, default=None)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    queue = _read_table(args.queue)
    min_reports = args.min_rows if args.min_reports < 0 else args.min_reports
    failures = _failures_for_queue(
        queue,
        reports_dir=args.reports_dir,
        min_rows=args.min_rows,
        expect_real=args.expect_real,
        expect_injected=args.expect_injected,
        min_reports=min_reports,
        apertures=tuple(args.apertures),
        require_wd1856=not args.skip_wd1856_check,
    )
    run_summary = {}
    if args.summary_json is not None and args.summary_json.exists():
        run_summary = json.loads(args.summary_json.read_text())
        n_leo_errors = int(run_summary.get("n_leo_errors", 0))
        if n_leo_errors > args.max_leo_errors:
            failures.append(f"n_leo_errors={n_leo_errors}; max allowed is {args.max_leo_errors}")
    summary = {
        "queue": str(args.queue),
        "reports_dir": str(args.reports_dir),
        "n_rows": int(len(queue)),
        "source_kind_counts": queue.get("source_kind", pd.Series("", index=queue.index)).astype(str).value_counts().to_dict(),
        "leo_class_counts": queue.get("leo_class", pd.Series("", index=queue.index)).fillna("").astype(str).value_counts().to_dict(),
        "passed": not failures,
        "failures": failures,
        "min_reports": int(min_reports),
        "apertures": list(args.apertures),
        "require_wd1856": not args.skip_wd1856_check,
    }
    if run_summary:
        summary["run_summary"] = run_summary
    if args.out_json is not None:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    if failures:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
