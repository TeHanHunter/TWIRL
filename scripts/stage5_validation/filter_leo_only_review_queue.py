#!/usr/bin/env python3
"""Create a review queue containing only full LEO-Vetter report rows."""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import pandas as pd


LABEL_COLUMNS = ("label", "label_source", "labeler", "notes")


def filter_leo_only(
    *,
    source_dir: Path,
    out_dir: Path,
    copy_reports: bool,
    clear_labels: bool,
) -> dict:
    queue_path = source_dir / "review_queue.csv"
    metrics_path = source_dir / "leo_metrics.csv"
    reports_dir = source_dir / "vet_reports"
    if not queue_path.exists():
        raise FileNotFoundError(f"missing review queue: {queue_path}")
    if not metrics_path.exists():
        raise FileNotFoundError(f"missing LEO metrics: {metrics_path}")
    if not reports_dir.exists():
        raise FileNotFoundError(f"missing vet report directory: {reports_dir}")

    queue = pd.read_csv(queue_path)
    metrics = pd.read_csv(metrics_path)
    for col in ("error", "plot_error"):
        if col not in metrics:
            metrics[col] = ""
    ok_metrics = metrics[
        metrics["error"].fillna("").astype(str).eq("")
        & metrics["plot_error"].fillna("").astype(str).eq("")
    ].copy()
    ok_ids = set(ok_metrics["review_id"].astype(str))
    out_queue = queue[queue["review_id"].astype(str).isin(ok_ids)].copy().reset_index(drop=True)
    if clear_labels:
        for col in LABEL_COLUMNS:
            if col in out_queue:
                out_queue[col] = ""

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "vet_reports").mkdir(parents=True, exist_ok=True)
    out_queue.to_csv(out_dir / "review_queue.csv", index=False)
    ok_metrics[
        ok_metrics["review_id"].astype(str).isin(set(out_queue["review_id"].astype(str)))
    ].to_csv(out_dir / "leo_metrics.csv", index=False)

    report_names = out_queue["leo_report_name"].fillna("").astype(str).tolist()
    missing_reports: list[str] = []
    copied_reports = 0
    with (out_dir / "leo_report_files.txt").open("w") as handle:
        for name in report_names:
            if not name:
                missing_reports.append(name)
                continue
            src = reports_dir / name
            if not src.exists():
                missing_reports.append(name)
                continue
            handle.write(f"vet_reports/{name}\n")
            if copy_reports:
                dst = out_dir / "vet_reports" / name
                if not dst.exists() or dst.stat().st_size != src.stat().st_size:
                    shutil.copy2(src, dst)
                    copied_reports += 1

    summary = {
        "source_dir": str(source_dir),
        "source_queue": str(queue_path),
        "source_vet_reports": str(reports_dir),
        "out_dir": str(out_dir),
        "rows_in": int(len(queue)),
        "rows_out": int(len(out_queue)),
        "dropped_non_leo_plot_rows": int(len(queue) - len(out_queue)),
        "copy_reports": bool(copy_reports),
        "copied_reports": int(copied_reports),
        "missing_reports": missing_reports,
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--copy-reports", action="store_true")
    ap.add_argument("--keep-labels", action="store_true")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    summary = filter_leo_only(
        source_dir=args.source_dir,
        out_dir=args.out_dir,
        copy_reports=args.copy_reports,
        clear_labels=not args.keep_labels,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    if summary["missing_reports"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
