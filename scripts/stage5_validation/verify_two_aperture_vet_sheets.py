#!/usr/bin/env python3
"""Verify TWIRL-native two-aperture vet-sheet outputs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


DEFAULT_QUEUE = Path("reports/stage5_validation/s56_mixed_teacher_queue_pdo/review_queue_1k.csv")
DEFAULT_OUT_DIR = Path("reports/stage5_validation/s56_mixed_teacher_queue_pdo/twirl_vet_sheets_adp015q_orcd")
DEFAULT_METRICS = Path("reports/stage5_validation/s56_mixed_teacher_queue_pdo/twirl_vet_metrics_adp015q_orcd.csv")
DEFAULT_SUMMARY = Path("reports/stage5_validation/s56_mixed_teacher_queue_pdo/twirl_two_aperture_vet_summary_adp015q_orcd.json")
DEFAULT_VERIFICATION = Path("reports/stage5_validation/s56_mixed_teacher_queue_pdo/twirl_two_aperture_vet_verification_adp015q_orcd.json")


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _status_counts(series: pd.Series) -> dict[str, int]:
    return {
        str(k): int(v)
        for k, v in series.fillna("").astype(str).value_counts().sort_index().items()
    }


def verify_twoap_outputs(
    *,
    queue_csv: Path,
    metrics_csv: Path,
    out_dir: Path,
    summary_json: Path | None,
    verification_json: Path | None,
    expected_rows: int | None,
    allow_reused: bool,
    require_pdf: bool,
) -> dict[str, Any]:
    failures: list[str] = []
    warnings: list[str] = []

    if not metrics_csv.exists():
        failures.append(f"missing metrics CSV: {metrics_csv}")
        metrics = pd.DataFrame()
    else:
        metrics = pd.read_csv(metrics_csv)

    queue = pd.read_csv(queue_csv) if queue_csv.exists() else pd.DataFrame()
    if not queue_csv.exists():
        failures.append(f"missing queue CSV: {queue_csv}")

    if expected_rows is None:
        expected_rows = int(len(queue)) if len(queue) else None
    if expected_rows is not None and len(metrics) != int(expected_rows):
        failures.append(f"metrics rows = {len(metrics)}; expected {int(expected_rows)}")
    if len(queue) and expected_rows is not None and int(expected_rows) > len(queue):
        failures.append(f"expected_rows={int(expected_rows)} exceeds queue rows={len(queue)}")

    required_columns = {
        "row_index",
        "review_id",
        "tic",
        "twirl_vet_status",
        "twirl_vet_sheet_name",
        "twirl_vet_sheet_pdf_name",
        "vet_branch",
        "anchor_aperture",
        "anchor_period_d",
        "anchor_duration_min",
    }
    missing_columns = sorted(required_columns - set(metrics.columns))
    if missing_columns:
        failures.append("missing metrics columns: " + ", ".join(missing_columns))

    status_counts = (
        _status_counts(metrics["twirl_vet_status"]) if "twirl_vet_status" in metrics else {}
    )
    allowed_status = {"ok", "reused"} if allow_reused else {"ok"}
    if "twirl_vet_status" in metrics:
        bad_status = ~metrics["twirl_vet_status"].fillna("").astype(str).isin(allowed_status)
        if bool(bad_status.any()):
            first = metrics.loc[bad_status, "twirl_vet_status"].fillna("").astype(str).iloc[0]
            failures.append(f"{int(bad_status.sum())} rows have bad twirl_vet_status; first={first!r}")

    if "row_index" in metrics and expected_rows is not None:
        row_index = pd.to_numeric(metrics["row_index"], errors="coerce")
        expected = set(range(int(expected_rows)))
        observed = set(int(v) for v in row_index[np.isfinite(row_index)].astype(int))
        missing = sorted(expected - observed)
        extra = sorted(observed - expected)
        if missing:
            failures.append(f"missing row_index values; first={missing[0]} n={len(missing)}")
        if extra:
            failures.append(f"unexpected row_index values; first={extra[0]} n={len(extra)}")

    if len(queue) and len(metrics) and {"row_index", "tic"}.issubset(metrics.columns):
        check = metrics.copy()
        check["_row_index_int"] = pd.to_numeric(check["row_index"], errors="coerce").astype("Int64")
        valid = check["_row_index_int"].notna() & check["_row_index_int"].between(0, len(queue) - 1)
        if bool(valid.any()):
            q_tic = pd.to_numeric(queue["tic"], errors="coerce").reset_index(drop=True)
            m_tic = pd.to_numeric(check.loc[valid, "tic"], errors="coerce")
            q_take = q_tic.iloc[check.loc[valid, "_row_index_int"].astype(int).to_numpy()].reset_index(drop=True)
            mismatch = m_tic.reset_index(drop=True).astype("Int64") != q_take.astype("Int64")
            if bool(mismatch.any()):
                failures.append(f"{int(mismatch.sum())} metrics rows have TIC mismatch against queue row_index")

    missing_png: list[str] = []
    missing_pdf: list[str] = []
    if "twirl_vet_sheet_name" in metrics:
        for name in metrics["twirl_vet_sheet_name"].fillna("").astype(str):
            if not name or not (out_dir / name).exists():
                missing_png.append(name)
    if require_pdf and "twirl_vet_sheet_pdf_name" in metrics:
        for name in metrics["twirl_vet_sheet_pdf_name"].fillna("").astype(str):
            if not name or not (out_dir / name).exists():
                missing_pdf.append(name)
    if missing_png:
        failures.append(f"missing PNG sheets: {len(missing_png)}; first={missing_png[0]!r}")
    if missing_pdf:
        failures.append(f"missing PDF sheets: {len(missing_pdf)}; first={missing_pdf[0]!r}")

    summary: dict[str, Any] | None = None
    if summary_json is not None:
        if summary_json.exists():
            summary = json.loads(summary_json.read_text())
            if expected_rows is not None and int(summary.get("n_rows", -1)) != int(expected_rows):
                failures.append(f"summary n_rows={summary.get('n_rows')}; expected {int(expected_rows)}")
            if str(summary.get("metrics_csv", "")) and Path(str(summary["metrics_csv"])).name != metrics_csv.name:
                warnings.append("summary metrics_csv basename differs from requested metrics CSV")
            summary_counts = summary.get("status_counts", {})
            if status_counts and summary_counts and {
                str(k): int(v) for k, v in summary_counts.items()
            } != status_counts:
                failures.append("summary status_counts do not match metrics CSV")
        else:
            failures.append(f"missing summary JSON: {summary_json}")

    payload = {
        "queue_csv": str(queue_csv),
        "metrics_csv": str(metrics_csv),
        "out_dir": str(out_dir),
        "summary_json": str(summary_json) if summary_json is not None else "",
        "expected_rows": int(expected_rows) if expected_rows is not None else None,
        "n_queue_rows": int(len(queue)),
        "n_metrics_rows": int(len(metrics)),
        "status_counts": status_counts,
        "allow_reused": bool(allow_reused),
        "require_pdf": bool(require_pdf),
        "n_missing_png": int(len(missing_png)),
        "n_missing_pdf": int(len(missing_pdf)),
        "warnings": warnings,
        "failures": failures,
        "passed": not failures,
    }
    if summary is not None:
        payload["summary_branch_name"] = summary.get("branch_name", "")
        payload["summary_apertures"] = summary.get("apertures", [])
        payload["summary_n_periods"] = summary.get("n_periods", None)

    if verification_json is not None:
        verification_json.parent.mkdir(parents=True, exist_ok=True)
        verification_json.write_text(json.dumps(payload, indent=2, sort_keys=True, default=_json_default) + "\n")
    return payload


def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--queue-csv", type=Path, default=DEFAULT_QUEUE)
    ap.add_argument("--metrics-csv", type=Path, default=DEFAULT_METRICS)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    ap.add_argument("--summary-json", type=Path, default=DEFAULT_SUMMARY)
    ap.add_argument("--verification-json", type=Path, default=DEFAULT_VERIFICATION)
    ap.add_argument("--expected-rows", type=int, default=None)
    ap.add_argument("--allow-reused", action="store_true")
    ap.add_argument("--no-pdf", action="store_true", help="Do not require PDF companions.")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    payload = verify_twoap_outputs(
        queue_csv=args.queue_csv,
        metrics_csv=args.metrics_csv,
        out_dir=args.out_dir,
        summary_json=args.summary_json,
        verification_json=args.verification_json,
        expected_rows=args.expected_rows,
        allow_reused=bool(args.allow_reused),
        require_pdf=not bool(args.no_pdf),
    )
    print(json.dumps(payload, indent=2, sort_keys=True, default=_json_default))
    return 0 if payload["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
