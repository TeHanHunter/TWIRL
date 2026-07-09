#!/usr/bin/env python3
"""Audit raw point counts for S56 recovery50 CNN tensor sizing."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.vetting.recovery50_cnn import (  # noqa: E402
    FOLD_VIEW_NAMES,
    TensorConfig,
    _phase_x_duration,
    _read_tensor_lc,
    _time_to_bjd,
    fold_view_ephemerides,
)
from twirl.vetting.adp_only import assert_adp_only_training_frame, validate_adp_only_apertures  # noqa: E402
from twirl.vetting.recovery50_teacher import json_default, read_table, _safe_float  # noqa: E402


DEFAULT_ROOT = REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue_2k"
DEFAULT_TRAINING_TABLE = DEFAULT_ROOT / "human_training_table_adp_only/human_vetting_training_table.csv"
DEFAULT_OUT_DIR = DEFAULT_ROOT / "cnn_raw_count_audit_adp_only"
DEFAULT_COMPACT_LC = (
    REPO_ROOT
    / "data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp_lc_export_pdo.h5"
)
DEFAULT_HLSP_ROOT = REPO_ROOT / "data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare"
DEFAULT_INJECTION_H5 = None


def _count_window(
    *,
    time_bjd: np.ndarray,
    quality: np.ndarray,
    period_d: float,
    t0_bjd: float,
    duration_min: float,
    window_durations: float,
) -> int:
    x = _phase_x_duration(time_bjd, period_d=period_d, t0_bjd=t0_bjd, duration_min=duration_min)
    good = (quality == 0) & np.isfinite(time_bjd) & np.isfinite(x) & (np.abs(x) <= window_durations)
    return int(good.sum())


def _quantiles(values: pd.Series) -> dict[str, float]:
    finite = pd.to_numeric(values, errors="coerce")
    finite = finite[np.isfinite(finite)]
    if finite.empty:
        return {}
    qs = [0.50, 0.90, 0.95, 0.99, 1.0]
    return {f"q{int(q * 100):02d}": float(finite.quantile(q)) for q in qs}


def audit_raw_counts(
    *,
    training_table: Path,
    out_dir: Path,
    compact_lc_h5: Path | None,
    hlsp_root: Path | None,
    injection_h5_override: Path | None,
    config: TensorConfig,
    max_rows: int | None,
    progress_every: int,
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = read_table(training_table).copy()
    assert_adp_only_training_frame(rows)
    validate_adp_only_apertures(config.apertures)
    if "row_id" not in rows:
        rows.insert(0, "row_id", np.arange(len(rows), dtype=int))
    if max_rows is not None:
        rows = rows.head(max_rows).copy()

    records: list[dict[str, Any]] = []
    status_counts: dict[str, int] = {}
    for processed, (_, row) in enumerate(rows.iterrows(), start=1):
        rec: dict[str, Any] = {
            "row_id": int(row.get("row_id", -1)),
            "review_id": str(row.get("review_id", "")),
            "tic": int(row.get("tic", -1)),
            "source_kind": str(row.get("source_kind", "")),
            "human_label": str(row.get("human_label", "")),
            "main_teacher_target": str(row.get("main_teacher_target", "")),
            "main_teacher_include": bool(row.get("main_teacher_include", False)),
            "harmonic_suspect": bool(row.get("harmonic_suspect", False)),
            "refold_factor": _safe_float(row.get("refold_factor")),
            "model_period_d": _safe_float(row.get("model_period_d")),
            "display_period_d": _safe_float(row.get("display_period_d")),
        }
        lc = _read_tensor_lc(
            row,
            compact_lc_h5=compact_lc_h5,
            hlsp_root=hlsp_root,
            injection_h5_override=injection_h5_override,
            apertures=config.apertures,
        )
        if lc is None:
            rec["status"] = "missing_light_curve"
            records.append(rec)
            status_counts[rec["status"]] = status_counts.get(rec["status"], 0) + 1
            continue
        time_bjd = _time_to_bjd(np.asarray(lc.time, dtype=np.float64))
        quality = np.asarray(lc.quality, dtype=np.int32)
        if not len(time_bjd):
            rec["status"] = "empty_light_curve"
            records.append(rec)
            status_counts[rec["status"]] = status_counts.get(rec["status"], 0) + 1
            continue
        rec["n_good_cadences"] = int(((quality == 0) & np.isfinite(time_bjd)).sum())
        for view_name, period_d, t0_bjd, duration_min in fold_view_ephemerides(row):
            rec[f"raw_folded_count_{view_name}"] = _count_window(
                time_bjd=time_bjd,
                quality=quality,
                period_d=period_d,
                t0_bjd=t0_bjd,
                duration_min=duration_min,
                window_durations=config.folded_window_durations,
            )
            rec[f"raw_context_count_{view_name}"] = _count_window(
                time_bjd=time_bjd,
                quality=quality,
                period_d=period_d,
                t0_bjd=t0_bjd,
                duration_min=duration_min,
                window_durations=config.context_window_durations,
            )
        rec["raw_folded_max_count"] = max(int(rec.get(f"raw_folded_count_{name}", 0)) for name in FOLD_VIEW_NAMES)
        rec["raw_context_max_count"] = max(int(rec.get(f"raw_context_count_{name}", 0)) for name in FOLD_VIEW_NAMES)
        rec["status"] = "ok"
        records.append(rec)
        status_counts[rec["status"]] = status_counts.get(rec["status"], 0) + 1
        if progress_every > 0 and (processed % progress_every == 0 or processed == len(rows)):
            print(f"[raw-count-audit] processed={processed:,}/{len(rows):,}", flush=True)

    table = pd.DataFrame(records)
    row_path = out_dir / "raw_count_rows.csv"
    table.to_csv(row_path, index=False)
    ok = table[table["status"].eq("ok")].copy()
    cap_candidates = [2048, 4096, 8192, 12288, 16384, 24576, 32768]
    summary = {
        "training_table": str(training_table),
        "out_dir": str(out_dir),
        "n_input_rows": int(len(rows)),
        "n_ok": int(len(ok)),
        "status_counts": {str(k): int(v) for k, v in sorted(status_counts.items())},
        "fold_view_names": list(FOLD_VIEW_NAMES),
        "folded_window_durations": float(config.folded_window_durations),
        "context_window_durations": float(config.context_window_durations),
        "n_good_cadences_quantiles": _quantiles(ok.get("n_good_cadences", pd.Series(dtype=float))),
        "raw_folded_max_count_quantiles": _quantiles(ok.get("raw_folded_max_count", pd.Series(dtype=float))),
        "raw_context_max_count_quantiles": _quantiles(ok.get("raw_context_max_count", pd.Series(dtype=float))),
        "truncated_rows_by_cap": {
            str(cap): {
                "raw_folded": int((ok.get("raw_folded_max_count", pd.Series(dtype=float)) > cap).sum()),
                "raw_context": int((ok.get("raw_context_max_count", pd.Series(dtype=float)) > cap).sum()),
            }
            for cap in cap_candidates
        },
        "outputs": {
            "rows": str(row_path),
            "summary": str(out_dir / "summary.json"),
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=json_default) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True, default=json_default))
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training-table", type=Path, default=DEFAULT_TRAINING_TABLE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--compact-lc-h5", type=Path, default=DEFAULT_COMPACT_LC)
    parser.add_argument("--hlsp-root", type=Path, default=DEFAULT_HLSP_ROOT)
    parser.add_argument("--injection-h5-override", type=Path, default=DEFAULT_INJECTION_H5)
    parser.add_argument("--apertures", default=",".join(TensorConfig().apertures))
    parser.add_argument("--folded-window-durations", type=float, default=4.0)
    parser.add_argument("--context-window-durations", type=float, default=12.0)
    parser.add_argument("--max-rows", type=int, default=None)
    parser.add_argument("--progress-every", type=int, default=100)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    apertures = tuple(part.strip() for part in args.apertures.split(",") if part.strip())
    cfg = TensorConfig(
        apertures=apertures,
        folded_window_durations=args.folded_window_durations,
        context_window_durations=args.context_window_durations,
    )
    audit_raw_counts(
        training_table=args.training_table,
        out_dir=args.out_dir,
        compact_lc_h5=args.compact_lc_h5,
        hlsp_root=args.hlsp_root,
        injection_h5_override=args.injection_h5_override,
        config=cfg,
        max_rows=args.max_rows,
        progress_every=args.progress_every,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
