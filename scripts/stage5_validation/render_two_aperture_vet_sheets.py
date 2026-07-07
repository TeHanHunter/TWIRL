#!/usr/bin/env python3
"""Render TWIRL-native two-aperture vet sheets for a review queue."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
for path in (SRC_ROOT, Path(__file__).resolve().parent):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from twirl.io.hlsp import read_hlsp  # noqa: E402
from twirl.io.compact_export import read_compact_lc_export, read_injected_lc_group  # noqa: E402
from twirl.search.bls import BLSConfig  # noqa: E402
from twirl.lightcurves.detrend_presets import (  # noqa: E402
    TWIRL_FS_V2_ADP015Q_BRANCH,
    compare_column_names,
)
from twirl.vetting.lightcurve_label_app import find_hlsp_path  # noqa: E402
from twirl.vetting.two_aperture import DEFAULT_TWO_APERTURE_APERTURES, render_two_aperture_sheet  # noqa: E402


DEFAULT_QUEUE = (
    REPO_ROOT / "reports/stage5_validation/s56_mixed_teacher_queue_pdo/review_queue_1k.csv"
)
DEFAULT_HLSP_ROOT = REPO_ROOT / "data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_adp015q_compare"
DEFAULT_OUT_DIR = REPO_ROOT / "reports/stage5_validation/s56_mixed_teacher_queue_pdo/twirl_vet_sheets"


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _clean(value: Any) -> Any:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass
    if isinstance(value, np.generic):
        return value.item()
    return value


def _row_id(row: pd.Series, fallback: int) -> str:
    value = _clean(row.get("review_id"))
    if value:
        return str(value).replace("/", "_").replace(":", "_")
    tic = _clean(row.get("tic")) or "unknown"
    return f"row{fallback:05d}_tic{tic}"


def _render_one(
    payload: tuple[
        dict[str, Any],
        int,
        str,
        str,
        str,
        str,
        str,
        tuple[str, ...],
        int,
        int,
        str,
    ]
) -> dict[str, Any]:
    (
        row,
        idx,
        hlsp_root_s,
        lc_export_h5_s,
        injection_h5_s,
        out_dir_s,
        branch_name,
        apertures,
        n_periods,
        n_peaks,
        overwrite_s,
    ) = payload
    hlsp_root = Path(hlsp_root_s)
    lc_export_h5 = Path(lc_export_h5_s) if lc_export_h5_s else None
    injection_h5 = Path(injection_h5_s) if injection_h5_s else None
    out_dir = Path(out_dir_s)
    overwrite = overwrite_s == "1"
    tic = int(float(row["tic"]))
    sector_value = row.get("sector")
    sector = int(float(sector_value)) if pd.notna(sector_value) else None
    base = _row_id(pd.Series(row), idx)
    out_path = out_dir / f"{base}_twirl_twoap_{branch_name}.png"
    record: dict[str, Any] = {
        "row_index": idx,
        "review_id": row.get("review_id", ""),
        "tic": tic,
        "sector": sector if sector is not None else "",
        "twirl_vet_sheet_name": out_path.name,
        "twirl_vet_sheet_pdf_name": out_path.with_suffix(".pdf").name,
        "twirl_vet_status": "",
    }
    if out_path.exists() and not overwrite:
        record["twirl_vet_status"] = "reused"
        return record
    if injection_h5 is not None:
        group_path = str(row.get("h5_group", "")).strip()
        record["twirl_vet_input"] = f"{injection_h5}:{group_path}" if group_path else str(injection_h5)
        if not group_path:
            record["twirl_vet_status"] = "missing_h5_group"
            return record
        lc = read_injected_lc_group(injection_h5, group_path=group_path, columns=apertures)
        if lc is None:
            record["twirl_vet_status"] = "missing_injection_lc"
            return record
    elif lc_export_h5 is not None:
        lc = read_compact_lc_export(lc_export_h5, tic=tic, columns=apertures)
        record["twirl_vet_input"] = str(lc_export_h5)
        if lc is None:
            record["twirl_vet_status"] = "missing_lc_export"
            return record
    else:
        path = find_hlsp_path(hlsp_root, tic, sector)
        record["twirl_vet_input"] = str(path) if path is not None else ""
        if path is None:
            record["twirl_vet_status"] = "missing_hlsp"
            return record
        lc = read_hlsp(path, columns=apertures)
        if lc is None:
            record["twirl_vet_status"] = "read_fail"
            return record
    missing = [ap for ap in apertures if ap not in lc.flux]
    if missing:
        record["twirl_vet_status"] = "missing_apertures:" + ",".join(missing)
        return record
    cfg = BLSConfig(
        apertures=apertures,
        n_periods=int(n_periods),
        n_peaks=int(n_peaks),
    )
    try:
        _, metrics = render_two_aperture_sheet(
            lc,
            out_path,
            branch_name=branch_name,
            cfg=cfg,
            apertures=apertures,
            anchor_aperture=apertures[0] if apertures else None,
            row_metadata=row,
        )
    except Exception as exc:
        record["twirl_vet_status"] = f"error:{type(exc).__name__}: {exc}"
        return record
    record.update(metrics)
    record["twirl_vet_status"] = "ok"
    return record


def render_queue(
    *,
    queue_csv: Path,
    hlsp_root: Path,
    lc_export_h5: Path | None,
    injection_h5: Path | None,
    out_dir: Path,
    metrics_csv: Path,
    summary_json: Path,
    branch_name: str,
    apertures: tuple[str, ...],
    limit: int,
    workers: int,
    n_periods: int,
    n_peaks: int,
    overwrite: bool,
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    queue = pd.read_csv(queue_csv)
    if limit > 0:
        queue = queue.head(limit).copy()
    payloads = [
        (
            row,
            idx,
            str(hlsp_root),
            str(lc_export_h5) if lc_export_h5 else "",
            str(injection_h5) if injection_h5 else "",
            str(out_dir),
            branch_name,
            apertures,
            int(n_periods),
            int(n_peaks),
            "1" if overwrite else "0",
        )
        for idx, row in enumerate(queue.to_dict("records"))
    ]
    rows: list[dict[str, Any]] = []
    workers = max(1, int(workers))
    if workers <= 1:
        for i, payload in enumerate(payloads, start=1):
            rows.append(_render_one(payload))
            if i % 25 == 0:
                print(f"[two-aperture] rendered {i:,}/{len(payloads):,}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for i, record in enumerate(ex.map(_render_one, payloads, chunksize=2), start=1):
                rows.append(record)
                if i % 25 == 0:
                    print(f"[two-aperture] rendered {i:,}/{len(payloads):,}", flush=True)
    metrics = pd.DataFrame(rows)
    metrics_csv.parent.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(metrics_csv, index=False)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "queue_csv": str(queue_csv),
        "hlsp_root": str(hlsp_root),
        "lc_export_h5": str(lc_export_h5) if lc_export_h5 else "",
        "injection_h5": str(injection_h5) if injection_h5 else "",
        "out_dir": str(out_dir),
        "metrics_csv": str(metrics_csv),
        "branch_name": branch_name,
        "apertures": list(apertures),
        "n_rows": int(len(metrics)),
        "status_counts": {
            str(k): int(v)
            for k, v in metrics.get("twirl_vet_status", pd.Series(dtype=str)).fillna("").astype(str).value_counts().sort_index().items()
        },
        "n_periods": int(n_periods),
        "n_peaks": int(n_peaks),
    }
    summary_json.parent.mkdir(parents=True, exist_ok=True)
    summary_json.write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n"
    )
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--queue-csv", type=Path, default=DEFAULT_QUEUE)
    ap.add_argument("--hlsp-root", type=Path, default=DEFAULT_HLSP_ROOT)
    ap.add_argument(
        "--lc-export-h5",
        type=Path,
        default=None,
        help=(
            "Optional compact LC HDF5 export. When supplied, target light curves "
            "are read from /targets/{tic:016d} instead of FITS."
        ),
    )
    ap.add_argument(
        "--injection-h5",
        type=Path,
        default=None,
        help=(
            "Optional injected-light-curve HDF5. When supplied, rows are read "
            "from their h5_group values instead of real FITS/compact exports."
        ),
    )
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    ap.add_argument("--metrics-csv", type=Path, default=None)
    ap.add_argument("--summary-json", type=Path, default=None)
    ap.add_argument("--branch-name", default=TWIRL_FS_V2_ADP015Q_BRANCH)
    ap.add_argument(
        "--apertures",
        nargs="+",
        default=list(DEFAULT_TWO_APERTURE_APERTURES),
        help=(
            "Two flux columns to render/search. Default uses the new "
            "twirl-fs-v2-adp015q compare columns."
        ),
    )
    ap.add_argument(
        "--aperture-tag",
        default=None,
        help=(
            "Shortcut for compare columns by tag, e.g. ADP015 -> "
            "DET_FLUX_ADP015_SML DET_FLUX_ADP015."
        ),
    )
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--n-periods", type=int, default=20_000)
    ap.add_argument("--n-peaks", type=int, default=10)
    ap.add_argument("--overwrite", action="store_true")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    metrics_csv = args.metrics_csv or (args.out_dir.parent / "twirl_vet_metrics.csv")
    summary_json = args.summary_json or (args.out_dir.parent / "twirl_two_aperture_vet_summary.json")
    apertures = tuple(args.apertures)
    if args.aperture_tag:
        col = compare_column_names(args.aperture_tag)
        apertures = (col["small"], col["primary"])
    summary = render_queue(
        queue_csv=args.queue_csv,
        hlsp_root=args.hlsp_root,
        lc_export_h5=args.lc_export_h5,
        injection_h5=args.injection_h5,
        out_dir=args.out_dir,
        metrics_csv=metrics_csv,
        summary_json=summary_json,
        branch_name=args.branch_name,
        apertures=apertures,
        limit=args.limit,
        workers=args.workers,
        n_periods=args.n_periods,
        n_peaks=args.n_peaks,
        overwrite=bool(args.overwrite),
    )
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
