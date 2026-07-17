#!/usr/bin/env python3
"""Merge chunked injected-light-curve recovery sweeps.

Chunked runs keep the expensive BLS pass restartable. Each chunk is produced by
``run_injection_recovery_mode_sweep.py`` and contains one sweep subdirectory
with ``injection_bls_recoveries.csv``. This script concatenates those recovery
rows, rebuilds the standard review queue, and writes the same recovery-mode
summary files as a non-chunked sweep.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
SRC_ROOT = REPO_ROOT / "src"
for path in (SCRIPT_DIR, SRC_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from build_s56_pretriage_review_queue import _finalize_queue  # noqa: E402
from run_injection_recovery_mode_sweep import DEFAULT_DURATIONS_MIN, SweepConfig  # noqa: E402
from summarize_injection_recovery_modes import summarize  # noqa: E402


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def _chunk_recovery_paths(chunk_root: Path, sweep_name: str) -> list[Path]:
    return sorted(chunk_root.glob(f"chunk_*/{sweep_name}/injection_bls_recoveries.csv"))


def merge_chunks(
    *,
    chunk_root: Path,
    out_dir: Path,
    sweep_name: str,
    apertures: tuple[str, ...],
    n_periods: int,
    survival_csv: Path | None,
    summary_aperture: str,
    expect_n: int,
) -> dict[str, Any]:
    paths = _chunk_recovery_paths(chunk_root, sweep_name)
    if not paths:
        raise FileNotFoundError(f"no chunk recovery CSVs found under {chunk_root} for {sweep_name}")

    frames: list[pd.DataFrame] = []
    source_rows: list[dict[str, Any]] = []
    for path in paths:
        frame = pd.read_csv(path)
        if "injection_id" not in frame:
            raise ValueError(f"missing injection_id in {path}")
        frames.append(frame)
        source_rows.append(
            {
                "chunk_csv": str(path),
                "n_rows": int(len(frame)),
                "first_injection_id": str(frame["injection_id"].iloc[0]) if len(frame) else "",
                "last_injection_id": str(frame["injection_id"].iloc[-1]) if len(frame) else "",
            }
        )

    merged = pd.concat(frames, ignore_index=True, sort=False)
    duplicated = merged["injection_id"].astype(str).duplicated(keep=False)
    if duplicated.any():
        dup_ids = sorted(merged.loc[duplicated, "injection_id"].astype(str).unique())
        raise ValueError(f"duplicate injection_id rows in chunks: {dup_ids[:10]}")
    merged = merged.sort_values("injection_id", kind="stable").reset_index(drop=True)
    if expect_n > 0 and len(merged) != expect_n:
        raise ValueError(f"expected {expect_n} merged rows, found {len(merged)}")

    sweep_dir = out_dir / sweep_name
    sweep_dir.mkdir(parents=True, exist_ok=True)
    merged_csv = sweep_dir / "injection_bls_recoveries.csv"
    merged.to_csv(merged_csv, index=False)

    queue = _finalize_queue(pd.DataFrame(), merged)
    pre_leo_csv = sweep_dir / "review_queue_pre_leo.csv"
    queue_csv = sweep_dir / "review_queue.csv"
    queue.to_csv(pre_leo_csv, index=False)
    queue.to_csv(queue_csv, index=False)

    summary = summarize(queue_csv, survival_csv, summary_aperture, sweep_dir / "recovery_mode_summary")
    sweep = SweepConfig(name=sweep_name, apertures=apertures, n_periods=n_periods)
    overview = pd.DataFrame(
        [
            {
                **asdict(sweep),
                "apertures": "+".join(apertures),
                "durations_min": "+".join(f"{value:g}" for value in DEFAULT_DURATIONS_MIN),
                "injection_id_filter_n": 0,
                "n_chunks": len(paths),
                "n_rows": summary["n_rows"],
                "strict_recovered": summary["strict_recovery_status_counts"].get("bls_recovered", 0),
                "top1_recovered": summary["recovery_mode_counts"].get("bls_top1_recovered", 0),
                "exact_topn": summary["recovery_mode_counts"].get("bls_topn_recovered", 0),
                "harmonic_topn": summary["recovery_mode_counts"].get("bls_topn_harmonic_match", 0),
                "unmatched": summary["recovery_mode_counts"].get("bls_peak_mismatch", 0),
                "queue_csv": str(queue_csv),
                "summary_json": str(sweep_dir / "recovery_mode_summary" / "summary.json"),
            }
        ]
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    overview.to_csv(out_dir / "sweep_overview.csv", index=False)
    pd.DataFrame(source_rows).to_csv(out_dir / "chunk_sources.csv", index=False)

    run_summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "chunk_root": str(chunk_root),
        "out_dir": str(out_dir),
        "sweep_name": sweep_name,
        "n_chunks": len(paths),
        "n_rows": int(len(queue)),
        "expect_n": int(expect_n),
        "merged_csv": str(merged_csv),
        "queue_csv": str(queue_csv),
        "summary": summary,
    }
    (out_dir / "merge_summary.json").write_text(
        json.dumps(run_summary, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )
    return run_summary


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--sweep-name", default="small_pair_200k")
    parser.add_argument("--apertures", nargs="+", default=["DET_FLUX_ADP_SML", "DET_FLUX_SML"])
    parser.add_argument("--n-periods", type=int, default=200_000)
    parser.add_argument("--survival-csv", type=Path, default=None)
    parser.add_argument("--summary-aperture", default="DET_FLUX_ADP_SML")
    parser.add_argument("--expect-n", type=int, default=0)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    summary = merge_chunks(
        chunk_root=args.chunk_root,
        out_dir=args.out_dir,
        sweep_name=args.sweep_name,
        apertures=tuple(args.apertures),
        n_periods=args.n_periods,
        survival_csv=args.survival_csv,
        summary_aperture=args.summary_aperture,
        expect_n=args.expect_n,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
