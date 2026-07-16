#!/usr/bin/env python3
"""Build training-compatible A2v1 candidate rows for teacher-v1 inference."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.teacher_candidates import (  # noqa: E402
    enrich_candidate_metadata,
    normalize_a2v1_peak_candidates,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--adp-peaks", type=Path, required=True)
    parser.add_argument("--compact-lc", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--sector", type=int)
    parser.add_argument("--small-peaks-per-tic", type=int, default=3)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--progress-every", type=int, default=500)
    args = parser.parse_args()
    peaks = (
        pd.read_parquet(args.adp_peaks)
        if args.adp_peaks.suffix.lower() == ".parquet"
        else pd.read_csv(args.adp_peaks, low_memory=False)
    )
    candidates = normalize_a2v1_peak_candidates(
        peaks,
        small_peaks_per_tic=args.small_peaks_per_tic,
        sector=args.sector,
    )
    candidates, summary = enrich_candidate_metadata(
        candidates,
        compact_lc_path=args.compact_lc,
        workers=args.workers,
        progress_every=args.progress_every,
    )
    if not summary["passed"]:
        raise RuntimeError(f"candidate metadata failed: {summary['metadata_status_counts']}")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    output = args.out_dir / "teacher_v1_scoring_candidates.parquet"
    try:
        candidates.to_parquet(output, compression="zstd", index=False)
    except (ImportError, ModuleNotFoundError, ValueError):
        output = output.with_suffix(".csv")
        candidates.to_csv(output, index=False)
    summary["outputs"] = {
        "candidate_table": str(output),
        "summary": str(args.out_dir / "summary.json"),
    }
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
