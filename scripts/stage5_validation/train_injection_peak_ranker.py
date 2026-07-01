#!/usr/bin/env python3
"""Train/evaluate a BLS peak ranker from injected peak-truth rows."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.search.peak_ranker import PeakRankerConfig, train_peak_ranker, write_peak_ranker_outputs


DEFAULT_PEAK_TABLE = (
    REPO_ROOT
    / "reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/"
    / "peak_training/s56_20k_injection_bls_peaks_chunked.csv"
)
DEFAULT_OUT_DIR = (
    REPO_ROOT
    / "reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/"
    / "peak_ranker"
)


def _read_table(path: Path):
    import pandas as pd

    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".json", ".jsonl"}:
        return pd.read_json(path, lines=suffix == ".jsonl")
    raise ValueError(f"unsupported table format: {path}")


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peak-table", type=Path, default=DEFAULT_PEAK_TABLE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--id-column", default="injection_id")
    parser.add_argument("--validation-fraction", type=float, default=0.20)
    parser.add_argument("--test-fraction", type=float, default=0.20)
    parser.add_argument("--random-state", type=int, default=56)
    parser.add_argument("--max-iter", type=int, default=500)
    parser.add_argument("--l2", type=float, default=1.0e-3)
    parser.add_argument("--no-scored-table", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    peaks = _read_table(args.peak_table)
    cfg = PeakRankerConfig(
        id_column=args.id_column,
        validation_fraction=args.validation_fraction,
        test_fraction=args.test_fraction,
        random_state=args.random_state,
        max_iter=args.max_iter,
        l2=args.l2,
    )
    result = train_peak_ranker(peaks, cfg)
    result.summary.update({"peak_table": str(args.peak_table), "out_dir": str(args.out_dir)})
    write_peak_ranker_outputs(result, args.out_dir, write_scored=not args.no_scored_table)
    all_summary = result.summary["splits"]["all"]
    model_top1 = all_summary["model_recall_at_k"][1]
    bls_top1 = all_summary["bls_rank_recall_at_k"][1]
    model_top5 = all_summary["model_recall_at_k"][5]
    print("[peak-ranker] complete")
    print(f"  rows: {len(result.scored_peaks):,}")
    print(f"  injections: {result.summary['n_injections']:,}")
    print(f"  rankable injections: {result.summary['n_rankable_injections']:,}")
    print(
        "  all recall@1 model/BLS: "
        f"{model_top1['n']}/{model_top1['denom']} vs {bls_top1['n']}/{bls_top1['denom']}"
    )
    print(f"  all recall@5 model: {model_top5['n']}/{model_top5['denom']}")
    print(f"  out: {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
