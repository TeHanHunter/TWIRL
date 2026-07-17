#!/usr/bin/env python3
"""Apply a trained injected-truth peak ranker to real or injected BLS peaks."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.search.peak_ranker import load_model, score_peak_table, select_ranked_ephemerides


DEFAULT_MODEL = (
    REPO_ROOT
    / "reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/"
    / "peak_ranker/peak_ranker_model.npz"
)
DEFAULT_PEAK_TABLE = (
    REPO_ROOT / "data_local/stage2/bls_first_pass_v2/sector_0056/candidates.parquet"
)
DEFAULT_OUT_DIR = (
    REPO_ROOT
    / "reports/stage5_validation/s56_ranker_selected_real_candidates"
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


def _write_table(df, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".parquet":
        df.to_parquet(path, index=False)
    else:
        df.to_csv(path, index=False)


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--peak-table", type=Path, default=DEFAULT_PEAK_TABLE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--id-column", default="tic")
    parser.add_argument("--score-prefix", default="ranker")
    parser.add_argument("--top-n", type=int, default=3)
    parser.add_argument("--candidate-only", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--write-scored", action=argparse.BooleanOptionalAction, default=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    model = load_model(args.model)
    peaks = _read_table(args.peak_table)
    scored = score_peak_table(
        peaks,
        model,
        score_prefix=args.score_prefix,
        candidate_only=bool(args.candidate_only),
        group_column=args.id_column,
    )
    selected = select_ranked_ephemerides(
        scored,
        id_column=args.id_column,
        score_prefix=args.score_prefix,
        top_n=args.top_n,
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.write_scored:
        _write_table(scored, args.out_dir / "scored_peaks.parquet")
    _write_table(selected, args.out_dir / "selected_ephemerides.csv")
    summary = {
        "model": str(args.model),
        "peak_table": str(args.peak_table),
        "out_dir": str(args.out_dir),
        "id_column": args.id_column,
        "top_n": int(args.top_n),
        "n_input_rows": int(len(peaks)),
        "n_scored_rows": int(len(scored)),
        "n_selected_rows": int(len(selected)),
        "n_selected_targets": int(selected[args.id_column].nunique()) if args.id_column in selected else 0,
        "score_columns": [col for col in scored.columns if col.startswith(f"{args.score_prefix}_p_")],
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print("[apply-peak-ranker] complete")
    print(f"  scored rows: {len(scored):,}")
    print(f"  selected rows: {len(selected):,}")
    print(f"  out: {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
