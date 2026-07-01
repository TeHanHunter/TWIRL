#!/usr/bin/env python3
"""Verify ranker-selected real-candidate ephemerides before LEO queue building."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SELECTED = (
    REPO_ROOT / "reports/stage5_validation/s56_ranker_selected_real_candidates/selected_ephemerides.csv"
)
KNOWN_APERTURES = {
    "DET_FLUX_SML",
    "DET_FLUX",
    "DET_FLUX_LAG",
    "DET_FLUX_ADP",
    "DET_FLUX_ADP_SML",
    "DET_FLUX_ADP_LAG",
}


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


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


def _finite(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df:
        return pd.Series(False, index=df.index)
    value = pd.to_numeric(df[col], errors="coerce")
    return np.isfinite(value)


def _positive(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df:
        return pd.Series(False, index=df.index)
    value = pd.to_numeric(df[col], errors="coerce")
    return np.isfinite(value) & value.gt(0)


def _blank_or_missing(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df:
        return pd.Series(True, index=df.index)
    return _as_text(df[col]).eq("")


def _score_column(df: pd.DataFrame, prefix: str) -> str:
    preferred = f"{prefix}_p_signal_peak"
    if preferred in df:
        return preferred
    score_cols = [col for col in df.columns if col.startswith(f"{prefix}_p_")]
    return score_cols[0] if score_cols else ""


def verify_selected_ephemerides(
    *,
    selected_ephemerides: Path,
    id_column: str = "tic",
    score_prefix: str = "ranker",
    top_n: int = 3,
    min_rows: int = 1,
    min_targets: int = 1,
) -> dict[str, Any]:
    df = _read_table(selected_ephemerides)
    failures: list[str] = []
    required = (
        id_column,
        "period_d",
        "t0_bjd",
        "duration_min",
        "ranker_selection_rank",
    )
    missing = [col for col in required if col not in df]
    if missing:
        failures.append("missing required columns: " + ", ".join(missing))
        return {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "selected_ephemerides": str(selected_ephemerides),
            "n_rows": int(len(df)),
            "n_targets": 0,
            "score_column": "",
            "passed": False,
            "failures": failures,
        }

    if len(df) < min_rows:
        failures.append(f"rows={len(df)}; expected at least {min_rows}")

    target = pd.to_numeric(df[id_column], errors="coerce")
    n_targets = int(target.dropna().nunique())
    if n_targets < min_targets:
        failures.append(f"unique {id_column}={n_targets}; expected at least {min_targets}")
    bad_target = int((~np.isfinite(target)).sum())
    if bad_target:
        failures.append(f"{bad_target} rows have non-finite {id_column}")

    for col in ("period_d", "duration_min"):
        bad = int((~_positive(df, col)).sum())
        if bad:
            failures.append(f"{bad} rows have non-positive/non-finite {col}")
    for col in ("t0_bjd",):
        bad = int((~_finite(df, col)).sum())
        if bad:
            failures.append(f"{bad} rows have non-finite {col}")

    rank = pd.to_numeric(df["ranker_selection_rank"], errors="coerce")
    bad_rank = ~(np.isfinite(rank) & rank.ge(1))
    if top_n > 0:
        bad_rank = bad_rank | rank.gt(top_n)
    if bool(bad_rank.any()):
        failures.append(f"{int(bad_rank.sum())} rows have invalid ranker_selection_rank")

    rows_per_target = df.loc[np.isfinite(target)].groupby(target.loc[np.isfinite(target)].astype(np.int64)).size()
    max_rows_per_target = int(rows_per_target.max()) if not rows_per_target.empty else 0
    if top_n > 0 and max_rows_per_target > top_n:
        failures.append(f"max rows per {id_column}={max_rows_per_target}; expected <= {top_n}")
    if top_n > 0 and not rows_per_target.empty:
        duplicate_ranks = int(
            df.assign(_target=target, _rank=rank)
            .dropna(subset=["_target", "_rank"])
            .duplicated(["_target", "_rank"])
            .sum()
        )
        if duplicate_ranks:
            failures.append(f"duplicate {id_column}/ranker_selection_rank rows: {duplicate_ranks}")

    aperture_col = "aperture" if "aperture" in df else ("rep_aperture" if "rep_aperture" in df else "")
    if aperture_col:
        aperture = _as_text(df[aperture_col])
        bad_aperture = ~aperture.isin(KNOWN_APERTURES)
        if bool(bad_aperture.any()):
            failures.append(f"{int(bad_aperture.sum())} rows have unknown {aperture_col}")

    score_col = _score_column(df, score_prefix)
    if score_col:
        score = pd.to_numeric(df[score_col], errors="coerce")
        bad_score = ~(np.isfinite(score) & score.ge(0) & score.le(1))
        if bool(bad_score.any()):
            failures.append(f"{int(bad_score.sum())} rows have invalid {score_col}")
    else:
        failures.append(f"missing {score_prefix}_p_* score column")

    for col in ("injection_id", "truth_period_d", "truth_t0_bjd", "truth_radius_rearth"):
        bad_truth = ~_blank_or_missing(df, col)
        if bool(bad_truth.any()):
            failures.append(f"{int(bad_truth.sum())} real rows have non-empty {col}")

    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "selected_ephemerides": str(selected_ephemerides),
        "id_column": id_column,
        "score_prefix": score_prefix,
        "score_column": score_col,
        "top_n": int(top_n),
        "n_rows": int(len(df)),
        "n_targets": n_targets,
        "max_rows_per_target": max_rows_per_target,
        "rank_counts": {
            str(k): int(v)
            for k, v in rank.dropna().astype(int).value_counts().sort_index().items()
        },
        "aperture_counts": (
            _as_text(df[aperture_col]).value_counts().sort_index().astype(int).to_dict()
            if aperture_col
            else {}
        ),
        "passed": not failures,
        "failures": failures,
    }
    return payload


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected-ephemerides", type=Path, default=DEFAULT_SELECTED)
    parser.add_argument("--out-json", type=Path, default=None)
    parser.add_argument("--id-column", default="tic")
    parser.add_argument("--score-prefix", default="ranker")
    parser.add_argument("--top-n", type=int, default=3)
    parser.add_argument("--min-rows", type=int, default=1)
    parser.add_argument("--min-targets", type=int, default=1)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    result = verify_selected_ephemerides(
        selected_ephemerides=args.selected_ephemerides,
        id_column=args.id_column,
        score_prefix=args.score_prefix,
        top_n=args.top_n,
        min_rows=args.min_rows,
        min_targets=args.min_targets,
    )
    text = json.dumps(result, indent=2, sort_keys=True, default=_json_default) + "\n"
    if args.out_json is not None:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(text)
    print(text, end="")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
