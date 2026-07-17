#!/usr/bin/env python3
"""Verify the real S56 BLS peak table before injected-ranker application."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_TABLE = REPO_ROOT / "data_local/stage2/bls_first_pass_v2/sector_0056/candidates.csv"
REQUIRED_COLUMNS = (
    "tic",
    "sector",
    "aperture",
    "peak_rank",
    "period_d",
    "t0_bjd",
    "duration_min",
    "sde",
    "status",
)


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


def _counts(series: pd.Series, *, limit: int = 30) -> dict[str, int]:
    return {
        str(k): int(v)
        for k, v in series.fillna("").astype(str).value_counts().sort_index().head(limit).items()
    }


def _numeric_summary(series: pd.Series) -> dict[str, float]:
    value = pd.to_numeric(series, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    if value.empty:
        return {}
    return {
        "min": float(value.min()),
        "p10": float(value.quantile(0.10)),
        "median": float(value.median()),
        "p90": float(value.quantile(0.90)),
        "max": float(value.max()),
    }


def verify_real_peak_table(
    *,
    peak_table: Path,
    min_rows: int = 1_000_000,
    min_tics: int = 18_000,
    min_peak_ranks: int = 10,
    min_apertures: int = 2,
    min_rows_per_tic_max: int = 10,
) -> dict[str, Any]:
    df = _read_table(peak_table)
    failures: list[str] = []
    missing = [col for col in REQUIRED_COLUMNS if col not in df]
    if missing:
        failures.append("missing required columns: " + ", ".join(missing))

    if len(df) < min_rows:
        failures.append(f"rows={len(df)}; expected at least {min_rows}")

    if "tic" in df:
        tic = pd.to_numeric(df["tic"], errors="coerce")
        n_tic = int(tic.dropna().nunique())
        rows_per_tic = df.loc[tic.notna()].groupby(tic.loc[tic.notna()].astype(np.int64)).size()
    else:
        n_tic = 0
        rows_per_tic = pd.Series(dtype=int)
    if n_tic < min_tics:
        failures.append(f"unique TICs={n_tic}; expected at least {min_tics}")

    max_rows_per_tic = int(rows_per_tic.max()) if not rows_per_tic.empty else 0
    if max_rows_per_tic < min_rows_per_tic_max:
        failures.append(f"max rows per TIC={max_rows_per_tic}; expected at least {min_rows_per_tic_max}")

    ok = pd.Series(True, index=df.index)
    if "status" in df:
        ok = df["status"].fillna("").astype(str).eq("ok")
    if "peak_rank" in df:
        peak_rank = pd.to_numeric(df["peak_rank"], errors="coerce")
        candidate = ok & peak_rank.gt(0)
        n_peak_ranks = int(peak_rank.loc[candidate].dropna().nunique())
    else:
        peak_rank = pd.Series(np.nan, index=df.index)
        candidate = ok
        n_peak_ranks = 0
    if n_peak_ranks < min_peak_ranks:
        failures.append(f"positive peak-rank count={n_peak_ranks}; expected at least {min_peak_ranks}")

    if "aperture" in df:
        n_apertures = int(df.loc[candidate, "aperture"].fillna("").astype(str).nunique())
    else:
        n_apertures = 0
    if n_apertures < min_apertures:
        failures.append(f"aperture count={n_apertures}; expected at least {min_apertures}")

    for col in ("tic", "period_d", "t0_bjd", "sde"):
        if col in df:
            bad = int((candidate & ~_finite(df, col)).sum())
            if bad:
                failures.append(f"{bad} candidate rows have non-finite {col}")
    if "duration_min" in df:
        bad_duration = int((candidate & ~_positive(df, "duration_min")).sum())
        if bad_duration:
            failures.append(f"{bad_duration} candidate rows have non-positive/non-finite duration_min")

    if "sector" in df:
        bad_sector = int((pd.to_numeric(df["sector"], errors="coerce") != 56).sum())
        if bad_sector:
            failures.append(f"{bad_sector} rows are not sector 56")

    targets_with_multiple_rows = int((rows_per_tic.ge(2)).sum()) if not rows_per_tic.empty else 0
    targets_with_min_peak_rows = int((rows_per_tic.ge(min_rows_per_tic_max)).sum()) if not rows_per_tic.empty else 0
    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "peak_table": str(peak_table),
        "n_rows": int(len(df)),
        "n_unique_tic": n_tic,
        "n_rankable_rows": int(candidate.sum()) if len(df) else 0,
        "n_positive_peak_ranks": n_peak_ranks,
        "n_apertures": n_apertures,
        "max_rows_per_tic": max_rows_per_tic,
        "n_tics_with_multiple_rows": targets_with_multiple_rows,
        "n_tics_with_min_peak_rows": targets_with_min_peak_rows,
        "rows_per_tic_summary": _numeric_summary(rows_per_tic),
        "status_counts": _counts(df["status"]) if "status" in df else {},
        "aperture_counts": _counts(df["aperture"]) if "aperture" in df else {},
        "peak_rank_counts": _counts(peak_rank.astype("Int64").astype(str)) if "peak_rank" in df else {},
        "passed": not failures,
        "failures": failures,
    }
    return payload


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peak-table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--out-json", type=Path, default=None)
    parser.add_argument("--min-rows", type=int, default=1_000_000)
    parser.add_argument("--min-tics", type=int, default=18_000)
    parser.add_argument("--min-peak-ranks", type=int, default=10)
    parser.add_argument("--min-apertures", type=int, default=2)
    parser.add_argument("--min-rows-per-tic-max", type=int, default=10)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    result = verify_real_peak_table(
        peak_table=args.peak_table,
        min_rows=args.min_rows,
        min_tics=args.min_tics,
        min_peak_ranks=args.min_peak_ranks,
        min_apertures=args.min_apertures,
        min_rows_per_tic_max=args.min_rows_per_tic_max,
    )
    text = json.dumps(result, indent=2, sort_keys=True, default=_json_default) + "\n"
    if args.out_json is not None:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(text)
    print(text, end="")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
