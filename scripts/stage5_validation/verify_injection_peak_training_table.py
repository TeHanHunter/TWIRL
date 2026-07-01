#!/usr/bin/env python3
"""Verify injected BLS peak table readiness before ranker training."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_TABLE = (
    REPO_ROOT
    / "reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/"
    / "peak_training/s56_20k_injection_bls_peaks.csv"
)
REQUIRED_COLUMNS = (
    "injection_id",
    "tic",
    "aperture",
    "status",
    "is_candidate_peak",
    "peak_rank",
    "is_injected_signal_peak",
    "match_kind",
    "truth_period_d",
    "truth_t0_bjd",
    "truth_duration_min",
    "period_d",
    "t0_bjd",
    "duration_min",
    "sde",
)
CADENCE_DIAGNOSTIC_COLUMNS = (
    "n_cad_total",
    "n_cad_quality",
    "n_cad_kept",
    "n_cad_edge_trimmed",
    "n_cad_sigma_clipped",
    "dropout_frac",
    "quality_dropout_frac",
)
WINDOW_OVERLAP_COLUMNS = (
    "transit_window_match",
    "transit_window_phase_period_d",
    "transit_window_delta_min",
    "transit_window_overlap_min",
    "transit_window_overlap_fraction",
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


def _as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    text = series.fillna("").astype(str).str.strip().str.lower()
    return text.isin({"true", "1", "yes", "y"})


def _finite(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df:
        return pd.Series(False, index=df.index)
    value = pd.to_numeric(df[col], errors="coerce")
    return np.isfinite(value)


def _nonnegative(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df:
        return pd.Series(False, index=df.index)
    value = pd.to_numeric(df[col], errors="coerce")
    return np.isfinite(value) & value.ge(0)


def _counts(series: pd.Series, *, limit: int = 40) -> dict[str, int]:
    return {
        str(k): int(v)
        for k, v in series.fillna("").astype(str).value_counts().sort_index().head(limit).items()
    }


def verify_injection_peak_table(
    *,
    peak_table: Path,
    min_injections: int = 1,
    min_candidate_rows: int = 1,
    min_apertures: int = 1,
    min_positive_peak_ranks: int = 1,
    require_cadence_diagnostics: bool = True,
    require_signal_peaks: bool = True,
) -> dict[str, Any]:
    df = _read_table(peak_table)
    failures: list[str] = []
    missing = [col for col in REQUIRED_COLUMNS if col not in df]
    if missing:
        failures.append("missing required columns: " + ", ".join(missing))

    missing_cadence = [col for col in CADENCE_DIAGNOSTIC_COLUMNS if col not in df]
    if require_cadence_diagnostics and missing_cadence:
        failures.append("missing cadence diagnostic columns: " + ", ".join(missing_cadence))

    if "injection_id" in df:
        n_injections = int(df["injection_id"].fillna("").astype(str).nunique())
    else:
        n_injections = 0
    if n_injections < min_injections:
        failures.append(f"unique injections={n_injections}; expected at least {min_injections}")

    candidate = pd.Series(False, index=df.index)
    if "is_candidate_peak" in df:
        candidate = _as_bool(df["is_candidate_peak"])
    elif "peak_rank" in df:
        candidate = pd.to_numeric(df["peak_rank"], errors="coerce").gt(0)
    n_candidate = int(candidate.sum())
    if n_candidate < min_candidate_rows:
        failures.append(f"candidate peak rows={n_candidate}; expected at least {min_candidate_rows}")

    if "aperture" in df:
        n_apertures = int(df.loc[candidate, "aperture"].fillna("").astype(str).nunique())
    else:
        n_apertures = 0
    if n_apertures < min_apertures:
        failures.append(f"aperture count={n_apertures}; expected at least {min_apertures}")

    if "peak_rank" in df:
        peak_rank = pd.to_numeric(df["peak_rank"], errors="coerce")
        n_positive_peak_ranks = int(peak_rank.loc[candidate & peak_rank.gt(0)].dropna().nunique())
    else:
        peak_rank = pd.Series(np.nan, index=df.index)
        n_positive_peak_ranks = 0
    if n_positive_peak_ranks < min_positive_peak_ranks:
        failures.append(
            f"positive peak-rank count={n_positive_peak_ranks}; expected at least {min_positive_peak_ranks}"
        )

    signal_peak = pd.Series(False, index=df.index)
    if "is_injected_signal_peak" in df:
        signal_peak = _as_bool(df["is_injected_signal_peak"])
    n_signal_rows = int((candidate & signal_peak).sum())
    if require_signal_peaks and n_signal_rows <= 0:
        failures.append("no candidate rows are labeled as injected signal peaks")

    for col in ("truth_period_d", "truth_t0_bjd", "truth_duration_min"):
        if col in df:
            bad = int((~_finite(df, col)).sum())
            if bad:
                failures.append(f"{bad} rows have non-finite {col}")
    for col in ("period_d", "t0_bjd", "duration_min", "sde"):
        if col in df:
            bad = int((candidate & ~_finite(df, col)).sum())
            if bad:
                failures.append(f"{bad} candidate rows have non-finite {col}")
    if require_cadence_diagnostics:
        for col in CADENCE_DIAGNOSTIC_COLUMNS:
            if col in df:
                bad = int((candidate & ~_nonnegative(df, col)).sum())
                if bad:
                    failures.append(f"{bad} candidate rows have negative/non-finite {col}")

    window_columns_present = [col for col in WINDOW_OVERLAP_COLUMNS if col in df]
    if window_columns_present:
        missing_window = [col for col in WINDOW_OVERLAP_COLUMNS if col not in df]
        if missing_window:
            failures.append("partial transit-window overlap schema; missing: " + ", ".join(missing_window))
        signal_window = pd.Series(False, index=df.index)
        if "transit_window_match" in df:
            signal_window = _as_bool(df["transit_window_match"])
        bad_signal_window = int((candidate & signal_peak & ~signal_window).sum())
        if bad_signal_window:
            failures.append(f"{bad_signal_window} signal peak rows lack transit-window overlap")
        for col in (
            "transit_window_phase_period_d",
            "transit_window_delta_min",
            "transit_window_overlap_min",
            "transit_window_overlap_fraction",
        ):
            if col in df:
                bad = int((candidate & ~_finite(df, col)).sum())
                if bad:
                    failures.append(f"{bad} candidate rows have non-finite {col}")

    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "peak_table": str(peak_table),
        "n_rows": int(len(df)),
        "n_injections": n_injections,
        "n_candidate_rows": n_candidate,
        "n_signal_peak_rows": n_signal_rows,
        "n_apertures": n_apertures,
        "n_positive_peak_ranks": n_positive_peak_ranks,
        "require_cadence_diagnostics": bool(require_cadence_diagnostics),
        "require_signal_peaks": bool(require_signal_peaks),
        "cadence_diagnostic_columns": [col for col in CADENCE_DIAGNOSTIC_COLUMNS if col in df],
        "transit_window_overlap_columns": window_columns_present,
        "status_counts": _counts(df["status"]) if "status" in df else {},
        "match_kind_counts": _counts(df["match_kind"]) if "match_kind" in df else {},
        "aperture_counts": _counts(df["aperture"]) if "aperture" in df else {},
        "search_branch_counts": _counts(df["search_branch"]) if "search_branch" in df else {},
        "peak_rank_counts": _counts(peak_rank.astype("Int64").astype(str)) if "peak_rank" in df else {},
        "passed": not failures,
        "failures": failures,
    }
    return payload


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peak-table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--out-json", type=Path, default=None)
    parser.add_argument("--min-injections", type=int, default=1)
    parser.add_argument("--min-candidate-rows", type=int, default=1)
    parser.add_argument("--min-apertures", type=int, default=1)
    parser.add_argument("--min-positive-peak-ranks", type=int, default=1)
    parser.add_argument("--require-cadence-diagnostics", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--require-signal-peaks", action=argparse.BooleanOptionalAction, default=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    result = verify_injection_peak_table(
        peak_table=args.peak_table,
        min_injections=args.min_injections,
        min_candidate_rows=args.min_candidate_rows,
        min_apertures=args.min_apertures,
        min_positive_peak_ranks=args.min_positive_peak_ranks,
        require_cadence_diagnostics=bool(args.require_cadence_diagnostics),
        require_signal_peaks=bool(args.require_signal_peaks),
    )
    text = json.dumps(result, indent=2, sort_keys=True, default=_json_default) + "\n"
    if args.out_json is not None:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(text)
    print(text, end="")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
