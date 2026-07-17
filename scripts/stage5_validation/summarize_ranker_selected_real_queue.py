#!/usr/bin/env python3
"""Summarize ranker-selected real-candidate ephemerides before LEO rendering."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


DEFAULT_SELECTED = (
    Path(__file__).resolve().parents[2]
    / "reports/stage5_validation/s56_ranker_selected_real_candidates_pdo/selected_ephemerides.csv"
)
DEFAULT_REVIEW_QUEUE = (
    Path(__file__).resolve().parents[2]
    / "reports/stage5_validation/s56_ranker_selected_real_review_queue_pdo/review_queue.csv"
)
DEFAULT_VERIFICATION = DEFAULT_REVIEW_QUEUE.with_name("verification.json")
DEFAULT_OUT_DIR = DEFAULT_REVIEW_QUEUE.parent / "ranker_selection_summary"

PERIOD_BINS = (0.0, 0.12, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, np.inf)
PERIOD_LABELS = ("<0.12", "0.12-0.25", "0.25-0.5", "0.5-1", "1-2", "2-5", "5-10", ">10")
TMAG_BINS = (-np.inf, 17.0, 18.0, 19.0, np.inf)
TMAG_LABELS = ("<17", "17-18", "18-19", ">19")
DURATION_BINS = (0.0, 3.0, 6.0, 10.0, 20.0, 30.0, np.inf)
DURATION_LABELS = ("<3", "3-6", "6-10", "10-20", "20-30", ">30")


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


def _score_column(df: pd.DataFrame) -> str:
    for col in ("ranker_p_signal_peak", "ranker_score"):
        if col in df:
            return col
    score_cols = [col for col in df.columns if col.startswith("ranker_p_")]
    if score_cols:
        return score_cols[0]
    return ""


def _int_counts(series: pd.Series) -> dict[str, int]:
    return {str(k): int(v) for k, v in series.fillna("").astype(str).value_counts().sort_index().items()}


def _numeric_summary(df: pd.DataFrame, cols: tuple[str, ...]) -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {}
    for col in cols:
        if col not in df:
            continue
        value = pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
        if value.empty:
            continue
        out[col] = {
            "min": float(value.min()),
            "p10": float(value.quantile(0.10)),
            "median": float(value.median()),
            "p90": float(value.quantile(0.90)),
            "max": float(value.max()),
        }
    return out


def _add_bins(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "period_d" in out:
        out["period_bin"] = pd.cut(
            pd.to_numeric(out["period_d"], errors="coerce"),
            PERIOD_BINS,
            labels=PERIOD_LABELS,
        ).astype(object)
    if "tmag" in out:
        out["tmag_bin"] = pd.cut(
            pd.to_numeric(out["tmag"], errors="coerce"),
            TMAG_BINS,
            labels=TMAG_LABELS,
        ).astype(object)
    if "duration_min" in out:
        out["duration_bin"] = pd.cut(
            pd.to_numeric(out["duration_min"], errors="coerce"),
            DURATION_BINS,
            labels=DURATION_LABELS,
        ).astype(object)
    return out


def _coverage_table(df: pd.DataFrame, col: str) -> pd.DataFrame:
    if col not in df:
        return pd.DataFrame(columns=[col, "n"])
    counts = df[col].fillna("").astype(str).value_counts().sort_index()
    return counts.rename_axis(col).reset_index(name="n")


def _top_rows(df: pd.DataFrame, score_col: str, n: int = 25) -> pd.DataFrame:
    keep = [
        col
        for col in (
            "review_id",
            "tic",
            "ranker_selection_rank",
            "aperture",
            "rep_aperture",
            "peak_rank",
            "period_d",
            "duration_min",
            "depth",
            "depth_snr",
            "sde",
            "sde_max",
            score_col,
            "tmag",
        )
        if col and col in df
    ]
    if not keep:
        return pd.DataFrame()
    sort_col = score_col if score_col and score_col in df else ("sde" if "sde" in df else keep[0])
    return df.sort_values(sort_col, ascending=False, na_position="last").loc[:, keep].head(n).copy()


def summarize_ranker_selection(
    *,
    selected_ephemerides: Path,
    review_queue: Path | None,
    verification_json: Path | None,
    out_dir: Path,
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    selected = _add_bins(_read_table(selected_ephemerides))
    review = _add_bins(_read_table(review_queue)) if review_queue is not None and review_queue.exists() else pd.DataFrame()
    score_col = _score_column(selected)
    if score_col and score_col in selected:
        selected[score_col] = pd.to_numeric(selected[score_col], errors="coerce")

    working = review if not review.empty else selected
    if "aperture" not in working and "rep_aperture" in working:
        working = working.assign(aperture=working["rep_aperture"])

    verification = (
        json.loads(verification_json.read_text())
        if verification_json is not None and verification_json.exists()
        else {}
    )
    if score_col:
        selected["_ranker_score"] = pd.to_numeric(selected[score_col], errors="coerce")
    else:
        selected["_ranker_score"] = np.nan

    rows_per_target = selected.groupby("tic").size() if "tic" in selected else pd.Series(dtype=int)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "selected_ephemerides": str(selected_ephemerides),
        "review_queue": str(review_queue) if review_queue is not None else "",
        "verification_json": str(verification_json) if verification_json is not None else "",
        "out_dir": str(out_dir),
        "score_column": score_col,
        "n_selected_rows": int(len(selected)),
        "n_selected_targets": int(selected["tic"].nunique()) if "tic" in selected else 0,
        "n_review_rows": int(len(review)) if not review.empty else 0,
        "n_review_targets": int(review["tic"].nunique()) if "tic" in review else 0,
        "rows_per_target": _int_counts(rows_per_target.astype(str)) if len(rows_per_target) else {},
        "selection_rank_counts": _int_counts(selected["ranker_selection_rank"]) if "ranker_selection_rank" in selected else {},
        "aperture_counts": _int_counts(working["aperture"]) if "aperture" in working else {},
        "period_bin_counts": _int_counts(working["period_bin"]) if "period_bin" in working else {},
        "tmag_bin_counts": _int_counts(working["tmag_bin"]) if "tmag_bin" in working else {},
        "duration_bin_counts": _int_counts(working["duration_bin"]) if "duration_bin" in working else {},
        "numeric_summary": _numeric_summary(
            working,
            ("period_d", "duration_min", "depth", "depth_snr", "sde", "sde_max", "tmag"),
        ),
        "score_summary": _numeric_summary(selected, ("_ranker_score",)),
        "verification_passed": verification.get("passed", None),
        "verification_failures": verification.get("failures", []),
    }

    _coverage_table(working, "aperture").to_csv(out_dir / "by_aperture.csv", index=False)
    _coverage_table(working, "period_bin").to_csv(out_dir / "by_period_bin.csv", index=False)
    _coverage_table(working, "tmag_bin").to_csv(out_dir / "by_tmag_bin.csv", index=False)
    _coverage_table(working, "duration_bin").to_csv(out_dir / "by_duration_bin.csv", index=False)
    _coverage_table(selected, "ranker_selection_rank").to_csv(out_dir / "by_selection_rank.csv", index=False)
    _top_rows(selected, score_col).to_csv(out_dir / "top_ranker_rows.csv", index=False)
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    (out_dir / "summary.md").write_text(_render_markdown(summary))
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return summary


def _render_counts(counts: dict[str, int], *, limit: int = 12) -> list[str]:
    if not counts:
        return ["- none"]
    return [f"- `{key}`: `{value}`" for key, value in list(counts.items())[:limit]]


def _render_markdown(summary: dict[str, Any]) -> str:
    lines = [
        "# Ranker-Selected Real Queue Summary",
        "",
        f"Selected rows: `{summary['n_selected_rows']}` across `{summary['n_selected_targets']}` TICs.",
        f"Review rows: `{summary['n_review_rows']}` across `{summary['n_review_targets']}` TICs.",
        f"Score column: `{summary['score_column'] or 'none'}`.",
        f"Verification passed: `{summary['verification_passed']}`.",
        "",
        "## Selection Ranks",
        "",
        *_render_counts(summary["selection_rank_counts"]),
        "",
        "## Apertures",
        "",
        *_render_counts(summary["aperture_counts"]),
        "",
        "## Period Bins",
        "",
        *_render_counts(summary["period_bin_counts"]),
        "",
        "## Tmag Bins",
        "",
        *_render_counts(summary["tmag_bin_counts"]),
    ]
    failures = summary.get("verification_failures") or []
    if failures:
        lines.extend(["", "## Verification Failures", ""])
        lines.extend(f"- {failure}" for failure in failures)
    lines.append("")
    return "\n".join(lines)


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected-ephemerides", type=Path, default=DEFAULT_SELECTED)
    parser.add_argument("--review-queue", type=Path, default=DEFAULT_REVIEW_QUEUE)
    parser.add_argument("--verification-json", type=Path, default=DEFAULT_VERIFICATION)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    summarize_ranker_selection(
        selected_ephemerides=args.selected_ephemerides,
        review_queue=args.review_queue,
        verification_json=args.verification_json,
        out_dir=args.out_dir,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
