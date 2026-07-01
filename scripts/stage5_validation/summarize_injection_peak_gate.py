#!/usr/bin/env python3
"""Summarize whether injected truth is recoverable in BLS peak lists.

This is the decision gate between search work and human-label scale-up:

* if injected truth is often present in top-N but not top-1, train/apply a peak
  ranker before LEO/human review;
* if injected truth is absent from top-N in important parameter regions,
  improve the search/detrending before asking for more labels;
* cells with low empirical top-N recovery are marked as below-threshold for the
  current BLS setup, not as human-label failures.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd


DEFAULT_PEAK_TABLE = (
    Path(__file__).resolve().parents[2]
    / "reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/"
    / "peak_training/s56_20k_injection_bls_peaks_chunked.csv"
)
DEFAULT_OUT_DIR = (
    Path(__file__).resolve().parents[2]
    / "reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/"
    / "peak_training_gate"
)
DEFAULT_PERIOD_BINS = (0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0)
DEFAULT_RADIUS_BINS = (0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 16.0, 24.0)
DEFAULT_TMAG_BINS = (0.0, 17.0, 18.0, 19.0, 99.0)


def _read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".json", ".jsonl"}:
        return pd.read_json(path, lines=suffix == ".jsonl")
    raise ValueError(f"unsupported table format: {path}")


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _parse_float_tuple(raw: str) -> tuple[float, ...]:
    return tuple(float(part.strip()) for part in str(raw).split(",") if part.strip())


def _label_intervals(bins: tuple[float, ...], *, fmt: str = "g") -> list[str]:
    labels = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        labels.append(f"[{lo:{fmt}},{hi:{fmt}})")
    return labels


def _first_existing(df: pd.DataFrame, candidates: tuple[str, ...]) -> str | None:
    for col in candidates:
        if col in df.columns:
            return col
    return None


def build_injection_level_table(peaks: pd.DataFrame, *, top_k: int = 20) -> pd.DataFrame:
    """Collapse a peak table to one row per injected light curve."""

    if "injection_id" not in peaks:
        raise KeyError("peak table must include injection_id")
    if "is_injected_signal_peak" not in peaks:
        raise KeyError("peak table must include is_injected_signal_peak")
    if "peak_rank" not in peaks:
        raise KeyError("peak table must include peak_rank")

    work = peaks.copy()
    if "is_candidate_peak" in work:
        work["_is_candidate_peak"] = work["is_candidate_peak"].fillna(False).astype(bool)
    else:
        work["_is_candidate_peak"] = True
    work["_is_signal_peak"] = work["is_injected_signal_peak"].fillna(False).astype(bool)
    work["_peak_rank"] = pd.to_numeric(work["peak_rank"], errors="coerce")
    if "exact_ephemeris_match" in work:
        work["_is_exact_peak"] = work["exact_ephemeris_match"].fillna(False).astype(bool)
    else:
        work["_is_exact_peak"] = False
    if "harmonic_ephemeris_match" in work:
        work["_is_harmonic_peak"] = work["harmonic_ephemeris_match"].fillna(False).astype(bool)
    else:
        work["_is_harmonic_peak"] = False
    work["_signal_top1"] = work["_is_signal_peak"] & work["_peak_rank"].le(1)
    work[f"_signal_top{top_k}"] = work["_is_signal_peak"] & work["_peak_rank"].le(top_k)
    work[f"_exact_top{top_k}"] = work["_is_exact_peak"] & work["_peak_rank"].le(top_k)
    work[f"_harmonic_top{top_k}"] = work["_is_harmonic_peak"] & work["_peak_rank"].le(top_k)

    truth_cols = [
        col
        for col in (
            "tic",
            "sector",
            "cam",
            "ccd",
            "tmag",
            "signal_family",
            "truth_period_d",
            "truth_t0_bjd",
            "truth_duration_min",
            "truth_depth",
            "truth_model_depth",
            "truth_sampled_model_depth",
            "truth_radius_rearth",
            "truth_radius_rwd",
            "truth_impact_b",
            "truth_n_in_transit",
            "truth_n_good_in_transit",
            "truth_grid_cell_id",
            "truth_grid_period_bin",
            "truth_grid_radius_bin",
        )
        if col in work.columns
    ]
    first = work.groupby("injection_id", dropna=False)[truth_cols].first() if truth_cols else pd.DataFrame()

    grouped = work.groupby("injection_id", dropna=False)
    out = pd.DataFrame(index=grouped.size().index)
    out["n_peak_rows"] = grouped.size()
    out["n_candidate_peak_rows"] = grouped["_is_candidate_peak"].sum().astype(int)
    out["top1_match"] = grouped["_signal_top1"].any()
    out[f"top{top_k}_match"] = grouped[f"_signal_top{top_k}"].any()
    out[f"top{top_k}_exact_match"] = grouped[f"_exact_top{top_k}"].any()
    out[f"top{top_k}_harmonic_match"] = grouped[f"_harmonic_top{top_k}"].any()
    signal_ranks = work.loc[work["_is_signal_peak"]].groupby("injection_id")["_peak_rank"].min()
    out["best_signal_peak_rank"] = out.index.to_series().map(signal_ranks)
    out["ranking_loss"] = out[f"top{top_k}_match"] & ~out["top1_match"]
    out[f"not_in_top{top_k}"] = ~out[f"top{top_k}_match"]
    if not first.empty:
        out = out.join(first)
    return out.reset_index()


def add_recovery_bins(
    injections: pd.DataFrame,
    *,
    period_bins: tuple[float, ...],
    radius_bins: tuple[float, ...],
    tmag_bins: tuple[float, ...],
) -> pd.DataFrame:
    out = injections.copy()
    period_col = _first_existing(out, ("truth_period_d", "period_d"))
    radius_col = _first_existing(out, ("truth_radius_rearth", "radius_rearth"))
    tmag_col = _first_existing(out, ("tmag", "truth_tmag"))
    if period_col:
        out["period_bin"] = pd.cut(
            pd.to_numeric(out[period_col], errors="coerce"),
            bins=list(period_bins),
            labels=_label_intervals(period_bins),
            include_lowest=True,
        ).astype(object)
    else:
        out["period_bin"] = ""
    if radius_col:
        out["radius_bin"] = pd.cut(
            pd.to_numeric(out[radius_col], errors="coerce"),
            bins=list(radius_bins),
            labels=_label_intervals(radius_bins),
            include_lowest=True,
        ).astype(object)
    else:
        out["radius_bin"] = ""
    if tmag_col:
        out["tmag_bin"] = pd.cut(
            pd.to_numeric(out[tmag_col], errors="coerce"),
            bins=list(tmag_bins),
            labels=_label_intervals(tmag_bins),
            include_lowest=True,
        ).astype(object)
    else:
        out["tmag_bin"] = ""
    return out


def build_cell_recovery_table(
    injections: pd.DataFrame,
    *,
    top_k: int,
    min_cell_count: int,
    recovery_threshold: float,
) -> pd.DataFrame:
    topk_col = f"top{top_k}_match"
    rows = []
    for keys, group in injections.groupby(["tmag_bin", "period_bin", "radius_bin"], dropna=False):
        tmag_bin, period_bin, radius_bin = keys
        n = int(len(group))
        top1_n = int(group["top1_match"].fillna(False).sum())
        topk_n = int(group[topk_col].fillna(False).sum())
        supported = n >= min_cell_count
        topk_frac = float(topk_n / n) if n else float("nan")
        if not supported:
            gate = "unsupported_cell"
        elif topk_frac >= recovery_threshold:
            gate = "above_empirical_threshold"
        else:
            gate = "below_empirical_threshold"
        rows.append(
            {
                "tmag_bin": str(tmag_bin),
                "period_bin": str(period_bin),
                "radius_bin": str(radius_bin),
                "n": n,
                "top1_n": top1_n,
                "top1_fraction": float(top1_n / n) if n else float("nan"),
                f"top{top_k}_n": topk_n,
                f"top{top_k}_fraction": topk_frac,
                "ranking_loss_n": int(group["ranking_loss"].fillna(False).sum()),
                "supported": bool(supported),
                "recovery_gate": gate,
            }
        )
    return pd.DataFrame(rows).sort_values(["tmag_bin", "period_bin", "radius_bin"]).reset_index(drop=True)


def annotate_injections_with_cells(injections: pd.DataFrame, cells: pd.DataFrame, *, top_k: int) -> pd.DataFrame:
    if cells.empty:
        return injections.copy()
    key_cols = ["tmag_bin", "period_bin", "radius_bin"]
    metric_cols = ["n", f"top{top_k}_fraction", "recovery_gate"]
    merged = injections.merge(
        cells[key_cols + metric_cols].rename(
            columns={
                "n": "recovery_cell_n",
                f"top{top_k}_fraction": f"recovery_cell_top{top_k}_fraction",
            }
        ),
        on=key_cols,
        how="left",
    )
    return merged


def _fraction(n: int, denom: int) -> float:
    return float(n / denom) if denom else float("nan")


def summarize_gate(
    injections: pd.DataFrame,
    cells: pd.DataFrame,
    *,
    top_k: int,
    min_topk_fraction_for_ranker: float,
    min_ranking_loss_fraction: float,
    ranker_summary: dict[str, Any] | None = None,
) -> dict[str, Any]:
    n = int(len(injections))
    topk_col = f"top{top_k}_match"
    top1_n = int(injections["top1_match"].fillna(False).sum())
    topk_n = int(injections[topk_col].fillna(False).sum())
    ranking_loss_n = int(injections["ranking_loss"].fillna(False).sum())
    below_col = f"not_in_top{top_k}"
    below_n = int(injections[below_col].fillna(False).sum())
    summary: dict[str, Any] = {
        "n_injections": n,
        "top_k": int(top_k),
        "top1": {"n": top1_n, "fraction": _fraction(top1_n, n)},
        f"top{top_k}": {"n": topk_n, "fraction": _fraction(topk_n, n)},
        "ranking_loss": {"n": ranking_loss_n, "fraction": _fraction(ranking_loss_n, n)},
        f"not_in_top{top_k}": {"n": below_n, "fraction": _fraction(below_n, n)},
        "cell_counts": {
            str(key): int(value)
            for key, value in cells["recovery_gate"].value_counts().sort_index().items()
        }
        if "recovery_gate" in cells
        else {},
    }
    topk_fraction = summary[f"top{top_k}"]["fraction"]
    ranking_loss_fraction = summary["ranking_loss"]["fraction"]
    if np.isfinite(topk_fraction) and topk_fraction < min_topk_fraction_for_ranker:
        recommendation = "improve_search_or_limit_to_empirical_recovery_region"
    elif np.isfinite(ranking_loss_fraction) and ranking_loss_fraction >= min_ranking_loss_fraction:
        recommendation = "train_peak_ranker_before_leo_human_review"
    else:
        recommendation = "bls_top_peak_is_adequate_for_current_recovered_region"
    summary["recommendation"] = recommendation
    if ranker_summary:
        summary["ranker_summary"] = ranker_summary
    return summary


def _extract_ranker_summary(path: Path | None) -> dict[str, Any] | None:
    if path is None:
        return None
    raw = json.loads(path.read_text())
    out: dict[str, Any] = {
        "path": str(path),
        "n_injections": raw.get("n_injections"),
        "n_rankable_injections": raw.get("n_rankable_injections"),
    }
    all_split = raw.get("splits", {}).get("all", {})
    for name in ("model_recall_at_k", "bls_rank_recall_at_k", "sde_recall_at_k"):
        if name in all_split:
            out[name] = all_split[name]
    return out


def write_markdown_report(summary: dict[str, Any], *, out_path: Path) -> None:
    top_k = int(summary["top_k"])
    top1 = summary["top1"]
    topk = summary[f"top{top_k}"]
    ranking_loss = summary["ranking_loss"]
    not_topk = summary[f"not_in_top{top_k}"]
    lines = [
        "# S56 Peak-Recovery Gate",
        "",
        f"- injections audited: `{summary['n_injections']:,}`",
        f"- top-1 truth recovery: `{top1['n']:,}/{summary['n_injections']:,}` "
        f"(`{top1['fraction']:.3f}`)",
        f"- top-{top_k} truth recovery: `{topk['n']:,}/{summary['n_injections']:,}` "
        f"(`{topk['fraction']:.3f}`)",
        f"- ranking losses: `{ranking_loss['n']:,}/{summary['n_injections']:,}` "
        f"(`{ranking_loss['fraction']:.3f}`)",
        f"- not in top-{top_k}: `{not_topk['n']:,}/{summary['n_injections']:,}` "
        f"(`{not_topk['fraction']:.3f}`)",
        f"- recommendation: `{summary['recommendation']}`",
        "",
        "Interpretation:",
        "",
        "- `ranking_loss` means BLS found the injected ephemeris in the top-N list, but not at rank 1.",
        "- `not_in_topN` means the current search did not surface the injected ephemeris in the audited peak list.",
        "- Low-recovery cells should be treated as below-threshold for this search setup, not as human-label failures.",
        "",
    ]
    if summary.get("cell_counts"):
        lines.append("Cell support/recovery gates:")
        for key, value in summary["cell_counts"].items():
            lines.append(f"- `{key}`: `{value}` cells")
        lines.append("")
    if summary.get("ranker_summary"):
        lines.append("Ranker summary was attached in the JSON output.")
        lines.append("")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peak-table", type=Path, default=DEFAULT_PEAK_TABLE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument("--period-bins", default=",".join(str(v) for v in DEFAULT_PERIOD_BINS))
    parser.add_argument("--radius-bins", default=",".join(str(v) for v in DEFAULT_RADIUS_BINS))
    parser.add_argument("--tmag-bins", default=",".join(str(v) for v in DEFAULT_TMAG_BINS))
    parser.add_argument("--min-cell-count", type=int, default=20)
    parser.add_argument("--recovery-threshold", type=float, default=0.50)
    parser.add_argument("--min-topk-fraction-for-ranker", type=float, default=0.50)
    parser.add_argument("--min-ranking-loss-fraction", type=float, default=0.05)
    parser.add_argument("--ranker-summary-json", type=Path, default=None)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    peaks = _read_table(args.peak_table)
    injections = build_injection_level_table(peaks, top_k=args.top_k)
    injections = add_recovery_bins(
        injections,
        period_bins=_parse_float_tuple(args.period_bins),
        radius_bins=_parse_float_tuple(args.radius_bins),
        tmag_bins=_parse_float_tuple(args.tmag_bins),
    )
    cells = build_cell_recovery_table(
        injections,
        top_k=args.top_k,
        min_cell_count=args.min_cell_count,
        recovery_threshold=args.recovery_threshold,
    )
    annotated = annotate_injections_with_cells(injections, cells, top_k=args.top_k)
    summary = summarize_gate(
        injections,
        cells,
        top_k=args.top_k,
        min_topk_fraction_for_ranker=args.min_topk_fraction_for_ranker,
        min_ranking_loss_fraction=args.min_ranking_loss_fraction,
        ranker_summary=_extract_ranker_summary(args.ranker_summary_json),
    )
    summary.update(
        {
            "peak_table": str(args.peak_table),
            "out_dir": str(args.out_dir),
            "period_bins": list(_parse_float_tuple(args.period_bins)),
            "radius_bins": list(_parse_float_tuple(args.radius_bins)),
            "tmag_bins": list(_parse_float_tuple(args.tmag_bins)),
            "min_cell_count": int(args.min_cell_count),
            "recovery_threshold": float(args.recovery_threshold),
        }
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    annotated.to_csv(args.out_dir / "injection_recovery_gate_by_injection.csv", index=False)
    cells.to_csv(args.out_dir / "period_radius_tmag_recovery_cells.csv", index=False)
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )
    write_markdown_report(summary, out_path=args.out_dir / "summary.md")

    print("[peak-gate] complete")
    print(f"  injections: {summary['n_injections']:,}")
    print(
        f"  top1/top{args.top_k}: "
        f"{summary['top1']['n']:,}/{summary[f'top{args.top_k}']['n']:,}"
    )
    print(f"  recommendation: {summary['recommendation']}")
    print(f"  out: {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
