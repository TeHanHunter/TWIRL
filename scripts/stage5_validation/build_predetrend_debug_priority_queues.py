#!/usr/bin/env python3
"""Build priority debug queues from pre-detrend recovery diagnostics.

The output is not a human-vetting queue. It is a compact target list for the
next search/detrending iteration:

* SNR-qualified BLS ranking losses, where the injected signal survived but BLS
  chose the wrong peak.
* Rows that become empirically recoverable only when a better aperture is used.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def _read_csv(path: Path, required: set[str]) -> pd.DataFrame:
    df = pd.read_csv(path)
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"{path} missing required columns: {missing}")
    return df


def _priority_reason(row: pd.Series) -> str:
    reasons: list[str] = []
    if row.get("failure_class") == "snr_qualified_bls_ranking_loss":
        reasons.append("search_ranking_loss")
    if bool(row.get("newly_recoverable_by_aperture", False)):
        reasons.append("newly_recoverable_by_best_aperture")
    if not reasons:
        reasons.append("context")
    return "+".join(reasons)


def _merge_debug_tables(failure_csv: Path, aperture_csv: Path) -> pd.DataFrame:
    failure = _read_csv(failure_csv, {"injection_id", "failure_class"})
    aperture = _read_csv(aperture_csv, {"injection_id", "best_aperture", "best_snr", "newly_recoverable_by_aperture"})
    keep_aperture = [
        col
        for col in (
            "injection_id",
            "current_snr",
            "current_depth_retention",
            "best_aperture",
            "best_snr",
            "best_depth_retention",
            "best_n_good_in_transit",
            "current_snr_ge_threshold",
            "best_snr_ge_threshold",
            "newly_recoverable_by_aperture",
            "snr_gain",
            "snr_gain_factor",
        )
        if col in aperture.columns
    ]
    merged = failure.merge(aperture.loc[:, keep_aperture], on="injection_id", how="left", suffixes=("", "_best"))
    if "best_aperture" not in merged:
        merged["best_aperture"] = ""
    merged["recommended_aperture"] = merged["best_aperture"].fillna("").astype(str)
    merged["priority_reason"] = merged.apply(_priority_reason, axis=1)
    return merged


def _compact_columns(df: pd.DataFrame) -> list[str]:
    preferred = [
        "priority_reason",
        "review_id",
        "injection_id",
        "tic",
        "tmag",
        "failure_class",
        "recommended_aperture",
        "rep_aperture",
        "current_snr",
        "best_snr",
        "snr_gain",
        "snr_gain_factor",
        "current_depth_retention",
        "best_depth_retention",
        "topn_recovery_status",
        "recovery_status",
        "recovery_mode",
        "period_d",
        "truth_period_d",
        "truth_t0_bjd",
        "truth_duration_min",
        "truth_model_depth",
        "truth_sampled_model_depth",
        "truth_radius_rearth",
        "truth_impact_b",
        "leo_class",
        "leo_report_name",
        "source_h5",
        "h5_group",
    ]
    cols = [col for col in preferred if col in df.columns]
    extras = [col for col in df.columns if col not in cols and col.startswith(("sde_", "rank_", "period_rel_err_", "topn_"))]
    return cols + extras


def build_priority_queues(
    failure_csv: Path,
    aperture_csv: Path,
    out_dir: Path,
    min_best_snr: float = 7.0,
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    merged = _merge_debug_tables(failure_csv, aperture_csv)
    ranking_loss = merged[merged["failure_class"].eq("snr_qualified_bls_ranking_loss")].copy()
    newly_aperture = merged[
        merged["newly_recoverable_by_aperture"].fillna(False).astype(bool)
        & pd.to_numeric(merged["best_snr"], errors="coerce").ge(min_best_snr)
    ].copy()

    priority = pd.concat([ranking_loss, newly_aperture], ignore_index=True)
    priority = priority.drop_duplicates("injection_id", keep="first").copy()
    if not priority.empty:
        priority["priority_reason"] = priority.apply(_priority_reason, axis=1)
        priority = priority.sort_values(
            ["priority_reason", "best_snr", "snr_gain"],
            ascending=[True, False, False],
            na_position="last",
        )

    compact_cols = _compact_columns(priority)
    priority.loc[:, compact_cols].to_csv(out_dir / "priority_debug_queue.csv", index=False)
    ranking_loss.loc[:, _compact_columns(ranking_loss)].to_csv(out_dir / "snr_qualified_bls_ranking_loss_queue.csv", index=False)
    newly_aperture.loc[:, _compact_columns(newly_aperture)].to_csv(out_dir / "newly_recoverable_by_best_aperture_queue.csv", index=False)

    ids = priority["injection_id"].dropna().astype(str).tolist()
    (out_dir / "priority_injection_ids.txt").write_text("\n".join(ids) + ("\n" if ids else ""))
    if "tic" in priority:
        tics = priority["tic"].dropna().astype(int).astype(str).drop_duplicates().tolist()
        (out_dir / "priority_tics.txt").write_text("\n".join(tics) + ("\n" if tics else ""))

    summary: dict[str, Any] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "failure_csv": str(failure_csv),
        "aperture_csv": str(aperture_csv),
        "min_best_snr": float(min_best_snr),
        "n_failure_rows": int(len(merged)),
        "n_snr_qualified_bls_ranking_loss": int(len(ranking_loss)),
        "n_newly_recoverable_by_best_aperture": int(len(newly_aperture)),
        "n_priority_unique": int(len(priority)),
        "priority_reason_counts": {
            str(key): int(value)
            for key, value in priority["priority_reason"].value_counts(dropna=False).sort_index().items()
        }
        if len(priority)
        else {},
        "recommended_aperture_counts": {
            str(key): int(value)
            for key, value in priority["recommended_aperture"].value_counts(dropna=False).sort_index().items()
        }
        if len(priority)
        else {},
        "outputs": {
            "priority_debug_queue": str(out_dir / "priority_debug_queue.csv"),
            "snr_qualified_bls_ranking_loss_queue": str(out_dir / "snr_qualified_bls_ranking_loss_queue.csv"),
            "newly_recoverable_by_best_aperture_queue": str(out_dir / "newly_recoverable_by_best_aperture_queue.csv"),
            "priority_injection_ids": str(out_dir / "priority_injection_ids.txt"),
            "priority_tics": str(out_dir / "priority_tics.txt"),
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--failure-csv", type=Path, required=True)
    parser.add_argument("--aperture-csv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--min-best-snr", type=float, default=7.0)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    build_priority_queues(args.failure_csv, args.aperture_csv, args.out_dir, min_best_snr=args.min_best_snr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
