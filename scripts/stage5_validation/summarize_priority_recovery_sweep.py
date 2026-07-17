#!/usr/bin/env python3
"""Summarize targeted priority recovery sweeps.

This consumes the output of ``run_s56_predetrend_priority_recovery_sweep_pdo.sh``
and joins each sweep queue back to the compact priority debug queue. The goal is
to answer whether small-aperture BLS actually recovers the rows that the
truth-window signal-survival diagnostic says should now be SNR-qualified.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


MODE_ORDER = (
    "bls_top1_recovered",
    "bls_topn_recovered",
    "bls_topn_harmonic_match",
    "bls_peak_mismatch",
)
MATCH_MODES = {"bls_top1_recovered", "bls_topn_recovered", "bls_topn_harmonic_match"}


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def _mode_series(df: pd.DataFrame) -> pd.Series:
    if "topn_recovery_status" in df:
        topn = df["topn_recovery_status"].fillna("").astype(str)
        strict = df.get("recovery_status", pd.Series("", index=df.index)).fillna("").astype(str)
        return topn.where(topn.ne(""), strict)
    return df.get("recovery_status", pd.Series("", index=df.index)).fillna("").astype(str)


def _ordered_counts(series: pd.Series) -> dict[str, int]:
    counts = series.value_counts(dropna=False).to_dict()
    out: dict[str, int] = {}
    for key in MODE_ORDER:
        if key in counts:
            out[key] = int(counts.pop(key))
    for key in sorted(counts):
        out[str(key)] = int(counts[key])
    return out


def _read_overview(sweep_dir: Path) -> pd.DataFrame:
    path = sweep_dir / "sweep_overview.csv"
    if not path.exists():
        raise FileNotFoundError(f"missing sweep overview: {path}")
    return pd.read_csv(path)


def _queue_path(row: pd.Series, sweep_dir: Path) -> Path:
    raw = str(row.get("queue_csv", ""))
    if raw:
        path = Path(raw)
        if path.exists():
            return path
        candidate = sweep_dir / path
        if candidate.exists():
            return candidate
    fallback = sweep_dir / str(row["name"]) / "review_queue.csv"
    if fallback.exists():
        return fallback
    raise FileNotFoundError(f"missing queue for sweep {row.get('name')}: {raw}")


def _summarize_joined(joined: pd.DataFrame) -> dict[str, Any]:
    mode = joined["sweep_recovery_mode"]
    matched = mode.isin(MATCH_MODES)
    out: dict[str, Any] = {
        "n": int(len(joined)),
        "n_any_match": int(matched.sum()),
        "frac_any_match": float(matched.mean()) if len(joined) else float("nan"),
        "recovery_mode_counts": _ordered_counts(mode),
    }
    for col in ("priority_reason", "recommended_aperture", "failure_class"):
        if col not in joined:
            continue
        table = pd.crosstab(joined[col].fillna("").astype(str), mode, dropna=False)
        out[f"by_{col}"] = {
            str(idx): {str(key): int(value) for key, value in row.items()}
            for idx, row in table.iterrows()
        }
    return out


def summarize_priority_sweep(
    sweep_dir: Path,
    priority_queue_csv: Path,
    out_dir: Path,
    baseline_sweep: str = "adp_priority_5k",
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    priority = pd.read_csv(priority_queue_csv)
    if "injection_id" not in priority:
        raise ValueError(f"priority queue missing injection_id: {priority_queue_csv}")
    overview = _read_overview(sweep_dir)

    per_sweep_rows: list[dict[str, Any]] = []
    joined_tables: dict[str, pd.DataFrame] = {}
    baseline_modes: pd.Series | None = None
    baseline_name_found = ""
    for _, sweep in overview.iterrows():
        name = str(sweep["name"])
        queue = pd.read_csv(_queue_path(sweep, sweep_dir))
        queue["sweep_recovery_mode"] = _mode_series(queue)
        keep = [
            col
            for col in (
                "injection_id",
                "sweep_recovery_mode",
                "recovery_status",
                "topn_recovery_status",
                "period_d",
                "sde_max",
                "rep_aperture",
                "n_topn_apertures_agree",
                "n_harmonic_apertures_agree",
            )
            if col in queue.columns
        ]
        joined = priority.merge(queue.loc[:, keep], on="injection_id", how="left", suffixes=("", "_sweep"))
        joined["sweep_name"] = name
        joined["sweep_any_match"] = joined["sweep_recovery_mode"].isin(MATCH_MODES)
        joined_tables[name] = joined
        if name == baseline_sweep:
            baseline_modes = joined.set_index("injection_id")["sweep_recovery_mode"]
            baseline_name_found = name

        summary = _summarize_joined(joined)
        per_sweep_rows.append(
            {
                "name": name,
                "apertures": sweep.get("apertures", ""),
                "n_periods": sweep.get("n_periods", ""),
                "durations_min": sweep.get("durations_min", ""),
                "n_rows": summary["n"],
                "n_any_match": summary["n_any_match"],
                "frac_any_match": summary["frac_any_match"],
                "top1": summary["recovery_mode_counts"].get("bls_top1_recovered", 0),
                "exact_topn": summary["recovery_mode_counts"].get("bls_topn_recovered", 0),
                "harmonic_topn": summary["recovery_mode_counts"].get("bls_topn_harmonic_match", 0),
                "unmatched": summary["recovery_mode_counts"].get("bls_peak_mismatch", 0),
            }
        )

    per_sweep = pd.DataFrame(per_sweep_rows)
    per_sweep.to_csv(out_dir / "priority_sweep_summary.csv", index=False)

    delta_rows: list[dict[str, Any]] = []
    if baseline_modes is not None:
        baseline_matched = baseline_modes.isin(MATCH_MODES)
        for name, joined in joined_tables.items():
            current = joined.set_index("injection_id")["sweep_recovery_mode"]
            current_matched = current.isin(MATCH_MODES)
            common = current.index.intersection(baseline_matched.index)
            newly = common[current_matched.loc[common] & ~baseline_matched.loc[common]]
            lost = common[~current_matched.loc[common] & baseline_matched.loc[common]]
            delta_rows.append(
                {
                    "name": name,
                    "baseline_sweep": baseline_name_found,
                    "new_matches_vs_baseline": int(len(newly)),
                    "lost_matches_vs_baseline": int(len(lost)),
                    "net_match_gain_vs_baseline": int(len(newly) - len(lost)),
                }
            )
            detail = joined[joined["injection_id"].isin(newly)].copy()
            detail.to_csv(out_dir / f"{name}_new_matches_vs_{baseline_name_found}.csv", index=False)
    deltas = pd.DataFrame(delta_rows)
    if not deltas.empty:
        deltas.to_csv(out_dir / "priority_sweep_deltas_vs_baseline.csv", index=False)

    for name, joined in joined_tables.items():
        joined.to_csv(out_dir / f"{name}_joined_priority_queue.csv", index=False)

    best = per_sweep.sort_values(["n_any_match", "top1"], ascending=[False, False]).head(1)
    summary_json: dict[str, Any] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sweep_dir": str(sweep_dir),
        "priority_queue_csv": str(priority_queue_csv),
        "baseline_sweep": baseline_sweep,
        "baseline_found": bool(baseline_modes is not None),
        "n_priority_rows": int(len(priority)),
        "n_sweeps": int(len(per_sweep)),
        "best_sweep": best.iloc[0].to_dict() if len(best) else {},
        "per_sweep": per_sweep.to_dict(orient="records"),
        "delta_vs_baseline": deltas.to_dict(orient="records") if not deltas.empty else [],
        "outputs": {
            "priority_sweep_summary": str(out_dir / "priority_sweep_summary.csv"),
            "priority_sweep_deltas_vs_baseline": str(out_dir / "priority_sweep_deltas_vs_baseline.csv"),
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary_json, indent=2, sort_keys=True, default=_json_default) + "\n")
    print(json.dumps(summary_json, indent=2, sort_keys=True, default=_json_default))
    return summary_json


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sweep-dir", type=Path, required=True)
    parser.add_argument("--priority-queue-csv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--baseline-sweep", default="adp_priority_5k")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    summarize_priority_sweep(
        args.sweep_dir,
        args.priority_queue_csv,
        args.out_dir,
        baseline_sweep=args.baseline_sweep,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
