#!/usr/bin/env python3
"""Build a small visual-audit queue from strict BLS-recovered injections."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))


DEFAULT_RECOVERY_CSV = (
    REPO_ROOT
    / "reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo"
    / "small_pair_200k/injection_bls_recoveries.csv"
)
DEFAULT_OUT_DIR = REPO_ROOT / "reports/stage5_validation/s56_bls_top1_recovery_vet_check"
DEFAULT_BRANCH = "current_adp"
DEFAULT_APERTURES = ("DET_FLUX_ADP_SML", "DET_FLUX_SML")


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _safe_id(value: Any) -> str:
    return str(value).replace("/", "_").replace(":", "_")


def _cut_with_fallback(values: pd.Series, edges: list[float], prefix: str) -> pd.Series:
    vals = pd.to_numeric(values, errors="coerce")
    labels = [f"{prefix}{i:02d}" for i in range(len(edges) - 1)]
    out = pd.cut(vals, bins=edges, labels=labels, include_lowest=True)
    return out.astype("object").fillna(f"{prefix}na")


def classify_recovery(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["strict_top1_recovered"] = out["recovery_status"].astype(str).eq("bls_recovered")
    out["topn_exact_recovered"] = out["topn_recovery_status"].astype(str).eq("bls_topn_recovered")
    out["topn_harmonic_recovered"] = out["topn_recovery_status"].astype(str).eq("bls_topn_harmonic_match")
    out["any_exact_or_harmonic_recovered"] = (
        out["strict_top1_recovered"] | out["topn_exact_recovered"] | out["topn_harmonic_recovered"]
    )
    return out


def select_visual_audit_rows(frame: pd.DataFrame, *, n: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    df = classify_recovery(frame)
    numeric_cols = [
        "truth_period_d",
        "truth_radius_rearth",
        "truth_sampled_model_depth",
        "truth_n_good_in_transit",
        "tmag",
        "sde_max",
    ]
    for col in numeric_cols:
        if col in df:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    keep = (
        df["strict_top1_recovered"]
        & np.isfinite(df["truth_period_d"])
        & np.isfinite(df["truth_radius_rearth"])
        & np.isfinite(df["tmag"])
        & df["source_h5"].astype(str).ne("")
        & df["h5_group"].astype(str).ne("")
    )
    if "truth_n_good_in_transit" in df:
        keep &= df["truth_n_good_in_transit"].fillna(0) >= 2
    candidates = df.loc[keep].copy()
    if len(candidates) < n:
        raise ValueError(f"only {len(candidates):,} strict top-1 candidates available; requested {n:,}")

    candidates["selection_period_bin"] = _cut_with_fallback(
        candidates["truth_period_d"],
        [0.08, 0.12, 0.20, 0.40, 0.80, 1.60, 3.20, 6.40, 13.10],
        "p",
    )
    candidates["selection_radius_bin"] = _cut_with_fallback(
        candidates["truth_radius_rearth"],
        [0.20, 0.50, 1.00, 2.00, 4.00, 8.00, 13.00, 18.10],
        "r",
    )
    candidates["selection_tmag_bin"] = _cut_with_fallback(
        candidates["tmag"],
        [0.0, 15.5, 16.5, 17.5, 18.5, 25.0],
        "t",
    )
    candidates["selection_cell"] = (
        candidates["selection_period_bin"].astype(str)
        + "_"
        + candidates["selection_radius_bin"].astype(str)
        + "_"
        + candidates["selection_tmag_bin"].astype(str)
    )

    selected_positions: list[int] = []
    groups = list(candidates.groupby("selection_cell", observed=True).indices.items())
    order = rng.permutation(len(groups))
    for group_index in order:
        _, indices = groups[int(group_index)]
        selected_positions.append(int(rng.choice(indices)))
        if len(selected_positions) == n:
            break

    if len(selected_positions) < n:
        remaining = np.setdiff1d(np.arange(len(candidates), dtype=int), np.asarray(selected_positions, dtype=int))
        fill = rng.choice(remaining, size=n - len(selected_positions), replace=False)
        selected_positions.extend(int(v) for v in fill)

    selected = candidates.iloc[selected_positions].sample(frac=1.0, random_state=seed).reset_index(drop=True)
    selected["selection_seed"] = int(seed)
    selected["selection_bucket"] = "strict_top1_bls_recovered_injection_visual_audit"
    selected["selection_weight"] = 1.0
    selected["visual_audit_row"] = np.arange(len(selected), dtype=int)
    selected["twirl_vet_sheet_name"] = selected["review_id"].map(
        lambda value: f"{_safe_id(value)}_twirl_twoap_{DEFAULT_BRANCH}.png"
    )
    selected["twirl_vet_sheet_pdf_name"] = selected["review_id"].map(
        lambda value: f"{_safe_id(value)}_twirl_twoap_{DEFAULT_BRANCH}.pdf"
    )
    return selected


def build_queue(*, recovery_csv: Path, out_dir: Path, n: int, seed: int) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    recovery = pd.read_csv(recovery_csv)
    queue = select_visual_audit_rows(recovery, n=n, seed=seed)
    queue_csv = out_dir / f"review_queue_{n}.csv"
    queue.to_csv(queue_csv, index=False)
    labels_csv = out_dir / "human_labels_vetted.csv"
    if not labels_csv.exists():
        labels_csv.write_text("review_id,tic,label,notes,updated_utc\n")

    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "recovery_csv": str(recovery_csv),
        "queue_csv": str(queue_csv),
        "n_rows": int(len(queue)),
        "seed": int(seed),
        "selection": "stratified random sample from strict top-1 BLS-recovered injections",
        "recovery_map_input_rows": int(len(recovery)),
        "strict_top1_recovered_rows": int(recovery["recovery_status"].astype(str).eq("bls_recovered").sum()),
        "topn_exact_rows": int(recovery["topn_recovery_status"].astype(str).eq("bls_topn_recovered").sum()),
        "topn_harmonic_rows": int(recovery["topn_recovery_status"].astype(str).eq("bls_topn_harmonic_match").sum()),
        "branch_for_render": DEFAULT_BRANCH,
        "apertures_for_render": list(DEFAULT_APERTURES),
        "period_bin_counts": {
            str(k): int(v) for k, v in queue["selection_period_bin"].value_counts().sort_index().items()
        },
        "radius_bin_counts": {
            str(k): int(v) for k, v in queue["selection_radius_bin"].value_counts().sort_index().items()
        },
        "tmag_bin_counts": {
            str(k): int(v) for k, v in queue["selection_tmag_bin"].value_counts().sort_index().items()
        },
        "depth_pct_range": [
            float(100.0 * queue["truth_sampled_model_depth"].min()),
            float(100.0 * queue["truth_sampled_model_depth"].max()),
        ],
        "period_d_range": [
            float(queue["truth_period_d"].min()),
            float(queue["truth_period_d"].max()),
        ],
        "radius_rearth_range": [
            float(queue["truth_radius_rearth"].min()),
            float(queue["truth_radius_rearth"].max()),
        ],
    }
    summary_json = out_dir / "summary.json"
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--recovery-csv", type=Path, default=DEFAULT_RECOVERY_CSV)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--seed", type=int, default=560100)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    summary = build_queue(
        recovery_csv=args.recovery_csv,
        out_dir=args.out_dir,
        n=int(args.n),
        seed=int(args.seed),
    )
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
