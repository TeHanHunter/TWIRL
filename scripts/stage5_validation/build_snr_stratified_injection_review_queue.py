#!/usr/bin/env python3
"""Build an SNR-stratified injected-candidate review queue.

The first 1k pre-detrend queue is useful, but a random/order-preserving pass is
dominated by low-SNR unmatched injections. This helper keeps the already
rendered LEO reports from a base queue, merges signal-survival and top-N BLS
diagnostics, and selects rows across recovery/SNR strata so human vetting sees
recoverable, ambiguous, harmonic, and low-SNR examples in controlled numbers.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
from typing import Any

import numpy as np
import pandas as pd


SNR_BINS = (-np.inf, 1.0, 3.0, 5.0, 7.0, 10.0, 20.0, np.inf)
SNR_LABELS = ("snr_lt1", "snr_1_3", "snr_3_5", "snr_5_7", "snr_7_10", "snr_10_20", "snr_gt20")
RARE_MODES = ("bls_top1_recovered", "bls_topn_recovered", "bls_topn_harmonic_match")


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def _read_survival(path: Path, aperture: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "aperture" in df:
        df = df[df["aperture"].astype(str).eq(aperture)].copy()
    keep = [
        col
        for col in (
            "injection_id",
            "delta_signal_multi_snr_mad",
            "delta_signal_snr_mad",
            "depth_retention_frac",
            "depth_retention_vs_model_frac",
            "measured_delta_depth",
            "n_good_ap_in_transit",
            "n_good_ap_oot",
        )
        if col in df.columns
    ]
    return df.loc[:, keep].drop_duplicates("injection_id", keep="first")


def _read_topn(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    keep = ["injection_id"]
    keep.extend(
        col
        for col in df.columns
        if (
            col.startswith("topn_")
            or col.startswith("top_")
            or col.startswith("harmonic_")
            or col.startswith("n_topn_")
            or col.startswith("n_harmonic_")
        )
    )
    keep = list(dict.fromkeys(col for col in keep if col in df.columns))
    return df.loc[:, keep].drop_duplicates("injection_id", keep="first")


def _recovery_mode(row: pd.Series) -> str:
    topn = str(row.get("topn_recovery_status", "") or "")
    if topn and topn != "nan":
        return topn
    status = str(row.get("recovery_status", "") or "")
    return "bls_top1_recovered" if status == "bls_recovered" else status


def _assign_stratum(row: pd.Series) -> str:
    mode = str(row["recovery_mode"])
    if mode in RARE_MODES:
        return mode
    snr = row.get("delta_signal_multi_snr_mad", np.nan)
    try:
        snr_value = float(snr)
    except (TypeError, ValueError):
        snr_value = np.nan
    if not np.isfinite(snr_value):
        return f"{mode}_snr_unknown"
    label = pd.cut([snr_value], SNR_BINS, labels=SNR_LABELS)[0]
    return f"{mode}_{label}"


def _target_counts(strata: pd.Series, n_rows: int, min_per_stratum: int, rare_cap: int) -> dict[str, int]:
    available = strata.value_counts().sort_index()
    targets = {stratum: 0 for stratum in available.index}
    remaining = int(n_rows)

    rare_first = [stratum for stratum in RARE_MODES if stratum in available.index]
    other_strata = [stratum for stratum in available.index if stratum not in set(rare_first)]
    for stratum in [*rare_first, *other_strata]:
        count = available[stratum]
        if remaining <= 0:
            break
        if stratum in RARE_MODES:
            take = min(int(count), int(rare_cap), remaining)
        else:
            take = min(int(count), int(min_per_stratum), remaining)
        targets[stratum] = take
        remaining -= take

    if remaining <= 0:
        return targets

    # Fill remaining slots roughly proportional to availability after the
    # minimum pass, so low-SNR controls are present but cannot consume the
    # entire queue before rare/high-value modes are included.
    while remaining > 0:
        candidates = {
            stratum: int(available[stratum]) - int(targets[stratum])
            for stratum in available.index
            if int(available[stratum]) > int(targets[stratum])
        }
        if not candidates:
            break
        total_left = sum(candidates.values())
        changed = False
        for stratum, count_left in candidates.items():
            if remaining <= 0:
                break
            take = max(1, int(round(remaining * count_left / total_left)))
            take = min(take, count_left, remaining)
            targets[stratum] += take
            remaining -= take
            changed = True
        if not changed:
            break
    return targets


def _select_rows(df: pd.DataFrame, n_rows: int, random_state: int, min_per_stratum: int, rare_cap: int) -> pd.DataFrame:
    rng = np.random.default_rng(int(random_state))
    targets = _target_counts(df["review_stratum"], n_rows, min_per_stratum, rare_cap)
    pieces: list[pd.DataFrame] = []
    for stratum, target in targets.items():
        if target <= 0:
            continue
        sub = df[df["review_stratum"].eq(stratum)].copy()
        if len(sub) > target:
            sub = sub.sample(n=target, random_state=int(rng.integers(0, 2**31 - 1)))
        pieces.append(sub)
    if not pieces:
        return df.head(0).copy()
    out = pd.concat(pieces, ignore_index=True, sort=False)
    out = out.sample(frac=1.0, random_state=random_state).reset_index(drop=True)
    return out


def _link_or_copy_reports(base_dir: Path, out_dir: Path, mode: str) -> str:
    src = base_dir / "vet_reports"
    dst = out_dir / "vet_reports"
    if dst.exists() or dst.is_symlink():
        return str(dst)
    if not src.exists():
        return ""
    if mode == "copy":
        shutil.copytree(src, dst)
    else:
        dst.symlink_to(src.resolve(), target_is_directory=True)
    return str(dst)


def build_queue(
    *,
    base_queue_csv: Path,
    topn_queue_csv: Path,
    survival_csv: Path,
    out_dir: Path,
    aperture: str,
    n_rows: int,
    random_state: int,
    min_per_stratum: int,
    rare_cap: int,
    reports_mode: str,
) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    base = pd.read_csv(base_queue_csv)
    topn = _read_topn(topn_queue_csv)
    survival = _read_survival(survival_csv, aperture)
    merged = base.merge(topn, on="injection_id", how="left")
    merged = merged.merge(survival, on="injection_id", how="left")
    merged["recovery_mode"] = merged.apply(_recovery_mode, axis=1)
    merged["review_stratum"] = merged.apply(_assign_stratum, axis=1)
    selected = _select_rows(merged, n_rows, random_state, min_per_stratum, rare_cap)

    selected.to_csv(out_dir / "review_queue.csv", index=False)
    selected["review_stratum"].value_counts().sort_index().rename_axis("review_stratum").reset_index(name="count").to_csv(
        out_dir / "stratum_counts.csv",
        index=False,
    )
    selected["recovery_mode"].value_counts().sort_index().rename_axis("recovery_mode").reset_index(name="count").to_csv(
        out_dir / "recovery_mode_counts.csv",
        index=False,
    )
    reports_path = _link_or_copy_reports(base_queue_csv.parent, out_dir, reports_mode)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "base_queue_csv": str(base_queue_csv),
        "topn_queue_csv": str(topn_queue_csv),
        "survival_csv": str(survival_csv),
        "aperture": aperture,
        "n_base_rows": int(len(base)),
        "n_selected_rows": int(len(selected)),
        "random_state": int(random_state),
        "min_per_stratum": int(min_per_stratum),
        "rare_cap": int(rare_cap),
        "reports_path": reports_path,
        "selected_recovery_modes": selected["recovery_mode"].value_counts().sort_index().to_dict(),
        "selected_strata": selected["review_stratum"].value_counts().sort_index().to_dict(),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    return selected


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-queue-csv", type=Path, required=True)
    parser.add_argument("--topn-queue-csv", type=Path, required=True)
    parser.add_argument("--survival-csv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--aperture", default="DET_FLUX_ADP")
    parser.add_argument("--n-rows", type=int, default=300)
    parser.add_argument("--random-state", type=int, default=56)
    parser.add_argument("--min-per-stratum", type=int, default=20)
    parser.add_argument("--rare-cap", type=int, default=80)
    parser.add_argument("--reports-mode", choices=("symlink", "copy", "none"), default="symlink")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    selected = build_queue(
        base_queue_csv=args.base_queue_csv,
        topn_queue_csv=args.topn_queue_csv,
        survival_csv=args.survival_csv,
        out_dir=args.out_dir,
        aperture=args.aperture,
        n_rows=args.n_rows,
        random_state=args.random_state,
        min_per_stratum=args.min_per_stratum,
        rare_cap=args.rare_cap,
        reports_mode=args.reports_mode,
    )
    print(f"[snr-queue] selected {len(selected):,} rows -> {args.out_dir / 'review_queue.csv'}")
    print(selected["review_stratum"].value_counts().sort_index().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
