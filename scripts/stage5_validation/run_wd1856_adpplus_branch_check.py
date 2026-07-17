#!/usr/bin/env python3
"""Check WD 1856 recovery for ADP/ADP+ two-aperture branches."""
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
SCRIPT_DIR = Path(__file__).resolve().parent
SRC_ROOT = REPO_ROOT / "src"
for path in (SCRIPT_DIR, SRC_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from twirl.io.hlsp import read_hlsp  # noqa: E402
from twirl.search.bls import BLSConfig, run_bls_on_lc  # noqa: E402
from twirl.vetting.adpplus import ADP_VET_APERTURES, branch_by_name, branched_light_curve  # noqa: E402
from twirl.vetting.lightcurve_label_app import find_hlsp_path  # noqa: E402


DEFAULT_OUT_DIR = Path("reports/stage5_validation/wd1856_adpplus_branch_check")
DEFAULT_HLSP_ROOT = Path("data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare")
DEFAULT_BRANCHES = "current_adp,adpplus_0p20,adpplus_0p12,adpplus_adaptive"
WD1856_TIC = 267574918
WD1856_PERIOD_D = 1.407939211
WD1856_T0_BJD = 2458779.375083
WD1856_DURATION_MIN = 8.0
HARMONIC_FACTORS = (0.25, 1.0 / 3.0, 0.5, 1.0, 2.0, 3.0, 4.0)


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _best_ephemeris_for_mask(lc, aperture: str, cfg: BLSConfig) -> tuple[float, float, float, str]:
    current = branch_by_name("current_adp")
    current_lc, _ = branched_light_curve(lc, aperture, current)
    res = run_bls_on_lc(current_lc, cfg, aperture=aperture)
    if not getattr(res, "peaks", None):
        return float("nan"), float("nan"), float("nan"), "no_current_peak"
    peak = res.peaks[0]
    return float(peak.period_d), float(peak.t0_bjd), float(peak.duration_min), "current_peak"


def _phase_error_min(t0_bjd: float, truth_t0_bjd: float, period_d: float) -> float:
    if not all(np.isfinite([t0_bjd, truth_t0_bjd, period_d])) or period_d <= 0:
        return float("nan")
    delta_d = ((t0_bjd - truth_t0_bjd + 0.5 * period_d) % period_d) - 0.5 * period_d
    return float(delta_d * 1440.0)


def _label_against_ephemeris(
    *,
    period_d: float,
    t0_bjd: float,
    duration_min: float,
    truth_period_d: float,
    truth_t0_bjd: float,
    truth_duration_min: float,
    period_tol: float = 0.02,
) -> dict[str, Any]:
    if np.isfinite(period_d) and np.isfinite(truth_period_d) and truth_period_d > 0:
        period_rel_err = float(abs(period_d - truth_period_d) / truth_period_d)
        period_ratio = float(period_d / truth_period_d)
        harmonic_factor = min(HARMONIC_FACTORS, key=lambda value: abs(period_ratio - value) / value)
        harmonic_period_rel_err = float(abs(period_ratio - harmonic_factor) / harmonic_factor)
    else:
        period_rel_err = float("nan")
        period_ratio = float("nan")
        harmonic_factor = float("nan")
        harmonic_period_rel_err = float("nan")
    t0_phase_err_min = _phase_error_min(t0_bjd, truth_t0_bjd, truth_period_d)
    phase_period = min(period_d, truth_period_d) if np.isfinite(period_d) and np.isfinite(truth_period_d) else float("nan")
    harmonic_t0_phase_err_min = _phase_error_min(t0_bjd, truth_t0_bjd, phase_period)
    tolerance_min = max(float(truth_duration_min), 20.0)
    exact = (
        np.isfinite(period_rel_err)
        and period_rel_err <= period_tol
        and np.isfinite(t0_phase_err_min)
        and abs(t0_phase_err_min) <= tolerance_min
    )
    is_exact_factor = np.isfinite(harmonic_factor) and abs(harmonic_factor - 1.0) < 1.0e-8
    harmonic = (
        not is_exact_factor
        and np.isfinite(harmonic_period_rel_err)
        and harmonic_period_rel_err <= period_tol
        and np.isfinite(harmonic_t0_phase_err_min)
        and abs(harmonic_t0_phase_err_min) <= tolerance_min
    )
    match_kind = "exact" if exact else ("harmonic" if harmonic else "mismatch")
    return {
        "period_rel_err": period_rel_err,
        "t0_phase_err_min": t0_phase_err_min,
        "period_ratio": period_ratio,
        "nearest_harmonic_factor": float(harmonic_factor),
        "harmonic_period_rel_err": harmonic_period_rel_err,
        "harmonic_t0_phase_err_min": harmonic_t0_phase_err_min,
        "exact_ephemeris_match": bool(exact),
        "harmonic_ephemeris_match": bool(harmonic),
        "is_injected_signal_peak": bool(exact or harmonic),
        "match_kind": match_kind,
        "t0_tolerance_min": tolerance_min,
    }


def _peak_rows(result, *, context: dict[str, Any], truth_period_d: float, truth_t0_bjd: float, truth_duration_min: float) -> list[dict[str, Any]]:
    if not getattr(result, "peaks", None):
        return [
            {
                **context,
                "status": result.status,
                "peak_rank": 0,
                "period_d": np.nan,
                "t0_bjd": np.nan,
                "duration_min": np.nan,
                "depth": np.nan,
                "depth_snr": np.nan,
                "sde": np.nan,
                "log_power": np.nan,
                "exact_ephemeris_match": False,
                "harmonic_ephemeris_match": False,
                "match_kind": "no_peak",
            }
        ]
    rows: list[dict[str, Any]] = []
    for peak in result.peaks:
        label = _label_against_ephemeris(
            period_d=float(peak.period_d),
            t0_bjd=float(peak.t0_bjd),
            duration_min=float(peak.duration_min),
            truth_period_d=float(truth_period_d),
            truth_t0_bjd=float(truth_t0_bjd),
            truth_duration_min=float(truth_duration_min),
        )
        rows.append(
            {
                **context,
                "status": result.status,
                "peak_rank": int(peak.peak_rank),
                "period_d": float(peak.period_d),
                "t0_bjd": float(peak.t0_bjd),
                "duration_min": float(peak.duration_min),
                "depth": float(peak.depth),
                "depth_snr": float(peak.depth_snr),
                "sde": float(peak.sde),
                "log_power": float(peak.log_power),
                **label,
            }
        )
    return rows


def run_check(args: argparse.Namespace) -> dict[str, Any]:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    path = find_hlsp_path(args.hlsp_root, args.tic, args.sector)
    if path is None:
        raise FileNotFoundError(f"no HLSP found for TIC {args.tic} under {args.hlsp_root}")
    lc = read_hlsp(path, columns=ADP_VET_APERTURES)
    if lc is None:
        raise RuntimeError(f"failed to read {path}")
    cfg = BLSConfig(apertures=ADP_VET_APERTURES, n_periods=args.n_periods, n_peaks=args.n_peaks)
    branches = [part.strip() for part in args.branches.split(",") if part.strip()]
    rows: list[dict[str, Any]] = []
    for aperture in ADP_VET_APERTURES:
        if aperture not in lc.flux:
            rows.append({"tic": args.tic, "aperture": aperture, "status": "missing_aperture"})
            continue
        for branch_name in branches:
            branch = branch_by_name(branch_name)
            modes = ["current"] if branch.is_identity else ["truth_masked", "iterative"]
            for mode in modes:
                if mode == "truth_masked":
                    period = args.truth_period_d
                    t0 = args.truth_t0_bjd
                    duration = args.truth_duration_min
                    mask_source = "truth"
                elif mode == "iterative":
                    period, t0, duration, mask_source = _best_ephemeris_for_mask(lc, aperture, cfg)
                else:
                    period = t0 = duration = np.nan
                    mask_source = "none"
                branch_lc, meta = branched_light_curve(
                    lc,
                    aperture,
                    branch,
                    period_d=period,
                    t0_bjd=t0,
                    duration_min=duration,
                )
                result = run_bls_on_lc(branch_lc, cfg, aperture=aperture)
                context = {
                    "tic": args.tic,
                    "sector": lc.sector,
                    "tmag": lc.tmag,
                    "aperture": aperture,
                    "branch": branch.name,
                    "mask_mode": mode,
                    "mask_source": mask_source,
                    "branch_status": meta.get("status", ""),
                    "branch_window_d": meta.get("window_d", np.nan),
                    "truth_period_d": args.truth_period_d,
                    "truth_t0_bjd": args.truth_t0_bjd,
                    "truth_duration_min": args.truth_duration_min,
                    "hlsp_path": str(path),
                }
                rows.extend(
                    _peak_rows(
                        result,
                        context=context,
                        truth_period_d=args.truth_period_d,
                        truth_t0_bjd=args.truth_t0_bjd,
                        truth_duration_min=args.truth_duration_min,
                    )
                )
    peaks = pd.DataFrame(rows)
    peaks_path = args.out_dir / "wd1856_branch_peaks.csv"
    peaks.to_csv(peaks_path, index=False)
    top1 = peaks[peaks.get("peak_rank", pd.Series(dtype=int)).fillna(0).astype(int).eq(1)].copy()
    summary_cols = [
        "aperture",
        "branch",
        "mask_mode",
        "status",
        "period_d",
        "t0_bjd",
        "duration_min",
        "sde",
        "depth",
        "period_rel_err",
        "t0_phase_err_min",
        "exact_ephemeris_match",
        "harmonic_ephemeris_match",
        "match_kind",
    ]
    summary_table = top1[[col for col in summary_cols if col in top1]].copy()
    summary_path = args.out_dir / "wd1856_top1_summary.csv"
    summary_table.to_csv(summary_path, index=False)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "hlsp_root": str(args.hlsp_root),
        "hlsp_path": str(path),
        "tic": int(args.tic),
        "sector": int(lc.sector),
        "truth_period_d": float(args.truth_period_d),
        "truth_t0_bjd": float(args.truth_t0_bjd),
        "truth_duration_min": float(args.truth_duration_min),
        "n_periods": int(args.n_periods),
        "n_peaks": int(args.n_peaks),
        "branches": branches,
        "apertures": list(ADP_VET_APERTURES),
        "outputs": {
            "peaks": str(peaks_path),
            "top1_summary": str(summary_path),
        },
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hlsp-root", type=Path, default=DEFAULT_HLSP_ROOT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--tic", type=int, default=WD1856_TIC)
    parser.add_argument("--sector", type=int, default=56)
    parser.add_argument("--truth-period-d", type=float, default=WD1856_PERIOD_D)
    parser.add_argument("--truth-t0-bjd", type=float, default=WD1856_T0_BJD)
    parser.add_argument("--truth-duration-min", type=float, default=WD1856_DURATION_MIN)
    parser.add_argument("--branches", default=DEFAULT_BRANCHES)
    parser.add_argument("--n-periods", type=int, default=50_000)
    parser.add_argument("--n-peaks", type=int, default=20)
    args = parser.parse_args(argv)
    summary = run_check(args)
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
