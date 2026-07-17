#!/usr/bin/env python3
"""Audit stronger ADP+ detrending branches against injected-signal BLS recovery."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import fields
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any

import h5py
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
SRC_ROOT = REPO_ROOT / "src"
for path in (SCRIPT_DIR, SRC_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from build_s56_pretriage_review_queue import _read_table  # noqa: E402
from twirl.io.hlsp import HLSPLightCurve, read_hlsp  # noqa: E402
from twirl.search.bls import BLSConfig, run_bls_on_lc  # noqa: E402
from twirl.vetting.adpplus import (  # noqa: E402
    ADP_VET_APERTURES,
    ADPPlusBranch,
    binned_trend_ptp,
    branch_by_name,
    branched_light_curve,
)
from twirl.vetting.lightcurve_label_app import find_hlsp_path  # noqa: E402


DEFAULT_INJECTION_CSV = (
    REPO_ROOT
    / "reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/"
    / "small_pair_200k/review_queue.csv"
)
DEFAULT_INJECTION_H5 = (
    REPO_ROOT
    / "data_local/stage3_injections/s56_twirlfs_v2_injection_training/"
    / "pdo_20k_predetrend_batman_periodradius_bright_balanced/injected_lightcurves.h5"
)
DEFAULT_REAL_CANDIDATES = (
    REPO_ROOT / "data_local/stage2/bls_first_pass_v2/sector_0056/vetted_per_tic_centroid.csv"
)
DEFAULT_HLSP_ROOT = REPO_ROOT / "data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare"
DEFAULT_OUT_DIR = REPO_ROOT / "reports/stage5_validation/s56_adpplus_bls_audit"
HARMONIC_FACTORS = (0.25, 1.0 / 3.0, 0.5, 1.0, 2.0, 3.0, 4.0)


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _finite_numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        return pd.Series(np.nan, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce").replace([np.inf, -np.inf], np.nan)


def _sample_balanced(frame: pd.DataFrame, n: int, columns: tuple[str, ...], seed: int) -> pd.DataFrame:
    if n <= 0 or frame.empty:
        return frame.head(0).copy()
    if len(frame) <= n:
        return frame.sample(frac=1.0, random_state=seed).reset_index(drop=True)
    work = frame.copy()
    for column in columns:
        if column not in work:
            work[column] = "missing"
        work[column] = work[column].astype(object).where(work[column].notna(), "missing").astype(str)
    groups = list(work.groupby(list(columns), dropna=False, sort=True))
    per_cell = max(1, int(np.ceil(n / max(len(groups), 1))))
    pieces: list[pd.DataFrame] = []
    for idx, (_, group) in enumerate(groups):
        pieces.append(group.sample(n=min(per_cell, len(group)), random_state=seed + idx))
    selected = pd.concat(pieces, ignore_index=False)
    if len(selected) > n:
        selected = selected.sample(n=n, random_state=seed + 10_000)
    elif len(selected) < n:
        remaining = work.loc[[i for i in work.index if i not in set(selected.index)]]
        if not remaining.empty:
            fill = remaining.sample(n=min(n - len(selected), len(remaining)), random_state=seed + 20_000)
            selected = pd.concat([selected, fill], ignore_index=False)
    return selected.sample(frac=1.0, random_state=seed + 30_000).reset_index(drop=True)


def select_injection_sample(path: Path, *, n: int, seed: int) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "source_kind" in df:
        df = df[df["source_kind"].fillna("").astype(str).eq("injection_recovery")].copy()
    for col in ("tmag", "truth_period_d", "truth_radius_rearth"):
        df[col] = _finite_numeric(df, col)
    df = df.dropna(subset=["injection_id", "h5_group", "tmag", "truth_period_d", "truth_radius_rearth"]).copy()
    df["tmag_bin"] = pd.cut(df["tmag"], [-np.inf, 17, 18, 19, np.inf], labels=["lt17", "17_18", "18_19", "gt19"])
    df["period_bin"] = pd.cut(
        df["truth_period_d"],
        [0.0, 0.12, 0.2, 0.35, 0.6, 1.0, 1.8, 3.2, 5.6, 10.0, np.inf],
        labels=["p000_012", "p012_020", "p020_035", "p035_060", "p060_100", "p100_180", "p180_320", "p320_560", "p560_1000", "p1000p"],
    )
    df["radius_bin"] = pd.cut(
        df["truth_radius_rearth"],
        [0.0, 0.5, 1, 2, 4, 8, 12, np.inf],
        labels=["r000_05", "r05_1", "r1_2", "r2_4", "r4_8", "r8_12", "r12p"],
    )
    if "topn_recovery_status" not in df:
        df["topn_recovery_status"] = df.get("recovery_status", "missing")
    return _sample_balanced(df, n, ("tmag_bin", "period_bin", "radius_bin", "topn_recovery_status"), seed)


def select_real_sample(path: Path, *, n: int, seed: int) -> pd.DataFrame:
    df = _read_table(path).copy()
    df["tic"] = _finite_numeric(df, "tic").astype("Int64")
    df = df.dropna(subset=["tic"]).copy()
    for col in ("period_d", "duration_min", "sde_max", "tmag", "n_apertures_agree"):
        df[col] = _finite_numeric(df, col)
    vet = df.get("vet_class", pd.Series("", index=df.index)).fillna("").astype(str)
    sde = df["sde_max"]
    duration = df["duration_min"]
    agree = df["n_apertures_agree"]
    q40 = float(sde.quantile(0.40)) if sde.notna().any() else 8.0
    q80_dur = float(duration.quantile(0.80)) if duration.notna().any() else 30.0
    buckets = pd.Series("real_flat_control", index=df.index, dtype=object)
    buckets.loc[vet.eq("planet_candidate") | sde.ge(sde.quantile(0.85)).fillna(False)] = "real_high_sde"
    buckets.loc[vet.isin(["pceb_grid_ceiling", "sub_roche_pceb_suspect"])] = "real_eb_pceb"
    buckets.loc[duration.ge(max(15.0, q80_dur)).fillna(False)] = "real_broad_trend"
    buckets.loc[agree.le(1).fillna(False)] = "real_aperture_disagreement"
    buckets.loc[sde.lt(q40).fillna(False)] = "real_low_sde_control"
    df["real_audit_bucket"] = buckets
    return _sample_balanced(df, n, ("real_audit_bucket",), seed)


def _peak_records(result: Any) -> list[dict[str, Any]]:
    if not getattr(result, "peaks", None):
        return [
            {
                "peak_rank": 0,
                "period_d": np.nan,
                "t0_bjd": np.nan,
                "duration_min": np.nan,
                "depth": np.nan,
                "depth_snr": np.nan,
                "sde": np.nan,
                "log_power": np.nan,
            }
        ]
    return [
        {
            "peak_rank": int(pk.peak_rank),
            "period_d": float(pk.period_d),
            "t0_bjd": float(pk.t0_bjd),
            "duration_min": float(pk.duration_min),
            "depth": float(pk.depth),
            "depth_snr": float(pk.depth_snr),
            "sde": float(pk.sde),
            "log_power": float(pk.log_power),
        }
        for pk in result.peaks
    ]


def _phase_error_min(t0_bjd: float, truth_t0_bjd: float, period_d: float) -> float:
    if not all(np.isfinite([t0_bjd, truth_t0_bjd, period_d])) or period_d <= 0:
        return float("nan")
    delta_d = ((t0_bjd - truth_t0_bjd + 0.5 * period_d) % period_d) - 0.5 * period_d
    return float(delta_d * 1440.0)


def _label_peak_against_truth(
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
    t0_tol_min = max(float(truth_duration_min), 20.0) if np.isfinite(truth_duration_min) else 20.0
    exact = (
        np.isfinite(period_rel_err)
        and period_rel_err <= period_tol
        and np.isfinite(t0_phase_err_min)
        and abs(t0_phase_err_min) <= t0_tol_min
    )
    is_exact_factor = np.isfinite(harmonic_factor) and abs(harmonic_factor - 1.0) < 1.0e-8
    harmonic = (
        not is_exact_factor
        and np.isfinite(harmonic_period_rel_err)
        and harmonic_period_rel_err <= period_tol
        and np.isfinite(harmonic_t0_phase_err_min)
        and abs(harmonic_t0_phase_err_min) <= t0_tol_min
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
        "t0_tolerance_min": float(t0_tol_min),
    }


def _align_flux_to_time(source_time: np.ndarray, source_flux: np.ndarray, target_time: np.ndarray) -> np.ndarray:
    if len(source_time) == len(target_time) and np.allclose(source_time, target_time, rtol=0, atol=1.0e-8, equal_nan=True):
        return np.asarray(source_flux, dtype=np.float64)
    good = np.isfinite(source_time) & np.isfinite(source_flux)
    if np.count_nonzero(good) < 2:
        raise ValueError("cannot align fallback flux with fewer than two finite cadences")
    order = np.argsort(source_time[good], kind="stable")
    return np.interp(
        np.asarray(target_time, dtype=np.float64),
        np.asarray(source_time[good], dtype=np.float64)[order],
        np.asarray(source_flux[good], dtype=np.float64)[order],
        left=np.nan,
        right=np.nan,
    )


def _fallback_hlsp_flux(hlsp_root: Path | None, *, tic: int, sector: int, aperture: str, target_time: np.ndarray) -> np.ndarray | None:
    if hlsp_root is None:
        return None
    path = find_hlsp_path(hlsp_root, tic, sector)
    if path is None:
        return None
    lc = read_hlsp(path, columns=(aperture,))
    if lc is None or aperture not in lc.flux:
        return None
    return _align_flux_to_time(lc.time, lc.flux[aperture], target_time)


def _h5_lc(path: Path, group_path: str, aperture: str, *, injected: bool, hlsp_root: Path | None = None) -> HLSPLightCurve:
    with h5py.File(path, "r") as h5:
        group = h5[group_path]
        time = np.asarray(group["time"], dtype=np.float64)
        quality = np.asarray(group["quality"], dtype=np.int32)
        orbitid = np.asarray(group["orbitid"], dtype=np.int16)
        tic = int(group.attrs["tic"])
        sector = int(group.attrs["sector"])
        dataset = f"{aperture}_{'injected' if injected else 'original'}"
        if dataset not in group:
            dataset = f"{aperture}_injected" if injected else aperture
        if dataset in group:
            flux = np.asarray(group[dataset], dtype=np.float64)
        else:
            original_name = f"{aperture}_original"
            if original_name in group:
                base_flux = np.asarray(group[original_name], dtype=np.float64)
            elif aperture in group:
                base_flux = np.asarray(group[aperture], dtype=np.float64)
            else:
                base_flux = _fallback_hlsp_flux(hlsp_root, tic=tic, sector=sector, aperture=aperture, target_time=time)
            if base_flux is None:
                raise KeyError(f"missing {dataset} in {path}:{group_path}")
            if injected:
                if "transit_model" not in group:
                    raise KeyError(f"missing transit_model for fallback injection in {path}:{group_path}")
                flux = np.asarray(base_flux, dtype=np.float64) * np.asarray(group["transit_model"], dtype=np.float64)
            else:
                flux = np.asarray(base_flux, dtype=np.float64)
        payload = {
            "tic": tic,
            "tmag": float(group.attrs["tessmag"]),
            "sector": sector,
            "cam": int(group.attrs["camera"]),
            "ccd": int(group.attrs["ccd"]),
            "ra": float("nan"),
            "dec": float("nan"),
            "time": time,
            "cadenceno": np.arange(len(time), dtype=np.int32),
            "orbitid": orbitid,
            "quality": quality,
            "flux": {aperture: flux},
            "path": Path(f"{path}:{group_path}"),
        }
    accepted = {field.name for field in fields(HLSPLightCurve)}
    return HLSPLightCurve(**{key: value for key, value in payload.items() if key in accepted})


def _truth_from_row(row: dict[str, Any]) -> dict[str, float]:
    return {
        "truth_period_d": float(row.get("truth_period_d", np.nan)),
        "truth_t0_bjd": float(row.get("truth_t0_bjd", np.nan)),
        "truth_duration_min": float(row.get("truth_duration_min", np.nan)),
        "truth_depth": float(row.get("truth_depth", np.nan)),
        "truth_radius_rearth": float(row.get("truth_radius_rearth", np.nan)),
        "truth_impact_b": float(row.get("truth_impact_b", np.nan)),
        "truth_inclination_deg": float(row.get("truth_inclination_deg", np.nan)),
    }


def _run_branch(
    lc: HLSPLightCurve,
    *,
    aperture: str,
    branch: ADPPlusBranch,
    mask_mode: str,
    cfg: BLSConfig,
    truth: dict[str, float] | None,
) -> tuple[Any, dict[str, Any], HLSPLightCurve]:
    if branch.is_identity:
        branch_lc, meta = branched_light_curve(lc, aperture, branch)
        return run_bls_on_lc(branch_lc, cfg, aperture=aperture), dict(meta), branch_lc
    if mask_mode == "truth_masked" and truth is not None:
        period = truth["truth_period_d"]
        t0 = truth["truth_t0_bjd"]
        duration = truth["truth_duration_min"]
    elif mask_mode == "iterative":
        current_lc, _ = branched_light_curve(lc, aperture, branch_by_name("current_adp"))
        current_res = run_bls_on_lc(current_lc, cfg, aperture=aperture)
        if getattr(current_res, "peaks", None):
            peak = current_res.peaks[0]
            period = float(peak.period_d)
            t0 = float(peak.t0_bjd)
            duration = float(peak.duration_min)
        else:
            period = t0 = duration = np.nan
    else:
        period = t0 = duration = np.nan
    branch_lc, meta = branched_light_curve(
        lc,
        aperture,
        branch,
        period_d=period,
        t0_bjd=t0,
        duration_min=duration,
    )
    return run_bls_on_lc(branch_lc, cfg, aperture=aperture), dict(meta), branch_lc


def _depth_retention(
    original_lc: HLSPLightCurve | None,
    injected_lc: HLSPLightCurve,
    aperture: str,
    truth: dict[str, float],
) -> float:
    if original_lc is None:
        return float("nan")
    period = truth["truth_period_d"]
    t0 = truth["truth_t0_bjd"]
    duration = truth["truth_duration_min"]
    if not all(np.isfinite([period, t0, duration])):
        return float("nan")
    from twirl.vetting.adpplus import transit_window_mask

    in_tr = transit_window_mask(
        injected_lc.time,
        period_d=period,
        t0_bjd=t0,
        duration_min=duration,
        width_factor=1.0,
        min_width_min=duration,
    )
    good = (injected_lc.quality == 0) & np.isfinite(injected_lc.flux[aperture]) & np.isfinite(original_lc.flux[aperture])
    if np.count_nonzero(good & in_tr) < 1 or np.count_nonzero(good & ~in_tr) < 20:
        return float("nan")
    delta = injected_lc.flux[aperture] - original_lc.flux[aperture]
    depth = -(float(np.nanmedian(delta[good & in_tr])) - float(np.nanmedian(delta[good & ~in_tr])))
    truth_depth = truth.get("truth_depth", np.nan)
    return depth / truth_depth if np.isfinite(truth_depth) and truth_depth > 0 else float("nan")


def _process_injection(payload: tuple[dict[str, Any], str, str, list[str], list[str], int, int]) -> list[dict[str, Any]]:
    row, h5_path_s, hlsp_root_s, branch_names, mask_modes, n_periods, n_peaks = payload
    h5_path = Path(h5_path_s)
    hlsp_root = Path(hlsp_root_s) if hlsp_root_s else None
    truth = _truth_from_row(row)
    cfg = BLSConfig(apertures=ADP_VET_APERTURES, n_periods=n_periods, n_peaks=n_peaks)
    records: list[dict[str, Any]] = []
    for aperture in ADP_VET_APERTURES:
        try:
            inj_lc = _h5_lc(h5_path, str(row["h5_group"]), aperture, injected=True, hlsp_root=hlsp_root)
            try:
                orig_base = _h5_lc(h5_path, str(row["h5_group"]), aperture, injected=False, hlsp_root=hlsp_root)
            except Exception:
                orig_base = None
        except Exception as exc:
            records.append({"source_kind": "injection_recovery", "injection_id": row.get("injection_id"), "aperture": aperture, "status": f"read_error:{exc}"})
            continue
        for branch_name in branch_names:
            branch = branch_by_name(branch_name)
            modes = ["current"] if branch.is_identity else mask_modes
            for mode in modes:
                try:
                    res, meta, branch_lc = _run_branch(inj_lc, aperture=aperture, branch=branch, mask_mode=mode, cfg=cfg, truth=truth)
                    orig_branch_lc = None
                    if orig_base is not None:
                        try:
                            if branch.is_identity:
                                orig_branch_lc, _ = branched_light_curve(orig_base, aperture, branch)
                            elif mode == "truth_masked":
                                orig_branch_lc, _ = branched_light_curve(
                                    orig_base,
                                    aperture,
                                    branch,
                                    period_d=truth["truth_period_d"],
                                    t0_bjd=truth["truth_t0_bjd"],
                                    duration_min=truth["truth_duration_min"],
                                )
                        except Exception:
                            orig_branch_lc = None
                    retention = _depth_retention(orig_branch_lc, branch_lc, aperture, truth)
                    for peak in _peak_records(res):
                        label = (
                            _label_peak_against_truth(
                                period_d=peak["period_d"],
                                t0_bjd=peak["t0_bjd"],
                                duration_min=peak["duration_min"],
                                truth_period_d=truth["truth_period_d"],
                                truth_t0_bjd=truth["truth_t0_bjd"],
                                truth_duration_min=truth["truth_duration_min"],
                            )
                            if peak["peak_rank"] > 0
                            else {"is_injected_signal_peak": False, "exact_ephemeris_match": False, "harmonic_ephemeris_match": False, "match_kind": "no_peak"}
                        )
                        records.append(
                            {
                                "source_kind": "injection_recovery",
                                "injection_id": row.get("injection_id"),
                                "tic": row.get("tic"),
                                "tmag": row.get("tmag"),
                                "aperture": aperture,
                                "branch": branch.name,
                                "mask_mode": mode,
                                "status": res.status,
                                "branch_status": meta.get("status", ""),
                                "branch_window_d": meta.get("window_d", np.nan),
                                "depth_retention_frac": retention,
                                "trend_ptp": binned_trend_ptp(branch_lc.time, branch_lc.flux[aperture], branch_lc.quality),
                                **truth,
                                **peak,
                                **label,
                            }
                        )
                except Exception as exc:
                    records.append({"source_kind": "injection_recovery", "injection_id": row.get("injection_id"), "aperture": aperture, "branch": branch.name, "mask_mode": mode, "status": f"error:{type(exc).__name__}: {exc}"})
    return records


def _process_real(payload: tuple[dict[str, Any], str, list[str], int, int]) -> list[dict[str, Any]]:
    row, hlsp_root_s, branch_names, n_periods, n_peaks = payload
    hlsp_root = Path(hlsp_root_s)
    tic = int(float(row["tic"]))
    sector = int(float(row["sector"])) if row.get("sector") not in (None, "") and pd.notna(row.get("sector")) else None
    path = find_hlsp_path(hlsp_root, tic, sector)
    if path is None:
        return [{"source_kind": "real_candidate", "tic": tic, "status": "missing_hlsp"}]
    lc = read_hlsp(path, columns=ADP_VET_APERTURES)
    if lc is None:
        return [{"source_kind": "real_candidate", "tic": tic, "status": "read_fail"}]
    cfg = BLSConfig(apertures=ADP_VET_APERTURES, n_periods=n_periods, n_peaks=n_peaks)
    records: list[dict[str, Any]] = []
    for aperture in ADP_VET_APERTURES:
        if aperture not in lc.flux:
            records.append({"source_kind": "real_candidate", "tic": tic, "aperture": aperture, "status": "missing_aperture"})
            continue
        for branch_name in branch_names:
            branch = branch_by_name(branch_name)
            mode = "current" if branch.is_identity else "iterative"
            try:
                res, meta, branch_lc = _run_branch(lc, aperture=aperture, branch=branch, mask_mode=mode, cfg=cfg, truth=None)
                for peak in _peak_records(res):
                    records.append(
                        {
                            "source_kind": "real_candidate",
                            "tic": tic,
                            "sector": row.get("sector", ""),
                            "tmag": row.get("tmag", np.nan),
                            "real_audit_bucket": row.get("real_audit_bucket", ""),
                            "aperture": aperture,
                            "branch": branch.name,
                            "mask_mode": mode,
                            "status": res.status,
                            "branch_status": meta.get("status", ""),
                            "branch_window_d": meta.get("window_d", np.nan),
                            "trend_ptp": binned_trend_ptp(branch_lc.time, branch_lc.flux[aperture], branch_lc.quality),
                            **peak,
                        }
                    )
            except Exception as exc:
                records.append({"source_kind": "real_candidate", "tic": tic, "aperture": aperture, "branch": branch.name, "mask_mode": mode, "status": f"error:{type(exc).__name__}: {exc}"})
    return records


def _run_payloads(payloads: list[Any], fn, workers: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    workers = max(1, int(workers))
    if workers <= 1:
        for idx, payload in enumerate(payloads, start=1):
            rows.extend(fn(payload))
            if idx % 25 == 0:
                print(f"[adpplus-audit] processed {idx:,}/{len(payloads):,}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for idx, batch in enumerate(ex.map(fn, payloads, chunksize=1), start=1):
                rows.extend(batch)
                if idx % 25 == 0:
                    print(f"[adpplus-audit] processed {idx:,}/{len(payloads):,}", flush=True)
    return rows


def summarize(inj: pd.DataFrame, real: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    columns = [
        "source_kind",
        "branch",
        "mask_mode",
        "aperture",
        "n_objects",
        "top1_exact_n",
        "top1_exact_frac",
        "topn_exact_n",
        "topn_exact_frac",
        "topn_exact_or_harmonic_n",
        "topn_exact_or_harmonic_frac",
        "median_depth_retention",
        "median_trend_ptp",
        "median_sde_rank1",
    ]
    if not inj.empty and "peak_rank" in inj:
        candidate = inj[inj["peak_rank"].fillna(0).astype(int) > 0].copy()
        for keys, group in candidate.groupby(["branch", "mask_mode", "aperture"], dropna=False):
            branch, mode, aperture = keys
            top1 = group[group["peak_rank"].astype(int).eq(1)]
            by_inj = group.groupby("injection_id", dropna=False)
            exact_any = by_inj["exact_ephemeris_match"].any()
            harmonic_any = by_inj["harmonic_ephemeris_match"].any()
            rows.append(
                {
                    "source_kind": "injection_recovery",
                    "branch": branch,
                    "mask_mode": mode,
                    "aperture": aperture,
                    "n_objects": int(group["injection_id"].nunique()),
                    "top1_exact_n": int(top1["exact_ephemeris_match"].fillna(False).astype(bool).sum()),
                    "top1_exact_frac": float(top1["exact_ephemeris_match"].fillna(False).astype(bool).mean()) if len(top1) else np.nan,
                    "topn_exact_n": int(exact_any.sum()),
                    "topn_exact_frac": float(exact_any.mean()) if len(exact_any) else np.nan,
                    "topn_exact_or_harmonic_n": int((exact_any | harmonic_any).sum()),
                    "topn_exact_or_harmonic_frac": float((exact_any | harmonic_any).mean()) if len(exact_any) else np.nan,
                    "median_depth_retention": float(np.nanmedian(pd.to_numeric(group["depth_retention_frac"], errors="coerce"))),
                    "median_trend_ptp": float(np.nanmedian(pd.to_numeric(group["trend_ptp"], errors="coerce"))),
                    "median_sde_rank1": float(np.nanmedian(pd.to_numeric(top1["sde"], errors="coerce"))) if len(top1) else np.nan,
                }
            )
    if not real.empty and "peak_rank" in real:
        candidate = real[real["peak_rank"].fillna(0).astype(int) > 0].copy()
        for keys, group in candidate.groupby(["branch", "mask_mode", "aperture"], dropna=False):
            branch, mode, aperture = keys
            top1 = group[group["peak_rank"].astype(int).eq(1)]
            rows.append(
                {
                    "source_kind": "real_candidate",
                    "branch": branch,
                    "mask_mode": mode,
                    "aperture": aperture,
                    "n_objects": int(group["tic"].nunique()),
                    "top1_exact_n": np.nan,
                    "top1_exact_frac": np.nan,
                    "topn_exact_n": np.nan,
                    "topn_exact_frac": np.nan,
                    "topn_exact_or_harmonic_n": np.nan,
                    "topn_exact_or_harmonic_frac": np.nan,
                    "median_depth_retention": np.nan,
                    "median_trend_ptp": float(np.nanmedian(pd.to_numeric(group["trend_ptp"], errors="coerce"))),
                    "median_sde_rank1": float(np.nanmedian(pd.to_numeric(top1["sde"], errors="coerce"))) if len(top1) else np.nan,
                }
            )
    out = pd.DataFrame(rows, columns=columns)
    if out.empty:
        return out
    return out.sort_values(["source_kind", "aperture", "branch", "mask_mode"])


def run_audit(args: argparse.Namespace) -> dict[str, Any]:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    branches = [part.strip() for part in args.branches.split(",") if part.strip()]
    for branch in branches:
        branch_by_name(branch)
    mask_modes = [part.strip() for part in args.mask_modes.split(",") if part.strip()]
    injected_sample = select_injection_sample(args.injection_csv, n=args.n_injected, seed=args.random_state)
    real_sample = select_real_sample(args.real_candidates, n=args.n_real, seed=args.random_state + 1)
    injected_sample.to_csv(args.out_dir / "audit_sample_injected.csv", index=False)
    real_sample.to_csv(args.out_dir / "audit_sample_real.csv", index=False)

    inj_payloads = [
        (row, str(args.injection_h5), str(args.hlsp_root), branches, mask_modes, int(args.n_periods), int(args.n_peaks))
        for row in injected_sample.to_dict("records")
    ]
    real_payloads = [
        (row, str(args.hlsp_root), branches, int(args.n_periods), int(args.n_peaks))
        for row in real_sample.to_dict("records")
    ]
    print(f"[adpplus-audit] injections={len(inj_payloads):,} real={len(real_payloads):,} branches={branches}", flush=True)
    injection_rows = _run_payloads(inj_payloads, _process_injection, args.workers) if inj_payloads else []
    real_rows = _run_payloads(real_payloads, _process_real, args.workers) if real_payloads else []
    inj_df = pd.DataFrame(injection_rows)
    real_df = pd.DataFrame(real_rows)
    inj_df.to_csv(args.out_dir / "injection_branch_peak_table.csv", index=False)
    real_df.to_csv(args.out_dir / "real_branch_peak_table.csv", index=False)
    summary_df = summarize(inj_df, real_df)
    summary_df.to_csv(args.out_dir / "branch_summary.csv", index=False)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "out_dir": str(args.out_dir),
        "injection_csv": str(args.injection_csv),
        "injection_h5": str(args.injection_h5),
        "real_candidates": str(args.real_candidates),
        "hlsp_root": str(args.hlsp_root),
        "n_injected_sample": int(len(injected_sample)),
        "n_real_sample": int(len(real_sample)),
        "branches": branches,
        "mask_modes": mask_modes,
        "apertures": list(ADP_VET_APERTURES),
        "n_periods": int(args.n_periods),
        "n_peaks": int(args.n_peaks),
        "outputs": {
            "audit_sample_injected": str(args.out_dir / "audit_sample_injected.csv"),
            "audit_sample_real": str(args.out_dir / "audit_sample_real.csv"),
            "injection_branch_peak_table": str(args.out_dir / "injection_branch_peak_table.csv"),
            "real_branch_peak_table": str(args.out_dir / "real_branch_peak_table.csv"),
            "branch_summary": str(args.out_dir / "branch_summary.csv"),
            "summary_json": str(args.out_dir / "summary.json"),
        },
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--injection-csv", type=Path, default=DEFAULT_INJECTION_CSV)
    ap.add_argument("--injection-h5", type=Path, default=DEFAULT_INJECTION_H5)
    ap.add_argument("--real-candidates", type=Path, default=DEFAULT_REAL_CANDIDATES)
    ap.add_argument("--hlsp-root", type=Path, default=DEFAULT_HLSP_ROOT)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    ap.add_argument("--n-injected", type=int, default=3000)
    ap.add_argument("--n-real", type=int, default=2000)
    ap.add_argument("--branches", default="current_adp,adpplus_0p20,adpplus_0p12,adpplus_adaptive")
    ap.add_argument("--mask-modes", default="truth_masked,iterative")
    ap.add_argument("--n-periods", type=int, default=50_000)
    ap.add_argument("--n-peaks", type=int, default=20)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--random-state", type=int, default=5611)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    summary = run_audit(args)
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
