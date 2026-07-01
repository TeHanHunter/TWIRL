#!/usr/bin/env python3
"""Build a peak-level BLS training table from injected light curves.

The existing injection-recovery review queues keep one summary row per
injection. That is the right format for recovery maps, but it is too collapsed
for training a search/ranking model: the ranker needs one row per BLS peak, with
an explicit truth label saying whether that peak corresponds to the injected
ephemeris.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, fields
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

from twirl.io.hlsp import BJDREFI, HLSPLightCurve  # noqa: E402
from twirl.search.bls import BLSConfig, run_bls_on_lc  # noqa: E402
from twirl.search.candidates import result_to_rows  # noqa: E402


DEFAULT_INJECTION_H5 = (
    REPO_ROOT
    / "data_local/stage3_injections/s56_twirlfs_v2_injection_training/"
    / "pdo_20k_predetrend_batman_periodradius_bright_balanced/injected_lightcurves.h5"
)
DEFAULT_OUT_TABLE = (
    REPO_ROOT
    / "reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/"
    / "peak_training/s56_20k_injection_bls_peaks.csv"
)
DEFAULT_APERTURES = ("DET_FLUX_ADP_SML", "DET_FLUX_SML")
DEFAULT_HARMONIC_FACTORS = (0.25, 1.0 / 3.0, 0.5, 1.0, 2.0, 3.0, 4.0)
PERIOD_RECOVERY_TOL = 0.02
MIN_WINDOW_OVERLAP_FRACTION = 0.50


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _as_float(value: Any, default: float = np.nan) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return float(default)
    return out if np.isfinite(out) else float(default)


def _phase_error_min(t0_bjd: float, truth_t0_bjd: float, period_d: float) -> float:
    if not all(np.isfinite([t0_bjd, truth_t0_bjd, period_d])) or period_d <= 0:
        return float("nan")
    delta_d = ((t0_bjd - truth_t0_bjd + 0.5 * period_d) % period_d) - 0.5 * period_d
    return float(delta_d * 1440.0)


def _nearest_harmonic_info(
    period_d: float,
    truth_period_d: float,
    *,
    harmonic_factors: tuple[float, ...] = DEFAULT_HARMONIC_FACTORS,
) -> dict[str, Any]:
    if not all(np.isfinite([period_d, truth_period_d])) or period_d <= 0 or truth_period_d <= 0:
        return {
            "period_ratio": float("nan"),
            "nearest_harmonic_factor": float("nan"),
            "harmonic_period_rel_err": float("nan"),
        }
    ratio = float(period_d / truth_period_d)
    factor = min(harmonic_factors, key=lambda value: abs(ratio - value) / value)
    return {
        "period_ratio": ratio,
        "nearest_harmonic_factor": float(factor),
        "harmonic_period_rel_err": float(abs(ratio - factor) / factor),
    }


def _transit_window_overlap_info(
    *,
    period_d: float,
    t0_bjd: float,
    duration_min: float,
    truth_period_d: float,
    truth_t0_bjd: float,
    truth_duration_min: float,
) -> dict[str, Any]:
    """Return phase-space overlap between the BLS and injected transit windows."""

    period = _as_float(period_d)
    truth_period = _as_float(truth_period_d)
    candidate_duration = _as_float(duration_min)
    truth_duration = _as_float(truth_duration_min)
    if not np.isfinite(candidate_duration) or candidate_duration <= 0:
        candidate_duration = truth_duration
    if not np.isfinite(truth_duration) or truth_duration <= 0:
        truth_duration = candidate_duration

    if (
        not np.isfinite(period)
        or period <= 0
        or not np.isfinite(truth_period)
        or truth_period <= 0
        or not np.isfinite(candidate_duration)
        or candidate_duration <= 0
        or not np.isfinite(truth_duration)
        or truth_duration <= 0
    ):
        return {
            "transit_window_phase_period_d": float("nan"),
            "transit_window_delta_min": float("nan"),
            "transit_window_overlap_min": float("nan"),
            "transit_window_overlap_fraction": float("nan"),
        }

    phase_period = min(period, truth_period)
    center_delta_min = _phase_error_min(t0_bjd, truth_t0_bjd, phase_period)
    if not np.isfinite(center_delta_min):
        overlap_min = float("nan")
        overlap_fraction = float("nan")
    else:
        half_width_sum = 0.5 * (candidate_duration + truth_duration)
        overlap_min = max(0.0, half_width_sum - abs(center_delta_min))
        overlap_fraction = overlap_min / min(candidate_duration, truth_duration)
        overlap_fraction = min(1.0, max(0.0, overlap_fraction))

    return {
        "transit_window_phase_period_d": float(phase_period),
        "transit_window_delta_min": float(center_delta_min),
        "transit_window_overlap_min": float(overlap_min),
        "transit_window_overlap_fraction": float(overlap_fraction),
    }


def label_peak_against_injection(
    *,
    period_d: float,
    t0_bjd: float,
    duration_min: float,
    truth_period_d: float,
    truth_t0_bjd: float,
    truth_duration_min: float,
    period_tol: float = PERIOD_RECOVERY_TOL,
    min_window_overlap_fraction: float = MIN_WINDOW_OVERLAP_FRACTION,
    harmonic_factors: tuple[float, ...] = DEFAULT_HARMONIC_FACTORS,
) -> dict[str, Any]:
    """Classify one BLS peak against injected truth.

    A positive label requires period or harmonic agreement and direct overlap
    between the predicted BLS transit window and injected transit window. This
    makes the training target match the ephemeris that downstream vetters need.
    """

    period = _as_float(period_d)
    t0 = _as_float(t0_bjd)
    truth_period = _as_float(truth_period_d)
    truth_t0 = _as_float(truth_t0_bjd)
    duration = _as_float(duration_min)
    truth_duration = _as_float(truth_duration_min, 20.0)
    t0_tol_min = max(float(truth_duration), 20.0) if np.isfinite(truth_duration) else 20.0
    window = _transit_window_overlap_info(
        period_d=period,
        t0_bjd=t0,
        duration_min=duration,
        truth_period_d=truth_period,
        truth_t0_bjd=truth_t0,
        truth_duration_min=truth_duration,
    )
    window_match = (
        np.isfinite(window["transit_window_overlap_fraction"])
        and float(window["transit_window_overlap_fraction"]) >= float(min_window_overlap_fraction)
    )

    if np.isfinite(period) and np.isfinite(truth_period) and truth_period > 0:
        period_rel_err = float(abs(period - truth_period) / truth_period)
    else:
        period_rel_err = float("nan")
    t0_phase_err_min = _phase_error_min(t0, truth_t0, truth_period)
    exact_ephemeris_match = (
        np.isfinite(period_rel_err)
        and period_rel_err <= period_tol
        and np.isfinite(t0_phase_err_min)
        and abs(t0_phase_err_min) <= t0_tol_min
        and window_match
    )

    harmonic = _nearest_harmonic_info(period, truth_period, harmonic_factors=harmonic_factors)
    harmonic_factor = _as_float(harmonic["nearest_harmonic_factor"])
    harmonic_t0_err = window["transit_window_delta_min"]
    is_exact_factor = np.isfinite(harmonic_factor) and abs(harmonic_factor - 1.0) < 1.0e-8
    harmonic_ephemeris_match = (
        not is_exact_factor
        and np.isfinite(harmonic["harmonic_period_rel_err"])
        and float(harmonic["harmonic_period_rel_err"]) <= period_tol
        and np.isfinite(harmonic_t0_err)
        and abs(harmonic_t0_err) <= t0_tol_min
        and window_match
    )
    if exact_ephemeris_match:
        match_kind = "exact"
    elif harmonic_ephemeris_match:
        match_kind = "harmonic"
    else:
        match_kind = "mismatch"

    return {
        "period_rel_err": period_rel_err,
        "t0_phase_err_min": t0_phase_err_min,
        "period_ratio": harmonic["period_ratio"],
        "nearest_harmonic_factor": harmonic["nearest_harmonic_factor"],
        "harmonic_period_rel_err": harmonic["harmonic_period_rel_err"],
        "harmonic_t0_phase_err_min": harmonic_t0_err,
        "exact_ephemeris_match": bool(exact_ephemeris_match),
        "harmonic_ephemeris_match": bool(harmonic_ephemeris_match),
        "is_injected_signal_peak": bool(exact_ephemeris_match or harmonic_ephemeris_match),
        "match_kind": match_kind,
        "t0_tolerance_min": float(t0_tol_min),
        "min_window_overlap_fraction": float(min_window_overlap_fraction),
        "transit_window_match": bool(window_match),
        **window,
    }


def _injection_apertures(attrs: dict[str, Any], requested: tuple[str, ...]) -> tuple[str, ...]:
    raw = attrs.get("apertures", "")
    if raw:
        try:
            available = tuple(str(item) for item in json.loads(str(raw)))
        except json.JSONDecodeError:
            available = tuple(part for part in str(raw).split("|") if part)
    else:
        available = (str(attrs.get("aperture", requested[0])),)
    selected = tuple(ap for ap in requested if ap in available)
    return selected or available


def _injection_lc_from_group(h5_path: Path, group_path: str, aperture: str) -> HLSPLightCurve:
    import h5py

    with h5py.File(h5_path, "r") as h5:
        group = h5[group_path]
        time = np.asarray(group["time"], dtype=np.float64)
        quality = np.asarray(group["quality"], dtype=np.int32)
        orbitid = np.asarray(group["orbitid"], dtype=np.int16)
        injected_name = f"{aperture}_injected"
        if injected_name in group:
            flux = np.asarray(group[injected_name], dtype=np.float64)
        elif "flux_injected" in group and str(group.attrs.get("aperture", "")) == aperture:
            flux = np.asarray(group["flux_injected"], dtype=np.float64)
        else:
            raise KeyError(f"missing injected aperture {aperture!r} in {h5_path}:{group_path}")
        tic = int(group.attrs["tic"])
        sector = int(group.attrs["sector"])
        cam = int(group.attrs["camera"])
        ccd = int(group.attrs["ccd"])
        tmag = float(group.attrs["tessmag"])
    payload = {
        "tic": tic,
        "tmag": tmag,
        "sector": sector,
        "cam": cam,
        "ccd": ccd,
        "ra": float("nan"),
        "dec": float("nan"),
        "time": time,
        "cadenceno": np.arange(len(time), dtype=np.int32),
        "orbitid": orbitid,
        "quality": quality,
        "flux": {aperture: flux},
        "path": Path(f"{h5_path}:{group_path}"),
    }
    accepted = {field.name for field in fields(HLSPLightCurve)}
    return HLSPLightCurve(**{key: value for key, value in payload.items() if key in accepted})


def _group_attrs(h5_path: Path, group_path: str) -> dict[str, Any]:
    import h5py

    with h5py.File(h5_path, "r") as h5:
        return dict(h5[group_path].attrs.items())


def _base_truth_record(attrs: dict[str, Any], group_path: str, h5_path: Path) -> dict[str, Any]:
    injection_id = str(attrs.get("injection_id", Path(group_path).name))
    return {
        "injection_id": injection_id,
        "source_h5": str(h5_path),
        "h5_group": group_path,
        "tic": int(attrs.get("tic", -1)),
        "sector": int(attrs.get("sector", -1)),
        "cam": int(attrs.get("camera", -1)),
        "ccd": int(attrs.get("ccd", -1)),
        "tmag": _as_float(attrs.get("tessmag")),
        "signal_family": str(attrs.get("signal_family", "")),
        "truth_period_d": _as_float(attrs.get("period_d")),
        "truth_t0_bjd": _as_float(attrs.get("t0_bjd")),
        "truth_duration_min": _as_float(attrs.get("duration_min")),
        "truth_depth": _as_float(attrs.get("depth")),
        "truth_model_depth": _as_float(attrs.get("model_depth", attrs.get("depth"))),
        "truth_sampled_model_depth": _as_float(attrs.get("sampled_model_depth", attrs.get("depth"))),
        "truth_radius_rearth": _as_float(attrs.get("radius_rearth")),
        "truth_radius_rwd": _as_float(attrs.get("radius_rwd")),
        "truth_impact_b": _as_float(attrs.get("impact_b")),
        "truth_a_over_rwd": _as_float(attrs.get("a_over_rwd")),
        "truth_inclination_deg": _as_float(attrs.get("inclination_deg")),
        "truth_n_in_transit": int(attrs.get("n_in_transit", -1)),
        "truth_n_good_in_transit": int(attrs.get("n_good_in_transit", -1)),
        "truth_injection_model": str(attrs.get("injection_model", "")),
        "truth_sampling_mode": str(attrs.get("sampling_mode", "")),
        "truth_grid_cell_id": str(attrs.get("grid_cell_id", "")),
        "truth_grid_period_bin": str(attrs.get("grid_period_bin", "")),
        "truth_grid_depth_bin": str(attrs.get("grid_depth_bin", "")),
        "truth_grid_radius_bin": str(attrs.get("grid_radius_bin", "")),
    }


def _branch_record(search_branch: str, cfg: BLSConfig) -> dict[str, Any]:
    return {
        "search_branch": str(search_branch),
        "bls_p_min_d": float(cfg.p_min_d),
        "bls_p_max_cap_d": float(cfg.p_max_cap_d),
        "bls_max_period_fraction": float(cfg.max_period_fraction),
        "bls_n_periods": int(cfg.n_periods),
        "bls_n_peaks": int(cfg.n_peaks),
        "bls_period_mask_frac": float(cfg.period_mask_frac),
        "bls_period_bin_edges": json.dumps(list(cfg.period_bin_edges)),
        "bls_max_peaks_per_period_bin": int(cfg.max_peaks_per_period_bin),
        "bls_sigma_clip": float(cfg.sigma_clip),
        "bls_min_cadences": int(cfg.min_cadences),
        "bls_durations_min": json.dumps(list(cfg.durations_min)),
    }


def peak_rows_for_injection(
    payload: tuple[str, str, tuple[str, ...], dict[str, Any], str, float, float, tuple[float, ...]]
) -> list[dict[str, Any]]:
    (
        h5_path_s,
        group_path,
        requested_apertures,
        cfg_kwargs,
        search_branch,
        period_tol,
        min_window_overlap_fraction,
        harmonic_factors,
    ) = payload
    h5_path = Path(h5_path_s)
    attrs = _group_attrs(h5_path, group_path)
    apertures = _injection_apertures(attrs, requested_apertures)
    truth = _base_truth_record(attrs, group_path, h5_path)
    cfg = BLSConfig(**cfg_kwargs)
    branch = _branch_record(search_branch, cfg)

    rows: list[dict[str, Any]] = []
    for aperture in apertures:
        try:
            lc = _injection_lc_from_group(h5_path, group_path, aperture)
            result = run_bls_on_lc(lc, cfg, aperture=aperture)
            peak_rows = [
                row
                for row in result_to_rows(result, run_id=f"s56-injection-peak-training:{search_branch}")
                if row.get("status") == "ok" and int(row.get("peak_rank", 0)) > 0
            ]
            if not peak_rows:
                rows.append(
                    {
                        **truth,
                        **branch,
                        "aperture": aperture,
                        "status": result.status,
                        "peak_rank": 0,
                        "is_candidate_peak": False,
                        "is_injected_signal_peak": False,
                        "match_kind": "no_peak",
                        "n_cad_total": result.n_cad_total,
                        "n_cad_quality": result.n_cad_quality,
                        "n_cad_kept": result.n_cad_kept,
                        "n_cad_edge_trimmed": result.n_cad_edge_trimmed,
                        "n_cad_sigma_clipped": result.n_cad_sigma_clipped,
                        "dropout_frac": result.dropout_frac,
                        "quality_dropout_frac": result.quality_dropout_frac,
                        "baseline_d": result.baseline_d,
                        "n_orbits": result.n_orbits,
                    }
                )
                continue
            for peak in peak_rows:
                period = _as_float(peak.get("period_d"))
                duration = _as_float(peak.get("duration_min"))
                label = label_peak_against_injection(
                    period_d=period,
                    t0_bjd=_as_float(peak.get("t0_bjd")),
                    duration_min=duration,
                    truth_period_d=truth["truth_period_d"],
                    truth_t0_bjd=truth["truth_t0_bjd"],
                    truth_duration_min=truth["truth_duration_min"],
                    period_tol=period_tol,
                    min_window_overlap_fraction=min_window_overlap_fraction,
                    harmonic_factors=harmonic_factors,
                )
                rows.append(
                    {
                        **truth,
                        **branch,
                        "aperture": aperture,
                        "status": result.status,
                        "is_candidate_peak": True,
                        "peak_rank": int(peak.get("peak_rank", 0)),
                        "period_d": period,
                        "t0_bjd": _as_float(peak.get("t0_bjd")),
                        "duration_min": duration,
                        "qtran": duration / 1440.0 / period if period > 0 and duration > 0 else np.nan,
                        "depth": _as_float(peak.get("depth")),
                        "depth_snr": _as_float(peak.get("depth_snr")),
                        "sde": _as_float(peak.get("sde")),
                        "log_power": _as_float(peak.get("log_power")),
                        "n_cad_total": result.n_cad_total,
                        "n_cad_quality": result.n_cad_quality,
                        "n_cad_kept": result.n_cad_kept,
                        "n_cad_edge_trimmed": result.n_cad_edge_trimmed,
                        "n_cad_sigma_clipped": result.n_cad_sigma_clipped,
                        "dropout_frac": result.dropout_frac,
                        "quality_dropout_frac": result.quality_dropout_frac,
                        "baseline_d": result.baseline_d,
                        "n_orbits": result.n_orbits,
                        **label,
                    }
                )
        except Exception as exc:
            rows.append(
                {
                    **truth,
                    **branch,
                    "aperture": aperture,
                    "status": f"error:{type(exc).__name__}: {exc}",
                    "peak_rank": 0,
                    "is_candidate_peak": False,
                    "is_injected_signal_peak": False,
                    "match_kind": "error",
                }
            )
    return rows


def _read_injection_keys(path: Path, n_injections: int, limit_keys: list[str] | None) -> list[str]:
    import h5py

    with h5py.File(path, "r") as h5:
        keys = sorted(h5["injections"].keys())
    if limit_keys is not None:
        wanted = set(limit_keys)
        keys = [key for key in keys if key in wanted]
    if n_injections > 0:
        keys = keys[:n_injections]
    return keys


def build_peak_training_table(
    *,
    injection_h5: Path,
    apertures: tuple[str, ...] = DEFAULT_APERTURES,
    n_injections: int = 0,
    workers: int = 1,
    n_periods: int = 200_000,
    n_peaks: int = 20,
    durations_min: tuple[float, ...] = BLSConfig.durations_min,
    p_min_d: float = BLSConfig.p_min_d,
    max_period_fraction: float = BLSConfig.max_period_fraction,
    p_max_cap_d: float = BLSConfig.p_max_cap_d,
    sigma_clip: float = BLSConfig.sigma_clip,
    min_cadences: int = BLSConfig.min_cadences,
    period_mask_frac: float = BLSConfig.period_mask_frac,
    period_bin_edges: tuple[float, ...] = BLSConfig.period_bin_edges,
    max_peaks_per_period_bin: int = BLSConfig.max_peaks_per_period_bin,
    period_tol: float = PERIOD_RECOVERY_TOL,
    min_window_overlap_fraction: float = MIN_WINDOW_OVERLAP_FRACTION,
    harmonic_factors: tuple[float, ...] = DEFAULT_HARMONIC_FACTORS,
    limit_keys: list[str] | None = None,
    search_branch: str = "standard",
) -> pd.DataFrame:
    injection_h5 = Path(injection_h5)
    if not injection_h5.exists():
        raise FileNotFoundError(f"missing injection HDF5: {injection_h5}")
    keys = _read_injection_keys(injection_h5, n_injections, limit_keys)
    cfg_kwargs = asdict(
        BLSConfig(
            apertures=apertures,
            n_periods=int(n_periods),
            n_peaks=int(n_peaks),
            durations_min=durations_min,
            p_min_d=float(p_min_d),
            max_period_fraction=float(max_period_fraction),
            p_max_cap_d=float(p_max_cap_d),
            sigma_clip=float(sigma_clip),
            min_cadences=int(min_cadences),
            period_mask_frac=float(period_mask_frac),
            period_bin_edges=tuple(float(value) for value in period_bin_edges),
            max_peaks_per_period_bin=int(max_peaks_per_period_bin),
        )
    )
    payloads = [
        (
            str(injection_h5),
            f"/injections/{key}",
            apertures,
            cfg_kwargs,
            str(search_branch),
            float(period_tol),
            float(min_window_overlap_fraction),
            harmonic_factors,
        )
        for key in keys
    ]
    print(
        f"[peak-training] running BLS for {len(payloads):,} injections, "
        f"{len(apertures)} requested apertures, top {n_peaks} peaks, branch={search_branch}",
        flush=True,
    )
    rows: list[dict[str, Any]] = []
    workers = max(1, int(workers))
    if workers <= 1:
        for idx, payload in enumerate(payloads, start=1):
            rows.extend(peak_rows_for_injection(payload))
            if idx % 50 == 0:
                print(f"  [peak-training] {idx:,}/{len(payloads):,}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for idx, batch in enumerate(ex.map(peak_rows_for_injection, payloads, chunksize=2), start=1):
                rows.extend(batch)
                if idx % 50 == 0:
                    print(f"  [peak-training] {idx:,}/{len(payloads):,}", flush=True)
    return pd.DataFrame(rows)


def summarize_peak_table(df: pd.DataFrame) -> dict[str, Any]:
    if df.empty:
        return {
            "n_peak_rows": 0,
            "n_injections": 0,
            "n_exact_match_injections": 0,
            "n_harmonic_match_injections": 0,
            "n_any_match_injections": 0,
        }
    peaks = df[df["is_candidate_peak"].fillna(False).astype(bool)].copy()
    grouped = df.groupby("injection_id", dropna=False)
    any_match = grouped["is_injected_signal_peak"].any()
    exact = grouped["exact_ephemeris_match"].any() if "exact_ephemeris_match" in df else pd.Series(dtype=bool)
    harmonic = grouped["harmonic_ephemeris_match"].any() if "harmonic_ephemeris_match" in df else pd.Series(dtype=bool)
    summary: dict[str, Any] = {
        "n_peak_rows": int(len(df)),
        "n_candidate_peak_rows": int(len(peaks)),
        "n_injections": int(df["injection_id"].nunique()),
        "n_any_match_injections": int(any_match.sum()),
        "any_match_fraction": float(any_match.mean()) if len(any_match) else float("nan"),
        "n_exact_match_injections": int(exact.sum()) if len(exact) else 0,
        "n_harmonic_match_injections": int(harmonic.sum()) if len(harmonic) else 0,
        "match_kind_counts": {
            str(key): int(value)
            for key, value in df["match_kind"].fillna("").astype(str).value_counts().sort_index().items()
        },
    }
    if "search_branch" in df:
        summary["search_branch_counts"] = {
            str(key): int(value)
            for key, value in df["search_branch"].fillna("").astype(str).value_counts().sort_index().items()
        }
    if len(peaks):
        matched = peaks[peaks["is_injected_signal_peak"].fillna(False).astype(bool)].copy()
        for k in (1, 2, 3, 5, 10, 20):
            eligible = grouped["peak_rank"].max() >= k
            within = (
                matched[matched["peak_rank"] <= k].groupby("injection_id")["is_injected_signal_peak"].any()
                if len(matched)
                else pd.Series(dtype=bool)
            )
            recovered = eligible.index.to_series().map(within).eq(True)
            denom = int(eligible.sum())
            summary[f"recall_at_{k}"] = {
                "n": int(recovered[eligible].sum()) if denom else 0,
                "denom": denom,
                "fraction": float(recovered[eligible].mean()) if denom else float("nan"),
            }
    return summary


def _write_table(df: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".parquet":
        df.to_parquet(path, index=False)
    else:
        df.to_csv(path, index=False)
    return path


def _parse_float_tuple(raw: str) -> tuple[float, ...]:
    return tuple(float(part.strip()) for part in str(raw).split(",") if part.strip())


def _parse_key_file(path: Path | None) -> list[str] | None:
    if path is None:
        return None
    return [line.strip() for line in Path(path).read_text().splitlines() if line.strip()]


def _build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--injection-h5", type=Path, default=DEFAULT_INJECTION_H5)
    ap.add_argument("--out-table", type=Path, default=DEFAULT_OUT_TABLE)
    ap.add_argument("--summary-json", type=Path, default=None)
    ap.add_argument("--apertures", nargs="+", default=list(DEFAULT_APERTURES))
    ap.add_argument("--n-injections", type=int, default=0, help="0 means all injections.")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--n-periods", type=int, default=200_000)
    ap.add_argument("--n-peaks", type=int, default=20)
    ap.add_argument("--durations-min", default=",".join(str(v) for v in BLSConfig.durations_min))
    ap.add_argument("--p-min-d", type=float, default=BLSConfig.p_min_d)
    ap.add_argument("--max-period-fraction", type=float, default=BLSConfig.max_period_fraction)
    ap.add_argument("--p-max-cap-d", type=float, default=BLSConfig.p_max_cap_d)
    ap.add_argument("--sigma-clip", type=float, default=BLSConfig.sigma_clip)
    ap.add_argument("--min-cadences", type=int, default=BLSConfig.min_cadences)
    ap.add_argument("--period-mask-frac", type=float, default=BLSConfig.period_mask_frac)
    ap.add_argument(
        "--period-bin-edges",
        default=",".join(str(v) for v in BLSConfig.period_bin_edges),
        help="Optional comma-separated period bin edges for peak-quota branches.",
    )
    ap.add_argument(
        "--max-peaks-per-period-bin",
        type=int,
        default=BLSConfig.max_peaks_per_period_bin,
        help="Optional maximum peaks from each period bin; 0 disables the quota.",
    )
    ap.add_argument("--period-tol", type=float, default=PERIOD_RECOVERY_TOL)
    ap.add_argument("--min-window-overlap-fraction", type=float, default=MIN_WINDOW_OVERLAP_FRACTION)
    ap.add_argument("--harmonic-factors", default=",".join(str(v) for v in DEFAULT_HARMONIC_FACTORS))
    ap.add_argument("--limit-keys-file", type=Path, default=None)
    ap.add_argument("--search-branch", default="standard")
    ap.add_argument("--overwrite", action="store_true")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    if args.out_table.exists() and not args.overwrite:
        print(f"[peak-training] output exists; pass --overwrite: {args.out_table}", file=sys.stderr)
        return 2
    summary_json = args.summary_json or args.out_table.with_name(args.out_table.stem + "_summary.json")
    df = build_peak_training_table(
        injection_h5=args.injection_h5,
        apertures=tuple(args.apertures),
        n_injections=args.n_injections,
        workers=args.workers,
        n_periods=args.n_periods,
        n_peaks=args.n_peaks,
        durations_min=_parse_float_tuple(args.durations_min),
        p_min_d=args.p_min_d,
        max_period_fraction=args.max_period_fraction,
        p_max_cap_d=args.p_max_cap_d,
        sigma_clip=args.sigma_clip,
        min_cadences=args.min_cadences,
        period_mask_frac=args.period_mask_frac,
        period_bin_edges=_parse_float_tuple(args.period_bin_edges),
        max_peaks_per_period_bin=args.max_peaks_per_period_bin,
        period_tol=args.period_tol,
        min_window_overlap_fraction=args.min_window_overlap_fraction,
        harmonic_factors=_parse_float_tuple(args.harmonic_factors),
        limit_keys=_parse_key_file(args.limit_keys_file),
        search_branch=args.search_branch,
    )
    _write_table(df, args.out_table)
    summary = summarize_peak_table(df)
    summary.update(
        {
            "injection_h5": str(args.injection_h5),
            "out_table": str(args.out_table),
            "apertures": list(args.apertures),
            "n_periods": int(args.n_periods),
            "n_peaks": int(args.n_peaks),
            "search_branch": str(args.search_branch),
            "p_min_d": float(args.p_min_d),
            "p_max_cap_d": float(args.p_max_cap_d),
            "max_period_fraction": float(args.max_period_fraction),
            "period_mask_frac": float(args.period_mask_frac),
            "period_bin_edges": list(_parse_float_tuple(args.period_bin_edges)),
            "max_peaks_per_period_bin": int(args.max_peaks_per_period_bin),
            "sigma_clip": float(args.sigma_clip),
            "min_cadences": int(args.min_cadences),
            "period_tol": float(args.period_tol),
            "min_window_overlap_fraction": float(args.min_window_overlap_fraction),
            "harmonic_factors": list(_parse_float_tuple(args.harmonic_factors)),
        }
    )
    summary_json.parent.mkdir(parents=True, exist_ok=True)
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    print("[peak-training] complete")
    print(f"  rows: {len(df):,}")
    print(f"  table: {args.out_table}")
    print(f"  summary: {summary_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
