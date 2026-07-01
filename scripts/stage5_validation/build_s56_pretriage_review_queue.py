#!/usr/bin/env python3
"""Build the S56 LEO-assisted review queue before human triage.

The review unit is a recovered candidate, not a raw light curve:

1. select a stratified real-candidate sample from the S56 BLS/vetter table,
2. run BLS on LC-level injected examples,
3. render WD-tuned LEO-Vetter PDF reports for every review row when possible,
4. write a browser-ready candidate CSV for ``run_lightcurve_vetting_app.py``.
"""
from __future__ import annotations

import argparse
import contextlib
from dataclasses import asdict, fields
from datetime import datetime, timezone
import io
import importlib.util
import json
import math
import os
from pathlib import Path
import signal
import sys
import tempfile
import warnings
from typing import Any

import numpy as np
import pandas as pd
from astropy.io import fits

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.io.hlsp import (  # noqa: E402
    BJDREFI,
    HLSPLightCurve,
    quality_mask,
    read_hlsp,
    tglc_mad_error,
)
from twirl.search.bls import BLSConfig, run_bls_on_lc  # noqa: E402
from twirl.search.candidates import result_to_rows  # noqa: E402
from twirl.vetting.lightcurve_label_app import find_hlsp_path  # noqa: E402

WD_1856_TIC = 267574918
DEFAULT_REAL_CANDIDATES = (
    REPO_ROOT / "data_local/stage2/bls_first_pass_v2/sector_0056/vetted_per_tic_centroid.parquet"
)
DEFAULT_HLSP_ROOT = (
    REPO_ROOT / "data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare"
)
DEFAULT_INJECTION_H5 = (
    REPO_ROOT
    / "data_local/stage3_injections/s56_twirlfs_v2_injection_training/injections_900/injected_lightcurves.h5"
)
DEFAULT_OUT_DIR = REPO_ROOT / "reports/stage5_validation/s56_pretriage_review_queue"
DEFAULT_STAR_CATALOG = (
    REPO_ROOT / "data_local/catalogs/twirl_master_catalog/twirl_wd_master_catalog_v0_ticmatched.fits"
)

LABEL_COLUMNS = ("label", "label_source", "labeler", "notes")
KNOWN_APERTURES = (
    "DET_FLUX_SML",
    "DET_FLUX",
    "DET_FLUX_LAG",
    "DET_FLUX_ADP",
    "DET_FLUX_ADP_SML",
    "DET_FLUX_ADP_LAG",
)
PERIOD_RECOVERY_TOL = 0.02
HARMONIC_PERIOD_FACTORS = (
    0.25,
    1.0 / 3.0,
    0.5,
    1.0,
    2.0,
    3.0,
    4.0,
    5.0,
    6.0,
    8.0,
    10.0,
    12.0,
)


def _has_parquet_engine() -> bool:
    return (
        importlib.util.find_spec("pyarrow") is not None
        or importlib.util.find_spec("fastparquet") is not None
    )


def _read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        try:
            return pd.read_parquet(path)
        except ImportError:
            csv_path = path.with_suffix(".csv")
            if csv_path.exists():
                return pd.read_csv(csv_path)
            raise
    if suffix == ".csv":
        return pd.read_csv(path)
    raise ValueError(f"unsupported table format: {path}")


def _write_table(df: pd.DataFrame, path: Path) -> Path:
    if path.suffix.lower() == ".parquet" and _has_parquet_engine():
        df.to_parquet(path, compression="zstd")
        return path
    csv_path = path.with_suffix(".csv")
    df.to_csv(csv_path, index=False)
    return csv_path


def _existing_table_path(path: Path) -> Path | None:
    if path.exists():
        return path
    csv_path = path.with_suffix(".csv")
    if csv_path.exists():
        return csv_path
    return None
REVIEW_COLUMNS = (
    "review_id",
    "source_kind",
    "source_bucket",
    "tic",
    "sector",
    "cam",
    "ccd",
    "tmag",
    "vet_class",
    "class_rank",
    "blind_rank",
    "period_d",
    "t0_bjd",
    "duration_min",
    "depth",
    "depth_snr",
    "sde_max",
    "rep_aperture",
    "n_apertures_agree",
    "apertures_agree",
    "centroid_status",
    "centroid_pass",
    "centroid_delta_pix",
    "centroid_z",
    "n_in_transit",
    "n_oot_band",
    "recovery_status",
    "injection_id",
    "signal_family",
    "truth_period_d",
    "truth_t0_bjd",
    "truth_duration_min",
    "truth_depth",
    "truth_target_depth",
    "truth_geometric_depth",
    "truth_model_depth",
    "truth_sampled_model_depth",
    "truth_radius_rearth",
    "truth_radius_rwd",
    "truth_impact_b",
    "truth_a_over_rwd",
    "truth_inclination_deg",
    "truth_duration_model",
    "truth_injection_model",
    "truth_sampling_mode",
    "truth_grid_cell_id",
    "truth_grid_period_bin",
    "truth_grid_depth_bin",
    "truth_grid_radius_bin",
    "truth_source_kind",
    "truth_source_bucket",
    "period_rel_err",
    "t0_phase_err_min",
    "topn_recovery_status",
    "n_topn_apertures_agree",
    "topn_apertures_agree",
    "n_harmonic_apertures_agree",
    "harmonic_apertures_agree",
    "star_source",
    "star_source_id",
    "star_atmosphere",
    "star_rad_rsun",
    "star_e_rad_rsun",
    "star_mass_msun",
    "star_e_mass_msun",
    "star_teff_k",
    "star_e_teff_k",
    "star_logg_cgs",
    "star_e_logg_cgs",
    "star_rho_g_cm3",
    "star_radius_source",
    "leo_class",
    "leo_report_name",
)


def _clean(value: Any) -> Any:
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass
    if isinstance(value, np.generic):
        return value.item()
    return value


def _phase_error_min(t0_bjd: float, truth_t0_bjd: float, period_d: float) -> float:
    if not all(np.isfinite([t0_bjd, truth_t0_bjd, period_d])) or period_d <= 0:
        return float("nan")
    delta_d = ((t0_bjd - truth_t0_bjd + 0.5 * period_d) % period_d) - 0.5 * period_d
    return float(delta_d * 1440.0)


def _period_ratio_info(period_d: float, truth_period_d: float) -> dict[str, Any]:
    if not all(np.isfinite([period_d, truth_period_d])) or period_d <= 0 or truth_period_d <= 0:
        return {
            "period_ratio": float("nan"),
            "harmonic_factor": float("nan"),
            "harmonic_rel_err": float("nan"),
            "harmonic_period_match": False,
        }
    ratio = float(period_d / truth_period_d)
    factor = min(HARMONIC_PERIOD_FACTORS, key=lambda value: abs(ratio - value) / value)
    rel_err = float(abs(ratio - factor) / factor)
    return {
        "period_ratio": ratio,
        "harmonic_factor": float(factor),
        "harmonic_rel_err": rel_err,
        "harmonic_period_match": bool(rel_err <= PERIOD_RECOVERY_TOL),
    }


def _harmonic_phase_error_min(t0_bjd: float, truth_t0_bjd: float, period_d: float, truth_period_d: float) -> float:
    if not all(np.isfinite([period_d, truth_period_d])) or period_d <= 0 or truth_period_d <= 0:
        return float("nan")
    return _phase_error_min(t0_bjd, truth_t0_bjd, min(period_d, truth_period_d))


def _topn_truth_match_summary(
    peak_rows: list[dict[str, Any]],
    *,
    truth_period_d: float,
    truth_t0_bjd: float,
    truth_duration_min: float,
) -> dict[str, Any]:
    """Summarize whether BLS found the truth ephemeris in top-N peaks.

    ``recovery_status`` remains a strict top-peak quantity. This helper records
    the extra information needed to separate true non-detections from BLS
    ranking/harmonic failures in injected samples.
    """
    t0_tol_min = max(float(truth_duration_min), 20.0) if np.isfinite(truth_duration_min) else 20.0
    empty = {
        "top_period_ratio": float("nan"),
        "top_harmonic_factor": float("nan"),
        "top_harmonic_rel_err": float("nan"),
        "top_harmonic_t0_phase_err_min": float("nan"),
        "top_harmonic_period_match": False,
        "top_harmonic_ephemeris_match": False,
        "topn_exact_recovered": False,
        "topn_exact_peak_rank": float("nan"),
        "topn_exact_sde": float("nan"),
        "topn_exact_period_rel_err": float("nan"),
        "topn_exact_t0_phase_err_min": float("nan"),
        "topn_harmonic_match": False,
        "topn_harmonic_peak_rank": float("nan"),
        "topn_harmonic_sde": float("nan"),
        "topn_harmonic_period_ratio": float("nan"),
        "topn_harmonic_factor": float("nan"),
        "topn_harmonic_rel_err": float("nan"),
        "topn_harmonic_t0_phase_err_min": float("nan"),
    }
    if not peak_rows:
        return empty

    ordered = sorted(
        peak_rows,
        key=lambda row: int(row.get("peak_rank", 10**9)) if str(row.get("peak_rank", "")).strip() else 10**9,
    )
    top = ordered[0]
    top_period = float(top.get("period_d", np.nan))
    top_t0 = float(top.get("t0_bjd", np.nan))
    top_info = _period_ratio_info(top_period, truth_period_d)
    top_harm_t0 = _harmonic_phase_error_min(top_t0, truth_t0_bjd, top_period, truth_period_d)
    out = {
        **empty,
        "top_period_ratio": top_info["period_ratio"],
        "top_harmonic_factor": top_info["harmonic_factor"],
        "top_harmonic_rel_err": top_info["harmonic_rel_err"],
        "top_harmonic_t0_phase_err_min": top_harm_t0,
        "top_harmonic_period_match": bool(top_info["harmonic_period_match"]),
        "top_harmonic_ephemeris_match": bool(
            top_info["harmonic_period_match"] and np.isfinite(top_harm_t0) and abs(top_harm_t0) <= t0_tol_min
        ),
    }

    exact_matches: list[tuple[dict[str, Any], float, float]] = []
    harmonic_matches: list[tuple[dict[str, Any], dict[str, Any], float]] = []
    for row in ordered:
        period = float(row.get("period_d", np.nan))
        t0 = float(row.get("t0_bjd", np.nan))
        period_rel_err = (
            abs(period - truth_period_d) / truth_period_d
            if np.isfinite(period) and np.isfinite(truth_period_d) and truth_period_d > 0
            else float("nan")
        )
        t0_phase_err_min = _phase_error_min(t0, truth_t0_bjd, truth_period_d)
        if np.isfinite(period_rel_err) and period_rel_err <= PERIOD_RECOVERY_TOL and abs(t0_phase_err_min) <= t0_tol_min:
            exact_matches.append((row, period_rel_err, t0_phase_err_min))

        info = _period_ratio_info(period, truth_period_d)
        harmonic_t0_err = _harmonic_phase_error_min(t0, truth_t0_bjd, period, truth_period_d)
        is_exact_factor = np.isfinite(info["harmonic_factor"]) and abs(float(info["harmonic_factor"]) - 1.0) < 1.0e-8
        if (
            bool(info["harmonic_period_match"])
            and not is_exact_factor
            and np.isfinite(harmonic_t0_err)
            and abs(harmonic_t0_err) <= t0_tol_min
        ):
            harmonic_matches.append((row, info, harmonic_t0_err))

    if exact_matches:
        row, period_rel_err, t0_phase_err_min = exact_matches[0]
        out.update(
            {
                "topn_exact_recovered": True,
                "topn_exact_peak_rank": row.get("peak_rank", ""),
                "topn_exact_sde": float(row.get("sde", np.nan)),
                "topn_exact_period_rel_err": period_rel_err,
                "topn_exact_t0_phase_err_min": t0_phase_err_min,
            }
        )
    if harmonic_matches:
        row, info, t0_phase_err_min = harmonic_matches[0]
        out.update(
            {
                "topn_harmonic_match": True,
                "topn_harmonic_peak_rank": row.get("peak_rank", ""),
                "topn_harmonic_sde": float(row.get("sde", np.nan)),
                "topn_harmonic_period_ratio": info["period_ratio"],
                "topn_harmonic_factor": info["harmonic_factor"],
                "topn_harmonic_rel_err": info["harmonic_rel_err"],
                "topn_harmonic_t0_phase_err_min": t0_phase_err_min,
            }
        )
    return out


def _finite_int_attr(attrs: dict[str, Any], key: str) -> int | None:
    value = attrs.get(key, None)
    try:
        out = int(value)
    except (TypeError, ValueError):
        return None
    return out if out >= 0 else None


def _injection_transit_counts(group: Any, attrs: dict[str, Any]) -> tuple[int, int]:
    """Return total and QUALITY=0 in-transit cadence counts for an injection.

    Early pre-detrend HDF5 products stored the boolean ``in_transit`` dataset
    but not the count attributes. Keep the review table backward-compatible so
    completeness diagnostics do not silently report zeros for those files.
    """
    n_in = _finite_int_attr(attrs, "n_in_transit")
    n_good = _finite_int_attr(attrs, "n_good_in_transit")
    if n_in is not None and n_good is not None:
        return n_in, n_good
    if "in_transit" not in group:
        fallback = n_good if n_good is not None else (n_in if n_in is not None else 0)
        return n_in if n_in is not None else fallback, fallback
    in_transit = np.asarray(group["in_transit"], dtype=bool)
    computed_in = int(np.count_nonzero(in_transit))
    if "quality" in group:
        quality = np.asarray(group["quality"])
        computed_good = int(np.count_nonzero(in_transit & (quality == 0)))
    else:
        computed_good = computed_in
    return n_in if n_in is not None else computed_in, n_good if n_good is not None else computed_good


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


def _injection_apertures_from_attrs(attrs: dict[str, Any], requested: tuple[str, ...]) -> tuple[str, ...]:
    raw = attrs.get("apertures", "")
    available: tuple[str, ...]
    if raw:
        try:
            decoded = json.loads(str(raw))
            available = tuple(str(v) for v in decoded)
        except json.JSONDecodeError:
            available = tuple(str(raw).split("|"))
    else:
        available = (str(attrs.get("aperture", requested[0])),)
    selected = tuple(ap for ap in requested if ap in available)
    return selected or available


def _best_injection_peak(
    *,
    h5_path: Path,
    group_path: str,
    apertures: tuple[str, ...],
    cfg_kwargs: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    per_aperture: dict[str, dict[str, Any]] = {}
    recovery_by_aperture: dict[str, str] = {}
    best_row: dict[str, Any] | None = None
    cfg = BLSConfig(apertures=apertures, **cfg_kwargs)

    import h5py

    with h5py.File(h5_path, "r") as h5:
        group = h5[group_path]
        attrs = dict(group.attrs.items())
        n_in_transit, n_good_in_transit = _injection_transit_counts(group, attrs)

    truth_period = float(attrs.get("period_d", np.nan))
    truth_t0 = float(attrs.get("t0_bjd", np.nan))
    truth_duration = float(attrs.get("duration_min", np.nan))

    for aperture in apertures:
        try:
            lc = _injection_lc_from_group(h5_path, group_path, aperture)
            result = run_bls_on_lc(lc, cfg, aperture=aperture)
            rows = result_to_rows(result, run_id="s56-injection-pretriage")
            peak_rows = [r for r in rows if r.get("status") == "ok" and int(r.get("peak_rank", 0)) > 0]
            row = max(peak_rows, key=lambda r: float(r.get("sde", -np.inf))) if peak_rows else rows[0]
        except Exception as exc:
            peak_rows = []
            row = {
                "status": f"error:{type(exc).__name__}: {exc}",
                "aperture": aperture,
                "sde": float("nan"),
                "peak_rank": "",
            }

        period = float(row.get("period_d", np.nan))
        t0_bjd = float(row.get("t0_bjd", np.nan))
        period_rel_err = (
            abs(period - truth_period) / truth_period
            if np.isfinite(period) and np.isfinite(truth_period) and truth_period > 0
            else float("nan")
        )
        t0_phase_err_min = _phase_error_min(t0_bjd, truth_t0, truth_period)
        recovered = (
            row.get("status") == "ok"
            and np.isfinite(period)
            and period_rel_err <= PERIOD_RECOVERY_TOL
            and abs(t0_phase_err_min) <= max(truth_duration, 20.0)
        )
        truth_match = _topn_truth_match_summary(
            peak_rows,
            truth_period_d=truth_period,
            truth_t0_bjd=truth_t0,
            truth_duration_min=truth_duration,
        )
        recovery_by_aperture[aperture] = "bls_recovered" if recovered else str(row.get("status", "failed"))
        per_aperture[aperture] = {
            "row": row,
            "period_rel_err": period_rel_err,
            "t0_phase_err_min": t0_phase_err_min,
            "recovered": recovered,
            "topn_recovered": bool(truth_match["topn_exact_recovered"]),
            "harmonic_matched": bool(truth_match["topn_harmonic_match"]),
            "truth_match": truth_match,
        }
        if best_row is None or float(row.get("sde", -np.inf)) > float(best_row.get("sde", -np.inf)):
            best_row = row

    if best_row is None:
        best_row = {"status": "no_apertures", "aperture": apertures[0], "sde": float("nan")}
    return attrs, per_aperture, best_row


def _run_one_injection_bls(payload: tuple[str, str, tuple[str, ...], dict[str, Any]]) -> dict[str, Any]:
    h5_path_s, group_path, requested_apertures, cfg_kwargs = payload
    h5_path = Path(h5_path_s)
    import h5py

    with h5py.File(h5_path, "r") as h5:
        group = h5[group_path]
        attrs = dict(group.attrs.items())
        n_in_transit, n_good_in_transit = _injection_transit_counts(group, attrs)

    apertures = _injection_apertures_from_attrs(attrs, requested_apertures)
    attrs, per_aperture, row = _best_injection_peak(
        h5_path=h5_path,
        group_path=group_path,
        apertures=apertures,
        cfg_kwargs=cfg_kwargs,
    )
    recovered_apertures = [ap for ap, rec in per_aperture.items() if bool(rec["recovered"])]
    if recovered_apertures:
        row = max(
            (per_aperture[ap]["row"] for ap in recovered_apertures),
            key=lambda r: float(r.get("sde", -np.inf)),
        )
    aperture = str(row.get("aperture", apertures[0]))

    period = float(row.get("period_d", np.nan))
    truth_period = float(attrs.get("period_d", np.nan))
    truth_t0 = float(attrs.get("t0_bjd", np.nan))
    t0_bjd = float(row.get("t0_bjd", np.nan))
    period_rel_err = (
        abs(period - truth_period) / truth_period
        if np.isfinite(period) and np.isfinite(truth_period) and truth_period > 0
        else float("nan")
    )
    t0_phase_err_min = _phase_error_min(t0_bjd, truth_t0, truth_period)
    recovery_status = "bls_recovered" if recovered_apertures else "bls_peak_mismatch"
    if row.get("status") != "ok" or not np.isfinite(period):
        recovery_status = f"bls_{row.get('status', 'failed')}"
    topn_recovered_apertures = [ap for ap, rec in per_aperture.items() if bool(rec.get("topn_recovered"))]
    harmonic_match_apertures = [ap for ap, rec in per_aperture.items() if bool(rec.get("harmonic_matched"))]
    if recovered_apertures:
        topn_recovery_status = "bls_top1_recovered"
    elif topn_recovered_apertures:
        topn_recovery_status = "bls_topn_recovered"
    elif harmonic_match_apertures:
        topn_recovery_status = "bls_topn_harmonic_match"
    else:
        topn_recovery_status = recovery_status

    injection_id = str(attrs.get("injection_id", Path(group_path).name))
    out = {
        "review_id": f"inj:{injection_id}",
        "source_kind": "injection_recovery",
        "source_bucket": f"injection_{attrs.get('signal_family', 'unknown')}_{recovery_status}",
        "tic": int(attrs.get("tic", row.get("tic", -1))),
        "sector": int(attrs.get("sector", row.get("sector", 56))),
        "cam": int(attrs.get("camera", row.get("cam", -1))),
        "ccd": int(attrs.get("ccd", row.get("ccd", -1))),
        "tmag": float(attrs.get("tessmag", row.get("tmag", np.nan))),
        "vet_class": f"injected_{attrs.get('signal_family', 'unknown')}",
        "class_rank": "",
        "blind_rank": "",
        "period_d": period,
        "t0_bjd": t0_bjd,
        "duration_min": float(row.get("duration_min", np.nan)),
        "depth": float(row.get("depth", np.nan)),
        "depth_snr": float(row.get("depth_snr", np.nan)),
        "sde_max": float(row.get("sde", np.nan)),
        "rep_aperture": aperture,
        "n_apertures_agree": len(recovered_apertures),
        "apertures_agree": ",".join(recovered_apertures),
        "centroid_status": "not_applicable_injection",
        "centroid_pass": "",
        "centroid_delta_pix": "",
        "centroid_z": "",
        "n_in_transit": n_good_in_transit,
        "n_oot_band": "",
        "recovery_status": recovery_status,
        "injection_id": injection_id,
        "signal_family": str(attrs.get("signal_family", "")),
        "truth_period_d": truth_period,
        "truth_t0_bjd": truth_t0,
        "truth_duration_min": float(attrs.get("duration_min", np.nan)),
        "truth_depth": float(attrs.get("depth", np.nan)),
        "truth_target_depth": _clean(attrs.get("target_depth", "")),
        "truth_geometric_depth": _clean(attrs.get("geometric_depth", "")),
        "truth_model_depth": _clean(attrs.get("model_depth", attrs.get("depth", ""))),
        "truth_sampled_model_depth": _clean(attrs.get("sampled_model_depth", "")),
        "truth_radius_rearth": _clean(attrs.get("radius_rearth", "")),
        "truth_radius_rwd": _clean(attrs.get("radius_rwd", "")),
        "truth_impact_b": _clean(attrs.get("impact_b", "")),
        "truth_a_over_rwd": _clean(attrs.get("a_over_rwd", "")),
        "truth_inclination_deg": _clean(attrs.get("inclination_deg", "")),
        "truth_duration_model": _clean(attrs.get("duration_model", "")),
        "truth_injection_model": _clean(attrs.get("injection_model", "")),
        "truth_sampling_mode": _clean(attrs.get("sampling_mode", "")),
        "truth_grid_cell_id": _clean(attrs.get("grid_cell_id", "")),
        "truth_grid_period_bin": _clean(attrs.get("grid_period_bin", "")),
        "truth_grid_depth_bin": _clean(attrs.get("grid_depth_bin", "")),
        "truth_grid_radius_bin": _clean(attrs.get("grid_radius_bin", "")),
        "truth_n_in_transit": n_in_transit,
        "truth_n_good_in_transit": n_good_in_transit,
        "truth_source_kind": "injection_recovery",
        "truth_source_bucket": f"injection_{attrs.get('signal_family', 'unknown')}_{recovery_status}",
        "period_rel_err": period_rel_err,
        "t0_phase_err_min": t0_phase_err_min,
        "topn_recovery_status": topn_recovery_status,
        "n_topn_apertures_agree": len(topn_recovered_apertures),
        "topn_apertures_agree": ",".join(topn_recovered_apertures),
        "n_harmonic_apertures_agree": len(harmonic_match_apertures),
        "harmonic_apertures_agree": ",".join(harmonic_match_apertures),
        "source_h5": str(h5_path),
        "h5_group": group_path,
        "bls_status": str(row.get("status", "")),
        "leo_class": "",
        "leo_report_name": "",
    }
    for ap, rec in per_aperture.items():
        ap_row = rec["row"]
        out[f"sde_{ap}"] = float(ap_row.get("sde", np.nan))
        out[f"rank_{ap}"] = ap_row.get("peak_rank", "")
        out[f"period_rel_err_{ap}"] = rec["period_rel_err"]
        out[f"t0_phase_err_min_{ap}"] = rec["t0_phase_err_min"]
        out[f"recovery_status_{ap}"] = "bls_recovered" if rec["recovered"] else str(ap_row.get("status", "failed"))
        match = rec["truth_match"]
        out[f"top_period_ratio_{ap}"] = match["top_period_ratio"]
        out[f"top_harmonic_factor_{ap}"] = match["top_harmonic_factor"]
        out[f"top_harmonic_rel_err_{ap}"] = match["top_harmonic_rel_err"]
        out[f"top_harmonic_t0_phase_err_min_{ap}"] = match["top_harmonic_t0_phase_err_min"]
        out[f"top_harmonic_period_match_{ap}"] = match["top_harmonic_period_match"]
        out[f"top_harmonic_ephemeris_match_{ap}"] = match["top_harmonic_ephemeris_match"]
        out[f"topn_exact_recovered_{ap}"] = match["topn_exact_recovered"]
        out[f"topn_exact_peak_rank_{ap}"] = match["topn_exact_peak_rank"]
        out[f"topn_exact_sde_{ap}"] = match["topn_exact_sde"]
        out[f"topn_exact_period_rel_err_{ap}"] = match["topn_exact_period_rel_err"]
        out[f"topn_exact_t0_phase_err_min_{ap}"] = match["topn_exact_t0_phase_err_min"]
        out[f"topn_harmonic_match_{ap}"] = match["topn_harmonic_match"]
        out[f"topn_harmonic_peak_rank_{ap}"] = match["topn_harmonic_peak_rank"]
        out[f"topn_harmonic_sde_{ap}"] = match["topn_harmonic_sde"]
        out[f"topn_harmonic_period_ratio_{ap}"] = match["topn_harmonic_period_ratio"]
        out[f"topn_harmonic_factor_{ap}"] = match["topn_harmonic_factor"]
        out[f"topn_harmonic_rel_err_{ap}"] = match["topn_harmonic_rel_err"]
        out[f"topn_harmonic_t0_phase_err_min_{ap}"] = match["topn_harmonic_t0_phase_err_min"]
    return out


def run_injection_bls(
    *,
    injection_h5: Path,
    apertures: tuple[str, ...],
    n_injections: int,
    workers: int,
    n_periods: int,
    limit_keys: list[str] | None = None,
    cfg_kwargs_extra: dict[str, Any] | None = None,
) -> pd.DataFrame:
    import h5py
    from concurrent.futures import ProcessPoolExecutor

    if not injection_h5.exists():
        raise FileNotFoundError(f"missing injection HDF5: {injection_h5}")
    with h5py.File(injection_h5, "r") as h5:
        keys = sorted(h5["injections"].keys())
    if limit_keys is not None:
        wanted = set(limit_keys)
        keys = [k for k in keys if k in wanted]
    if n_injections > 0:
        keys = keys[:n_injections]
    cfg_kwargs = {"n_periods": int(n_periods), "n_peaks": 5}
    if cfg_kwargs_extra:
        cfg_kwargs.update(cfg_kwargs_extra)
    payloads = [(str(injection_h5), f"/injections/{key}", apertures, cfg_kwargs) for key in keys]
    print(
        f"[review] running BLS on {len(payloads):,} injected light curves "
        f"across {len(apertures)} apertures",
        flush=True,
    )
    if workers <= 1:
        rows = []
        for idx, payload in enumerate(payloads, 1):
            rows.append(_run_one_injection_bls(payload))
            if idx % 50 == 0:
                print(f"  [review] injection BLS {idx:,}/{len(payloads):,}", flush=True)
    else:
        rows = []
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for idx, row in enumerate(ex.map(_run_one_injection_bls, payloads, chunksize=4), 1):
                rows.append(row)
                if idx % 50 == 0:
                    print(f"  [review] injection BLS {idx:,}/{len(payloads):,}", flush=True)
    return pd.DataFrame(rows)


def _take_bucket(
    df: pd.DataFrame,
    *,
    name: str,
    mask: pd.Series,
    n: int,
    used: set[int],
    random_state: int,
    random: bool = False,
) -> pd.DataFrame:
    sub = df[mask & ~df["tic"].astype(int).isin(used)].copy()
    if sub.empty or n <= 0:
        return sub.head(0)
    if random:
        sub = sub.sample(n=min(n, len(sub)), random_state=random_state)
    else:
        sub = sub.sort_values(["class_rank", "sde_max"], ascending=[True, False]).head(n)
    sub["source_bucket"] = name
    used.update(sub["tic"].astype(int).tolist())
    return sub


def _attach_real_review_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["review_id"] = out["tic"].map(lambda tic: f"real:{int(tic)}")
    out["source_kind"] = "real_candidate"
    out["recovery_status"] = "real_candidate"
    out["injection_id"] = ""
    out["signal_family"] = ""
    out["truth_period_d"] = ""
    out["truth_t0_bjd"] = ""
    out["truth_duration_min"] = ""
    out["truth_depth"] = ""
    out["truth_target_depth"] = ""
    out["truth_geometric_depth"] = ""
    out["truth_model_depth"] = ""
    out["truth_sampled_model_depth"] = ""
    out["truth_radius_rearth"] = ""
    out["truth_radius_rwd"] = ""
    out["truth_impact_b"] = ""
    out["truth_a_over_rwd"] = ""
    out["truth_inclination_deg"] = ""
    out["truth_duration_model"] = ""
    out["truth_injection_model"] = ""
    out["truth_sampling_mode"] = ""
    out["truth_grid_cell_id"] = ""
    out["truth_grid_period_bin"] = ""
    out["truth_grid_depth_bin"] = ""
    out["truth_grid_radius_bin"] = ""
    out["truth_source_kind"] = "real_candidate"
    out["truth_source_bucket"] = out["source_bucket"]
    out["period_rel_err"] = ""
    out["t0_phase_err_min"] = ""
    out["leo_class"] = ""
    out["leo_report_name"] = ""
    return out


def _attach_ranker_real_review_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize ranker-selected BLS peaks into browser review rows."""

    out = df.copy()
    out["tic"] = out["tic"].astype(int)
    if "ranker_selection_rank" not in out:
        out["ranker_selection_rank"] = 1
    out["ranker_selection_rank"] = pd.to_numeric(
        out["ranker_selection_rank"],
        errors="coerce",
    ).fillna(1).astype(int)
    if "sde_max" not in out and "sde" in out:
        out["sde_max"] = out["sde"]
    if "rep_aperture" not in out and "aperture" in out:
        out["rep_aperture"] = out["aperture"]
    if "vet_class" not in out:
        out["vet_class"] = "ranker_selected"
    out["source_bucket"] = "real_ranker_selected"
    if "class_rank" not in out:
        out["class_rank"] = out["ranker_selection_rank"]
    if "blind_rank" not in out:
        out["blind_rank"] = np.arange(1, len(out) + 1)
    if "n_apertures_agree" not in out:
        out["n_apertures_agree"] = 1
    if "apertures_agree" not in out:
        out["apertures_agree"] = out.get("rep_aperture", "")
    for col, value in (
        ("centroid_status", "not_run"),
        ("centroid_pass", ""),
        ("centroid_delta_pix", ""),
        ("centroid_z", ""),
        ("n_in_transit", ""),
        ("n_oot_band", ""),
    ):
        if col not in out:
            out[col] = value
    out = _attach_real_review_columns(out)
    out["review_id"] = out.apply(
        lambda row: f"real:{int(row['tic'])}:ranker:{int(row['ranker_selection_rank'])}",
        axis=1,
    )
    out["recovery_status"] = "real_ranker_selected"
    out["truth_source_bucket"] = out["source_bucket"]
    return out


def select_real_candidates(path: Path, n_real: int, random_state: int) -> pd.DataFrame:
    if n_real <= 0:
        return pd.DataFrame()
    df = _read_table(path).copy()
    df["tic"] = df["tic"].astype(int)
    used: set[int] = set()
    pieces: list[pd.DataFrame] = []

    wd = _take_bucket(
        df,
        name="wd1856_benchmark",
        mask=df["tic"].eq(WD_1856_TIC),
        n=1,
        used=used,
        random_state=random_state,
    )
    pieces.append(wd)
    remaining = max(n_real - sum(len(p) for p in pieces), 0)
    quotas = [
        ("planet_centroid_pass", df["vet_class"].eq("planet_candidate") & df.get("centroid_pass", False).eq(True), 0.25),
        ("planet_candidate_top", df["vet_class"].eq("planet_candidate"), 0.25),
        ("sub_roche_pceb_suspect", df["vet_class"].eq("sub_roche_pceb_suspect"), 0.20),
        ("pceb_grid_ceiling", df["vet_class"].eq("pceb_grid_ceiling"), 0.20),
        ("high_sde_all_classes", df["sde_max"].notna(), 0.10),
    ]
    for name, mask, frac in quotas:
        take_n = max(1, int(round(n_real * frac)))
        pieces.append(
            _take_bucket(
                df,
                name=name,
                mask=mask,
                n=min(take_n, remaining),
                used=used,
                random_state=random_state,
            )
        )
        remaining = max(n_real - sum(len(p) for p in pieces), 0)
    if remaining > 0:
        pieces.append(
            _take_bucket(
                df,
                name="random_background",
                mask=df["sde_max"].notna(),
                n=remaining,
                used=used,
                random_state=random_state,
                random=True,
            )
        )
    out = pd.concat([p for p in pieces if len(p)], ignore_index=True)
    return _attach_real_review_columns(out.head(n_real).copy())


def select_real_candidates_stratified(path: Path, n_real: int, random_state: int) -> pd.DataFrame:
    """Select a broad, blinded real-candidate sample across BLS/vetter space."""
    if n_real <= 0:
        return pd.DataFrame()
    df = _read_table(path).copy()
    df["tic"] = df["tic"].astype(int)
    if "sde_max" in df:
        df["sde_max"] = pd.to_numeric(df["sde_max"], errors="coerce")
    if "period_d" in df:
        df["period_d"] = pd.to_numeric(df["period_d"], errors="coerce")
    if "depth" in df:
        df["depth"] = pd.to_numeric(df["depth"], errors="coerce")
    aperture_agree = (
        pd.to_numeric(df["n_apertures_agree"], errors="coerce")
        if "n_apertures_agree" in df
        else pd.Series(np.nan, index=df.index)
    )

    rng = np.random.default_rng(random_state)
    used: set[int] = set()
    pieces: list[pd.DataFrame] = []
    finite_sde = df["sde_max"].replace([np.inf, -np.inf], np.nan).dropna()
    q40 = float(finite_sde.quantile(0.40)) if len(finite_sde) else -np.inf
    q70 = float(finite_sde.quantile(0.70)) if len(finite_sde) else -np.inf
    q90 = float(finite_sde.quantile(0.90)) if len(finite_sde) else -np.inf

    quotas = [
        ("real_top_sde", df["sde_max"].ge(q90), 0.20, False),
        ("real_planet_candidate", df["vet_class"].eq("planet_candidate"), 0.18, False),
        ("real_short_period", df["period_d"].between(0.08, 0.35, inclusive="both"), 0.12, False),
        ("real_pceb_like", df["vet_class"].isin(["sub_roche_pceb_suspect", "pceb_grid_ceiling"]), 0.15, False),
        ("real_single_aperture", aperture_agree.eq(1), 0.10, True),
        ("real_mid_sde", df["sde_max"].between(q40, q70, inclusive="both"), 0.15, True),
        ("real_low_sde_control", df["sde_max"].lt(q40) & df["sde_max"].notna(), 0.10, True),
    ]
    for name, mask, frac, random_pick in quotas:
        take_n = max(1, int(round(n_real * frac)))
        piece = _take_bucket(
            df,
            name=name,
            mask=mask,
            n=take_n,
            used=used,
            random_state=random_state,
            random=random_pick,
        )
        pieces.append(piece)

    remaining = max(n_real - sum(len(p) for p in pieces), 0)
    if remaining > 0:
        pieces.append(
            _take_bucket(
                df,
                name="real_stratified_fill",
                mask=df["sde_max"].notna(),
                n=remaining,
                used=used,
                random_state=random_state,
                random=True,
            )
        )

    out = pd.concat([p for p in pieces if len(p)], ignore_index=True)
    if len(out) > n_real:
        out = out.sample(n=n_real, random_state=random_state).reset_index(drop=True)
    else:
        out = out.reset_index(drop=True)
    order = rng.permutation(len(out))
    out = out.iloc[order].reset_index(drop=True)
    return _attach_real_review_columns(out)


def select_ranker_real_candidates(path: Path, n_real: int, random_state: int) -> pd.DataFrame:
    """Select ranker-chosen real ephemerides for LEO/human review."""

    if n_real <= 0:
        return pd.DataFrame()
    df = _read_table(path).copy()
    if df.empty:
        return _attach_ranker_real_review_columns(df)
    if "tic" not in df:
        raise KeyError("ranker-selected real candidate table must include tic")
    for required in ("period_d", "t0_bjd", "duration_min"):
        if required not in df:
            raise KeyError(f"ranker-selected real candidate table must include {required}")
    score_cols = [col for col in df.columns if col.startswith("ranker_p_")]
    if "ranker_p_signal_peak" in df.columns:
        df["_ranker_sort_score"] = pd.to_numeric(df["ranker_p_signal_peak"], errors="coerce")
    elif score_cols:
        df["_ranker_sort_score"] = pd.to_numeric(df[score_cols[0]], errors="coerce")
    elif "sde" in df:
        df["_ranker_sort_score"] = pd.to_numeric(df["sde"], errors="coerce")
    else:
        df["_ranker_sort_score"] = np.nan
    if "ranker_selection_rank" in df:
        df["_ranker_selection_rank_sort"] = pd.to_numeric(df["ranker_selection_rank"], errors="coerce")
    else:
        df["_ranker_selection_rank_sort"] = 1
    df = df.sort_values(
        ["_ranker_selection_rank_sort", "_ranker_sort_score"],
        ascending=[True, False],
        kind="stable",
    )
    out = df.head(n_real).copy()
    out = out.drop(columns=[c for c in ("_ranker_sort_score", "_ranker_selection_rank_sort") if c in out])
    if len(out) > 1:
        out = out.sample(frac=1.0, random_state=random_state).reset_index(drop=True)
    return _attach_ranker_real_review_columns(out)


_RHO_WD_G_CM3 = 3.85e5
_G_CGS = 6.67430e-8
_MSUN_G = 1.98847e33
_RSUN_CM = 6.957e10
_CANONICAL_WD_STAR = {
    "rad": 0.013,
    "e_rad": 0.002,
    "mass": 0.6,
    "e_mass": 0.05,
    "Teff": 10000.0,
    "e_Teff": 1000.0,
    "rho": _RHO_WD_G_CM3,
    "u1": 0.05,
    "u2": 0.05,
    "Rs": 0.013,
}


def _as_finite_float(value: Any, default: float = float("nan")) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if np.isfinite(out) else default


def _radius_rsun_from_mass_logg(mass_msun: float, logg_cgs: float) -> float:
    mass = _as_finite_float(mass_msun)
    logg = _as_finite_float(logg_cgs)
    if not np.isfinite(mass) or not np.isfinite(logg) or mass <= 0:
        return float("nan")
    return float(np.sqrt(_G_CGS * mass * _MSUN_G / (10.0 ** logg)) / _RSUN_CM)


def _radius_uncertainty_rsun(
    radius_rsun: float,
    mass_msun: float,
    e_mass_msun: float,
    e_logg_cgs: float,
) -> float:
    radius = _as_finite_float(radius_rsun)
    mass = _as_finite_float(mass_msun)
    e_mass = _as_finite_float(e_mass_msun)
    e_logg = _as_finite_float(e_logg_cgs)
    if radius <= 0 or mass <= 0 or not np.isfinite(radius):
        return float("nan")
    terms = []
    if np.isfinite(e_mass) and e_mass > 0:
        terms.append((e_mass / mass) ** 2)
    if np.isfinite(e_logg) and e_logg > 0:
        terms.append((np.log(10.0) * e_logg) ** 2)
    if not terms:
        return float("nan")
    return float(0.5 * radius * np.sqrt(sum(terms)))


def _density_g_cm3(mass_msun: float, radius_rsun: float) -> float:
    mass = _as_finite_float(mass_msun)
    radius = _as_finite_float(radius_rsun)
    if mass <= 0 or radius <= 0 or not np.isfinite(mass) or not np.isfinite(radius):
        return float("nan")
    volume = 4.0 / 3.0 * np.pi * (radius * _RSUN_CM) ** 3
    return float(mass * _MSUN_G / volume)


def _valid_wd_star_record(record: dict[str, Any]) -> bool:
    radius = _as_finite_float(record.get("star_rad_rsun"))
    mass = _as_finite_float(record.get("star_mass_msun"))
    teff = _as_finite_float(record.get("star_teff_k"))
    logg = _as_finite_float(record.get("star_logg_cgs"))
    return (
        0.002 <= radius <= 0.08
        and 0.15 <= mass <= 1.4
        and 2_000.0 <= teff <= 120_000.0
        and 6.0 <= logg <= 10.0
    )


def _canonical_star_record(tic: int, *, source_id: Any = "") -> dict[str, Any]:
    return {
        "star_source": "canonical_fallback",
        "star_source_id": str(_clean(source_id)) if _clean(source_id) != "" else "",
        "star_atmosphere": "",
        "star_rad_rsun": _CANONICAL_WD_STAR["rad"],
        "star_e_rad_rsun": _CANONICAL_WD_STAR["e_rad"],
        "star_mass_msun": _CANONICAL_WD_STAR["mass"],
        "star_e_mass_msun": _CANONICAL_WD_STAR["e_mass"],
        "star_teff_k": _CANONICAL_WD_STAR["Teff"],
        "star_e_teff_k": _CANONICAL_WD_STAR["e_Teff"],
        "star_logg_cgs": "",
        "star_e_logg_cgs": "",
        "star_rho_g_cm3": _CANONICAL_WD_STAR["rho"],
        "star_radius_source": "canonical",
    }


def _star_record_from_catalog_row(
    row: Any,
    *,
    tic: int,
    source_id: Any,
    atmosphere: str,
) -> dict[str, Any] | None:
    suffix = str(atmosphere)
    mass = _as_finite_float(row[f"mass_{suffix}"])
    e_mass = _as_finite_float(row[f"emass_{suffix}"])
    teff = _as_finite_float(row[f"teff_{suffix}"])
    e_teff = _as_finite_float(row[f"eteff_{suffix}"])
    logg = _as_finite_float(row[f"logg_{suffix}"])
    e_logg = _as_finite_float(row[f"elogg_{suffix}"])
    radius = _radius_rsun_from_mass_logg(mass, logg)
    e_radius = _radius_uncertainty_rsun(radius, mass, e_mass, e_logg)
    if not np.isfinite(e_radius) or e_radius <= 0:
        e_radius = 0.2 * radius if np.isfinite(radius) and radius > 0 else float("nan")
    try:
        clean_source_id = str(int(source_id))
    except (TypeError, ValueError, OverflowError):
        clean_source_id = str(_clean(source_id)) if _clean(source_id) != "" else ""
    record = {
        "star_source": f"GF21_{suffix}",
        "star_source_id": clean_source_id,
        "star_atmosphere": suffix,
        "star_rad_rsun": radius,
        "star_e_rad_rsun": e_radius,
        "star_mass_msun": mass,
        "star_e_mass_msun": e_mass if np.isfinite(e_mass) else 0.1 * mass,
        "star_teff_k": teff,
        "star_e_teff_k": e_teff if np.isfinite(e_teff) else 1000.0,
        "star_logg_cgs": logg,
        "star_e_logg_cgs": e_logg if np.isfinite(e_logg) else "",
        "star_rho_g_cm3": _density_g_cm3(mass, radius),
        "star_radius_source": "sqrt_GM_over_g",
    }
    return record if _valid_wd_star_record(record) else None


def load_wd_star_catalog(
    catalog_path: Path | None,
    *,
    tics: set[int],
    tic_column: str = "tic_id",
    atmosphere_priority: tuple[str, ...] = ("H", "He", "mixed"),
) -> dict[int, dict[str, Any]]:
    """Load TIC-keyed WD host parameters from the local GF21/TWIRL catalog.

    The Gentile-Fusillo catalog stores atmosphere-fit mass and logg, but not a
    direct radius column. Radius is derived from ``g = GM/R^2`` and the chosen
    atmosphere solution. Rows without a valid atmosphere solution are left for
    canonical fallback during annotation.
    """

    if catalog_path is None:
        return {}
    catalog_path = Path(catalog_path)
    if not catalog_path.exists():
        raise FileNotFoundError(f"missing WD star catalog: {catalog_path}")
    wanted = {int(tic) for tic in tics if np.isfinite(float(tic))}
    if not wanted:
        return {}

    required = {"source_id", tic_column}
    for suffix in atmosphere_priority:
        required.update(
            {
                f"teff_{suffix}",
                f"eteff_{suffix}",
                f"logg_{suffix}",
                f"elogg_{suffix}",
                f"mass_{suffix}",
                f"emass_{suffix}",
            }
        )

    records: dict[int, dict[str, Any]] = {}
    with fits.open(catalog_path, memmap=True) as hdul:
        data = hdul[1].data
        names = set(data.columns.names)
        missing = sorted(required - names)
        if missing:
            raise KeyError(f"WD star catalog missing columns: {missing}")
        tic_values = np.asarray(data[tic_column], dtype=np.int64)
        idx = np.flatnonzero(np.isin(tic_values, list(wanted)))
        for row_idx in idx:
            tic = int(tic_values[row_idx])
            if tic in records:
                continue
            row = data[row_idx]
            source_id = row["source_id"]
            for suffix in atmosphere_priority:
                record = _star_record_from_catalog_row(
                    row,
                    tic=tic,
                    source_id=source_id,
                    atmosphere=suffix,
                )
                if record is not None:
                    records[tic] = record
                    break
    return records


def annotate_star_parameters(
    queue: pd.DataFrame,
    *,
    star_records: dict[int, dict[str, Any]],
) -> pd.DataFrame:
    if queue.empty:
        return queue
    out = queue.copy()
    for col in REVIEW_COLUMNS:
        if col.startswith("star_") and col not in out:
            out[col] = ""
    for idx, row in out.iterrows():
        tic = int(row["tic"])
        record = star_records.get(tic, _canonical_star_record(tic))
        for key, value in record.items():
            out.loc[idx, key] = _clean(value)
    return out


def star_source_counts(queue: pd.DataFrame) -> dict[str, int]:
    if "star_source" not in queue:
        return {}
    return {
        str(key): int(value)
        for key, value in queue["star_source"].fillna("").astype(str).value_counts().sort_index().items()
    }


def _wd_star_for(row_or_tic: Any) -> dict[str, Any]:
    if isinstance(row_or_tic, pd.Series):
        tic = int(row_or_tic.get("tic", -1))
        radius = _as_finite_float(row_or_tic.get("star_rad_rsun"))
        mass = _as_finite_float(row_or_tic.get("star_mass_msun"))
        teff = _as_finite_float(row_or_tic.get("star_teff_k"))
        if np.isfinite(radius) and radius > 0 and np.isfinite(mass) and mass > 0 and np.isfinite(teff):
            star = {
                "rad": radius,
                "e_rad": max(_as_finite_float(row_or_tic.get("star_e_rad_rsun"), 0.2 * radius), 1.0e-5),
                "mass": mass,
                "e_mass": max(_as_finite_float(row_or_tic.get("star_e_mass_msun"), 0.1 * mass), 1.0e-5),
                "Teff": teff,
                "e_Teff": max(_as_finite_float(row_or_tic.get("star_e_teff_k"), 1000.0), 1.0),
                "rho": _as_finite_float(row_or_tic.get("star_rho_g_cm3"), _density_g_cm3(mass, radius)),
                "u1": _CANONICAL_WD_STAR["u1"],
                "u2": _CANONICAL_WD_STAR["u2"],
                "Rs": radius,
                "id": tic,
                "tic": tic,
                "source_id": _clean(row_or_tic.get("star_source_id", "")),
                "star_source": _clean(row_or_tic.get("star_source", "")),
                "logg": _clean(row_or_tic.get("star_logg_cgs", "")),
                "e_logg": _clean(row_or_tic.get("star_e_logg_cgs", "")),
            }
            return star
    else:
        tic = int(row_or_tic)
    star = dict(_CANONICAL_WD_STAR)
    star["id"] = int(tic)
    star["tic"] = int(tic)
    star["star_source"] = "canonical_fallback"
    return star


def _leo_label(fa: bool, fp: bool) -> str:
    if (not fa) and (not fp):
        return "PC"
    if fa:
        return "FA"
    return "FP"


def _aperture_for_row(row: pd.Series, fallback: str) -> str:
    rep = str(row.get("rep_aperture", "")).strip()
    if rep in KNOWN_APERTURES:
        return rep
    return fallback


def _mad_error_for_arrays(time: np.ndarray, flux: np.ndarray, quality: np.ndarray) -> np.ndarray:
    good = (quality == 0) & np.isfinite(time) & np.isfinite(flux)
    err = np.full_like(flux, np.nan, dtype=np.float64)
    if np.any(good):
        med = np.nanmedian(flux[good])
        mad = np.nanmedian(np.abs(flux[good] - med))
        sigma = float(1.4826 * mad)
        err[good] = sigma if np.isfinite(sigma) and sigma > 0 else np.nanstd(flux[good])
    return err


def _write_fallback_leo_report(
    report_path: Path,
    *,
    tlc: Any,
    label: str,
    aperture: str,
    plot_error: str,
) -> None:
    """Write a minimal triage PDF when LEO metrics succeed but plotting fails."""
    import matplotlib.pyplot as plt

    report_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(3, 1, figsize=(10, 7), constrained_layout=True)
    rel_time = np.asarray(tlc.time, dtype=float) - BJDREFI
    axes[0].plot(rel_time, tlc.raw, ".", ms=2, color="0.45", label="raw")
    axes[0].plot(rel_time, tlc.flux, ".", ms=2, color="tab:blue", label="detrended")
    axes[0].set_xlabel("BJD - 2457000")
    axes[0].set_ylabel("Flux")
    axes[0].legend(loc="best", fontsize=8)

    phase_hours = np.asarray(tlc.phase, dtype=float) * float(tlc.per) * 24.0
    axes[1].plot(phase_hours, tlc.flux, ".", ms=2, color="tab:blue")
    axes[1].axvspan(-0.5 * float(tlc.dur) * 24.0, 0.5 * float(tlc.dur) * 24.0, color="tab:red", alpha=0.15)
    axes[1].set_xlim(-4.0 * float(tlc.dur) * 24.0, 4.0 * float(tlc.dur) * 24.0)
    axes[1].set_xlabel("Hours from midtransit")
    axes[1].set_ylabel("Detrended flux")

    lines = [
        f"TIC {tlc.tic}",
        f"WD-tuned LEO class: {label}",
        f"Aperture: {aperture}",
        f"Period: {tlc.per:.6f} d",
        f"Epoch: {tlc.epo:.6f} BJD",
        f"Duration: {tlc.dur * 1440.0:.2f} min",
        f"MES: {tlc.metrics.get('MES', np.nan):.3g}",
        f"N_transit: {tlc.metrics.get('N_transit', np.nan)}",
        f"n_in: {tlc.metrics.get('n_in', np.nan)}",
        "",
        "LEO plot_summary failed after metric computation:",
        plot_error,
    ]
    axes[2].axis("off")
    axes[2].text(0.0, 1.0, "\n".join(lines), va="top", ha="left", fontsize=9, family="monospace")
    fig.savefig(report_path, bbox_inches="tight", dpi=150)
    plt.close(fig)


def _finite_float(value: Any, default: float = float("nan")) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return float(default)
    return out if np.isfinite(out) else float(default)


def _phasefold_for_plot(time: np.ndarray, period: float, epoch: float) -> np.ndarray:
    phase = np.mod(time - epoch + 0.5 * period, period) / period - 0.5
    return np.asarray(phase, dtype=np.float64)


def _refresh_plot_windows(tlc: Any, period: float, epoch: float, duration: float) -> None:
    tlc.per = float(period)
    tlc.epo = float(epoch)
    tlc.dur = float(duration)
    tlc.qtran = float(duration / period)
    tlc.phase = _phasefold_for_plot(np.asarray(tlc.time, dtype=np.float64), tlc.per, tlc.epo)
    tlc.in_tran = np.abs(tlc.phase) < 0.5 * tlc.qtran
    tlc.near_tran = np.abs(tlc.phase) < tlc.qtran
    tlc.fit_tran = np.abs(tlc.phase) < 2.0 * tlc.qtran
    tlc.before_tran = (tlc.phase < -0.5 * tlc.qtran) & (tlc.phase > -1.5 * tlc.qtran)
    tlc.after_tran = (tlc.phase > 0.5 * tlc.qtran) & (tlc.phase < 1.5 * tlc.qtran)
    phase2 = np.mod(np.asarray(tlc.time, dtype=np.float64) - tlc.epo, 2.0 * tlc.per) / tlc.per
    phase2[phase2 > 1.0] -= 2.0
    tlc.odd_tran = np.abs(phase2) < 0.5 * tlc.qtran
    tlc.even_tran = np.abs(phase2) > 1.0 - 0.5 * tlc.qtran
    tlc.epochs = np.round((np.asarray(tlc.time, dtype=np.float64) - tlc.epo) / tlc.per)
    tlc.tran_epochs = np.unique(tlc.epochs[tlc.in_tran])
    tlc.N_transit = len(tlc.tran_epochs)
    tlc.n_in = int(np.count_nonzero(tlc.in_tran))
    tlc.n_before = int(np.count_nonzero(tlc.before_tran))
    tlc.n_after = int(np.count_nonzero(tlc.after_tran))


def _model_grid_is_safe(period: float, duration: float, *, max_points: int = 50_000) -> bool:
    if not all(np.isfinite([period, duration])) or period <= 0 or duration <= 0:
        return False
    points = 100.0 * period / duration
    return np.isfinite(points) and 10 <= points <= max_points


def _leo_plot_window_has_finite_data(
    time: np.ndarray,
    flux: np.ndarray,
    period: float,
    epoch: float,
    duration: float,
) -> bool:
    """Mirror LEO's individual-transit plot window preflight.

    LEO's ``plot_summary`` prefers fitted transit/trapezoid ephemerides over
    the BLS ephemeris when their AIC values are finite. At 200 s cadence a
    fitted duration can become narrower than the sampled cadence grid, causing
    the individual-transit panel to call ``nanmax`` on an empty in-window
    array. Treat those fitted models as unsafe for plotting only.
    """

    if not all(np.isfinite([period, epoch, duration])) or period <= 0 or duration <= 0:
        return False
    time = np.asarray(time, dtype=np.float64)
    flux = np.asarray(flux, dtype=np.float64)
    phase = np.mod(time - epoch, period) / period
    phase[phase > 0.5] -= 1.0
    near_tran = np.abs(phase) < duration / period
    finite = np.isfinite(phase) & np.isfinite(flux)
    return bool(np.any(finite & near_tran) and np.any(finite & ~near_tran))


def _leo_plot_copy(tlc: Any) -> tuple[Any, str]:
    """Return a LEO object copy whose metrics are safe for plot_summary only."""
    import copy

    plot_tlc = copy.copy(tlc)
    plot_tlc.metrics = dict(tlc.metrics)
    metrics = plot_tlc.metrics
    notes: list[str] = []

    period = _finite_float(metrics.get("per", getattr(tlc, "per", np.nan)), getattr(tlc, "per", np.nan))
    epoch = _finite_float(metrics.get("epo", getattr(tlc, "epo", np.nan)), getattr(tlc, "epo", np.nan))
    duration = _finite_float(metrics.get("dur", getattr(tlc, "dur", np.nan)), getattr(tlc, "dur", np.nan))
    if not np.isfinite(period) or period <= 0:
        period = max(_finite_float(getattr(tlc, "per", np.nan), 1.0), 1.0e-3)
        notes.append("period")
    if not np.isfinite(epoch):
        epoch = _finite_float(getattr(tlc, "epo", np.nan), float(np.nanmedian(tlc.time)))
        notes.append("epoch")
    if not np.isfinite(duration) or duration <= 0:
        duration = max(_finite_float(getattr(tlc, "dur", np.nan), 0.01 * period), 1.0e-5)
        notes.append("duration")
    duration = float(np.clip(duration, 2.0e-4 * period, 0.20 * period))

    # LEO's individual-transit plot requires at least one near-transit and one
    # out-of-transit cadence. Widen only the display window when the BLS/fit
    # ephemeris is too narrow for the cadence grid.
    for _ in range(12):
        _refresh_plot_windows(plot_tlc, period, epoch, duration)
        if np.any(plot_tlc.near_tran) and np.any(~plot_tlc.near_tran):
            break
        duration = min(duration * 2.0, 0.20 * period)
    metrics["per"] = period
    metrics["epo"] = epoch
    metrics["dur"] = duration
    metrics["qtran"] = duration / period
    metrics["N_transit"] = plot_tlc.N_transit
    metrics["n_in"] = plot_tlc.n_in
    metrics["n_before"] = plot_tlc.n_before
    metrics["n_after"] = plot_tlc.n_after

    zpt = _finite_float(metrics.get("zpt", getattr(tlc, "zpt", np.nan)), float(np.nanmedian(plot_tlc.flux)))
    dep = _finite_float(metrics.get("dep", getattr(tlc, "dep", np.nan)), 0.0)
    dep = float(np.clip(dep if dep > 0 else abs(dep), 1.0e-8, 0.999))
    plot_tlc.zpt = zpt
    plot_tlc.dep = dep
    metrics["zpt"] = zpt
    metrics["dep"] = dep

    rp = float(np.sqrt(dep))
    a_rs = max(1.0 / max(np.pi * metrics["qtran"], 1.0e-6), 1.5)
    defaults = {
        "transit_per": period,
        "transit_epo": epoch,
        "transit_dur": duration,
        "transit_RpRs": rp,
        "transit_RpRs_err": 0.0,
        "transit_aRs": a_rs,
        "transit_aRs_err": 0.0,
        "transit_b": 0.0,
        "transit_b_err": 0.0,
        "transit_u1": 0.05,
        "transit_u2": 0.05,
        "transit_zpt": zpt,
        "transit_odd_RpRs": rp,
        "transit_odd_epo": epoch,
        "transit_even_RpRs": rp,
        "transit_even_epo": epoch,
        "trap_per": period,
        "trap_epo": epoch,
        "trap_dur": duration,
        "trap_dep": dep,
        "trap_qtran": metrics["qtran"],
        "trap_qin": 0.5,
        "trap_zpt": zpt,
        "trap_odd_dep": dep,
        "trap_odd_epo": epoch,
        "trap_even_dep": dep,
        "trap_even_epo": epoch,
        "trap_sig_dep": np.nan,
        "transit_sig_dep": np.nan,
        "sig_dep": np.nan,
        "Rp": metrics.get("Rp", np.nan),
        "Rp_err": metrics.get("Rp_err", np.nan),
        "a": metrics.get("a", np.nan),
        "a_err": metrics.get("a_err", np.nan),
        "Seff": metrics.get("Seff", np.nan),
        "Seff_err": metrics.get("Seff_err", np.nan),
    }
    for key, default in defaults.items():
        value = _finite_float(metrics.get(key, np.nan), float(default))
        metrics[key] = value

    transit_safe = (
        _model_grid_is_safe(metrics["transit_per"], metrics["transit_dur"])
        and _leo_plot_window_has_finite_data(
            np.asarray(plot_tlc.time, dtype=np.float64),
            np.asarray(plot_tlc.flux, dtype=np.float64),
            metrics["transit_per"],
            metrics["transit_epo"],
            metrics["transit_dur"],
        )
        and metrics["transit_RpRs"] > 0
        and metrics["transit_aRs"] > metrics["transit_b"]
        and metrics["transit_RpRs"] < 20
    )
    trap_safe = (
        _model_grid_is_safe(metrics["trap_per"], metrics["trap_dur"])
        and _leo_plot_window_has_finite_data(
            np.asarray(plot_tlc.time, dtype=np.float64),
            np.asarray(plot_tlc.flux, dtype=np.float64),
            metrics["trap_per"],
            metrics["trap_epo"],
            metrics["trap_dur"],
        )
        and 0 < metrics["trap_qtran"] < 0.5
        and 0 <= metrics["trap_dep"] <= 1.0
    )
    if not transit_safe:
        metrics["transit_aic"] = np.nan
        notes.append("transit_model_disabled")
    if not trap_safe:
        metrics["trap_aic"] = np.nan
        notes.append("trap_model_disabled")

    phs_sec = _finite_float(metrics.get("phs_sec", np.nan), np.nan)
    dep_sec = _finite_float(metrics.get("dep_sec", np.nan), np.nan)
    sig_sec = _finite_float(metrics.get("sig_sec", np.nan), np.nan)
    if not all(np.isfinite([phs_sec, dep_sec, sig_sec])):
        metrics["phs_sec"] = np.nan
        metrics["dep_sec"] = np.nan
        metrics["sig_sec"] = np.nan
        notes.append("secondary_disabled")

    return plot_tlc, ",".join(dict.fromkeys(notes))


def _leo_cache_dirs() -> None:
    pid = os.getpid()
    mpl_cache = Path(tempfile.gettempdir()) / f"twirl_mplconfig_{pid}"
    mpl_cache.mkdir(parents=True, exist_ok=True)
    xdg_cache = Path(tempfile.gettempdir()) / f"twirl_cache_{pid}"
    xdg_cache.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(mpl_cache)
    os.environ["XDG_CACHE_HOME"] = str(xdg_cache)


def _leo_arrays_for_row(row: pd.Series, hlsp_root: Path, injection_h5: Path, aperture: str):
    source_kind = str(row.get("source_kind", ""))
    if source_kind == "injection_recovery":
        import h5py

        with h5py.File(injection_h5, "r") as h5:
            group = h5[str(row["h5_group"])]
            time = np.asarray(group["time"], dtype=np.float64) + BJDREFI
            injected_name = f"{aperture}_injected"
            original_name = f"{aperture}_original"
            if injected_name in group:
                flux = np.asarray(group[injected_name], dtype=np.float64)
            else:
                flux = np.asarray(group["flux_injected"], dtype=np.float64)
            raw = np.asarray(group[original_name], dtype=np.float64) if original_name in group else (
                np.asarray(group["flux_original"], dtype=np.float64)
                if "flux_original" in group
                else flux.copy()
            )
            quality = np.asarray(group["quality"], dtype=np.int32)
        err = _mad_error_for_arrays(time, flux, quality)
        good = (quality == 0) & np.isfinite(time) & np.isfinite(flux) & np.isfinite(err)
        return time[good], raw[good], flux[good], err[good]

    tic = int(row["tic"])
    sector = int(row.get("sector", 56))
    path = find_hlsp_path(hlsp_root, tic, sector)
    if path is None:
        raise FileNotFoundError(f"no HLSP for TIC {tic}")
    lc = read_hlsp(path, columns=("SAP_FLUX", aperture))
    if lc is None:
        raise RuntimeError(f"read_hlsp failed: {path}")
    sigma = tglc_mad_error(lc, aperture=aperture)
    good = quality_mask(lc, aperture) & np.isfinite(lc.flux.get("SAP_FLUX", np.nan)) & np.isfinite(sigma)
    return (
        lc.time[good] + BJDREFI,
        lc.flux["SAP_FLUX"][good].astype(float),
        lc.flux[aperture][good].astype(float),
        sigma[good].astype(float),
    )


class _LeoTimeoutError(TimeoutError):
    pass


class _leo_timeout:
    def __init__(self, seconds: int):
        self.seconds = int(seconds)
        self._old_handler = None

    def __enter__(self):
        if self.seconds <= 0:
            return self

        def _raise_timeout(signum, frame):
            raise _LeoTimeoutError(f"LEO row exceeded {self.seconds} s timeout")

        self._old_handler = signal.signal(signal.SIGALRM, _raise_timeout)
        signal.setitimer(signal.ITIMER_REAL, float(self.seconds))
        return self

    def __exit__(self, exc_type, exc, tb):
        if self.seconds > 0:
            signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, self._old_handler)
        return False


def _render_one_leo_report(payload: tuple[int, dict[str, Any], str, str, str, str, int, bool]) -> tuple[int, dict[str, Any], dict[str, Any]]:
    idx, row_dict, out_dir_s, hlsp_root_s, injection_h5_s, aperture, timeout_s, overwrite = payload
    _leo_cache_dirs()

    from leo_vetter.main import TCELightCurve
    from leo_vetter.plots import plot_summary
    from leo_vetter.wd_thresholds import check_thresholds_wd

    row = pd.Series(row_dict)
    out_dir = Path(out_dir_s)
    hlsp_root = Path(hlsp_root_s)
    injection_h5 = Path(injection_h5_s)
    reports_dir = out_dir / "vet_reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    rec: dict[str, Any] = {
        "review_id": row.get("review_id", ""),
        "tic": int(row.get("tic", -1)),
        "source_kind": row.get("source_kind", ""),
        "aperture_used": "",
        "error": "",
        "plot_error": "",
        "plot_sanitized": "",
    }
    queue_updates: dict[str, Any] = {}
    try:
        with _leo_timeout(timeout_s):
            row_aperture = _aperture_for_row(row, aperture)
            rec["aperture_used"] = row_aperture
            period = float(row["period_d"])
            t0 = float(row["t0_bjd"])
            duration_d = float(row["duration_min"]) / 1440.0
            if not all(np.isfinite([period, t0, duration_d])) or period <= 0 or duration_d <= 0:
                raise ValueError("invalid ephemeris")
            time, raw, flux, err = _leo_arrays_for_row(row, hlsp_root, injection_h5, row_aperture)
            if len(time) < 200:
                raise ValueError(f"too few LEO cadences: {len(time)}")
            tic = int(row["tic"])
            tlc = TCELightCurve(tic, time, raw, flux, err, period, t0, duration_d)
            star = _wd_star_for(row)
            sink = io.StringIO()
            with warnings.catch_warnings(), contextlib.redirect_stdout(sink):
                warnings.simplefilter("ignore", RuntimeWarning)
                tlc.compute_flux_metrics(star, cap_b=False)
                tlc.metrics["Rs"] = float(star["Rs"])
                fa = bool(check_thresholds_wd(tlc.metrics, "FA"))
                fp = bool(check_thresholds_wd(tlc.metrics, "FP"))
                label = _leo_label(fa, fp)
            report_name = (
                f"{label}_row{idx:05d}"
                f"_tic{tic:010d}_{row_aperture}_P{period:08.4f}d.pdf"
            )
            report_path = reports_dir / report_name
            if overwrite or not report_path.exists():
                with warnings.catch_warnings(), contextlib.redirect_stdout(sink):
                    warnings.simplefilter("ignore", RuntimeWarning)
                    plot_tlc, plot_sanitized = _leo_plot_copy(tlc)
                    rec["plot_sanitized"] = plot_sanitized
                    try:
                        plot_summary(plot_tlc, star, save_fig=True, save_file=str(report_path))
                    except Exception as plot_exc:
                        rec["plot_error"] = f"{type(plot_exc).__name__}: {plot_exc}"
                        _write_fallback_leo_report(
                            report_path,
                            tlc=tlc,
                            label=label,
                            aperture=row_aperture,
                            plot_error=rec["plot_error"],
                        )
            rec.update(dict(tlc.metrics))
            rec.update(
                {
                    "leo_FA": fa,
                    "leo_FP": fp,
                    "leo_PC": (not fa) and (not fp),
                    "leo_class": label,
                    "leo_report_name": report_name,
                    "leo_report_path": str(report_path),
                    "star_source": star.get("star_source", ""),
                    "star_source_id": star.get("source_id", ""),
                    "star_rad_rsun": star.get("rad", ""),
                    "star_mass_msun": star.get("mass", ""),
                    "star_teff_k": star.get("Teff", ""),
                    "star_rho_g_cm3": star.get("rho", ""),
                }
            )
            queue_updates = {"leo_class": label, "leo_report_name": report_name}
    except Exception as exc:
        rec["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        import matplotlib.pyplot as plt

        plt.close("all")
    return idx, rec, queue_updates


def render_leo_reports(
    queue: pd.DataFrame,
    *,
    out_dir: Path,
    hlsp_root: Path,
    injection_h5: Path,
    aperture: str,
    max_reports: int,
    timeout_s: int,
    overwrite: bool,
    workers: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    reports_dir = out_dir / "vet_reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    queue = queue.copy().reset_index(drop=True)
    n = len(queue) if max_reports <= 0 else min(max_reports, len(queue))
    workers = max(1, int(workers))
    print(
        f"[review] rendering LEO reports for {n:,} / {len(queue):,} rows "
        f"with {workers} worker(s)",
        flush=True,
    )
    payloads = [
        (
            idx,
            row.to_dict(),
            str(out_dir),
            str(hlsp_root),
            str(injection_h5),
            aperture,
            int(timeout_s),
            bool(overwrite),
        )
        for idx, (_, row) in enumerate(queue.head(n).iterrows(), start=1)
    ]
    rows_by_idx: dict[int, dict[str, Any]] = {}
    updates_by_idx: dict[int, dict[str, Any]] = {}
    if workers <= 1:
        iterator = (_render_one_leo_report(payload) for payload in payloads)
        for done, (idx, rec, updates) in enumerate(iterator, start=1):
            rows_by_idx[idx] = rec
            updates_by_idx[idx] = updates
            if done % 25 == 0:
                print(f"  [review] LEO {done:,}/{n:,}", flush=True)
    else:
        from concurrent.futures import ProcessPoolExecutor, as_completed

        with ProcessPoolExecutor(max_workers=workers) as ex:
            futures = [ex.submit(_render_one_leo_report, payload) for payload in payloads]
            for done, fut in enumerate(as_completed(futures), start=1):
                idx, rec, updates = fut.result()
                rows_by_idx[idx] = rec
                updates_by_idx[idx] = updates
                if done % 25 == 0:
                    print(f"  [review] LEO {done:,}/{n:,}", flush=True)
    for idx, updates in updates_by_idx.items():
        for col, value in updates.items():
            queue.loc[idx - 1, col] = value
    rows = [rows_by_idx[idx] for idx in sorted(rows_by_idx)]
    return queue, pd.DataFrame(rows)


def _finalize_queue(real: pd.DataFrame, injected: pd.DataFrame) -> pd.DataFrame:
    queue = pd.concat([real, injected], ignore_index=True, sort=False)
    if queue.empty:
        return queue
    for col in REVIEW_COLUMNS:
        if col not in queue.columns:
            queue[col] = ""
    for col in LABEL_COLUMNS:
        queue[col] = ""
    queue = queue.loc[:, list(REVIEW_COLUMNS) + list(LABEL_COLUMNS) + [c for c in queue.columns if c not in REVIEW_COLUMNS and c not in LABEL_COLUMNS]]
    if hasattr(queue, "map"):
        return queue.map(_clean)
    return queue.applymap(_clean)


def _blind_review_metadata(queue: pd.DataFrame) -> pd.DataFrame:
    """Hide provenance in columns displayed by the browser app.

    The truth/source columns remain available for post-label scoring, but the
    visible bucket no longer tells the reviewer whether a row is real or
    injected.
    """
    out = queue.copy()
    if "truth_source_bucket" not in out:
        out["truth_source_bucket"] = out.get("source_bucket", "")
    if "truth_source_kind" not in out:
        out["truth_source_kind"] = out.get("source_kind", "")
    if "source_bucket" in out:
        out["source_bucket"] = "review_candidate"
    if "vet_class" in out:
        out["vet_class"] = out["vet_class"].replace(
            to_replace=r"^injected_.*$",
            value="review_candidate",
            regex=True,
        )
    return out


def _build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--real-candidates", type=Path, default=DEFAULT_REAL_CANDIDATES)
    ap.add_argument("--hlsp-root", type=Path, default=DEFAULT_HLSP_ROOT)
    ap.add_argument("--injection-h5", type=Path, default=DEFAULT_INJECTION_H5)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    ap.add_argument(
        "--star-catalog",
        type=Path,
        default=None,
        help=(
            "Optional TWIRL/GF21 WD catalog with TIC matches and atmosphere-fit "
            "mass/logg/Teff columns. When provided, LEO receives per-target "
            "WD parameters with canonical fallback."
        ),
    )
    ap.add_argument(
        "--star-atmosphere-priority",
        default="H,He,mixed",
        help="Comma-separated atmosphere-fit priority for --star-catalog.",
    )
    ap.add_argument("--n-real", type=int, default=100)
    ap.add_argument("--n-injections", type=int, default=900)
    ap.add_argument(
        "--real-selection",
        choices=("legacy", "stratified_blind", "ranker_selected"),
        default="legacy",
        help=(
            "Real-candidate selection recipe. Use stratified_blind for the 10k mixed review set, "
            "or ranker_selected for ephemerides selected by the injected-truth peak ranker."
        ),
    )
    ap.add_argument(
        "--blind-review-metadata",
        action="store_true",
        help="Hide source-bucket and injected-class provenance in the browser-facing queue.",
    )
    ap.add_argument(
        "--shuffle-review-rows",
        action="store_true",
        help="Deterministically shuffle the final review queue before LEO rendering and CSV output.",
    )
    ap.add_argument(
        "--apertures",
        nargs="+",
        default=list(KNOWN_APERTURES),
        help="Apertures to use for injection BLS recovery. Real rows use their rep_aperture for LEO.",
    )
    ap.add_argument(
        "--aperture",
        default=None,
        help="Deprecated single-aperture alias and fallback LEO aperture.",
    )
    ap.add_argument("--workers", type=int, default=max(1, min(8, (os.cpu_count() or 2) // 2)))
    ap.add_argument("--n-periods", type=int, default=200_000)
    ap.add_argument("--random-state", type=int, default=56)
    ap.add_argument("--skip-leo", action="store_true")
    ap.add_argument("--reuse-injection-bls", action="store_true",
                    help="Reuse <out-dir>/injection_bls_recoveries.parquet instead of rerunning BLS.")
    ap.add_argument("--max-leo-reports", type=int, default=0,
                    help="Render at most N LEO reports. 0 means all rows.")
    ap.add_argument("--leo-timeout-s", type=int, default=300,
                    help="Per-row LEO timeout in seconds. Use 0 to disable.")
    ap.add_argument("--leo-workers", type=int, default=1,
                    help="Parallel LEO report workers. Keep at 1 for debugging.")
    ap.add_argument("--overwrite", action="store_true")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    apertures = tuple([args.aperture] if args.aperture else args.apertures)
    fallback_aperture = apertures[0]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.n_real > 0 and not args.real_candidates.exists():
        print(f"[review] missing real candidates: {args.real_candidates}", file=sys.stderr)
        return 2

    if args.real_selection == "ranker_selected":
        real = select_ranker_real_candidates(args.real_candidates, args.n_real, args.random_state)
    elif args.real_selection == "stratified_blind":
        real = select_real_candidates_stratified(args.real_candidates, args.n_real, args.random_state)
    else:
        real = select_real_candidates(args.real_candidates, args.n_real, args.random_state)
    injected = pd.DataFrame()
    if args.n_injections > 0:
        injection_bls_path = args.out_dir / "injection_bls_recoveries.parquet"
        existing_injection_bls_path = _existing_table_path(injection_bls_path)
        if args.reuse_injection_bls and existing_injection_bls_path is not None:
            injection_bls_path = existing_injection_bls_path
            injected = _read_table(injection_bls_path)
            if args.n_injections > 0:
                injected = injected.head(args.n_injections).copy()
            print(f"[review] reused {len(injected):,} injected BLS rows from {injection_bls_path}", flush=True)
        else:
            injected = run_injection_bls(
                injection_h5=args.injection_h5,
                apertures=apertures,
                n_injections=args.n_injections,
                workers=args.workers,
                n_periods=args.n_periods,
            )
            injection_bls_path = _write_table(injected, injection_bls_path)

    queue = _finalize_queue(real, injected)
    if args.shuffle_review_rows and len(queue) > 1:
        queue = queue.sample(frac=1.0, random_state=args.random_state).reset_index(drop=True)

    star_records: dict[int, dict[str, Any]] = {}
    atmosphere_priority = tuple(
        part.strip() for part in str(args.star_atmosphere_priority).split(",") if part.strip()
    )
    if args.star_catalog is not None:
        tics = {int(tic) for tic in queue["tic"].dropna().astype(int).tolist()}
        star_records = load_wd_star_catalog(
            args.star_catalog,
            tics=tics,
            atmosphere_priority=atmosphere_priority,
        )
        queue = annotate_star_parameters(queue, star_records=star_records)
        print(
            "[review] star catalog annotations: "
            f"{sum(1 for row in queue['star_source'].astype(str) if row != 'canonical_fallback'):,} "
            f"catalog-backed / {len(queue):,} rows",
            flush=True,
        )

    pre_leo_csv = args.out_dir / "review_queue_pre_leo.csv"
    queue.to_csv(pre_leo_csv, index=False)

    metrics = pd.DataFrame()
    leo_metrics_path = args.out_dir / "leo_metrics.parquet"
    if not args.skip_leo:
        queue, metrics = render_leo_reports(
            queue,
            out_dir=args.out_dir,
            hlsp_root=args.hlsp_root,
            injection_h5=args.injection_h5,
            aperture=fallback_aperture,
            max_reports=args.max_leo_reports,
            timeout_s=args.leo_timeout_s,
            overwrite=args.overwrite,
            workers=args.leo_workers,
        )
        leo_metrics_path = _write_table(metrics, leo_metrics_path)

    if args.blind_review_metadata:
        queue = _blind_review_metadata(queue)

    queue_csv = args.out_dir / "review_queue.csv"
    queue.to_csv(queue_csv, index=False)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "light_curve_product": "S56 TWIRL-FS v2 compare",
        "aperture": fallback_aperture,
        "apertures": list(apertures),
        "hlsp_root": str(args.hlsp_root),
        "real_candidates": str(args.real_candidates),
        "injection_h5": str(args.injection_h5),
        "n_real": int(len(real)),
        "n_injections": int(len(injected)),
        "n_review_rows": int(len(queue)),
        "real_selection": args.real_selection,
        "blind_review_metadata": bool(args.blind_review_metadata),
        "shuffle_review_rows": bool(args.shuffle_review_rows),
        "n_leo_reports_attempted": int(len(metrics)),
        "n_leo_errors": int(metrics["error"].astype(str).ne("").sum()) if "error" in metrics else 0,
        "n_leo_plot_errors": int(metrics["plot_error"].astype(str).ne("").sum()) if "plot_error" in metrics else 0,
        "leo_timeout_s": int(args.leo_timeout_s),
        "leo_workers": int(args.leo_workers),
        "star_catalog": str(args.star_catalog) if args.star_catalog is not None else "",
        "star_atmosphere_priority": list(atmosphere_priority),
        "star_source_counts": star_source_counts(queue),
        "bls_config": asdict(BLSConfig(apertures=apertures, n_periods=args.n_periods, n_peaks=5)),
        "outputs": {
            "review_queue_csv": str(queue_csv),
            "pre_leo_csv": str(pre_leo_csv),
            "leo_metrics": str(leo_metrics_path),
            "vet_reports": str(args.out_dir / "vet_reports"),
        },
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print("[review] complete")
    print(f"  rows: {len(queue):,}  real={len(real):,} injected={len(injected):,}")
    print(f"  queue: {queue_csv}")
    print(f"  reports: {args.out_dir / 'vet_reports'}")
    print(f"  summary: {args.out_dir / 'summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
