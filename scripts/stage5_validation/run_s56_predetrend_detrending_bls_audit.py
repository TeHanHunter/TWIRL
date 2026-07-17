#!/usr/bin/env python3
"""Audit raw-flux detrending strengths with injected-signal BLS recovery.

This is the production-branch comparison for S56 pre-detrend injections. It
starts from stored original/injected raw aperture flux, applies candidate
``FluxDetrendConfig`` methods, reruns BLS on the detrended injected curve, and
labels each BLS peak against the injected ephemeris.

This intentionally differs from the ADP+ vetting diagnostic: no candidate-window
masking is used for the deployable branches, and the detrending starts from raw
aperture flux rather than from an already detrended ADP light curve.
"""
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
SRC_ROOT = REPO_ROOT / "src"
SCRIPT_ROOT = REPO_ROOT / "scripts/stage5_validation"
for path in (SRC_ROOT, SCRIPT_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from twirl.io.hlsp import HLSPLightCurve  # noqa: E402
from twirl.lightcurves.detrend_presets import TWIRL_FS_V2_ADP015Q_BRANCH  # noqa: E402
from twirl.search.bls import BLSConfig, run_bls_on_lc  # noqa: E402
from build_injection_peak_training_table import (  # noqa: E402
    DEFAULT_HARMONIC_FACTORS,
    MIN_WINDOW_OVERLAP_FRACTION,
    PERIOD_RECOVERY_TOL,
    label_peak_against_injection,
)
from sweep_predetrend_detrending_methods import (  # noqa: E402
    METHOD_SPECS,
    binned_trend_metrics,
    detrend_method,
    robust_sigma,
)


DEFAULT_INJECTION_H5 = (
    REPO_ROOT
    / "data_local/stage3_injections/s56_twirlfs_v2_injection_training/"
    / "pdo_20k_predetrend_batman_periodradius_bright_balanced/injected_lightcurves.h5"
)
DEFAULT_OUT_DIR = REPO_ROOT / "reports/stage5_validation/s56_predetrend_detrending_bls_audit"
DEFAULT_METHODS = (
    "current_adp03q",
    TWIRL_FS_V2_ADP015Q_BRANCH,
    "spline_q02_gap02",
    "spline_q025_gap02",
    "spline_q05_gap02",
    "spline_q08_gap02",
    "current_uniform08_gap05",
)
DEFAULT_RAW_APERTURES = ("Small", "Primary")


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


def _as_int(value: Any, default: int = -1) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return int(default)


def _read_key_metadata(injection_h5: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    with h5py.File(injection_h5, "r") as h5:
        for key in sorted(h5["injections"].keys()):
            group = h5["injections"][key]
            attrs = group.attrs
            rows.append(
                {
                    "h5_key": str(key),
                    "h5_group": f"/injections/{key}",
                    "injection_id": str(attrs.get("injection_id", key)),
                    "tic": _as_int(attrs.get("tic")),
                    "tmag": _as_float(attrs.get("tessmag")),
                    "period_d": _as_float(attrs.get("period_d")),
                    "radius_rearth": _as_float(attrs.get("radius_rearth")),
                    "depth": _as_float(attrs.get("depth")),
                    "duration_min": _as_float(attrs.get("duration_min")),
                    "n_good_in_transit": _as_int(attrs.get("n_good_in_transit"), 0),
                }
            )
    return pd.DataFrame(rows)


def _default_manifest_path(injection_h5: Path) -> Path:
    return Path(injection_h5).with_name("injection_manifest.csv")


def _read_manifest_metadata(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "h5_group" not in df:
        if "injection_id" not in df:
            raise KeyError(f"manifest lacks h5_group and injection_id: {path}")
        df["h5_group"] = "/injections/" + df["injection_id"].astype(str)
    out = pd.DataFrame(
        {
            "h5_key": df["h5_group"].astype(str).str.rsplit("/", n=1).str[-1],
            "h5_group": df["h5_group"].astype(str),
            "injection_id": df["injection_id"].astype(str)
            if "injection_id" in df
            else df["h5_group"].astype(str).str.rsplit("/", n=1).str[-1],
            "tic": pd.to_numeric(df.get("tic", np.nan), errors="coerce"),
            "tmag": pd.to_numeric(df.get("tessmag", df.get("tmag", np.nan)), errors="coerce"),
            "period_d": pd.to_numeric(df.get("period_d", np.nan), errors="coerce"),
            "radius_rearth": pd.to_numeric(df.get("radius_rearth", np.nan), errors="coerce"),
            "depth": pd.to_numeric(df.get("depth", np.nan), errors="coerce"),
            "duration_min": pd.to_numeric(df.get("duration_min", np.nan), errors="coerce"),
            "n_good_in_transit": pd.to_numeric(df.get("n_good_in_transit", np.nan), errors="coerce"),
        }
    )
    return out


def _parse_key_file(path: Path | None) -> set[str] | None:
    if path is None:
        return None
    keys = set()
    for line in Path(path).read_text().splitlines():
        item = line.strip()
        if not item or item.startswith("#"):
            continue
        keys.add(item.split(",", 1)[0].strip())
    return keys


def _balanced_sample(meta: pd.DataFrame, *, n: int, seed: int) -> pd.DataFrame:
    if n <= 0 or len(meta) <= n:
        return meta.sort_values("h5_key").reset_index(drop=True)
    rng = np.random.default_rng(seed)
    df = meta.copy()
    df["_rand"] = rng.random(len(df))
    df["tmag_bin"] = pd.cut(
        pd.to_numeric(df["tmag"], errors="coerce"),
        bins=[-np.inf, 15.0, 17.0, 18.5, np.inf],
        labels=["tmag_lt15", "tmag_15_17", "tmag_17_18p5", "tmag_ge18p5"],
    ).astype(str)
    df["period_bin"] = pd.cut(
        np.log10(pd.to_numeric(df["period_d"], errors="coerce").clip(lower=1.0e-4)),
        bins=np.linspace(np.log10(0.08), np.log10(13.0), 6),
        include_lowest=True,
    ).astype(str)
    df["radius_bin"] = pd.cut(
        pd.to_numeric(df["radius_rearth"], errors="coerce"),
        bins=[-np.inf, 1.0, 2.0, 4.0, 8.0, 13.0, np.inf],
        labels=["r_lt1", "r_1_2", "r_2_4", "r_4_8", "r_8_13", "r_ge13"],
    ).astype(str)
    groups = list(df.groupby(["tmag_bin", "period_bin", "radius_bin"], dropna=False))
    quota = max(1, int(np.ceil(n / max(1, len(groups)))))
    chosen_parts: list[pd.DataFrame] = []
    for _, group in groups:
        chosen_parts.append(group.sort_values("_rand").head(quota))
    chosen = pd.concat(chosen_parts, ignore_index=True) if chosen_parts else df.head(0)
    if len(chosen) < n:
        remaining = df.loc[~df["h5_key"].isin(set(chosen["h5_key"]))].sort_values("_rand")
        chosen = pd.concat([chosen, remaining.head(n - len(chosen))], ignore_index=True)
    return chosen.sort_values("_rand").head(n).drop(columns=["_rand"]).reset_index(drop=True)


def select_sample(
    *,
    injection_h5: Path,
    injection_manifest: Path | None,
    n_injections: int,
    limit_keys_file: Path | None,
    random_state: int,
) -> pd.DataFrame:
    manifest_path = injection_manifest or _default_manifest_path(injection_h5)
    if manifest_path.exists():
        meta = _read_manifest_metadata(manifest_path)
    else:
        meta = _read_key_metadata(injection_h5)
    wanted = _parse_key_file(limit_keys_file)
    if wanted is not None:
        meta = meta[meta["h5_key"].isin(wanted) | meta["injection_id"].isin(wanted)].copy()
    if meta.empty:
        raise RuntimeError("no injections selected")
    return _balanced_sample(meta, n=n_injections, seed=random_state)


def _truth_record(group: h5py.Group, group_path: str, injection_h5: Path) -> dict[str, Any]:
    attrs = group.attrs
    return {
        "injection_id": str(attrs.get("injection_id", Path(group_path).name)),
        "source_h5": str(injection_h5),
        "h5_group": str(group_path),
        "tic": _as_int(attrs.get("tic")),
        "sector": _as_int(attrs.get("sector")),
        "cam": _as_int(attrs.get("camera")),
        "ccd": _as_int(attrs.get("ccd")),
        "tmag": _as_float(attrs.get("tessmag")),
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
        "truth_n_in_transit": _as_int(attrs.get("n_in_transit")),
        "truth_n_good_in_transit": _as_int(attrs.get("n_good_in_transit")),
        "truth_grid_cell_id": str(attrs.get("grid_cell_id", "")),
        "truth_grid_period_bin": str(attrs.get("grid_period_bin", "")),
        "truth_grid_depth_bin": str(attrs.get("grid_depth_bin", "")),
        "truth_grid_radius_bin": str(attrs.get("grid_radius_bin", "")),
    }


def _make_lc(
    *,
    group: h5py.Group,
    group_path: str,
    raw_aperture: str,
    method_name: str,
    flux: np.ndarray,
) -> HLSPLightCurve:
    attrs = group.attrs
    payload = {
        "tic": _as_int(attrs.get("tic")),
        "tmag": _as_float(attrs.get("tessmag")),
        "sector": _as_int(attrs.get("sector")),
        "cam": _as_int(attrs.get("camera")),
        "ccd": _as_int(attrs.get("ccd")),
        "ra": float("nan"),
        "dec": float("nan"),
        "time": np.asarray(group["time"], dtype=np.float64),
        "cadenceno": np.asarray(group["cadenceno"], dtype=np.int32)
        if "cadenceno" in group
        else np.arange(len(group["time"]), dtype=np.int32),
        "orbitid": np.asarray(group["orbitid"], dtype=np.int16),
        "quality": np.asarray(group["quality"], dtype=np.int32),
        "flux": {raw_aperture: np.asarray(flux, dtype=np.float64)},
        "path": Path(f"{group.file.filename}:{group_path}:{raw_aperture}:{method_name}"),
    }
    accepted = {field.name for field in fields(HLSPLightCurve)}
    return HLSPLightCurve(**{key: value for key, value in payload.items() if key in accepted})


def _depth_retention_and_snr(
    *,
    det_original: np.ndarray,
    det_injected: np.ndarray,
    quality: np.ndarray,
    in_transit: np.ndarray,
    truth_depth: float,
) -> dict[str, float]:
    finite = np.isfinite(det_original) & np.isfinite(det_injected)
    good = finite & (quality == 0)
    in_good = good & in_transit
    oot_good = good & ~in_transit
    if np.count_nonzero(in_good) < 1 or np.count_nonzero(oot_good) < 20:
        return {
            "depth_retention_frac": float("nan"),
            "delta_signal_multi_snr_mad": float("nan"),
            "delta_depth_abs": float("nan"),
        }
    delta = det_injected - det_original
    delta_in = float(np.nanmedian(delta[in_good]))
    delta_oot = float(np.nanmedian(delta[oot_good]))
    delta_depth = -(delta_in - delta_oot)
    retention = delta_depth / truth_depth if np.isfinite(truth_depth) and truth_depth > 0 else float("nan")
    sigma = robust_sigma(det_original[oot_good])
    per_cad_snr = delta_depth / sigma if np.isfinite(sigma) and sigma > 0 else float("nan")
    multi_snr = per_cad_snr * np.sqrt(np.count_nonzero(in_good)) if np.isfinite(per_cad_snr) else float("nan")
    return {
        "depth_retention_frac": retention,
        "delta_signal_multi_snr_mad": multi_snr,
        "delta_depth_abs": delta_depth,
    }


def _peak_records(result: Any) -> list[dict[str, Any]]:
    peaks = getattr(result, "peaks", None) or []
    if not peaks:
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
            "peak_rank": int(peak.peak_rank),
            "period_d": float(peak.period_d),
            "t0_bjd": float(peak.t0_bjd),
            "duration_min": float(peak.duration_min),
            "depth": float(peak.depth),
            "depth_snr": float(peak.depth_snr),
            "sde": float(peak.sde),
            "log_power": float(peak.log_power),
        }
        for peak in peaks
    ]


def process_payload(payload: tuple[str, str, tuple[str, ...], tuple[str, ...], dict[str, Any], float, float, tuple[float, ...]]) -> list[dict[str, Any]]:
    (
        injection_h5_s,
        group_path,
        method_names,
        raw_apertures,
        cfg_kwargs,
        period_tol,
        min_window_overlap_fraction,
        harmonic_factors,
    ) = payload
    injection_h5 = Path(injection_h5_s)
    cfg = BLSConfig(**cfg_kwargs)
    rows: list[dict[str, Any]] = []
    with h5py.File(injection_h5, "r") as h5:
        group = h5[group_path]
        truth = _truth_record(group, group_path, injection_h5)
        time = np.asarray(group["time"], dtype=np.float64)
        quality = np.asarray(group["quality"], dtype=np.int32)
        in_transit = np.asarray(group["in_transit"], dtype=bool) if "in_transit" in group else np.zeros(len(time), dtype=bool)
        for raw_aperture in raw_apertures:
            original_name = f"RAW_FLUX_{raw_aperture}_original"
            injected_name = f"RAW_FLUX_{raw_aperture}_injected"
            if original_name not in group or injected_name not in group:
                for method_name in method_names:
                    rows.append(
                        {
                            **truth,
                            "raw_aperture": raw_aperture,
                            "method": method_name,
                            "method_kind": "",
                            "status": f"missing_raw_aperture:{raw_aperture}",
                            "peak_rank": 0,
                            "is_candidate_peak": False,
                            "exact_ephemeris_match": False,
                            "harmonic_ephemeris_match": False,
                            "is_injected_signal_peak": False,
                            "match_kind": "missing_raw",
                        }
                    )
                continue
            original_raw = np.asarray(group[original_name], dtype=np.float64)
            injected_raw = np.asarray(group[injected_name], dtype=np.float64)
            for method_name in method_names:
                method = METHOD_SPECS[method_name]
                try:
                    det_original, _, meta = detrend_method(
                        time,
                        original_raw,
                        quality,
                        method,
                        in_transit=in_transit,
                    )
                    det_injected, _, _ = detrend_method(
                        time,
                        injected_raw,
                        quality,
                        method,
                        in_transit=in_transit,
                    )
                    lc = _make_lc(
                        group=group,
                        group_path=group_path,
                        raw_aperture=raw_aperture,
                        method_name=method_name,
                        flux=det_injected,
                    )
                    result = run_bls_on_lc(lc, cfg, aperture=raw_aperture)
                    depth_stats = _depth_retention_and_snr(
                        det_original=det_original,
                        det_injected=det_injected,
                        quality=quality,
                        in_transit=in_transit,
                        truth_depth=truth["truth_depth"],
                    )
                    trend = binned_trend_metrics(time, det_original, quality, bin_d=0.5)
                    for peak in _peak_records(result):
                        if int(peak["peak_rank"]) > 0:
                            label = label_peak_against_injection(
                                period_d=_as_float(peak["period_d"]),
                                t0_bjd=_as_float(peak["t0_bjd"]),
                                duration_min=_as_float(peak["duration_min"]),
                                truth_period_d=truth["truth_period_d"],
                                truth_t0_bjd=truth["truth_t0_bjd"],
                                truth_duration_min=truth["truth_duration_min"],
                                period_tol=period_tol,
                                min_window_overlap_fraction=min_window_overlap_fraction,
                                harmonic_factors=harmonic_factors,
                            )
                        else:
                            label = {
                                "exact_ephemeris_match": False,
                                "harmonic_ephemeris_match": False,
                                "is_injected_signal_peak": False,
                                "match_kind": "no_peak",
                            }
                        rows.append(
                            {
                                **truth,
                                "raw_aperture": raw_aperture,
                                "method": method_name,
                                "method_kind": meta.get("kind", ""),
                                "fit_mask": meta.get("fit_mask", method.fit_mask),
                                "branch_status": meta.get("status", ""),
                                "bkspace_d": meta.get("bkspace_d", np.nan),
                                "gap_split_d": meta.get("gap_split_d", np.nan),
                                "knot_strategy": meta.get("knot_strategy", ""),
                                "window_d": meta.get("window_d", np.nan),
                                "status": result.status,
                                "is_candidate_peak": int(peak["peak_rank"]) > 0,
                                "n_cad_total": result.n_cad_total,
                                "n_cad_quality": result.n_cad_quality,
                                "n_cad_kept": result.n_cad_kept,
                                "n_cad_sigma_clipped": result.n_cad_sigma_clipped,
                                "dropout_frac": result.dropout_frac,
                                "quality_dropout_frac": result.quality_dropout_frac,
                                "baseline_d": result.baseline_d,
                                "n_orbits": result.n_orbits,
                                "det_trend_ptp": trend["trend_ptp"],
                                "det_trend_sigma": trend["trend_sigma"],
                                **depth_stats,
                                **peak,
                                **label,
                            }
                        )
                except Exception as exc:
                    rows.append(
                        {
                            **truth,
                            "raw_aperture": raw_aperture,
                            "method": method_name,
                            "method_kind": method.kind,
                            "status": f"error:{type(exc).__name__}: {exc}",
                            "peak_rank": 0,
                            "is_candidate_peak": False,
                            "exact_ephemeris_match": False,
                            "harmonic_ephemeris_match": False,
                            "is_injected_signal_peak": False,
                            "match_kind": "error",
                        }
                    )
    return rows


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    for (method, raw_aperture), group in df.groupby(["method", "raw_aperture"], dropna=False):
        objects = group.groupby("injection_id", dropna=False)
        n_objects = int(objects.ngroups)
        top1 = group[group["peak_rank"].fillna(0).astype(int).eq(1)].copy()
        top1_exact_ids = set(top1.loc[top1["exact_ephemeris_match"].fillna(False).astype(bool), "injection_id"])
        exact_any = objects["exact_ephemeris_match"].any() if "exact_ephemeris_match" in group else pd.Series(dtype=bool)
        harmonic_any = objects["harmonic_ephemeris_match"].any() if "harmonic_ephemeris_match" in group else pd.Series(dtype=bool)
        status_counts = {
            str(key): int(value)
            for key, value in group["status"].fillna("").astype(str).value_counts().sort_index().items()
        }
        per_injection = group.drop_duplicates(["injection_id", "method", "raw_aperture"])
        rows.append(
            {
                "method": method,
                "raw_aperture": raw_aperture,
                "n_injections": n_objects,
                "top1_exact_n": int(len(top1_exact_ids)),
                "top1_exact_frac": float(len(top1_exact_ids) / n_objects) if n_objects else np.nan,
                "topn_exact_n": int(exact_any.sum()) if len(exact_any) else 0,
                "topn_exact_frac": float(exact_any.sum() / n_objects) if n_objects else np.nan,
                "topn_exact_or_harmonic_n": int((exact_any | harmonic_any).sum()) if len(exact_any) else 0,
                "topn_exact_or_harmonic_frac": float((exact_any | harmonic_any).sum() / n_objects) if n_objects else np.nan,
                "median_depth_retention": float(np.nanmedian(pd.to_numeric(per_injection["depth_retention_frac"], errors="coerce"))),
                "median_delta_signal_multi_snr_mad": float(np.nanmedian(pd.to_numeric(per_injection["delta_signal_multi_snr_mad"], errors="coerce"))),
                "median_det_trend_ptp": float(np.nanmedian(pd.to_numeric(per_injection["det_trend_ptp"], errors="coerce"))),
                "median_sde_rank1": float(np.nanmedian(pd.to_numeric(top1["sde"], errors="coerce"))) if len(top1) else np.nan,
                "ok_peak_rows": int(group["status"].fillna("").astype(str).eq("ok").sum()),
                "status_counts_json": json.dumps(status_counts, sort_keys=True),
            }
        )
    out = pd.DataFrame(rows)
    return out.sort_values(
        ["top1_exact_frac", "topn_exact_or_harmonic_frac", "median_det_trend_ptp"],
        ascending=[False, False, True],
    )


def _frame_to_markdown(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    cols = [
        "method",
        "raw_aperture",
        "n_injections",
        "top1_exact_frac",
        "topn_exact_or_harmonic_frac",
        "median_depth_retention",
        "median_det_trend_ptp",
        "median_sde_rank1",
    ]
    view = df[[col for col in cols if col in df.columns]].copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            view[col] = view[col].map(lambda value: f"{value:.4g}" if np.isfinite(value) else "")
    headers = [str(c) for c in view.columns]
    body = [[str(value) for value in row] for row in view.itertuples(index=False, name=None)]
    widths = [max(len(headers[i]), *(len(row[i]) for row in body)) for i in range(len(headers))]
    lines = [
        "| " + " | ".join(headers[i].ljust(widths[i]) for i in range(len(headers))) + " |",
        "| " + " | ".join("-" * widths[i] for i in range(len(headers))) + " |",
    ]
    lines.extend("| " + " | ".join(row[i].ljust(widths[i]) for i in range(len(headers))) + " |" for row in body)
    return "\n".join(lines)


def write_report(summary_df: pd.DataFrame, out_path: Path) -> None:
    lines = [
        "# Raw-Flux Detrending-Strength BLS Audit",
        "",
        "This audit starts from stored raw aperture flux, applies each detrending method, reruns BLS, and labels peaks against injected BATMAN truth.",
        "",
        "## Branch Summary",
        "",
        _frame_to_markdown(summary_df),
        "",
        "## Notes",
        "",
        "- `top1_exact_frac` is the strict deployable recovery rate: BLS rank 1 matches the injected period and transit window.",
        "- `topn_exact_or_harmonic_frac` includes exact and accepted harmonic matches among all retained peaks.",
        "- `median_depth_retention` measures injected-minus-original depth after detrending; values near 1 preserve the signal.",
        "- This script does not use transit-window masking for deployable branches.",
        "",
    ]
    out_path.write_text("\n".join(lines), encoding="utf-8")


def run_audit(args: argparse.Namespace) -> dict[str, Any]:
    methods = tuple(args.methods)
    for method in methods:
        if method not in METHOD_SPECS:
            raise KeyError(f"unknown detrending method: {method}")
    raw_apertures = tuple(args.raw_apertures)
    sample = select_sample(
        injection_h5=args.injection_h5,
        injection_manifest=args.injection_manifest,
        n_injections=args.n_injections,
        limit_keys_file=args.limit_keys_file,
        random_state=args.random_state,
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    sample.to_csv(args.out_dir / "audit_sample_injected.csv", index=False)
    cfg = BLSConfig(
        apertures=raw_apertures,
        n_periods=int(args.n_periods),
        n_peaks=int(args.n_peaks),
        durations_min=tuple(float(part) for part in args.durations_min.split(",") if part.strip()),
        p_min_d=float(args.p_min_d),
        max_period_fraction=float(args.max_period_fraction),
        p_max_cap_d=float(args.p_max_cap_d),
        sigma_clip=float(args.sigma_clip),
        min_cadences=int(args.min_cadences),
        period_mask_frac=float(args.period_mask_frac),
    )
    cfg_kwargs = {
        "apertures": cfg.apertures,
        "p_min_d": cfg.p_min_d,
        "max_period_fraction": cfg.max_period_fraction,
        "p_max_cap_d": cfg.p_max_cap_d,
        "durations_min": cfg.durations_min,
        "n_periods": cfg.n_periods,
        "n_peaks": cfg.n_peaks,
        "period_mask_frac": cfg.period_mask_frac,
        "period_bin_edges": cfg.period_bin_edges,
        "max_peaks_per_period_bin": cfg.max_peaks_per_period_bin,
        "min_cadences": cfg.min_cadences,
        "sigma_clip": cfg.sigma_clip,
        "orbit_edge_trim_d": cfg.orbit_edge_trim_d,
    }
    harmonic_factors = tuple(float(part) for part in args.harmonic_factors.split(",") if part.strip())
    payloads = [
        (
            str(args.injection_h5),
            str(row["h5_group"]),
            methods,
            raw_apertures,
            cfg_kwargs,
            float(args.period_tol),
            float(args.min_window_overlap_fraction),
            harmonic_factors,
        )
        for row in sample.to_dict("records")
    ]
    print(
        f"[predetrend-bls-audit] injections={len(payloads):,} methods={len(methods)} "
        f"raw_apertures={len(raw_apertures)} n_periods={args.n_periods}",
        flush=True,
    )
    rows: list[dict[str, Any]] = []
    workers = max(1, int(args.workers))
    if workers <= 1:
        for idx, payload in enumerate(payloads, start=1):
            rows.extend(process_payload(payload))
            if idx % 25 == 0:
                print(f"[predetrend-bls-audit] processed {idx:,}/{len(payloads):,}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for idx, batch in enumerate(ex.map(process_payload, payloads, chunksize=1), start=1):
                rows.extend(batch)
                if idx % 25 == 0:
                    print(f"[predetrend-bls-audit] processed {idx:,}/{len(payloads):,}", flush=True)
    peak_df = pd.DataFrame(rows)
    peak_table = args.out_dir / "detrending_bls_peak_table.csv"
    peak_df.to_csv(peak_table, index=False)
    summary_df = summarize(peak_df)
    summary_table = args.out_dir / "detrending_bls_branch_summary.csv"
    summary_df.to_csv(summary_table, index=False)
    report_md = args.out_dir / "summary.md"
    write_report(summary_df, report_md)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "injection_h5": str(args.injection_h5),
        "injection_manifest": str(args.injection_manifest or _default_manifest_path(args.injection_h5)),
        "out_dir": str(args.out_dir),
        "n_injections": int(len(sample)),
        "methods": list(methods),
        "raw_apertures": list(raw_apertures),
        "n_periods": int(args.n_periods),
        "n_peaks": int(args.n_peaks),
        "workers": int(args.workers),
        "outputs": {
            "audit_sample_injected": str(args.out_dir / "audit_sample_injected.csv"),
            "peak_table": str(peak_table),
            "branch_summary": str(summary_table),
            "summary_md": str(report_md),
            "summary_json": str(args.out_dir / "summary.json"),
        },
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--injection-h5", type=Path, default=DEFAULT_INJECTION_H5)
    ap.add_argument(
        "--injection-manifest",
        type=Path,
        default=None,
        help="Optional injection_manifest.csv. Defaults to the HDF5 sibling file when present.",
    )
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    ap.add_argument("--n-injections", type=int, default=3000, help="0 means all injections.")
    ap.add_argument("--method", dest="methods", action="append", default=None)
    ap.add_argument("--raw-aperture", dest="raw_apertures", action="append", default=None)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--n-periods", type=int, default=50_000)
    ap.add_argument("--n-peaks", type=int, default=20)
    ap.add_argument("--durations-min", default=",".join(str(v) for v in BLSConfig.durations_min))
    ap.add_argument("--p-min-d", type=float, default=BLSConfig.p_min_d)
    ap.add_argument("--max-period-fraction", type=float, default=BLSConfig.max_period_fraction)
    ap.add_argument("--p-max-cap-d", type=float, default=BLSConfig.p_max_cap_d)
    ap.add_argument("--sigma-clip", type=float, default=BLSConfig.sigma_clip)
    ap.add_argument("--min-cadences", type=int, default=BLSConfig.min_cadences)
    ap.add_argument("--period-mask-frac", type=float, default=BLSConfig.period_mask_frac)
    ap.add_argument("--period-tol", type=float, default=PERIOD_RECOVERY_TOL)
    ap.add_argument("--min-window-overlap-fraction", type=float, default=MIN_WINDOW_OVERLAP_FRACTION)
    ap.add_argument("--harmonic-factors", default=",".join(str(v) for v in DEFAULT_HARMONIC_FACTORS))
    ap.add_argument("--limit-keys-file", type=Path, default=None)
    ap.add_argument("--random-state", type=int, default=5612)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    args.methods = tuple(args.methods) if args.methods else DEFAULT_METHODS
    args.raw_apertures = tuple(args.raw_apertures) if args.raw_apertures else DEFAULT_RAW_APERTURES
    summary = run_audit(args)
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
