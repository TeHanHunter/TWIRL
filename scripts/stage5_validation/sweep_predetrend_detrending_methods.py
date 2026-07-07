#!/usr/bin/env python3
"""Compare detrending methods on stored pre-detrend injected raw flux.

This diagnostic reuses the raw original/injected aperture-flux arrays in the
pre-detrend injection HDF5. For each method it detrends the original and
injected raw light curves, then measures two quantities:

1. injected-signal survival at the truth in-transit cadences;
2. residual low-frequency trend in the original detrended curve.

The goal is to find methods that sit on the useful Pareto front: high injected
signal SNR / depth retention while still suppressing most slow trends.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any

import h5py
import numpy as np
import pandas as pd
from scipy.ndimage import median_filter, percentile_filter, uniform_filter1d
from scipy.signal import savgol_filter

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.lightcurves.flux_detrend import FluxDetrendConfig, flux_space_detrend_result  # noqa: E402
from twirl.lightcurves.detrend_presets import (  # noqa: E402
    TWIRL_FS_V2_ADP015Q_BRANCH,
    adp015q_config,
)
from twirl.plotting.style import apply_twirl_style, get_ordered_palette  # noqa: E402


RAW_APERTURES = ("Small", "Primary", "Large")
DEFAULT_METHODS = (
    "constant",
    "poly2_gap05",
    "poly3_gap05",
    "median015_gap05",
    "median020_gap05",
    "median025_gap05",
    "median03_gap05",
    "median04_gap05",
    "median05_gap05",
    "median06_gap05",
    "median10_gap05",
    "median20_gap05",
    "savgol03_p2_gap05",
    "savgol06_p2_gap05",
    "savgol10_p2_gap05",
    "savgol20_p2_gap05",
    "pctl75_03_gap05",
    "pctl60_06_gap05",
    "pctl75_06_gap05",
    "pctl90_06_gap05",
    "pctl75_10_gap05",
    "pctl75_20_gap05",
    "current_adp03q",
    TWIRL_FS_V2_ADP015Q_BRANCH,
    "spline_q02_gap02",
    "spline_q025_gap02",
    "spline_q05_gap02",
    "spline_q08_gap02",
    "current_uniform08_gap05",
    "spline_q12_gap05",
    "spline_uniform12_gap05",
    "spline_uniform20_gap05",
    "oracle_adp03q",
    "oracle_q05_gap02",
)


@dataclass(frozen=True)
class MethodSpec:
    name: str
    kind: str
    cfg: FluxDetrendConfig | None = None
    degree: int = 2
    gap_split_d: float = 0.5
    window_d: float = 1.0
    polyorder: int = 2
    percentile: float = 75.0
    smooth_d: float = 0.0
    fit_mask: str = "quality"


METHOD_SPECS: dict[str, MethodSpec] = {
    "constant": MethodSpec(name="constant", kind="constant"),
    "poly2_gap05": MethodSpec(name="poly2_gap05", kind="poly", degree=2, gap_split_d=0.5),
    "poly3_gap05": MethodSpec(name="poly3_gap05", kind="poly", degree=3, gap_split_d=0.5),
    "median015_gap05": MethodSpec(
        name="median015_gap05",
        kind="median_filter",
        window_d=0.15,
        gap_split_d=0.5,
    ),
    "median020_gap05": MethodSpec(
        name="median020_gap05",
        kind="median_filter",
        window_d=0.20,
        gap_split_d=0.5,
    ),
    "median025_gap05": MethodSpec(
        name="median025_gap05",
        kind="median_filter",
        window_d=0.25,
        gap_split_d=0.5,
    ),
    "median03_gap05": MethodSpec(
        name="median03_gap05",
        kind="median_filter",
        window_d=0.3,
        gap_split_d=0.5,
    ),
    "median04_gap05": MethodSpec(
        name="median04_gap05",
        kind="median_filter",
        window_d=0.4,
        gap_split_d=0.5,
    ),
    "median05_gap05": MethodSpec(
        name="median05_gap05",
        kind="median_filter",
        window_d=0.5,
        gap_split_d=0.5,
    ),
    "median06_gap05": MethodSpec(
        name="median06_gap05",
        kind="median_filter",
        window_d=0.6,
        gap_split_d=0.5,
    ),
    "median10_gap05": MethodSpec(
        name="median10_gap05",
        kind="median_filter",
        window_d=1.0,
        gap_split_d=0.5,
    ),
    "median20_gap05": MethodSpec(
        name="median20_gap05",
        kind="median_filter",
        window_d=2.0,
        gap_split_d=0.5,
    ),
    "savgol03_p2_gap05": MethodSpec(
        name="savgol03_p2_gap05",
        kind="savgol",
        window_d=0.3,
        polyorder=2,
        gap_split_d=0.5,
    ),
    "savgol06_p2_gap05": MethodSpec(
        name="savgol06_p2_gap05",
        kind="savgol",
        window_d=0.6,
        polyorder=2,
        gap_split_d=0.5,
    ),
    "savgol10_p2_gap05": MethodSpec(
        name="savgol10_p2_gap05",
        kind="savgol",
        window_d=1.0,
        polyorder=2,
        gap_split_d=0.5,
    ),
    "savgol20_p2_gap05": MethodSpec(
        name="savgol20_p2_gap05",
        kind="savgol",
        window_d=2.0,
        polyorder=2,
        gap_split_d=0.5,
    ),
    "pctl75_03_gap05": MethodSpec(
        name="pctl75_03_gap05",
        kind="percentile_filter",
        window_d=0.3,
        percentile=75.0,
        gap_split_d=0.5,
    ),
    "pctl60_06_gap05": MethodSpec(
        name="pctl60_06_gap05",
        kind="percentile_filter",
        window_d=0.6,
        percentile=60.0,
        gap_split_d=0.5,
    ),
    "pctl75_06_gap05": MethodSpec(
        name="pctl75_06_gap05",
        kind="percentile_filter",
        window_d=0.6,
        percentile=75.0,
        gap_split_d=0.5,
    ),
    "pctl90_06_gap05": MethodSpec(
        name="pctl90_06_gap05",
        kind="percentile_filter",
        window_d=0.6,
        percentile=90.0,
        gap_split_d=0.5,
    ),
    "pctl75_10_gap05": MethodSpec(
        name="pctl75_10_gap05",
        kind="percentile_filter",
        window_d=1.0,
        percentile=75.0,
        gap_split_d=0.5,
    ),
    "pctl75_20_gap05": MethodSpec(
        name="pctl75_20_gap05",
        kind="percentile_filter",
        window_d=2.0,
        percentile=75.0,
        gap_split_d=0.5,
    ),
    "current_adp03q": MethodSpec(
        name="current_adp03q",
        kind="spline",
        cfg=FluxDetrendConfig(
            bkspace_d=0.3,
            output_mode="subtractive",
            scale_strategy="auto",
            gap_split_d=0.2,
            knot_strategy="quantile",
        ),
    ),
    "spline_q015_gap02": MethodSpec(
        name="spline_q015_gap02",
        kind="spline",
        cfg=FluxDetrendConfig(
            bkspace_d=0.15,
            output_mode="subtractive",
            scale_strategy="auto",
            gap_split_d=0.2,
            knot_strategy="quantile",
        ),
    ),
    TWIRL_FS_V2_ADP015Q_BRANCH: MethodSpec(
        name=TWIRL_FS_V2_ADP015Q_BRANCH,
        kind="spline",
        cfg=adp015q_config(),
    ),
    "spline_q02_gap02": MethodSpec(
        name="spline_q02_gap02",
        kind="spline",
        cfg=FluxDetrendConfig(
            bkspace_d=0.2,
            output_mode="subtractive",
            scale_strategy="auto",
            gap_split_d=0.2,
            knot_strategy="quantile",
        ),
    ),
    "spline_q025_gap02": MethodSpec(
        name="spline_q025_gap02",
        kind="spline",
        cfg=FluxDetrendConfig(
            bkspace_d=0.25,
            output_mode="subtractive",
            scale_strategy="auto",
            gap_split_d=0.2,
            knot_strategy="quantile",
        ),
    ),
    "spline_q05_gap02": MethodSpec(
        name="spline_q05_gap02",
        kind="spline",
        cfg=FluxDetrendConfig(
            bkspace_d=0.5,
            output_mode="subtractive",
            scale_strategy="auto",
            gap_split_d=0.2,
            knot_strategy="quantile",
        ),
    ),
    "spline_q08_gap02": MethodSpec(
        name="spline_q08_gap02",
        kind="spline",
        cfg=FluxDetrendConfig(
            bkspace_d=0.8,
            output_mode="subtractive",
            scale_strategy="auto",
            gap_split_d=0.2,
            knot_strategy="quantile",
        ),
    ),
    "current_uniform08_gap05": MethodSpec(
        name="current_uniform08_gap05",
        kind="spline",
        cfg=FluxDetrendConfig(
            bkspace_d=0.8,
            output_mode="subtractive",
            scale_strategy="auto",
            gap_split_d=0.5,
            knot_strategy="uniform",
        ),
    ),
    "spline_q12_gap05": MethodSpec(
        name="spline_q12_gap05",
        kind="spline",
        cfg=FluxDetrendConfig(
            bkspace_d=1.2,
            output_mode="subtractive",
            scale_strategy="auto",
            gap_split_d=0.5,
            knot_strategy="quantile",
        ),
    ),
    "spline_uniform12_gap05": MethodSpec(
        name="spline_uniform12_gap05",
        kind="spline",
        cfg=FluxDetrendConfig(
            bkspace_d=1.2,
            output_mode="subtractive",
            scale_strategy="auto",
            gap_split_d=0.5,
            knot_strategy="uniform",
        ),
    ),
    "spline_uniform20_gap05": MethodSpec(
        name="spline_uniform20_gap05",
        kind="spline",
        cfg=FluxDetrendConfig(
            bkspace_d=2.0,
            output_mode="subtractive",
            scale_strategy="auto",
            gap_split_d=0.5,
            knot_strategy="uniform",
        ),
    ),
    "oracle_adp03q": MethodSpec(
        name="oracle_adp03q",
        kind="spline",
        fit_mask="oracle_nontransit",
        cfg=FluxDetrendConfig(
            bkspace_d=0.3,
            output_mode="subtractive",
            scale_strategy="auto",
            gap_split_d=0.2,
            knot_strategy="quantile",
        ),
    ),
    "oracle_q05_gap02": MethodSpec(
        name="oracle_q05_gap02",
        kind="spline",
        fit_mask="oracle_nontransit",
        cfg=FluxDetrendConfig(
            bkspace_d=0.5,
            output_mode="subtractive",
            scale_strategy="auto",
            gap_split_d=0.2,
            knot_strategy="quantile",
        ),
    ),
}


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def robust_sigma(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return float("nan")
    med = float(np.nanmedian(finite))
    mad = float(np.nanmedian(np.abs(finite - med)))
    sigma = 1.4826 * mad
    if np.isfinite(sigma) and sigma > 0:
        return sigma
    sigma = float(np.nanstd(finite))
    return sigma if np.isfinite(sigma) else float("nan")


def robust_abs_scale(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return float("nan")
    scale = float(np.nanpercentile(np.abs(finite), 68.0))
    if np.isfinite(scale) and scale > 0:
        return scale
    scale = float(np.nanmedian(np.abs(finite)))
    return scale if np.isfinite(scale) and scale > 0 else float("nan")


def load_injection_ids(path: Path | None) -> tuple[str, ...] | None:
    if path is None:
        return None
    out: list[str] = []
    seen: set[str] = set()
    for line in path.read_text().splitlines():
        item = line.strip()
        if not item or item.startswith("#"):
            continue
        item = item.split(",", 1)[0].strip()
        if item and item not in seen:
            out.append(item)
            seen.add(item)
    return tuple(out)


def time_gap_segments(time: np.ndarray, gap_split_d: float) -> list[np.ndarray]:
    finite = np.isfinite(time)
    if finite.sum() <= 1 or not np.isfinite(gap_split_d) or gap_split_d <= 0:
        return [np.arange(len(time), dtype=int)]
    order = np.argsort(np.where(finite, time, np.inf), kind="stable")
    finite_order = order[finite[order]]
    cuts = np.where(np.diff(time[finite_order]) > gap_split_d)[0] + 1
    if cuts.size == 0:
        return [np.arange(len(time), dtype=int)]
    starts = np.concatenate(([0], cuts))
    ends = np.concatenate((cuts, [finite_order.size]))
    return [np.sort(finite_order[start:end]) for start, end in zip(starts, ends)]


def polynomial_cotrend(
    time: np.ndarray,
    flux: np.ndarray,
    quality: np.ndarray,
    *,
    degree: int,
    gap_split_d: float,
) -> tuple[np.ndarray, str]:
    cotrend = np.full(len(time), np.nan, dtype=np.float64)
    statuses: list[str] = []
    for idx in time_gap_segments(time, gap_split_d):
        fit = np.isfinite(time[idx]) & np.isfinite(flux[idx]) & (quality[idx] == 0)
        if np.count_nonzero(fit) < max(20, degree + 3):
            med = np.nanmedian(flux[idx][fit]) if np.any(fit) else np.nan
            cotrend[idx] = med
            statuses.append("constant_insufficient_fit")
            continue
        t_fit = time[idx][fit]
        f_fit = flux[idx][fit]
        center = float(np.nanmedian(t_fit))
        scale_t = float(np.nanmax(np.abs(t_fit - center)))
        if not np.isfinite(scale_t) or scale_t <= 0:
            scale_t = 1.0
        x_fit = (t_fit - center) / scale_t
        x_all = (time[idx] - center) / scale_t
        deg = min(degree, max(0, np.count_nonzero(fit) - 1))
        try:
            coeffs = np.polyfit(x_fit, f_fit, deg=deg)
            cotrend[idx] = np.polyval(coeffs, np.clip(x_all, np.nanmin(x_fit), np.nanmax(x_fit)))
            statuses.append(f"poly{deg}")
        except Exception:
            med = np.nanmedian(f_fit)
            cotrend[idx] = med
            statuses.append("constant_poly_failed")
    return cotrend, ";".join(statuses)


def effective_quality(
    quality: np.ndarray,
    method: MethodSpec,
    *,
    in_transit: np.ndarray | None,
) -> np.ndarray:
    out = np.asarray(quality, dtype=np.int32).copy()
    if method.fit_mask == "quality":
        return out
    if method.fit_mask == "oracle_nontransit":
        if in_transit is None:
            raise ValueError(f"{method.name} requires in_transit mask")
        out[np.asarray(in_transit, dtype=bool)] = 1
        return out
    raise ValueError(f"unknown fit mask mode for {method.name}: {method.fit_mask!r}")


def odd_window_from_days(time: np.ndarray, window_d: float, *, min_window: int) -> int:
    finite_time = np.asarray(time, dtype=np.float64)
    finite_time = finite_time[np.isfinite(finite_time)]
    if finite_time.size < 2:
        return min_window | 1
    dt = np.diff(np.sort(finite_time))
    dt = dt[np.isfinite(dt) & (dt > 0)]
    cadence_d = float(np.nanmedian(dt)) if dt.size else np.nan
    if not np.isfinite(cadence_d) or cadence_d <= 0:
        return min_window | 1
    width = int(np.ceil(window_d / cadence_d))
    width = max(min_window, width)
    if width % 2 == 0:
        width += 1
    return width


def interpolate_fit_values(time: np.ndarray, flux: np.ndarray, fit: np.ndarray) -> np.ndarray:
    out = np.full(len(flux), np.nan, dtype=np.float64)
    finite_time = np.isfinite(time)
    good = finite_time & np.isfinite(flux) & fit
    if np.count_nonzero(good) == 0:
        fill = float(np.nanmedian(flux[np.isfinite(flux)])) if np.isfinite(flux).any() else np.nan
        out[finite_time] = fill
        return out
    if np.count_nonzero(good) == 1:
        out[finite_time] = float(flux[good][0])
        return out
    order_good = np.argsort(time[good], kind="stable")
    t_good = time[good][order_good]
    f_good = flux[good][order_good]
    order_all = np.argsort(time[finite_time], kind="stable")
    t_all = time[finite_time][order_all]
    interp = np.interp(t_all, t_good, f_good, left=f_good[0], right=f_good[-1])
    finite_idx = np.where(finite_time)[0][order_all]
    out[finite_idx] = interp
    return out


def filter_cotrend(
    time: np.ndarray,
    flux: np.ndarray,
    quality: np.ndarray,
    method: MethodSpec,
) -> tuple[np.ndarray, str]:
    cotrend = np.full(len(flux), np.nan, dtype=np.float64)
    statuses: list[str] = []
    for idx in time_gap_segments(time, method.gap_split_d):
        segment_time = time[idx]
        finite_segment = np.isfinite(segment_time)
        if np.count_nonzero(finite_segment) < 3:
            med = float(np.nanmedian(flux[idx][np.isfinite(flux[idx])])) if np.isfinite(flux[idx]).any() else np.nan
            cotrend[idx] = med
            statuses.append("constant_short_segment")
            continue
        order = np.argsort(np.where(finite_segment, segment_time, np.inf), kind="stable")
        ordered_idx = idx[order]
        finite_ordered_idx = ordered_idx[np.isfinite(time[ordered_idx])]
        t = time[finite_ordered_idx]
        f = flux[finite_ordered_idx]
        q = quality[finite_ordered_idx]
        fit = np.isfinite(t) & np.isfinite(f) & (q == 0)
        if np.count_nonzero(fit) < 20:
            med = float(np.nanmedian(f[fit])) if np.any(fit) else np.nan
            cotrend[finite_ordered_idx] = med
            statuses.append("constant_insufficient_fit")
            continue

        filled = interpolate_fit_values(t, f, fit)
        if method.kind == "median_filter":
            width = odd_window_from_days(t, method.window_d, min_window=5)
            width = min(width, len(filled) if len(filled) % 2 == 1 else max(1, len(filled) - 1))
            if width < 3:
                trend = np.full_like(filled, np.nanmedian(f[fit]), dtype=np.float64)
                status = "constant_short_filter"
            else:
                trend = median_filter(filled, size=width, mode="nearest")
                if method.smooth_d > 0:
                    smooth_width = odd_window_from_days(t, method.smooth_d, min_window=3)
                    trend = uniform_filter1d(trend, size=smooth_width, mode="nearest")
                status = f"median_filter_window{width}"
        elif method.kind == "percentile_filter":
            width = odd_window_from_days(t, method.window_d, min_window=5)
            width = min(width, len(filled) if len(filled) % 2 == 1 else max(1, len(filled) - 1))
            if width < 3:
                trend = np.full_like(filled, np.nanmedian(f[fit]), dtype=np.float64)
                status = "constant_short_percentile"
            else:
                trend = percentile_filter(
                    filled,
                    percentile=float(method.percentile),
                    size=width,
                    mode="nearest",
                )
                status = f"percentile_filter_window{width}_p{method.percentile:g}"
        elif method.kind == "savgol":
            min_width = method.polyorder + 3
            if min_width % 2 == 0:
                min_width += 1
            width = odd_window_from_days(t, method.window_d, min_window=min_width)
            width = min(width, len(filled) if len(filled) % 2 == 1 else max(1, len(filled) - 1))
            if width <= method.polyorder or width < 3:
                trend = np.full_like(filled, np.nanmedian(f[fit]), dtype=np.float64)
                status = "constant_short_savgol"
            else:
                polyorder = min(method.polyorder, width - 1)
                trend = savgol_filter(filled, window_length=width, polyorder=polyorder, mode="interp")
                status = f"savgol_window{width}_poly{polyorder}"
        else:
            raise ValueError(f"filter_cotrend cannot handle {method.kind}")

        cotrend[finite_ordered_idx] = trend
        statuses.append(status)
    if np.isnan(cotrend).any():
        med = float(np.nanmedian(cotrend[np.isfinite(cotrend)])) if np.isfinite(cotrend).any() else np.nan
        cotrend[~np.isfinite(cotrend)] = med
    return cotrend, ";".join(statuses)


def relative_subtractive(
    flux: np.ndarray,
    cotrend: np.ndarray,
    quality: np.ndarray,
) -> tuple[np.ndarray, float, str]:
    fit = np.isfinite(flux) & np.isfinite(cotrend) & (quality == 0)
    scale = robust_abs_scale(flux[fit])
    scale_source = "robust_abs"
    med_abs = abs(float(np.nanmedian(flux[fit]))) if np.any(fit) else np.nan
    sigma = robust_sigma(flux[fit])
    if np.isfinite(med_abs) and np.isfinite(sigma) and med_abs >= max(3.0 * sigma, 1.0e-12):
        scale = med_abs
        scale_source = "median_abs"
    out = np.full_like(flux, np.nan, dtype=np.float64)
    finite = np.isfinite(flux) & np.isfinite(cotrend) & np.isfinite(scale) & (scale > 0)
    out[finite] = 1.0 + (flux[finite] - cotrend[finite]) / scale
    return out, scale, scale_source


def detrend_method(
    time: np.ndarray,
    flux: np.ndarray,
    quality: np.ndarray,
    method: MethodSpec,
    *,
    in_transit: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    fit_quality = effective_quality(quality, method, in_transit=in_transit)
    if method.kind == "spline":
        assert method.cfg is not None
        result = flux_space_detrend_result(time, flux, quality=fit_quality, cfg=method.cfg)
        return (
            np.asarray(result.det_flux, dtype=np.float64),
            np.asarray(result.cotrend, dtype=np.float64),
            {
                "status": result.cotrend_status,
                "scale": result.scale,
                "scale_source": result.scale_source,
                "kind": method.kind,
                "bkspace_d": method.cfg.bkspace_d,
                "gap_split_d": method.cfg.gap_split_d,
                "knot_strategy": method.cfg.knot_strategy,
                "fit_mask": method.fit_mask,
            },
        )
    if method.kind == "poly":
        cotrend, status = polynomial_cotrend(
            time,
            flux,
            fit_quality,
            degree=method.degree,
            gap_split_d=method.gap_split_d,
        )
        det, scale, scale_source = relative_subtractive(flux, cotrend, fit_quality)
        return det, cotrend, {
            "status": status,
            "scale": scale,
            "scale_source": scale_source,
            "kind": method.kind,
            "degree": method.degree,
            "gap_split_d": method.gap_split_d,
            "fit_mask": method.fit_mask,
        }
    if method.kind in {"median_filter", "percentile_filter", "savgol"}:
        cotrend, status = filter_cotrend(time, flux, fit_quality, method)
        det, scale, scale_source = relative_subtractive(flux, cotrend, fit_quality)
        return det, cotrend, {
            "status": status,
            "scale": scale,
            "scale_source": scale_source,
            "kind": method.kind,
            "window_d": method.window_d,
            "percentile": method.percentile if method.kind == "percentile_filter" else np.nan,
            "polyorder": method.polyorder if method.kind == "savgol" else np.nan,
            "gap_split_d": method.gap_split_d,
            "fit_mask": method.fit_mask,
        }
    if method.kind == "constant":
        fit = np.isfinite(flux) & (fit_quality == 0)
        med = float(np.nanmedian(flux[fit])) if np.any(fit) else np.nan
        cotrend = np.full_like(flux, med, dtype=np.float64)
        det, scale, scale_source = relative_subtractive(flux, cotrend, fit_quality)
        return det, cotrend, {
            "status": "constant",
            "scale": scale,
            "scale_source": scale_source,
            "kind": method.kind,
            "fit_mask": method.fit_mask,
        }
    raise ValueError(f"unknown method kind: {method.kind}")


def binned_trend_metrics(
    time: np.ndarray,
    values: np.ndarray,
    quality: np.ndarray,
    *,
    bin_d: float,
) -> dict[str, float]:
    good = np.isfinite(time) & np.isfinite(values) & (quality == 0)
    if np.count_nonzero(good) < 50:
        return {"trend_ptp": np.nan, "trend_sigma": np.nan, "n_trend_bins": 0}
    t = time[good]
    y = values[good]
    lo = float(np.nanmin(t))
    hi = float(np.nanmax(t))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return {"trend_ptp": np.nan, "trend_sigma": np.nan, "n_trend_bins": 0}
    edges = np.arange(lo, hi + bin_d, bin_d)
    if edges.size < 3:
        edges = np.linspace(lo, hi, 4)
    medians: list[float] = []
    for a, b in zip(edges[:-1], edges[1:]):
        mask = (t >= a) & (t < b)
        if np.count_nonzero(mask) >= 10:
            medians.append(float(np.nanmedian(y[mask])))
    med = np.asarray(medians, dtype=np.float64)
    med = med[np.isfinite(med)]
    if med.size < 3:
        return {"trend_ptp": np.nan, "trend_sigma": np.nan, "n_trend_bins": int(med.size)}
    centered = med - np.nanmedian(med)
    return {
        "trend_ptp": float(np.nanpercentile(centered, 90) - np.nanpercentile(centered, 10)),
        "trend_sigma": robust_sigma(centered),
        "n_trend_bins": int(med.size),
    }


def normalized_raw_flux(raw: np.ndarray, quality: np.ndarray) -> np.ndarray:
    fit = np.isfinite(raw) & (quality == 0)
    med = float(np.nanmedian(raw[fit])) if np.any(fit) else np.nan
    scale = robust_abs_scale(raw[fit])
    out = np.full_like(raw, np.nan, dtype=np.float64)
    ok = np.isfinite(raw) & np.isfinite(med) & np.isfinite(scale) & (scale > 0)
    out[ok] = 1.0 + (raw[ok] - med) / scale
    return out


def float_attr(attrs: Any, key: str) -> float:
    try:
        return float(attrs.get(key, np.nan))
    except (TypeError, ValueError):
        return float("nan")


def measure_one(
    *,
    group: h5py.Group,
    injection_id: str,
    raw_aperture: str,
    method: MethodSpec,
    trend_bin_d: float,
) -> dict[str, Any] | None:
    original_name = f"RAW_FLUX_{raw_aperture}_original"
    injected_name = f"RAW_FLUX_{raw_aperture}_injected"
    if original_name not in group or injected_name not in group:
        return None
    time = np.asarray(group["time"], dtype=np.float64)
    quality = np.asarray(group["quality"], dtype=np.int32)
    in_transit = np.asarray(group["in_transit"], dtype=bool)
    original_raw = np.asarray(group[original_name], dtype=np.float64)
    injected_raw = np.asarray(group[injected_name], dtype=np.float64)
    det_original, cotrend, meta = detrend_method(
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

    finite = np.isfinite(det_original) & np.isfinite(det_injected)
    good = finite & (quality == 0)
    in_good = good & in_transit
    oot_good = good & ~in_transit
    if np.count_nonzero(in_good) < 1 or np.count_nonzero(oot_good) < 20:
        return None
    delta = det_injected - det_original
    delta_in = float(np.nanmedian(delta[in_good]))
    delta_oot = float(np.nanmedian(delta[oot_good]))
    delta_depth = -(delta_in - delta_oot)
    oot_sigma = robust_sigma(det_original[oot_good])
    per_cad_snr = delta_depth / oot_sigma if np.isfinite(oot_sigma) and oot_sigma > 0 else np.nan
    multi_snr = per_cad_snr * np.sqrt(np.count_nonzero(in_good)) if np.isfinite(per_cad_snr) else np.nan

    truth_sampled = float_attr(group.attrs, "sampled_model_depth")
    truth_model = float_attr(group.attrs, "model_depth")
    truth_depth = truth_sampled if np.isfinite(truth_sampled) and truth_sampled > 0 else truth_model
    retention = delta_depth / truth_depth if np.isfinite(truth_depth) and truth_depth > 0 else np.nan

    raw_norm = normalized_raw_flux(original_raw, quality)
    raw_trend = binned_trend_metrics(time, raw_norm, quality, bin_d=trend_bin_d)
    det_trend = binned_trend_metrics(time, det_original, quality, bin_d=trend_bin_d)
    raw_ptp = raw_trend["trend_ptp"]
    det_ptp = det_trend["trend_ptp"]
    trend_reduction = 1.0 - det_ptp / raw_ptp if np.isfinite(raw_ptp) and raw_ptp > 0 and np.isfinite(det_ptp) else np.nan

    return {
        "injection_id": injection_id,
        "tic": int(group.attrs.get("tic", -1)),
        "tmag": float_attr(group.attrs, "tessmag"),
        "raw_aperture": raw_aperture,
        "method": method.name,
        "method_kind": meta.get("kind", ""),
        "fit_mask": meta.get("fit_mask", method.fit_mask),
        "status": meta.get("status", ""),
        "scale": meta.get("scale", np.nan),
        "scale_source": meta.get("scale_source", ""),
        "bkspace_d": meta.get("bkspace_d", np.nan),
        "window_d": meta.get("window_d", np.nan),
        "polyorder": meta.get("polyorder", np.nan),
        "gap_split_d": meta.get("gap_split_d", np.nan),
        "knot_strategy": meta.get("knot_strategy", ""),
        "truth_period_d": float_attr(group.attrs, "period_d"),
        "truth_duration_min": float_attr(group.attrs, "duration_min"),
        "truth_depth": float_attr(group.attrs, "depth"),
        "truth_model_depth": truth_model,
        "truth_sampled_model_depth": truth_sampled,
        "truth_radius_rearth": float_attr(group.attrs, "radius_rearth"),
        "truth_impact_b": float_attr(group.attrs, "impact_b"),
        "n_good_in_transit": int(np.count_nonzero(in_good)),
        "n_good_oot": int(np.count_nonzero(oot_good)),
        "delta_in_median": delta_in,
        "delta_oot_median": delta_oot,
        "delta_depth_abs": delta_depth,
        "depth_retention_frac": retention,
        "oot_sigma_mad": oot_sigma,
        "delta_signal_snr_mad": per_cad_snr,
        "delta_signal_multi_snr_mad": multi_snr,
        "raw_trend_ptp": raw_ptp,
        "raw_trend_sigma": raw_trend["trend_sigma"],
        "det_trend_ptp": det_ptp,
        "det_trend_sigma": det_trend["trend_sigma"],
        "trend_reduction_frac": trend_reduction,
        "n_trend_bins": det_trend["n_trend_bins"],
        "median_det_original": float(np.nanmedian(det_original[oot_good])),
    }


def summarize(df: pd.DataFrame, *, snr_threshold: float) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (method, raw_aperture), group in df.groupby(["method", "raw_aperture"], dropna=False):
        finite_ret = pd.to_numeric(group["depth_retention_frac"], errors="coerce").replace([np.inf, -np.inf], np.nan)
        finite_snr = pd.to_numeric(group["delta_signal_multi_snr_mad"], errors="coerce").replace([np.inf, -np.inf], np.nan)
        finite_trend = pd.to_numeric(group["det_trend_ptp"], errors="coerce").replace([np.inf, -np.inf], np.nan)
        finite_reduction = pd.to_numeric(group["trend_reduction_frac"], errors="coerce").replace([np.inf, -np.inf], np.nan)
        rows.append(
            {
                "method": method,
                "raw_aperture": raw_aperture,
                "n": int(len(group)),
                "median_depth_retention": float(np.nanmedian(finite_ret)),
                "p16_depth_retention": float(np.nanpercentile(finite_ret.dropna(), 16)) if finite_ret.notna().any() else np.nan,
                "p84_depth_retention": float(np.nanpercentile(finite_ret.dropna(), 84)) if finite_ret.notna().any() else np.nan,
                "median_multi_snr": float(np.nanmedian(finite_snr)),
                "n_snr_ge_threshold": int((finite_snr >= snr_threshold).sum()),
                "frac_snr_ge_threshold": float((finite_snr >= snr_threshold).mean()),
                "median_det_trend_ptp": float(np.nanmedian(finite_trend)),
                "median_trend_reduction_frac": float(np.nanmedian(finite_reduction)),
                "median_tmag": float(np.nanmedian(pd.to_numeric(group["tmag"], errors="coerce"))),
            }
        )
    out = pd.DataFrame(rows)
    return out.sort_values(
        ["n_snr_ge_threshold", "median_trend_reduction_frac", "median_multi_snr"],
        ascending=[False, False, False],
    )


def frame_to_markdown(df: pd.DataFrame, *, floatfmt: str = ".4g") -> str:
    if df.empty:
        return "_No rows._"
    headers = [str(c) for c in df.columns]
    body: list[list[str]] = []
    for _, row in df.iterrows():
        values: list[str] = []
        for value in row:
            if isinstance(value, float):
                values.append(format(value, floatfmt) if np.isfinite(value) else "")
            else:
                values.append(str(value))
        body.append(values)
    widths = [max(len(headers[i]), *(len(r[i]) for r in body)) for i in range(len(headers))]
    lines = [
        "| " + " | ".join(headers[i].ljust(widths[i]) for i in range(len(headers))) + " |",
        "| " + " | ".join("-" * widths[i] for i in range(len(headers))) + " |",
    ]
    lines.extend("| " + " | ".join(row[i].ljust(widths[i]) for i in range(len(headers))) + " |" for row in body)
    return "\n".join(lines)


def write_report(summary_df: pd.DataFrame, out_path: Path, *, snr_threshold: float) -> None:
    lines = [
        "# Pre-Detrend Detrending Method Sweep",
        "",
        f"Empirical SNR threshold: `{snr_threshold:g}`.",
        "",
        "## Best Rows",
        "",
        frame_to_markdown(summary_df.head(20), floatfmt=".4g"),
        "",
        "## Interpretation Aid",
        "",
        "- `n_snr_ge_threshold` is the number of rows whose injected-minus-original truth-window signal reaches the SNR threshold after detrending.",
        "- `median_depth_retention` near `1` means the method preserves the injected BATMAN signal depth.",
        "- `median_det_trend_ptp` is the robust 90-10% range of binned original detrended flux; lower is better.",
        "- `median_trend_reduction_frac` compares that trend range to the raw normalized flux; higher is better.",
        "",
    ]
    out_path.write_text("\n".join(lines), encoding="utf-8")


def plot_pareto(summary_df: pd.DataFrame, out_dir: Path) -> dict[str, str]:
    paths: dict[str, str] = {}
    if summary_df.empty:
        return paths
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    apply_twirl_style("column")
    apertures = list(dict.fromkeys(summary_df["raw_aperture"].astype(str)))
    colors = dict(zip(apertures, get_ordered_palette(len(apertures), "viridis")))
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for ap, sub in summary_df.groupby("raw_aperture", dropna=False):
        ax.scatter(
            sub["median_det_trend_ptp"],
            sub["n_snr_ge_threshold"],
            s=45,
            label=str(ap),
            color=colors[str(ap)],
            edgecolors="black",
            linewidths=0.4,
        )
        for _, row in sub.iterrows():
            if row["n_snr_ge_threshold"] >= summary_df["n_snr_ge_threshold"].max() - 2:
                ax.annotate(
                    str(row["method"]),
                    (row["median_det_trend_ptp"], row["n_snr_ge_threshold"]),
                    xytext=(3, 3),
                    textcoords="offset points",
                    fontsize=6,
                )
    ax.set_xlabel("Residual trend 90-10% range")
    ax.set_ylabel("Rows with empirical signal SNR >= threshold")
    ax.legend(title="Raw aperture", fontsize=7, title_fontsize=8)
    fig.tight_layout()
    png = out_dir / "detrending_method_pareto.png"
    pdf = out_dir / "detrending_method_pareto.pdf"
    fig.savefig(png, dpi=180)
    fig.savefig(pdf)
    plt.close(fig)
    paths["pareto_png"] = str(png)
    paths["pareto_pdf"] = str(pdf)
    return paths


def run_sweep(
    *,
    injection_h5: Path,
    out_dir: Path,
    methods: tuple[str, ...],
    raw_apertures: tuple[str, ...],
    injection_ids: tuple[str, ...] | None,
    limit: int,
    trend_bin_d: float,
    snr_threshold: float,
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    method_specs = tuple(METHOD_SPECS[name] for name in methods)
    rows: list[dict[str, Any]] = []
    with h5py.File(injection_h5, "r") as h5:
        keys = sorted(h5["injections"].keys())
        if injection_ids is not None:
            wanted = set(injection_ids)
            keys = [key for key in keys if key in wanted]
        if limit > 0:
            keys = keys[:limit]
        for idx, key in enumerate(keys, 1):
            group = h5["injections"][key]
            if idx == 1 or idx % 25 == 0:
                print(f"[detrend-sweep] {idx}/{len(keys)} {key}", flush=True)
            for raw_ap in raw_apertures:
                for method in method_specs:
                    measured = measure_one(
                        group=group,
                        injection_id=key,
                        raw_aperture=raw_ap,
                        method=method,
                        trend_bin_d=trend_bin_d,
                    )
                    if measured is not None:
                        rows.append(measured)
    detail = pd.DataFrame(rows)
    detail.to_csv(out_dir / "detrending_method_detail.csv", index=False)
    summary_df = summarize(detail, snr_threshold=snr_threshold) if not detail.empty else pd.DataFrame()
    summary_df.to_csv(out_dir / "detrending_method_summary.csv", index=False)
    write_report(summary_df, out_dir / "summary.md", snr_threshold=snr_threshold)
    plot_paths = plot_pareto(summary_df, out_dir)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "injection_h5": str(injection_h5),
        "out_dir": str(out_dir),
        "n_injection_ids": len(keys),
        "n_detail_rows": int(len(detail)),
        "methods": list(methods),
        "raw_apertures": list(raw_apertures),
        "trend_bin_d": float(trend_bin_d),
        "snr_threshold": float(snr_threshold),
        "best_rows": summary_df.head(10).to_dict(orient="records") if not summary_df.empty else [],
        "plot_paths": plot_paths,
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--injection-h5", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--injection-id-file", type=Path, default=None)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--method", action="append", default=[])
    parser.add_argument("--raw-aperture", action="append", default=[])
    parser.add_argument("--trend-bin-d", type=float, default=0.5)
    parser.add_argument("--snr-threshold", type=float, default=7.0)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    methods = tuple(args.method) if args.method else DEFAULT_METHODS
    unknown = sorted(set(methods) - set(METHOD_SPECS))
    if unknown:
        raise SystemExit(f"unknown methods: {unknown}; valid={sorted(METHOD_SPECS)}")
    raw_apertures = tuple(args.raw_aperture) if args.raw_aperture else RAW_APERTURES
    unknown_ap = sorted(set(raw_apertures) - set(RAW_APERTURES))
    if unknown_ap:
        raise SystemExit(f"unknown raw apertures: {unknown_ap}; valid={RAW_APERTURES}")
    run_sweep(
        injection_h5=args.injection_h5,
        out_dir=args.out_dir,
        methods=methods,
        raw_apertures=raw_apertures,
        injection_ids=load_injection_ids(args.injection_id_file),
        limit=args.limit,
        trend_bin_d=args.trend_bin_d,
        snr_threshold=args.snr_threshold,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
