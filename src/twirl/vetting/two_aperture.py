"""Two-aperture TWIRL vet sheets for S56 human triage."""
from __future__ import annotations

from dataclasses import asdict
import os
from pathlib import Path
import tempfile
from typing import Any, Mapping

import numpy as np

_MPL_CACHE = Path(tempfile.gettempdir()) / "twirl_mplconfig"
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))

from twirl.io.hlsp import BJDREFI, HLSPLightCurve, quality_mask
from twirl.lightcurves.detrend_presets import (
    ADP015Q_COLUMN_TAG,
    TWIRL_FS_V2_ADP015Q_BRANCH,
    compare_column_names,
)
from twirl.plotting.style import apply_twirl_style
from twirl.search.bls import BLSConfig, run_bls_on_lc
from twirl.vetting.adpplus import (
    ADP_VET_APERTURES,
    ADPPlusBranch,
    binned_trend_ptp,
    branch_by_name,
    branched_light_curve,
)


_ADP015_COLUMNS = compare_column_names(ADP015Q_COLUMN_TAG)
DEFAULT_TWO_APERTURE_BRANCH = TWIRL_FS_V2_ADP015Q_BRANCH
DEFAULT_TWO_APERTURE_APERTURES: tuple[str, str] = (
    _ADP015_COLUMNS["small"],
    _ADP015_COLUMNS["primary"],
)
TWO_APERTURE_VET_SHEET_VERSION = "S56-ADP-HV2"

_FOLD_RAW_MARKER_SIZE = 5.2
_FOLD_RAW_ALPHA = 0.34
_FOLD_BIN_MARKER_SIZE = 4.6
_FOLD_BIN_LINE_WIDTH = 0.78


def _aperture_prefix(aperture: str) -> str:
    return aperture.replace("DET_FLUX_", "").lower()


def _phase_min(time_btjd: np.ndarray, *, period_d: float, t0_bjd: float) -> np.ndarray:
    t0_d = float(t0_bjd) - float(BJDREFI)
    phase_d = ((time_btjd - t0_d + 0.5 * period_d) % period_d) - 0.5 * period_d
    return phase_d * 1440.0


def _phase_cycle(time_btjd: np.ndarray, *, period_d: float, t0_bjd: float) -> np.ndarray:
    t0_d = float(t0_bjd) - float(BJDREFI)
    phase_d = ((time_btjd - t0_d + 0.5 * period_d) % period_d) - 0.5 * period_d
    return phase_d / float(period_d)


def _bin_xy(x: np.ndarray, y: np.ndarray, *, n_bins: int = 80) -> tuple[np.ndarray, np.ndarray]:
    bx, by, _, _ = _bin_stats(x, y, n_bins=n_bins)
    return bx, by


def _bin_stats(
    x: np.ndarray,
    y: np.ndarray,
    *,
    n_bins: int = 80,
    x_min: float | None = None,
    x_max: float | None = None,
    min_count: int = 2,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Median-bin a folded light curve and return point errors."""

    finite = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(finite) < 3:
        empty = np.array([], dtype=float)
        return empty, empty, empty, empty.astype(int)
    xg = x[finite]
    yg = y[finite]
    lo = float(np.nanmin(xg)) if x_min is None else float(x_min)
    hi = float(np.nanmax(xg)) if x_max is None else float(x_max)
    edges = np.linspace(lo, hi, n_bins + 1)
    if not np.all(np.isfinite(edges)) or np.unique(edges).size < 3:
        empty = np.array([], dtype=float)
        return empty, empty, empty, empty.astype(int)
    idx = np.digitize(xg, edges) - 1
    in_range = (idx >= 0) & (idx < n_bins)
    idx = idx[in_range]
    yg = yg[in_range]
    centers = 0.5 * (edges[:-1] + edges[1:])
    if idx.size == 0:
        empty = np.array([], dtype=float)
        return empty, empty, empty, empty.astype(int)
    order = np.argsort(idx, kind="stable")
    idx_s = idx[order]
    y_s = yg[order]
    unique_idx, starts, counts = np.unique(idx_s, return_index=True, return_counts=True)
    keep = counts >= int(min_count)
    if not np.any(keep):
        empty = np.array([], dtype=float)
        return empty, empty, empty, empty.astype(int)
    out_idx = unique_idx[keep]
    out_counts = counts[keep].astype(int)
    meds = np.full(out_idx.size, np.nan, dtype=float)
    errs = np.full(out_idx.size, np.nan, dtype=float)
    for j, (start, count) in enumerate(zip(starts[keep], out_counts)):
        vals = y_s[start : start + count]
        med = float(np.nanmedian(vals))
        meds[j] = med
        mad = 1.4826 * float(np.nanmedian(np.abs(vals - med)))
        errs[j] = mad / np.sqrt(count) if np.isfinite(mad) and count > 0 else np.nan
    ok = np.isfinite(meds)
    return centers[out_idx][ok], meds[ok], errs[ok], out_counts[ok]


def _fold_window_min(duration_min: float | None) -> float:
    if duration_min is None or not np.isfinite(duration_min) or duration_min <= 0:
        return 45.0
    return float(np.clip(max(18.0, 4.0 * float(duration_min)), 18.0, 60.0))


def _fold_bin_count(window_min: float, duration_min: float | None) -> int:
    if duration_min is None or not np.isfinite(duration_min) or duration_min <= 0:
        bin_width = 2.0
    else:
        bin_width = max(1.0, min(3.0, float(duration_min) / 3.0))
    return int(np.clip(np.ceil(2.0 * float(window_min) / bin_width), 18, 72))


def _fold_window_phase(period_d: float, duration_min: float | None) -> tuple[float, float]:
    window_min = _fold_window_min(duration_min)
    if not np.isfinite(period_d) or period_d <= 0:
        return 0.03, window_min
    window = window_min / (float(period_d) * 1440.0)
    return float(np.clip(window, 0.0015, 0.25)), window_min


def _event_numbers(time_btjd: np.ndarray, *, period_d: float, t0_bjd: float) -> np.ndarray:
    t0_d = float(t0_bjd) - float(BJDREFI)
    return np.rint((np.asarray(time_btjd, dtype=float) - t0_d) / float(period_d)).astype(int)


def _event_times_btjd(time_btjd: np.ndarray, *, period_d: float, t0_bjd: float) -> np.ndarray:
    finite = np.asarray(time_btjd, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0 or not np.isfinite(period_d) or period_d <= 0 or not np.isfinite(t0_bjd):
        return np.array([], dtype=float)
    t0_d = float(t0_bjd) - float(BJDREFI)
    n_lo = int(np.ceil((float(np.nanmin(finite)) - t0_d) / float(period_d)))
    n_hi = int(np.floor((float(np.nanmax(finite)) - t0_d) / float(period_d)))
    if n_hi < n_lo:
        return np.array([], dtype=float)
    return t0_d + np.arange(n_lo, n_hi + 1, dtype=float) * float(period_d)


def _mark_transit_windows(
    ax: Any,
    time_btjd: np.ndarray,
    *,
    period_d: float,
    t0_bjd: float,
    duration_min: float | None,
) -> None:
    events = _event_times_btjd(time_btjd, period_d=period_d, t0_bjd=t0_bjd)
    if events.size == 0:
        return
    event_no = _event_numbers(events, period_d=period_d, t0_bjd=t0_bjd)
    ymin, ymax = ax.get_ylim()
    yr = ymax - ymin
    y0 = ymin + 0.005 * yr
    y1 = ymin + 0.160 * yr
    even = (event_no % 2) == 0
    odd = ~even
    if np.any(even):
        ax.vlines(events[even], y0, y1, color="#1f4e79", lw=2.0, alpha=0.92)
    if np.any(odd):
        ax.vlines(events[odd], y0, y1, color="#b83227", lw=2.0, alpha=0.92)


def _norm_flux(lc: HLSPLightCurve, aperture: str) -> np.ndarray:
    flux = np.asarray(lc.flux[aperture], dtype=float)
    keep = quality_mask(lc, aperture)
    med = float(np.nanmedian(flux[keep])) if np.any(keep) else float(np.nanmedian(flux[np.isfinite(flux)]))
    if not np.isfinite(med) or med == 0:
        med = 1.0
    return flux / med


def _run_branch_bls(
    lc: HLSPLightCurve,
    *,
    aperture: str,
    branch: ADPPlusBranch,
    cfg: BLSConfig,
    period_d: float | None = None,
    t0_bjd: float | None = None,
    duration_min: float | None = None,
    return_periodogram: bool = True,
) -> tuple[Any, dict[str, np.ndarray] | None, HLSPLightCurve, dict[str, Any]]:
    branch_lc, meta = branched_light_curve(
        lc,
        aperture,
        branch,
        period_d=period_d,
        t0_bjd=t0_bjd,
        duration_min=duration_min,
    )
    result = run_bls_on_lc(branch_lc, cfg, aperture=aperture, return_periodogram=return_periodogram)
    if return_periodogram:
        res, spectrum = result
        return res, spectrum, branch_lc, dict(meta)
    return result, None, branch_lc, dict(meta)


def _best_peak_dict(result: Any) -> dict[str, float]:
    if not getattr(result, "peaks", None):
        return {
            "period_d": float("nan"),
            "t0_bjd": float("nan"),
            "duration_min": float("nan"),
            "depth": float("nan"),
            "depth_snr": float("nan"),
            "sde": float("nan"),
        }
    return asdict(result.peaks[0])


def _finite_peak_dict(peak: dict[str, Any] | None) -> dict[str, float] | None:
    if peak is None:
        return None
    try:
        period_d = float(peak.get("period_d", np.nan))
        t0_bjd = float(peak.get("t0_bjd", np.nan))
        duration_min = float(peak.get("duration_min", np.nan))
    except (TypeError, ValueError):
        return None
    if not (np.isfinite(period_d) and np.isfinite(t0_bjd) and np.isfinite(duration_min)):
        return None
    if period_d <= 0 or duration_min <= 0:
        return None
    return {
        "period_d": period_d,
        "t0_bjd": t0_bjd,
        "duration_min": duration_min,
        "depth": float(peak.get("depth", np.nan)),
        "depth_snr": float(peak.get("depth_snr", np.nan)),
        "sde": float(peak.get("sde", np.nan)),
    }


def _row_ephemeris_peak(row_metadata: dict[str, Any] | None) -> dict[str, float] | None:
    if not row_metadata:
        return None
    peak = _finite_peak_dict(
        {
            "period_d": row_metadata.get("period_d"),
            "t0_bjd": row_metadata.get("t0_bjd"),
            "duration_min": row_metadata.get("duration_min"),
            "depth": row_metadata.get("depth", np.nan),
            "depth_snr": row_metadata.get("depth_snr", np.nan),
            "sde": row_metadata.get("sde_max", row_metadata.get("sde", np.nan)),
        }
    )
    return peak


def _depth_at_ephemeris(
    lc: HLSPLightCurve,
    aperture: str,
    *,
    period_d: float,
    t0_bjd: float,
    duration_min: float,
) -> tuple[float, float, int]:
    y = _norm_flux(lc, aperture)
    keep = quality_mask(lc, aperture)
    phase = _phase_min(lc.time, period_d=period_d, t0_bjd=t0_bjd)
    half = 0.5 * float(duration_min)
    in_tr = keep & np.isfinite(y) & (np.abs(phase) <= half)
    oot = keep & np.isfinite(y) & (np.abs(phase) > max(2.0 * half, 20.0)) & (np.abs(phase) < 90.0)
    if np.count_nonzero(in_tr) < 1 or np.count_nonzero(oot) < 10:
        return float("nan"), float("nan"), int(np.count_nonzero(in_tr))
    oot_med = float(np.nanmedian(y[oot]))
    in_med = float(np.nanmedian(y[in_tr]))
    sigma = 1.4826 * float(np.nanmedian(np.abs(y[oot] - oot_med)))
    depth = oot_med - in_med
    snr = depth / sigma * np.sqrt(np.count_nonzero(in_tr)) if sigma > 0 else float("nan")
    return float(depth), float(snr), int(np.count_nonzero(in_tr))


def _even_odd_depth_metrics(
    lc: HLSPLightCurve,
    aperture: str,
    *,
    period_d: float,
    t0_bjd: float,
    duration_min: float,
) -> dict[str, float | int]:
    y = _norm_flux(lc, aperture)
    keep = quality_mask(lc, aperture)
    phase = _phase_min(lc.time, period_d=period_d, t0_bjd=t0_bjd)
    half = 0.5 * float(duration_min)
    in_tr = keep & np.isfinite(y) & (np.abs(phase) <= half)
    oot = keep & np.isfinite(y) & (np.abs(phase) > max(2.0 * half, 20.0)) & (np.abs(phase) < 90.0)
    event_no = _event_numbers(lc.time, period_d=period_d, t0_bjd=t0_bjd)
    if np.count_nonzero(oot) < 10:
        return {
            "even_depth": np.nan,
            "odd_depth": np.nan,
            "even_n_in": int(np.count_nonzero(in_tr & ((event_no % 2) == 0))),
            "odd_n_in": int(np.count_nonzero(in_tr & ((event_no % 2) != 0))),
            "even_odd_depth_delta": np.nan,
            "even_odd_sigma_delta": np.nan,
        }
    oot_med = float(np.nanmedian(y[oot]))
    sigma = 1.4826 * float(np.nanmedian(np.abs(y[oot] - oot_med)))
    values: dict[str, float | int] = {}
    errors: dict[str, float] = {}
    for name, parity in (("even", 0), ("odd", 1)):
        sel = in_tr & ((event_no % 2) == parity)
        n_in = int(np.count_nonzero(sel))
        values[f"{name}_n_in"] = n_in
        if n_in < 1:
            values[f"{name}_depth"] = np.nan
            errors[name] = np.nan
            continue
        depth = oot_med - float(np.nanmedian(y[sel]))
        values[f"{name}_depth"] = float(depth)
        errors[name] = sigma / np.sqrt(n_in) if sigma > 0 else np.nan
    even_depth = float(values.get("even_depth", np.nan))
    odd_depth = float(values.get("odd_depth", np.nan))
    delta = odd_depth - even_depth if np.isfinite(even_depth) and np.isfinite(odd_depth) else np.nan
    denom = np.hypot(errors.get("even", np.nan), errors.get("odd", np.nan))
    values["even_odd_depth_delta"] = float(delta) if np.isfinite(delta) else np.nan
    values["even_odd_sigma_delta"] = float(abs(delta) / denom) if np.isfinite(delta) and denom > 0 else np.nan
    return values


def _format_period(period_d: float) -> str:
    return f"{period_d:.6g} d" if np.isfinite(period_d) else "nan"


def _bold_period(period_d: float) -> str:
    if not np.isfinite(period_d):
        return "P=nan"
    return rf"$P=\mathbf{{{float(period_d):.6g}}}\,\mathrm{{d}}$"


def _duration_hours_label(duration_min: float | None) -> str:
    if duration_min is None or not np.isfinite(duration_min):
        return r"$t_\mathrm{dur}=\mathrm{nan}$"
    return rf"$t_\mathrm{{dur}}={float(duration_min) / 60.0:.3g}\,\mathrm{{hr}}$"


def _duration_days_label(duration_min: float | None) -> str:
    if duration_min is None or not np.isfinite(duration_min):
        return r"$t_\mathrm{dur}=\mathrm{nan}$"
    return rf"$t_\mathrm{{dur}}={float(duration_min) / 1440.0:.3g}\,\mathrm{{d}}$"


def _duration_days_text(duration_min: float | None) -> str:
    if duration_min is None or not np.isfinite(duration_min):
        return "nan d"
    return f"{float(duration_min) / 1440.0:.6g} d"


def _set_fold_ylim(
    ax: Any,
    y: np.ndarray,
    sel: np.ndarray,
    *,
    depth: float | None = None,
    shared_ylim: tuple[float, float] | None = None,
) -> None:
    if shared_ylim is not None:
        ax.set_ylim(*shared_ylim)
        return
    if np.count_nonzero(sel) < 3:
        return
    vals = np.asarray(y[sel], dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size < 3:
        return
    med = float(np.nanmedian(vals))
    mad = 1.4826 * float(np.nanmedian(np.abs(vals - med)))
    lo, hi = np.nanpercentile(vals, [2.0, 98.0])
    if np.isfinite(mad) and mad > 0:
        lo = max(float(lo), med - 6.0 * mad)
        hi = min(float(hi), med + 6.0 * mad)
    if depth is not None and np.isfinite(depth) and depth > 0:
        lo = min(float(lo), 1.0 - 1.25 * float(depth))
        hi = max(float(hi), 1.0 + 0.25 * float(depth))
    span = max(float(hi - lo), 0.02)
    ax.set_ylim(float(lo) - 0.12 * span, float(hi) + 0.12 * span)


def _shared_light_curve_ylim(
    lcs: Mapping[str, HLSPLightCurve],
    peaks: Mapping[str, Mapping[str, Any]],
) -> tuple[float, float]:
    """Return one robust normalized-flux scale for every LC panel on a sheet."""

    lows: list[float] = []
    highs: list[float] = []
    for aperture, blc in lcs.items():
        if aperture not in blc.flux:
            continue
        y = _norm_flux(blc, aperture)
        keep = quality_mask(blc, aperture) & np.isfinite(y)
        if np.count_nonzero(keep) < 3:
            continue
        lo, hi = np.nanpercentile(y[keep], [1.0, 99.5])
        peak = peaks.get(aperture, {})
        period_d = float(peak.get("period_d", np.nan))
        t0_bjd = float(peak.get("t0_bjd", np.nan))
        duration_min = float(peak.get("duration_min", np.nan))
        if (
            np.isfinite(period_d)
            and period_d > 0
            and np.isfinite(t0_bjd)
            and np.isfinite(duration_min)
            and duration_min > 0
        ):
            phase = _phase_cycle(blc.time, period_d=period_d, t0_bjd=t0_bjd)
            half_phase = 0.5 * duration_min / (period_d * 1440.0)
            primary = keep & (np.abs(phase) <= half_phase)
            secondary = keep & (np.abs(np.abs(phase) - 0.5) <= half_phase)
            event = primary | secondary
            n_event = int(np.count_nonzero(event))
            if n_event:
                event_qlo = 0.0 if n_event < 4 else 0.5
                event_lo, event_hi = np.nanpercentile(y[event], [event_qlo, 99.5])
                lo = min(float(lo), float(event_lo))
                hi = max(float(hi), float(event_hi))
            depth = float(peak.get("depth", np.nan))
            if np.isfinite(depth) and depth > 0:
                lo = min(float(lo), 1.0 - 1.30 * depth)
                hi = max(float(hi), 1.0 + 0.30 * depth)
        lows.append(float(lo))
        highs.append(float(hi))

    if not lows or not highs:
        return 0.95, 1.05
    lo = float(np.nanmin(lows))
    hi = float(np.nanmax(highs))
    span = max(hi - lo, 0.05)
    pad = 0.05 * span
    return lo - pad, hi + pad


def _annotate_fold_period(
    ax: Any,
    *,
    period_d: float,
    duration_min: float | None,
    sde: float | None = None,
) -> None:
    lines = [f"P={_format_period(period_d)}"]
    if duration_min is not None and np.isfinite(duration_min):
        lines.append(f"dur={_duration_days_text(duration_min)}")
    if sde is not None and np.isfinite(sde):
        lines.append(f"SDE={sde:.2g}")
    ax.text(
        0.02,
        0.96,
        "\n".join(lines),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        bbox={"facecolor": "white", "edgecolor": "0.85", "alpha": 0.86, "pad": 1.8},
    )


def _plot_full_lc(
    ax: Any,
    blc: HLSPLightCurve,
    aperture: str,
    *,
    branch_name: str,
    peak: dict[str, float],
    shared_ylim: tuple[float, float] | None = None,
) -> None:
    y = _norm_flux(blc, aperture)
    keep = quality_mask(blc, aperture)
    ax.scatter(blc.time[keep], y[keep], s=3.0, c="0.12", alpha=0.58, linewidths=0, rasterized=True)
    if keep.any() and shared_ylim is None:
        finite_keep = keep & np.isfinite(y)
        lo, hi = np.nanpercentile(y[finite_keep], [1.0, 99.5])
        period_d = float(peak.get("period_d", np.nan))
        t0_bjd = float(peak.get("t0_bjd", np.nan))
        duration_min = float(peak.get("duration_min", np.nan))
        if np.isfinite(period_d) and np.isfinite(t0_bjd) and np.isfinite(duration_min) and duration_min > 0:
            phase_min = _phase_min(blc.time, period_d=period_d, t0_bjd=t0_bjd)
            in_tr = finite_keep & (np.abs(phase_min) <= 0.5 * duration_min)
            n_in = int(np.count_nonzero(in_tr))
            if n_in > 0:
                qlo = 0.0 if n_in < 4 else 0.5
                in_lo, in_hi = np.nanpercentile(y[in_tr], [qlo, 99.5])
                lo = min(float(lo), float(in_lo))
                hi = max(float(hi), float(in_hi))
            depth = float(peak.get("depth", np.nan))
            if np.isfinite(depth) and depth > 0:
                lo = min(float(lo), 1.0 - 1.30 * depth)
                hi = max(float(hi), 1.0 + 0.30 * depth)
        pad = 0.05 * max(float(hi - lo), 0.05)
        ax.set_ylim(float(lo) - pad, float(hi) + pad)
    if keep.any():
        ax.set_xlim(float(np.nanmin(blc.time[keep])), float(np.nanmax(blc.time[keep])))
    if shared_ylim is not None:
        ax.set_ylim(*shared_ylim)
    if np.isfinite(peak.get("period_d", np.nan)) and np.isfinite(peak.get("t0_bjd", np.nan)):
        _mark_transit_windows(
            ax,
            blc.time[keep],
            period_d=float(peak["period_d"]),
            t0_bjd=float(peak["t0_bjd"]),
            duration_min=peak.get("duration_min", np.nan),
    )
    ax.set_title(f"{aperture} full LC; {_bold_period(float(peak.get('period_d', np.nan)))}", loc="left", fontsize=10.5)
    ax.set_ylabel("relative flux")


def _plot_full_phase_fold(
    ax: Any,
    blc: HLSPLightCurve,
    aperture: str,
    *,
    period_d: float,
    t0_bjd: float,
    duration_min: float,
    depth: float | None = None,
    title: str,
    shared_ylim: tuple[float, float] | None = None,
) -> None:
    y = _norm_flux(blc, aperture)
    keep = quality_mask(blc, aperture)
    phase = _phase_cycle(blc.time, period_d=period_d, t0_bjd=t0_bjd)
    sel = keep & np.isfinite(y) & np.isfinite(phase)
    ax.scatter(phase[sel], y[sel], s=2.6, c="0.42", alpha=0.28, linewidths=0, rasterized=True)
    zoom_window, zoom_window_min = _fold_window_phase(period_d, duration_min)
    zoom_bins = _fold_bin_count(zoom_window_min, duration_min)
    zoom_bin_width = 2.0 * float(zoom_window) / max(int(zoom_bins), 1)
    n_full_bins = int(np.ceil(1.0 / zoom_bin_width)) if zoom_bin_width > 0 else 120
    n_full_bins = max(24, n_full_bins)
    marker_size = float(np.clip(115.0 / np.sqrt(n_full_bins), 1.4, 3.0))
    bx, by, be, _ = _bin_stats(phase[sel], y[sel], n_bins=n_full_bins, x_min=-0.5, x_max=0.5, min_count=1)
    ax.scatter(
        bx,
        by,
        s=marker_size**2,
        c="#1f4e79",
        alpha=0.82,
        linewidths=0,
        rasterized=True,
    )
    half = 0.5 * float(duration_min) / (float(period_d) * 1440.0)
    if np.isfinite(half) and half > 0:
        ax.axvspan(-half, half, color="#d24b32", alpha=0.14, lw=0)
    ax.axvline(0.0, color="#9b2f25", lw=0.9, alpha=0.78)
    ax.set_xlim(-0.5, 0.5)
    _set_fold_ylim(ax, y, sel, depth=depth, shared_ylim=shared_ylim)
    ylo, yhi = ax.get_ylim()
    ypad = 0.04 * max(float(yhi - ylo), 1.0e-6)
    box_y0 = float(ylo) + ypad
    box_y1 = float(yhi) - ypad
    ax.plot(
        [-zoom_window, zoom_window, zoom_window, -zoom_window, -zoom_window],
        [box_y0, box_y0, box_y1, box_y1, box_y0],
        color="0.08",
        lw=0.85,
        alpha=0.92,
        solid_capstyle="butt",
        clip_on=True,
    )
    ax.set_title(title, loc="left", fontsize=10.2)
    ax.set_xlabel("orbital phase")
    ax.set_ylabel("relative flux")


def _plot_harmonic_fold(
    ax: Any,
    blc: HLSPLightCurve,
    aperture: str,
    *,
    period_d: float,
    t0_bjd: float,
    duration_min: float,
    depth: float | None,
    title: str,
    show_ylabel: bool,
    shared_ylim: tuple[float, float] | None = None,
) -> None:
    """Plot one raw-plus-binned full-cycle harmonic view."""

    y = _norm_flux(blc, aperture)
    keep = quality_mask(blc, aperture)
    phase = _phase_cycle(blc.time, period_d=period_d, t0_bjd=t0_bjd)
    phase = np.where(phase < -0.25, phase + 1.0, phase)
    sel = keep & np.isfinite(y) & np.isfinite(phase)
    ax.scatter(phase[sel], y[sel], s=2.2, c="0.45", alpha=0.22, linewidths=0, rasterized=True)
    zoom_window, zoom_window_min = _fold_window_phase(period_d, duration_min)
    zoom_bins = _fold_bin_count(zoom_window_min, duration_min)
    zoom_bin_width = 2.0 * float(zoom_window) / max(int(zoom_bins), 1)
    n_bins = int(np.ceil(1.0 / zoom_bin_width)) if zoom_bin_width > 0 else 100
    n_bins = int(np.clip(n_bins, 36, 180))
    bx, by, _, _ = _bin_stats(
        phase[sel],
        y[sel],
        n_bins=n_bins,
        x_min=-0.25,
        x_max=0.75,
        min_count=1,
    )
    ax.scatter(bx, by, s=7.0, c="#1f4e79", alpha=0.88, linewidths=0, rasterized=True)
    half = 0.5 * float(duration_min) / (float(period_d) * 1440.0)
    if np.isfinite(half) and half > 0:
        ax.axvspan(-half, half, color="#d24b32", alpha=0.13, lw=0)
    ax.axvline(0.0, color="#9b2f25", lw=0.95, alpha=0.84)
    ax.axvline(0.5, color="#1f4e79", lw=0.90, ls="--", alpha=0.82)
    ax.set_xlim(-0.25, 0.75)
    ax.set_xticks([-0.25, 0.0, 0.25, 0.5, 0.75])
    _set_fold_ylim(ax, y, sel, depth=depth, shared_ylim=shared_ylim)
    ax.set_title(title, loc="left", fontsize=8.4)
    ax.set_xlabel("phase", fontsize=8.0)
    ax.set_ylabel("relative flux" if show_ylabel else "", fontsize=8.0)
    ax.tick_params(labelsize=7.2)


def _harmonic_factor_label(factor: float) -> str:
    labels = {0.25: "P/4", 0.5: "P/2", 1.0: "P", 2.0: "2P", 4.0: "4P"}
    for value, label in labels.items():
        if np.isclose(float(factor), value):
            return label
    return f"{float(factor):g}P"


def _plot_periodogram(
    ax: Any,
    spec: dict[str, np.ndarray] | None,
    *,
    bls_peak: dict[str, float],
    review_peak: dict[str, float] | None = None,
) -> None:
    if spec is None:
        ax.text(0.5, 0.5, "No periodogram", ha="center", va="center")
        ax.set_axis_off()
        return
    ax.semilogx(spec["period"], spec["sde"], color="0.25", lw=0.55, rasterized=True)
    review_period = float((review_peak or bls_peak).get("period_d", np.nan))
    bls_period = float(bls_peak.get("period_d", np.nan))
    if np.isfinite(review_period):
        ax.axvline(
            review_period,
            color="#b83227",
            lw=1.15,
            ls="-",
            label=f"review P={review_period:.4g} d",
        )
    if np.isfinite(bls_period) and (
        not np.isfinite(review_period) or not np.isclose(bls_period, review_period, rtol=1.0e-4, atol=0.0)
    ):
        ax.axvline(
            bls_period,
            color="#1f4e79",
            lw=1.0,
            ls="--",
            label=f"aperture BLS max={bls_period:.4g} d",
        )
    ax.set_xlim(float(np.nanmin(spec["period"])), float(np.nanmax(spec["period"])))
    ax.set_ylabel("SDE")
    ax.set_xlabel("Period (d)")
    if ax.get_legend_handles_labels()[0]:
        ax.legend(loc="upper left", fontsize=6.6, framealpha=0.88, handlelength=1.8)
    ax.tick_params(labelsize=8.5)


def _plot_folded_bins(
    ax: Any,
    blc: HLSPLightCurve,
    aperture: str,
    *,
    period_d: float,
    t0_bjd: float,
    duration_min: float,
    depth: float | None = None,
    sde: float | None = None,
    title: str,
    show_duration_days: bool = False,
    phase_shift: float = 0.0,
    shared_ylim: tuple[float, float] | None = None,
) -> None:
    y = _norm_flux(blc, aperture)
    keep = quality_mask(blc, aperture)
    shifted_t0 = float(t0_bjd) + float(phase_shift) * float(period_d)
    phase = _phase_cycle(blc.time, period_d=period_d, t0_bjd=shifted_t0)
    window, window_min = _fold_window_phase(period_d, duration_min)
    sel = keep & np.isfinite(y) & np.isfinite(phase) & (np.abs(phase) <= window)
    ax.scatter(
        phase[sel],
        y[sel],
        s=_FOLD_RAW_MARKER_SIZE,
        c="0.55",
        alpha=_FOLD_RAW_ALPHA,
        linewidths=0,
        rasterized=True,
    )
    n_bins = _fold_bin_count(window_min, duration_min)
    bx, by, be, bn = _bin_stats(phase[sel], y[sel], n_bins=n_bins, x_min=-window, x_max=window, min_count=1)
    good_err = np.isfinite(be) & (be > 0)
    yerr = np.where(good_err, be, np.nan)
    ax.errorbar(
        bx,
        by,
        yerr=yerr,
        fmt="o",
        ms=_FOLD_BIN_MARKER_SIZE,
        lw=_FOLD_BIN_LINE_WIDTH,
        color="#1f4e79",
        ecolor="#1f4e79",
        alpha=0.92,
        capsize=0,
        linestyle="none",
    )
    half = 0.5 * float(duration_min) / (float(period_d) * 1440.0)
    ax.axvspan(-half, half, color="#d24b32", alpha=0.13, lw=0)
    ax.axvline(0.0, color="#9b2f25", lw=0.9, alpha=0.80)
    ax.set_xlim(-window, window)
    _set_fold_ylim(ax, y, sel, depth=depth, shared_ylim=shared_ylim)
    if phase_shift:
        title_bits = [title, r"$\phi=0.5$"]
    else:
        title_bits = [title]
    if show_duration_days:
        title_bits.append(_duration_days_label(duration_min))
    ax.set_title("; ".join(title_bits), loc="left", fontsize=10.2)
    ax.set_xlabel("orbital phase")
    ax.set_ylabel("relative flux")


def _plot_even_odd_bins(
    ax: Any,
    blc: HLSPLightCurve,
    aperture: str,
    *,
    period_d: float,
    t0_bjd: float,
    duration_min: float,
    title: str,
    metrics: dict[str, Any],
    prefix: str,
    metric_scope: str = "own",
    shared_ylim: tuple[float, float] | None = None,
) -> None:
    y = _norm_flux(blc, aperture)
    keep = quality_mask(blc, aperture)
    phase = _phase_cycle(blc.time, period_d=period_d, t0_bjd=t0_bjd)
    event_no = _event_numbers(blc.time, period_d=period_d, t0_bjd=t0_bjd)
    window, window_min = _fold_window_phase(period_d, duration_min)
    sel = keep & np.isfinite(y) & np.isfinite(phase) & (np.abs(phase) <= window)
    ax.scatter(
        phase[sel],
        y[sel],
        s=_FOLD_RAW_MARKER_SIZE,
        c="0.55",
        alpha=_FOLD_RAW_ALPHA,
        linewidths=0,
        rasterized=True,
    )
    colors = {"even": "#1f4e79", "odd": "#c05a27"}
    markers = {"even": "o", "odd": "s"}
    n_bins = _fold_bin_count(window_min, duration_min)
    for name, parity in (("even", 0), ("odd", 1)):
        part = sel & ((event_no % 2) == parity)
        bx, by, be, _ = _bin_stats(phase[part], y[part], n_bins=n_bins, x_min=-window, x_max=window, min_count=1)
        ax.errorbar(
            bx,
            by,
            yerr=np.where(np.isfinite(be) & (be > 0), be, np.nan),
            fmt=markers[name],
            ms=_FOLD_BIN_MARKER_SIZE,
            lw=_FOLD_BIN_LINE_WIDTH,
            color=colors[name],
            ecolor=colors[name],
            alpha=0.92,
            capsize=0,
            linestyle="none",
        )
    half = 0.5 * float(duration_min) / (float(period_d) * 1440.0)
    ax.axvspan(-half, half, color="#d24b32", alpha=0.11, lw=0)
    ax.axvline(0.0, color="#9b2f25", lw=0.9, alpha=0.75)
    ax.set_xlim(-window, window)
    _set_fold_ylim(ax, y, sel, shared_ylim=shared_ylim)
    ax.set_title(title, loc="left", fontsize=10.2)
    sigma_delta = metrics.get(f"{prefix}_{metric_scope}_even_odd_sigma_delta", np.nan)
    if np.isfinite(sigma_delta):
        ax.text(
            0.98,
            0.04,
            f"|odd-even|={float(sigma_delta):.2g} sigma",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8.5,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.80, "pad": 1.5},
        )
    ax.set_xlabel("orbital phase")
    ax.set_ylabel("relative flux")


def build_two_aperture_metrics(
    lc: HLSPLightCurve,
    *,
    branch: ADPPlusBranch,
    cfg: BLSConfig,
    apertures: tuple[str, ...] = DEFAULT_TWO_APERTURE_APERTURES,
    anchor_aperture: str | None = None,
    anchor_peak_override: dict[str, Any] | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Run BLS for both apertures and return metrics plus plotting payload."""

    results: dict[str, Any] = {}
    spectra: dict[str, dict[str, np.ndarray]] = {}
    lcs: dict[str, HLSPLightCurve] = {}
    branch_meta: dict[str, Any] = {}
    for ap in apertures:
        if ap not in lc.flux:
            continue
        if branch.is_identity:
            res, spec, branch_lc, meta = _run_branch_bls(lc, aperture=ap, branch=branch, cfg=cfg)
        else:
            current = branch_by_name("current_adp")
            current_res, _, _, _ = _run_branch_bls(lc, aperture=ap, branch=current, cfg=cfg, return_periodogram=False)
            peak = _best_peak_dict(current_res)
            res, spec, branch_lc, meta = _run_branch_bls(
                lc,
                aperture=ap,
                branch=branch,
                cfg=cfg,
                period_d=peak["period_d"],
                t0_bjd=peak["t0_bjd"],
                duration_min=peak["duration_min"],
            )
        results[ap] = res
        if spec is not None:
            spectra[ap] = spec
        lcs[ap] = branch_lc
        branch_meta[ap] = meta

    if anchor_aperture and anchor_aperture in results:
        anchor_ap = anchor_aperture
    else:
        anchor_ap = apertures[0] if apertures and apertures[0] in results else next(iter(results), "")
    override_peak = _finite_peak_dict(anchor_peak_override)
    anchor_peak = override_peak or (_best_peak_dict(results[anchor_ap]) if anchor_ap else {})
    metrics: dict[str, Any] = {
        "tic": lc.tic,
        "sector": lc.sector,
        "cam": lc.cam,
        "ccd": lc.ccd,
        "tmag": lc.tmag,
        "vet_branch": branch.name,
        "anchor_aperture": anchor_ap,
        "anchor_period_d": anchor_peak.get("period_d", np.nan),
        "anchor_t0_bjd": anchor_peak.get("t0_bjd", np.nan),
        "anchor_duration_min": anchor_peak.get("duration_min", np.nan),
        "anchor_sde": anchor_peak.get("sde", np.nan),
        "anchor_source": "row_metadata" if override_peak else "bls_anchor_aperture",
    }
    for ap, res in results.items():
        peak = _best_peak_dict(res)
        prefix = _aperture_prefix(ap)
        metrics[f"{prefix}_status"] = res.status
        for key, value in peak.items():
            metrics[f"{prefix}_{key}"] = value
        if np.isfinite(peak.get("period_d", np.nan)) and np.isfinite(peak.get("t0_bjd", np.nan)) and np.isfinite(peak.get("duration_min", np.nan)):
            for key, value in _even_odd_depth_metrics(
                lcs[ap],
                ap,
                period_d=float(peak["period_d"]),
                t0_bjd=float(peak["t0_bjd"]),
                duration_min=float(peak["duration_min"]),
            ).items():
                metrics[f"{prefix}_own_{key}"] = value
        if np.isfinite(metrics["anchor_period_d"]):
            depth, snr, n_in = _depth_at_ephemeris(
                lcs[ap],
                ap,
                period_d=float(metrics["anchor_period_d"]),
                t0_bjd=float(metrics["anchor_t0_bjd"]),
                duration_min=float(metrics["anchor_duration_min"]),
            )
            metrics[f"{prefix}_anchor_depth"] = depth
            metrics[f"{prefix}_anchor_snr"] = snr
            metrics[f"{prefix}_anchor_n_in"] = n_in
            if np.isfinite(metrics["anchor_t0_bjd"]) and np.isfinite(metrics["anchor_duration_min"]):
                for key, value in _even_odd_depth_metrics(
                    lcs[ap],
                    ap,
                    period_d=float(metrics["anchor_period_d"]),
                    t0_bjd=float(metrics["anchor_t0_bjd"]),
                    duration_min=float(metrics["anchor_duration_min"]),
                ).items():
                    metrics[f"{prefix}_anchor_{key}"] = value
        metrics[f"{prefix}_trend_ptp"] = binned_trend_ptp(lcs[ap].time, lcs[ap].flux[ap], lcs[ap].quality)
        metrics[f"{prefix}_branch_status"] = branch_meta[ap].get("status", "")
        metrics[f"{prefix}_branch_window_d"] = branch_meta[ap].get("window_d", np.nan)

    small_prefix = _aperture_prefix(apertures[0]) if len(apertures) >= 1 else ""
    primary_prefix = _aperture_prefix(apertures[1]) if len(apertures) >= 2 else ""
    small_p = metrics.get(f"{small_prefix}_period_d", np.nan)
    primary_p = metrics.get(f"{primary_prefix}_period_d", np.nan)
    if np.isfinite(small_p) and np.isfinite(primary_p) and small_p > 0:
        metrics["aperture_period_rel_delta"] = abs(float(primary_p) - float(small_p)) / float(small_p)
    else:
        metrics["aperture_period_rel_delta"] = np.nan
    small_depth = metrics.get(f"{small_prefix}_anchor_depth", np.nan)
    primary_depth = metrics.get(f"{primary_prefix}_anchor_depth", np.nan)
    metrics["aperture_depth_ratio_primary_over_small"] = (
        float(primary_depth) / float(small_depth)
        if np.isfinite(small_depth) and abs(float(small_depth)) > 0
        else np.nan
    )
    metrics["aperture_disagreement_flag"] = bool(
        np.isfinite(metrics["aperture_period_rel_delta"]) and metrics["aperture_period_rel_delta"] > 0.02
    )
    payload = {"results": results, "spectra": spectra, "lcs": lcs, "metrics": metrics}
    return metrics, payload


def measure_two_aperture_candidate_metadata(
    lc: HLSPLightCurve,
    *,
    anchor_peak: dict[str, Any],
    own_peaks: dict[str, dict[str, Any]],
    apertures: tuple[str, str],
) -> dict[str, Any]:
    """Measure scalar CNN metadata at fixed, already-searched ephemerides.

    Unlike :func:`build_two_aperture_metrics`, this helper never reruns BLS.
    It is used for large teacher-inference pools whose peak table has already
    been generated from the active ADP pair.
    """

    anchor = _finite_peak_dict(anchor_peak)
    if anchor is None:
        raise ValueError("anchor_peak has no finite ephemeris")
    metrics: dict[str, Any] = {
        "tic": int(lc.tic),
        "sector": int(lc.sector),
        "cam": int(lc.cam),
        "ccd": int(lc.ccd),
        "tmag": float(lc.tmag),
        "anchor_aperture": str(apertures[0]),
        "anchor_period_d": float(anchor["period_d"]),
        "anchor_t0_bjd": float(anchor["t0_bjd"]),
        "anchor_duration_min": float(anchor["duration_min"]),
        "anchor_sde": float(anchor.get("sde", np.nan)),
        "anchor_source": "precomputed_adp_peak",
    }
    for aperture in apertures:
        if aperture not in lc.flux:
            continue
        prefix = _aperture_prefix(aperture)
        own = _finite_peak_dict(own_peaks.get(aperture))
        if own is not None:
            for key, value in own.items():
                metrics[f"{prefix}_{key}"] = value
            for key, value in _even_odd_depth_metrics(
                lc,
                aperture,
                period_d=float(own["period_d"]),
                t0_bjd=float(own["t0_bjd"]),
                duration_min=float(own["duration_min"]),
            ).items():
                metrics[f"{prefix}_own_{key}"] = value
        depth, snr, n_in = _depth_at_ephemeris(
            lc,
            aperture,
            period_d=float(anchor["period_d"]),
            t0_bjd=float(anchor["t0_bjd"]),
            duration_min=float(anchor["duration_min"]),
        )
        metrics[f"{prefix}_anchor_depth"] = depth
        metrics[f"{prefix}_anchor_snr"] = snr
        metrics[f"{prefix}_anchor_n_in"] = n_in
        for key, value in _even_odd_depth_metrics(
            lc,
            aperture,
            period_d=float(anchor["period_d"]),
            t0_bjd=float(anchor["t0_bjd"]),
            duration_min=float(anchor["duration_min"]),
        ).items():
            metrics[f"{prefix}_anchor_{key}"] = value
        metrics[f"{prefix}_trend_ptp"] = binned_trend_ptp(
            lc.time, lc.flux[aperture], lc.quality
        )

    small_prefix = _aperture_prefix(apertures[0])
    primary_prefix = _aperture_prefix(apertures[1])
    small_period = float(metrics.get(f"{small_prefix}_period_d", np.nan))
    primary_period = float(metrics.get(f"{primary_prefix}_period_d", np.nan))
    metrics["aperture_period_rel_delta"] = (
        abs(primary_period - small_period) / small_period
        if np.isfinite(small_period) and small_period > 0 and np.isfinite(primary_period)
        else np.nan
    )
    small_depth = float(metrics.get(f"{small_prefix}_anchor_depth", np.nan))
    primary_depth = float(metrics.get(f"{primary_prefix}_anchor_depth", np.nan))
    metrics["aperture_depth_ratio_primary_over_small"] = (
        primary_depth / small_depth
        if np.isfinite(primary_depth) and np.isfinite(small_depth) and abs(small_depth) > 0
        else np.nan
    )
    metrics["aperture_disagreement_flag"] = bool(
        np.isfinite(metrics["aperture_period_rel_delta"])
        and metrics["aperture_period_rel_delta"] > 0.02
    )
    return metrics


def render_two_aperture_sheet(
    lc: HLSPLightCurve,
    out_path: Path,
    *,
    branch_name: str = DEFAULT_TWO_APERTURE_BRANCH,
    cfg: BLSConfig | None = None,
    apertures: tuple[str, ...] = DEFAULT_TWO_APERTURE_APERTURES,
    anchor_aperture: str | None = None,
    row_metadata: dict[str, Any] | None = None,
    use_row_ephemeris: bool = False,
    write_pdf: bool = True,
    harmonic_factors: tuple[float, ...] = (),
) -> tuple[Path, dict[str, Any]]:
    """Render a PNG/PDF two-aperture vet sheet and return metrics."""

    import matplotlib

    mpl_cache = Path(tempfile.gettempdir()) / "twirl_mplconfig"
    mpl_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_cache))
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    apply_twirl_style("full_page")
    branch = branch_by_name(branch_name)
    cfg = cfg or BLSConfig(apertures=apertures, n_periods=20_000, n_peaks=10)
    metrics, payload = build_two_aperture_metrics(
        lc,
        branch=branch,
        cfg=cfg,
        apertures=apertures,
        anchor_aperture=anchor_aperture,
        anchor_peak_override=_row_ephemeris_peak(row_metadata) if use_row_ephemeris else None,
    )
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    harmonic_factors = tuple(float(value) for value in harmonic_factors)
    if any(not np.isfinite(value) or value <= 0 for value in harmonic_factors):
        raise ValueError("harmonic factors must be finite and positive")
    if harmonic_factors:
        fig = plt.figure(figsize=(13.6, 14.4))
        outer = fig.add_gridspec(2, 1, height_ratios=[4.9, 1.75], hspace=0.22)
        gs = outer[0].subgridspec(
            5,
            4,
            height_ratios=[1.02, 0.82, 0.92, 0.92, 0.28],
            width_ratios=[0.82, 1.18, 0.82, 1.18],
            hspace=0.52,
            wspace=0.24,
        )
        harmonic_gs = outer[1].subgridspec(2, len(harmonic_factors), hspace=0.42, wspace=0.24)
    else:
        fig = plt.figure(figsize=(13.6, 10.6))
        gs = fig.add_gridspec(
            5,
            4,
            height_ratios=[1.02, 0.82, 0.92, 0.92, 0.28],
            width_ratios=[0.82, 1.18, 0.82, 1.18],
            hspace=0.52,
            wspace=0.24,
        )
        harmonic_gs = None
    apertures = [ap for ap in apertures if ap in payload["lcs"]]
    anchor_period = metrics.get("anchor_period_d", np.nan)
    anchor_t0 = metrics.get("anchor_t0_bjd", np.nan)
    anchor_dur = metrics.get("anchor_duration_min", np.nan)

    plot_peaks: dict[str, dict[str, float]] = {}
    for ap in apertures[:2]:
        peak = _best_peak_dict(payload["results"][ap])
        prefix = _aperture_prefix(ap)
        anchor_depth = metrics.get(f"{prefix}_anchor_depth", np.nan)
        if use_row_ephemeris and np.isfinite(anchor_period) and np.isfinite(anchor_t0) and np.isfinite(anchor_dur):
            peak = {
                "period_d": float(anchor_period),
                "t0_bjd": float(anchor_t0),
                "duration_min": float(anchor_dur),
                "depth": float(anchor_depth) if np.isfinite(anchor_depth) else np.nan,
                "depth_snr": metrics.get(f"{prefix}_anchor_snr", np.nan),
                "sde": metrics.get("anchor_sde", np.nan),
            }
        plot_peaks[ap] = peak
    shared_ylim = _shared_light_curve_ylim(
        {ap: payload["lcs"][ap] for ap in apertures[:2]},
        plot_peaks,
    )
    metrics["vet_sheet_version"] = TWO_APERTURE_VET_SHEET_VERSION
    metrics["shared_flux_ymin"] = float(shared_ylim[0])
    metrics["shared_flux_ymax"] = float(shared_ylim[1])
    metrics["matched_review_odd_even_sampling"] = True

    for col, ap in enumerate(apertures):
        if col > 1:
            break
        offset = 0 if col == 0 else 2
        blc: HLSPLightCurve = payload["lcs"][ap]
        peak = _best_peak_dict(payload["results"][ap])
        prefix = _aperture_prefix(ap)
        plot_peak = plot_peaks[ap]
        plot_title = "Review fold" if use_row_ephemeris else "Own fold"
        metric_scope = "anchor" if use_row_ephemeris else "own"

        ax = fig.add_subplot(gs[0, offset : offset + 2])
        _plot_full_lc(
            ax,
            blc,
            ap,
            branch_name=branch.name,
            peak=plot_peak,
            shared_ylim=shared_ylim,
        )

        ax = fig.add_subplot(gs[1, offset : offset + 2])
        if np.isfinite(plot_peak["period_d"]):
            _plot_full_phase_fold(
                ax,
                blc,
                ap,
                period_d=float(plot_peak["period_d"]),
                t0_bjd=float(plot_peak["t0_bjd"]),
                duration_min=float(plot_peak["duration_min"]),
                depth=float(plot_peak["depth"]) if np.isfinite(plot_peak.get("depth", np.nan)) else None,
                title="Full phase fold",
                shared_ylim=shared_ylim,
            )
        else:
            ax.text(0.5, 0.5, "No finite review ephemeris", ha="center", va="center")
            ax.set_axis_off()

        ax = fig.add_subplot(gs[2, offset])
        _plot_periodogram(
            ax,
            payload["spectra"].get(ap),
            bls_peak=peak,
            review_peak=plot_peak,
        )

        ax = fig.add_subplot(gs[2, offset + 1])
        if np.isfinite(plot_peak["period_d"]):
            _plot_folded_bins(
                ax,
                blc,
                ap,
                period_d=float(plot_peak["period_d"]),
                t0_bjd=float(plot_peak["t0_bjd"]),
                duration_min=float(plot_peak["duration_min"]),
                depth=float(plot_peak["depth"]) if np.isfinite(plot_peak.get("depth", np.nan)) else None,
                sde=float(plot_peak["sde"]) if np.isfinite(plot_peak.get("sde", np.nan)) else None,
                title=plot_title,
                show_duration_days=True,
                shared_ylim=shared_ylim,
            )
            ax.set_xlabel("")
        else:
            ax.text(0.5, 0.5, "No finite review ephemeris", ha="center", va="center")
            ax.set_axis_off()

        ax = fig.add_subplot(gs[3, offset])
        if np.isfinite(plot_peak["period_d"]):
            _plot_folded_bins(
                ax,
                blc,
                ap,
                period_d=float(plot_peak["period_d"]),
                t0_bjd=float(plot_peak["t0_bjd"]),
                duration_min=float(plot_peak["duration_min"]),
                depth=float(plot_peak["depth"]) if np.isfinite(plot_peak.get("depth", np.nan)) else None,
                sde=float(plot_peak["sde"]) if np.isfinite(plot_peak.get("sde", np.nan)) else None,
                title="Secondary fold",
                phase_shift=0.5,
                shared_ylim=shared_ylim,
            )
        else:
            ax.text(0.5, 0.5, "No finite review ephemeris", ha="center", va="center")
            ax.set_axis_off()

        ax = fig.add_subplot(gs[3, offset + 1])
        if np.isfinite(plot_peak["period_d"]):
            _plot_even_odd_bins(
                ax,
                blc,
                ap,
                period_d=float(plot_peak["period_d"]),
                t0_bjd=float(plot_peak["t0_bjd"]),
                duration_min=float(plot_peak["duration_min"]),
                title="Odd/even fold",
                metrics=metrics,
                prefix=prefix,
                metric_scope=metric_scope,
                shared_ylim=shared_ylim,
            )
        else:
            ax.text(0.5, 0.5, "No finite review ephemeris", ha="center", va="center")
            ax.set_axis_off()

        if harmonic_gs is not None:
            for factor_index, factor in enumerate(harmonic_factors):
                ax = fig.add_subplot(harmonic_gs[col, factor_index])
                harmonic_period = float(plot_peak["period_d"]) * float(factor)
                if np.isfinite(harmonic_period) and harmonic_period > 0:
                    _plot_harmonic_fold(
                        ax,
                        blc,
                        ap,
                        period_d=harmonic_period,
                        t0_bjd=float(plot_peak["t0_bjd"]),
                        duration_min=float(plot_peak["duration_min"]),
                        depth=float(plot_peak["depth"]) if np.isfinite(plot_peak.get("depth", np.nan)) else None,
                        title=f"{ap}: {_harmonic_factor_label(factor)} ({harmonic_period:.5g} d)",
                        show_ylabel=factor_index == 0,
                        shared_ylim=shared_ylim,
                    )
                else:
                    ax.text(0.5, 0.5, "No finite review ephemeris", ha="center", va="center")
                    ax.set_axis_off()

    ax_txt = fig.add_subplot(gs[4, :])
    ax_txt.axis("off")
    lines = [
        f"TIC {lc.tic}  S{lc.sector} cam{lc.cam}/ccd{lc.ccd}  T={lc.tmag:.2f}",
        (
            f"anchor={metrics.get('anchor_aperture', '')}  "
            f"P={metrics.get('anchor_period_d', np.nan):.6g} d  "
            f"T0={metrics.get('anchor_t0_bjd', np.nan):.6f}  "
            f"dur={_duration_days_text(metrics.get('anchor_duration_min', np.nan))}  "
            f"SDE={metrics.get('anchor_sde', np.nan):.3g}"
        ),
        (
            f"aperture period delta={metrics.get('aperture_period_rel_delta', np.nan):.3g}; "
            f"primary/small depth ratio={metrics.get('aperture_depth_ratio_primary_over_small', np.nan):.3g}; "
            f"aperture_disagreement={metrics.get('aperture_disagreement_flag')}"
        ),
    ]
    for ap in apertures[:2]:
        prefix = _aperture_prefix(ap)
        sigma_delta = metrics.get(f"{prefix}_own_even_odd_sigma_delta", np.nan)
        n_even = metrics.get(f"{prefix}_own_even_n_in", 0)
        n_odd = metrics.get(f"{prefix}_own_odd_n_in", 0)
        if np.isfinite(sigma_delta):
            lines.append(f"{ap} even/odd: |depth delta|={float(sigma_delta):.2g} sigma; n_even={int(n_even)}, n_odd={int(n_odd)}")
    ax_txt.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", family="monospace", fontsize=9.4)
    fig.subplots_adjust(top=0.980, bottom=0.035, left=0.060, right=0.992)
    fig.savefig(out_path, dpi=160)
    if write_pdf:
        try:
            fig.savefig(out_path.with_suffix(".pdf"))
        except Exception:
            pass
    plt.close(fig)
    metrics["twirl_vet_sheet_name"] = out_path.name
    metrics["twirl_vet_sheet_pdf_name"] = out_path.with_suffix(".pdf").name if write_pdf else ""
    metrics["harmonic_factors"] = ",".join(f"{value:g}" for value in harmonic_factors)
    return out_path, metrics
