"""Second-pass ADP detrending utilities for TWIRL vetting.

The ADP columns in the S56 compare product are already relative-flux light
curves.  The helpers here apply a conservative high-pass correction on top of
those columns without redefining the official light-curve product.
"""
from __future__ import annotations

from dataclasses import dataclass, fields

import numpy as np
from scipy.ndimage import median_filter

from twirl.io.hlsp import BJDREFI, HLSPLightCurve
from twirl.lightcurves.detrend_presets import TWIRL_FS_V2_ADP015Q_BRANCH


ADP_VET_APERTURES: tuple[str, str] = ("DET_FLUX_ADP_SML", "DET_FLUX_ADP")


@dataclass(frozen=True)
class ADPPlusBranch:
    name: str
    window_d: float | None
    gap_split_d: float = 0.2
    min_window_d: float = 0.12
    max_window_d: float = 0.50
    duration_window_factor: float = 12.0
    mask_duration_factor: float = 3.0
    mask_min_width_min: float = 20.0

    @property
    def is_identity(self) -> bool:
        return self.window_d is None and self.name in {
            "adp",
            "current",
            "current_adp",
            "identity",
            TWIRL_FS_V2_ADP015Q_BRANCH,
        }


DEFAULT_ADPPLUS_BRANCHES: tuple[ADPPlusBranch, ...] = (
    ADPPlusBranch("identity", None),
    ADPPlusBranch("current_adp", None),
    ADPPlusBranch(TWIRL_FS_V2_ADP015Q_BRANCH, None),
    ADPPlusBranch("adpplus_0p20", 0.20),
    ADPPlusBranch("adpplus_0p12", 0.12),
    ADPPlusBranch("adpplus_adaptive", None),
)


def branch_by_name(name: str) -> ADPPlusBranch:
    for branch in DEFAULT_ADPPLUS_BRANCHES:
        if branch.name == name:
            return branch
    raise KeyError(f"unknown ADP+ branch: {name}")


def robust_sigma(values: np.ndarray) -> float:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return float("nan")
    med = float(np.nanmedian(vals))
    mad = float(np.nanmedian(np.abs(vals - med)))
    sigma = 1.4826 * mad
    if np.isfinite(sigma) and sigma > 0:
        return float(sigma)
    std = float(np.nanstd(vals))
    return std if np.isfinite(std) and std > 0 else float("nan")


def transit_window_mask(
    time_d: np.ndarray,
    *,
    period_d: float | None,
    t0_bjd: float | None,
    duration_min: float | None,
    width_factor: float = 3.0,
    min_width_min: float = 20.0,
) -> np.ndarray:
    """Return a widened periodic transit mask in BTJD coordinates."""

    time = np.asarray(time_d, dtype=float)
    out = np.zeros(time.shape, dtype=bool)
    if period_d is None or t0_bjd is None or duration_min is None:
        return out
    if not np.isfinite(period_d) or not np.isfinite(t0_bjd) or not np.isfinite(duration_min):
        return out
    if period_d <= 0 or duration_min <= 0:
        return out
    width_min = max(float(min_width_min), float(width_factor) * float(duration_min))
    half_d = 0.5 * width_min / 1440.0
    t0_d = float(t0_bjd) - float(BJDREFI)
    phase_d = ((time - t0_d + 0.5 * float(period_d)) % float(period_d)) - 0.5 * float(period_d)
    return np.isfinite(time) & (np.abs(phase_d) <= half_d)


def _time_gap_segments(time: np.ndarray, gap_split_d: float) -> list[np.ndarray]:
    finite = np.isfinite(time)
    if finite.sum() <= 1 or not np.isfinite(gap_split_d) or gap_split_d <= 0:
        return [np.arange(len(time), dtype=int)]
    order = np.argsort(np.where(finite, time, np.inf), kind="stable")
    finite_order = order[finite[order]]
    cuts = np.where(np.diff(time[finite_order]) > gap_split_d)[0] + 1
    if cuts.size == 0:
        return [np.sort(finite_order)]
    starts = np.concatenate(([0], cuts))
    ends = np.concatenate((cuts, [finite_order.size]))
    return [np.sort(finite_order[start:end]) for start, end in zip(starts, ends)]


def _odd_window_from_days(time: np.ndarray, window_d: float) -> int:
    t = np.asarray(time, dtype=float)
    t = t[np.isfinite(t)]
    if t.size < 2:
        return 5
    dt = np.diff(np.sort(t))
    dt = dt[np.isfinite(dt) & (dt > 0)]
    cadence_d = float(np.nanmedian(dt)) if dt.size else np.nan
    if not np.isfinite(cadence_d) or cadence_d <= 0:
        return 5
    width = max(5, int(np.ceil(float(window_d) / cadence_d)))
    if width % 2 == 0:
        width += 1
    return width


def _interp_fit_values(time: np.ndarray, flux: np.ndarray, fit: np.ndarray) -> np.ndarray:
    out = np.full(len(flux), np.nan, dtype=float)
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


def branch_window_d(branch: ADPPlusBranch, duration_min: float | None) -> float:
    if branch.is_identity:
        return float("nan")
    if branch.name != "adpplus_adaptive" and branch.window_d is not None:
        return float(branch.window_d)
    duration = float(duration_min) if duration_min is not None and np.isfinite(duration_min) else 10.0
    window = branch.duration_window_factor * max(duration, 3.0) / 1440.0
    return float(np.clip(window, branch.min_window_d, branch.max_window_d))


def apply_adpplus_branch(
    time_d: np.ndarray,
    flux: np.ndarray,
    quality: np.ndarray,
    branch: ADPPlusBranch,
    *,
    period_d: float | None = None,
    t0_bjd: float | None = None,
    duration_min: float | None = None,
) -> tuple[np.ndarray, dict[str, float | str | int]]:
    """Return the branch flux and small metadata dict.

    For ADP+ branches, the baseline is fitted after masking the supplied
    ephemeris.  If no ephemeris is supplied, the branch is still usable as a
    blind high-pass, but the metadata records that no transit mask was used.
    """

    time = np.asarray(time_d, dtype=float)
    y = np.asarray(flux, dtype=float)
    q = np.asarray(quality)
    if branch.is_identity:
        return y.copy(), {
            "branch": branch.name,
            "status": "identity",
            "window_d": float("nan"),
            "n_masked": 0,
        }

    mask = transit_window_mask(
        time,
        period_d=period_d,
        t0_bjd=t0_bjd,
        duration_min=duration_min,
        width_factor=branch.mask_duration_factor,
        min_width_min=branch.mask_min_width_min,
    )
    window_d = branch_window_d(branch, duration_min)
    out = np.full_like(y, np.nan, dtype=float)
    statuses: list[str] = []
    for idx in _time_gap_segments(time, branch.gap_split_d):
        segment_fit = np.isfinite(time[idx]) & np.isfinite(y[idx]) & (q[idx] == 0) & ~mask[idx]
        if np.count_nonzero(segment_fit) < 20:
            med = float(np.nanmedian(y[idx][np.isfinite(y[idx])])) if np.isfinite(y[idx]).any() else np.nan
            baseline = np.full(idx.size, med, dtype=float)
            statuses.append("constant_insufficient_fit")
        else:
            t_seg = time[idx]
            y_seg = y[idx]
            filled = _interp_fit_values(t_seg, y_seg, segment_fit)
            width = _odd_window_from_days(t_seg, window_d)
            width = min(width, len(filled) if len(filled) % 2 else max(1, len(filled) - 1))
            if width < 3:
                med = float(np.nanmedian(y_seg[segment_fit]))
                baseline = np.full(idx.size, med, dtype=float)
                statuses.append("constant_short_filter")
            else:
                baseline = median_filter(filled, size=width, mode="nearest")
                statuses.append(f"median_window{width}")
        det = 1.0 + (y[idx] - baseline)
        good = np.isfinite(det) & (q[idx] == 0) & ~mask[idx]
        if np.count_nonzero(good) >= 10:
            det = det + (1.0 - float(np.nanmedian(det[good])))
        out[idx] = det
    return out, {
        "branch": branch.name,
        "status": ";".join(statuses),
        "window_d": float(window_d),
        "n_masked": int(mask.sum()),
    }


def branched_light_curve(
    lc: HLSPLightCurve,
    aperture: str,
    branch: ADPPlusBranch,
    *,
    period_d: float | None = None,
    t0_bjd: float | None = None,
    duration_min: float | None = None,
) -> tuple[HLSPLightCurve, dict[str, float | str | int]]:
    """Return an HLSPLightCurve containing one aperture after branch detrending."""

    if aperture not in lc.flux:
        raise KeyError(f"{aperture} not loaded in {lc.path}")
    flux, meta = apply_adpplus_branch(
        lc.time,
        lc.flux[aperture],
        lc.quality,
        branch,
        period_d=period_d,
        t0_bjd=t0_bjd,
        duration_min=duration_min,
    )
    payload = {
        "tic": lc.tic,
        "tmag": lc.tmag,
        "sector": lc.sector,
        "cam": lc.cam,
        "ccd": lc.ccd,
        "ra": getattr(lc, "ra", float("nan")),
        "dec": getattr(lc, "dec", float("nan")),
        "time": lc.time,
        "cadenceno": lc.cadenceno,
        "orbitid": lc.orbitid,
        "quality": lc.quality,
        "flux": {aperture: flux},
        "path": lc.path,
    }
    accepted = {field.name for field in fields(HLSPLightCurve)}
    out = HLSPLightCurve(**{key: value for key, value in payload.items() if key in accepted})
    return out, meta


def binned_trend_ptp(
    time_d: np.ndarray,
    flux: np.ndarray,
    quality: np.ndarray,
    *,
    bin_d: float = 0.25,
) -> float:
    good = np.isfinite(time_d) & np.isfinite(flux) & (quality == 0)
    if np.count_nonzero(good) < 50:
        return float("nan")
    t = np.asarray(time_d, dtype=float)[good]
    y = np.asarray(flux, dtype=float)[good]
    lo = float(np.nanmin(t))
    hi = float(np.nanmax(t))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return float("nan")
    edges = np.arange(lo, hi + bin_d, bin_d)
    if edges.size < 4:
        edges = np.linspace(lo, hi, 4)
    medians: list[float] = []
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (t >= a) & (t < b)
        if np.count_nonzero(sel) >= 10:
            medians.append(float(np.nanmedian(y[sel])))
    vals = np.asarray(medians, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size < 3:
        return float("nan")
    vals = vals - float(np.nanmedian(vals))
    return float(np.nanpercentile(vals, 90) - np.nanpercentile(vals, 10))
