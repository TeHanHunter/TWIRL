"""Flux-space BSpline cotrending for TWIRL.

Uses the same knot-spaced BSpline cotrend family as QLP `lctools
detrend`, but operates on linear flux instead of magnitude. QLP's
historical knot spacing is 0.3 d; TWIRL's default is deliberately longer
after injection sweeps showed 0.3 d absorbs multi-hour transit/eclipse
signals in faint zero-crossing light curves.

Key behavioral differences from QLP:

* **Negative flux is fine.** We do not clip flux <= 0; the BSpline
  fits the actual flux distribution including legitimate
  background-subtraction excursions. For T >= 19 hosts this restores
  the ~50% of cadences QLP silently dropped.
* **Subtractive residual, not ratio.** Detrended output is
  ``1 + (flux - spline) / scale`` instead of ``flux / spline``.
  Faint WDs frequently have splines passing through zero
  (background-subtracted flux near zero); the ratio form blows up
  there, while the subtractive form stays finite.
* **Robust positive normalization for faint targets.** The default
  ``scale_strategy="auto"`` uses ``abs(median(flux_fit))`` only when
  that median is well above the robust flux scatter. Otherwise it falls
  back to a positive robust absolute-flux scale. This avoids both
  near-zero blow-ups and sign inversions when the background-subtracted
  good-cadence median is negative.
* **Quality mask gates fit weights, not the input.** We pass all
  finite-flux cadences to the BSpline but down-weight QUALITY != 0
  cadences (or exclude them via fit weights). This way the smoothed
  curve is anchored on good cadences but the detrended output covers
  the full time grid.
* **Per-cadence error from MAD.** Detrended flux uncertainty is
  ``1.4826 * MAD`` of the detrended series, broadcast across all
  cadences. Same convention as ``twirl.io.hlsp.tglc_mad_error``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import LSQUnivariateSpline


@dataclass
class FluxDetrendConfig:
    bkspace_d: float = 0.8      # BSpline knot spacing in days; chosen to preserve
                                # 5 min-6 hr injected events in the S56 faint sweep
    k: int = 3                  # cubic
    sigma_clip: float = 5.0     # outlier rejection threshold for iterative fit
    max_iter: int = 5
    edge_pad_d: float = 0.5     # do not place knots within this distance of either
                                # data edge (avoids spline ringing at orbit
                                # boundaries that the QLP recipe also avoids)
    output_mode: str = "subtractive"  # "subtractive" or "divisive"
    scale_strategy: str = "auto"      # "median", "median_abs", or "auto"
    min_scale_abs: float = 1.0e-12
    min_scale_snr: float = 3.0        # auto mode trusts median only above this S/N


@dataclass
class FluxDetrendResult:
    det_flux: np.ndarray
    det_flux_err: np.ndarray
    cotrend: np.ndarray
    scale: float
    scale_source: str
    fit_count: int
    output_mode: str
    scale_strategy: str


def _build_knots(time: np.ndarray, bkspace_d: float, edge_pad_d: float) -> np.ndarray:
    """Interior knots for LSQUnivariateSpline, every ``bkspace_d`` days,
    padded ``edge_pad_d`` from the time range endpoints.
    """
    t0, t1 = time.min(), time.max()
    interior_start = t0 + edge_pad_d
    interior_end = t1 - edge_pad_d
    if interior_end <= interior_start:
        return np.array([], dtype=float)
    n_knots = max(1, int(np.floor((interior_end - interior_start) / bkspace_d)))
    return np.linspace(interior_start, interior_end, n_knots + 1)[1:-1]


def flux_space_detrend(
    time: np.ndarray,
    flux: np.ndarray,
    quality: np.ndarray | None = None,
    flux_err: np.ndarray | None = None,
    cfg: FluxDetrendConfig | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    result = flux_space_detrend_result(time, flux, quality=quality, flux_err=flux_err, cfg=cfg)
    return result.det_flux, result.det_flux_err, result.cotrend


def flux_space_detrend_result(
    time: np.ndarray,
    flux: np.ndarray,
    quality: np.ndarray | None = None,
    flux_err: np.ndarray | None = None,
    cfg: FluxDetrendConfig | None = None,
) -> FluxDetrendResult:
    """Fit a BSpline to ``flux(time)`` and return the relative-flux
    light curve, per-cadence error, cotrend curve, and scale diagnostics.

    Parameters
    ----------
    time : array, BTJD
        Sorted, finite. Length N.
    flux : array, linear units (e-/s or relative)
        May contain NaN and negative values. Length N.
    quality : array of int, optional
        TESS QUALITY flags; cadences with quality != 0 are excluded
        from the BSpline fit weighting but kept in the output.
    flux_err : array, optional
        Per-cadence input uncertainty for the fit. If None, all
        finite cadences get equal weight.
    cfg : FluxDetrendConfig

    Result fields
    -------------
    det_flux : array, same length as ``flux``
        Detrended relative flux: ``1 + (flux - spline(time)) /
        scale``. Median ~ 1 after downstream recentering, centered on
        a global positive scale rather than on the local spline value,
        so cadences where the spline crosses zero stay finite. NaN at
        cadences where input ``flux`` was NaN.
    det_flux_err : array, same length
        MAD-RMS of the detrended series, broadcast. Constant across
        good cadences, NaN elsewhere.
    cotrend : array, same length
        The fitted cotrend curve evaluated at every input time, for
        diagnostic plots.
    scale : float
        Denominator used by subtractive modes. NaN for divisive output.
    scale_source : str
        Diagnostic label for how ``scale`` was selected.

    Notes
    -----
    Iterative outlier rejection (``sigma_clip``, ``max_iter`` from
    ``cfg``) re-fits after masking |f - spline| > sigma_clip *
    1.4826 * MAD(f - spline). Converges typically in 2-3 iterations.
    """
    cfg = cfg or FluxDetrendConfig()
    _validate_cfg(cfg)
    time = np.asarray(time, dtype=np.float64)
    flux = np.asarray(flux, dtype=np.float64)
    n = len(flux)

    fit_mask = np.isfinite(time) & np.isfinite(flux)
    if quality is not None:
        fit_mask &= np.asarray(quality) == 0
    if flux_err is not None:
        fit_mask &= np.isfinite(flux_err) & (np.asarray(flux_err) > 0)

    if fit_mask.sum() < 50:
        # Insufficient data -- return the requested relative form against a
        # constant cotrend, with the same scale guard as the full path.
        med = np.nanmedian(flux[fit_mask]) if fit_mask.any() else np.nan
        cotrend = np.full(n, med, dtype=np.float64)
        det, scale, scale_source = _relative_flux(
            flux=flux,
            cotrend=cotrend,
            fit_flux=flux[fit_mask],
            cfg=cfg,
        )
        return FluxDetrendResult(
            det_flux=det,
            det_flux_err=_broadcast_mad(det),
            cotrend=cotrend,
            scale=scale,
            scale_source=scale_source,
            fit_count=int(fit_mask.sum()),
            output_mode=cfg.output_mode,
            scale_strategy=cfg.scale_strategy,
        )

    t_fit = time[fit_mask]
    f_fit = flux[fit_mask]
    w_fit = (
        1.0 / np.asarray(flux_err)[fit_mask]
        if flux_err is not None
        else np.ones_like(f_fit)
    )

    # Sort by time (required by LSQUnivariateSpline).
    order = np.argsort(t_fit)
    t_fit, f_fit, w_fit = t_fit[order], f_fit[order], w_fit[order]

    # LSQUnivariateSpline needs strictly increasing knots inside the data
    # range AND strictly increasing x — dedupe pathological identical times.
    uniq_mask = np.concatenate(([True], np.diff(t_fit) > 0))
    t_fit, f_fit, w_fit = t_fit[uniq_mask], f_fit[uniq_mask], w_fit[uniq_mask]

    knots = _build_knots(t_fit, cfg.bkspace_d, cfg.edge_pad_d)
    if len(knots) < 2:
        # Time baseline shorter than 3 * bkspace; fall back to a polynomial.
        spline_at = _safe_poly(t_fit, f_fit, w_fit, time, flux)
        det_flux, scale, scale_source = _relative_flux(
            flux=flux,
            cotrend=spline_at,
            fit_flux=f_fit,
            cfg=cfg,
        )
        return FluxDetrendResult(
            det_flux=det_flux,
            det_flux_err=_broadcast_mad(det_flux),
            cotrend=spline_at,
            scale=scale,
            scale_source=scale_source,
            fit_count=int(fit_mask.sum()),
            output_mode=cfg.output_mode,
            scale_strategy=cfg.scale_strategy,
        )

    good = np.ones_like(f_fit, dtype=bool)
    spl = None  # last successful spline; None means every iter raised
    for _ in range(cfg.max_iter):
        try:
            new_spl = LSQUnivariateSpline(
                t_fit[good], f_fit[good], t=knots, k=cfg.k, w=w_fit[good]
            )
        except Exception:
            # Knot configuration invalid (e.g. surviving t_fit[good] no longer
            # spans the knot range after sigma-clip). Keep the last successful
            # spl (if any) and stop iterating.
            break
        spl = new_spl
        resid = f_fit - spl(t_fit)
        sigma = _robust_sigma(resid)
        if not np.isfinite(sigma) or sigma <= 0:
            break
        new_good = np.abs(resid) < cfg.sigma_clip * sigma
        if np.array_equal(new_good, good):
            break
        good = new_good

    if spl is None:
        # Never produced a spline — polynomial fallback.
        spline_at = _safe_poly(t_fit, f_fit, w_fit, time, flux)
    else:
        spline_at = spl(time)

    det_flux, scale, scale_source = _relative_flux(
        flux=flux,
        cotrend=spline_at,
        fit_flux=f_fit,
        cfg=cfg,
    )
    det_flux_err = _broadcast_mad(det_flux)
    return FluxDetrendResult(
        det_flux=det_flux,
        det_flux_err=det_flux_err,
        cotrend=spline_at,
        scale=scale,
        scale_source=scale_source,
        fit_count=int(fit_mask.sum()),
        output_mode=cfg.output_mode,
        scale_strategy=cfg.scale_strategy,
    )


def _validate_cfg(cfg: FluxDetrendConfig) -> None:
    if cfg.output_mode not in {"subtractive", "divisive"}:
        raise ValueError(f"unknown output_mode: {cfg.output_mode!r}")
    if cfg.scale_strategy not in {"median", "median_abs", "auto"}:
        raise ValueError(f"unknown scale_strategy: {cfg.scale_strategy!r}")
    if cfg.min_scale_abs <= 0:
        raise ValueError("min_scale_abs must be positive")
    if cfg.min_scale_snr <= 0:
        raise ValueError("min_scale_snr must be positive")


def _relative_flux(
    *,
    flux: np.ndarray,
    cotrend: np.ndarray,
    fit_flux: np.ndarray,
    cfg: FluxDetrendConfig,
) -> tuple[np.ndarray, float, str]:
    """Convert raw flux + cotrend to relative output for one configured mode."""
    det = np.full_like(flux, np.nan, dtype=np.float64)
    finite = np.isfinite(flux)

    if cfg.output_mode == "divisive":
        ok = finite & np.isfinite(cotrend) & (np.abs(cotrend) > cfg.min_scale_abs)
        with np.errstate(divide="ignore", invalid="ignore"):
            det[ok] = flux[ok] / cotrend[ok]
        det[~np.isfinite(det)] = np.nan
        return det, np.nan, "local_cotrend"

    scale, scale_source = _choose_subtractive_scale(fit_flux, cfg)
    if not np.isfinite(scale) or scale == 0:
        return det, scale, scale_source
    with np.errstate(divide="ignore", invalid="ignore"):
        det[finite] = 1.0 + (flux[finite] - cotrend[finite]) / scale
    det[~np.isfinite(det)] = np.nan
    return det, float(scale), scale_source


def _choose_subtractive_scale(
    fit_flux: np.ndarray,
    cfg: FluxDetrendConfig,
) -> tuple[float, str]:
    finite = np.asarray(fit_flux, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return np.nan, "no_finite_fit_flux"

    med = float(np.nanmedian(finite))
    if cfg.scale_strategy == "median":
        if np.isfinite(med) and abs(med) > cfg.min_scale_abs:
            return med, "median"
        return np.nan, "median_unusable"

    med_abs = abs(med)
    if cfg.scale_strategy == "median_abs":
        if np.isfinite(med_abs) and med_abs > cfg.min_scale_abs:
            return med_abs, "median_abs"
        return np.nan, "median_abs_unusable"

    sigma = _robust_sigma(finite)
    robust_abs = _robust_abs_scale(finite)
    floor = np.nanmax([cfg.min_scale_abs, sigma, robust_abs])
    if not np.isfinite(floor) or floor <= 0:
        return np.nan, "auto_unusable"
    snr_floor = cfg.min_scale_abs
    if np.isfinite(sigma) and sigma > 0:
        snr_floor = max(snr_floor, cfg.min_scale_snr * sigma)
    if np.isfinite(med_abs) and med_abs >= snr_floor:
        return med_abs, "auto_median_abs"
    return float(floor), "auto_robust_abs"


def _robust_sigma(x: np.ndarray) -> float:
    finite = np.asarray(x, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return np.nan
    med = np.nanmedian(finite)
    mad = np.nanmedian(np.abs(finite - med))
    sigma = 1.4826 * float(mad)
    if np.isfinite(sigma) and sigma > 0:
        return sigma
    sigma = float(np.nanstd(finite))
    return sigma if np.isfinite(sigma) else np.nan


def _robust_abs_scale(x: np.ndarray) -> float:
    finite = np.asarray(x, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return np.nan
    abs_flux = np.abs(finite)
    scale = float(np.nanpercentile(abs_flux, 68.0))
    if np.isfinite(scale) and scale > 0:
        return scale
    scale = float(np.nanmedian(abs_flux))
    return scale if np.isfinite(scale) else np.nan


def _safe_poly(
    t_fit: np.ndarray,
    f_fit: np.ndarray,
    w_fit: np.ndarray,
    time: np.ndarray,
    flux: np.ndarray,
) -> np.ndarray:
    """Polynomial cotrend fallback. Tries deg-2 weighted; degrades to
    deg-1 unweighted; ultimately returns the median, so the caller
    always gets a finite ``spline_at`` of the right shape.
    """
    for deg, use_w in ((2, True), (1, True), (1, False), (0, False)):
        try:
            if use_w:
                coeffs = np.polyfit(t_fit, f_fit, deg=deg, w=w_fit)
            else:
                coeffs = np.polyfit(t_fit, f_fit, deg=deg)
            return np.polyval(coeffs, time)
        except Exception:
            continue
    med = np.nanmedian(flux)
    return np.full(len(time), med if np.isfinite(med) else 1.0)


def _broadcast_mad(x: np.ndarray) -> np.ndarray:
    """1.4826 * MAD of the finite values, broadcast to the input shape.
    NaN where the input is NaN. Matches ``twirl.io.hlsp.tglc_mad_error``.
    """
    finite = np.isfinite(x)
    if not finite.any():
        return np.full_like(x, np.nan)
    med = np.nanmedian(x[finite])
    mad = np.nanmedian(np.abs(x[finite] - med))
    sigma = 1.4826 * float(mad)
    out = np.full_like(x, np.nan, dtype=np.float64)
    out[finite] = sigma
    return out
