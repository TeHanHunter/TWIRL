"""Flux-space BSpline cotrending for TWIRL.

Mirrors QLP `lctools detrend`'s knot-spaced cotrend, but operates on
linear flux instead of magnitude. The QLP recipe is a BSpline with
knots every ``bkspace`` days plus iterative outlier sigma-clipping;
this implementation preserves that structure so the cadence-by-cadence
output is comparable to QLP for the bright regime where both produce
sensible results.

Key behavioral differences from QLP:

* **Negative flux is fine.** We do not clip flux <= 0; the BSpline
  fits the actual flux distribution including legitimate
  background-subtraction excursions. For T >= 19 hosts this restores
  the ~50% of cadences QLP silently dropped.
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
    bkspace_d: float = 0.3      # BSpline knot spacing in days, matches QLP default
    k: int = 3                  # cubic
    sigma_clip: float = 5.0     # outlier rejection threshold for iterative fit
    max_iter: int = 5
    edge_pad_d: float = 0.5     # do not place knots within this distance of either
                                # data edge (avoids spline ringing at orbit
                                # boundaries that the QLP recipe also avoids)


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
    """Fit a BSpline to ``flux(time)`` and return the relative-flux
    light curve and per-cadence error.

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

    Returns
    -------
    det_flux : array, same length as ``flux``
        Detrended relative flux (flux / spline(time)). Median ~ 1.
        NaN at cadences where input ``flux`` was NaN.
    det_flux_err : array, same length
        MAD-RMS of the detrended series, broadcast. Constant across
        good cadences, NaN elsewhere.
    spline : array, same length
        The fitted cotrend curve evaluated at every input time, for
        diagnostic plots.

    Notes
    -----
    Iterative outlier rejection (``sigma_clip``, ``max_iter`` from
    ``cfg``) re-fits after masking |f - spline| > sigma_clip *
    1.4826 * MAD(f - spline). Converges typically in 2-3 iterations.
    """
    cfg = cfg or FluxDetrendConfig()
    time = np.asarray(time, dtype=np.float64)
    flux = np.asarray(flux, dtype=np.float64)
    n = len(flux)

    fit_mask = np.isfinite(time) & np.isfinite(flux)
    if quality is not None:
        fit_mask &= np.asarray(quality) == 0
    if flux_err is not None:
        fit_mask &= np.isfinite(flux_err) & (np.asarray(flux_err) > 0)

    if fit_mask.sum() < 50:
        # Insufficient data — return median-normalized flux without detrend.
        med = np.nanmedian(flux[fit_mask]) if fit_mask.any() else np.nan
        det = flux / med if np.isfinite(med) and med != 0 else np.full(n, np.nan)
        return det, _broadcast_mad(det), np.full(n, med)

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
        det_flux = np.where(np.isfinite(flux), flux / spline_at, np.nan)
        return det_flux, _broadcast_mad(det_flux), spline_at

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
        mad = np.nanmedian(np.abs(resid - np.nanmedian(resid)))
        sigma = 1.4826 * mad if mad > 0 else np.std(resid)
        new_good = np.abs(resid) < cfg.sigma_clip * sigma
        if np.array_equal(new_good, good):
            break
        good = new_good

    if spl is None:
        # Never produced a spline — polynomial fallback.
        spline_at = _safe_poly(t_fit, f_fit, w_fit, time, flux)
    else:
        spline_at = spl(time)

    det_flux = np.where(np.isfinite(flux), flux / spline_at, np.nan)
    det_flux_err = _broadcast_mad(det_flux)
    return det_flux, det_flux_err, spline_at


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
