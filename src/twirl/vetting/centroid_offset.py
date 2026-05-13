"""TGLC-centroid on-target test — the TWIRL-canonical pixel-vetting check.

For a transit candidate with parameters (period, t0, duration), compare
the median per-cadence centroid during in-transit cadences to the median
during nearby out-of-transit cadences. An on-target signal leaves the
flux-weighted centroid unchanged; a background-EB contaminator pulls the
centroid toward the target (i.e. away from the contaminator) during
in-transit cadences, by an amount that scales with the contaminator's
distance and apparent depth.

Inputs are the columns already in our QLP-style HLSP FITS (read via
:func:`twirl.io.hlsp.read_hlsp`):

    TIME, SAP_X, SAP_Y, QUALITY, DET_FLUX  (and a (P, t0, dur) from BLS)

So the test runs on **every** TWIRL HLSP — no SPOC TPF dependency.

Why not LEO-Vetter's pixel ``offset`` test?
  LEO's :func:`leo_vetter.pixel.pixel_vetting` does PRF-fit difference
  imaging on SPOC Target Pixel Files. That is the gold standard but
  requires the target to be on the SPOC 2-min cadence list, which most
  faint WDs (T >= 16) are NOT. The TGLC flux-weighted centroid is less
  sensitive to sub-pixel offsets than a PRF fit but works for 100% of
  targets at the FFI cadence. Use this for production survey vetting;
  reserve LEO's offset for the small number of WDs with SPOC TPFs and
  any final candidate going to follow-up.

Edge cases this function handles explicitly:

* **Noisy centroid baseline**: if ``sigma_oot`` exceeds
  ``max_useful_sigma_pix`` (default 1.0 pix), the centroid is too noisy
  for the test to discriminate; we report ``status='uninformative'``
  rather than 'on_target' or 'off_target'. This happens for some
  faint or scattered-light-affected TICs.
* **Few in-transit cadences**: if there are < ``min_in_transit`` good
  in-transit cadences, status='too_few_in_transit'.
* **NaN/invalid centroids**: filtered before the comparison.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class CentroidResult:
    delta_x_pix: float       # median(x_in_transit) - median(x_oot)
    delta_y_pix: float
    delta_pix: float         # sqrt(dx^2 + dy^2)
    sigma_x_oot_pix: float   # std of x in OOT band (per-cadence noise)
    sigma_y_oot_pix: float
    sigma_oot_pix: float     # sqrt(sigma_x^2 + sigma_y^2); the "noise floor"
                              # for |Δ| in the same units
    n_in_transit: int
    n_oot: int
    z_score: float           # |Δ| / sigma_oot; values > 3 are "off-target"
    status: str              # 'on_target' | 'off_target' | 'uninformative' |
                              # 'too_few_in_transit' | 'no_data'
    centroid_pass: bool      # True iff status == 'on_target'


def centroid_in_out_shift(
    time_bjd: np.ndarray,
    sap_x: np.ndarray,
    sap_y: np.ndarray,
    quality: np.ndarray,
    period_d: float,
    t0_bjd: float,
    duration_min: float,
    n_sigma_threshold: float = 3.0,
    oot_window_factor: float = 3.0,
    min_in_transit: int = 5,
    min_oot: int = 20,
    max_useful_sigma_pix: float = 1.0,
) -> CentroidResult:
    """Compute the in-transit vs out-of-transit centroid shift.

    Parameters
    ----------
    time_bjd : array
        Per-cadence BJD (not BTJD).
    sap_x, sap_y : array
        Per-cadence flux-weighted centroid in pixels (TGLC ``SAP_X``,
        ``SAP_Y``).
    quality : array
        TESS QUALITY flags; cadences with quality != 0 are excluded.
    period_d, t0_bjd : float
        BLS-fitted period and epoch.
    duration_min : float
        BLS-fitted transit duration (minutes).
    n_sigma_threshold : float
        |Δ| > n_sigma_threshold * sigma_oot → off-target. Default 3.
    oot_window_factor : float
        Out-of-transit comparison band is
        ``[duration, oot_window_factor * duration]`` from transit center
        (i.e. immediately adjacent in phase, not the full LC). Default 3.
    min_in_transit, min_oot : int
        Minimum cadence counts for the test to be informative.
    max_useful_sigma_pix : float
        If ``sigma_oot_pix`` exceeds this, the centroid is too noisy
        to discriminate; status = 'uninformative'. Default 1.0 pix.

    Returns
    -------
    :class:`CentroidResult`
    """
    time = np.asarray(time_bjd, dtype=float)
    x = np.asarray(sap_x, dtype=float)
    y = np.asarray(sap_y, dtype=float)
    q = np.asarray(quality)
    keep = (q == 0) & np.isfinite(time) & np.isfinite(x) & np.isfinite(y)

    if keep.sum() < (min_in_transit + min_oot):
        return CentroidResult(
            np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
            int(keep.sum()), 0, np.nan, "no_data", False,
        )

    t, x, y = time[keep], x[keep], y[keep]
    dur_d = duration_min / 1440.0
    phase = ((t - t0_bjd + 0.5 * period_d) % period_d) - 0.5 * period_d

    in_transit = np.abs(phase) < 0.5 * dur_d
    oot = (np.abs(phase) > dur_d) & (np.abs(phase) < oot_window_factor * dur_d)

    n_in = int(in_transit.sum())
    n_oot = int(oot.sum())

    if n_in < min_in_transit:
        return CentroidResult(
            np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
            n_in, n_oot, np.nan, "too_few_in_transit", False,
        )
    if n_oot < min_oot:
        return CentroidResult(
            np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
            n_in, n_oot, np.nan, "no_data", False,
        )

    dx = float(np.median(x[in_transit]) - np.median(x[oot]))
    dy = float(np.median(y[in_transit]) - np.median(y[oot]))
    delta = float(np.hypot(dx, dy))
    sx = float(np.std(x[oot]))
    sy = float(np.std(y[oot]))
    sigma_oot = float(np.hypot(sx, sy))

    if not np.isfinite(sigma_oot) or sigma_oot == 0:
        return CentroidResult(
            dx, dy, delta, sx, sy, sigma_oot, n_in, n_oot,
            np.nan, "no_data", False,
        )

    if sigma_oot > max_useful_sigma_pix:
        # Centroid scatter is too large for this test to be discriminating.
        return CentroidResult(
            dx, dy, delta, sx, sy, sigma_oot, n_in, n_oot,
            delta / sigma_oot, "uninformative", False,
        )

    z = delta / sigma_oot
    if z > n_sigma_threshold:
        return CentroidResult(
            dx, dy, delta, sx, sy, sigma_oot, n_in, n_oot,
            z, "off_target", False,
        )
    return CentroidResult(
        dx, dy, delta, sx, sy, sigma_oot, n_in, n_oot,
        z, "on_target", True,
    )
