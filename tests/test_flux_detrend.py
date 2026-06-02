from __future__ import annotations

import numpy as np

from twirl.lightcurves.flux_detrend import (
    FluxDetrendConfig,
    flux_space_detrend_result,
)


def _binned_amp(time: np.ndarray, flux: np.ndarray, bin_d: float = 0.25) -> float:
    vals = []
    edges = np.arange(float(np.nanmin(time)), float(np.nanmax(time)) + bin_d, bin_d)
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (time >= lo) & (time < hi) & np.isfinite(flux)
        if m.sum() >= 10:
            vals.append(float(np.nanmedian(flux[m])))
    vals = np.asarray(vals)
    return float(np.nanpercentile(vals, 95) - np.nanpercentile(vals, 5))


def test_gap_split_prevents_inter_orbit_polynomial_fallback() -> None:
    rng = np.random.default_rng(42)
    t1 = np.linspace(0.0, 4.0, 900)
    t2 = np.linspace(10.0, 14.0, 950)
    trend1 = 900.0 * np.exp(-t1 / 0.7) - 250.0 + 30.0 * np.sin(2.0 * t1)
    trend2 = (
        700.0 * np.exp(-(t2 - 10.0) / 0.8)
        - 150.0
        + 20.0 * np.sin(2.0 * (t2 - 10.0))
    )
    flux = np.concatenate([trend1, trend2]) + rng.normal(0.0, 80.0, t1.size + t2.size)
    time = np.concatenate([t1, t2])
    quality = np.zeros_like(time, dtype=int)
    flux_err = np.full_like(time, 80.0)

    no_split = flux_space_detrend_result(
        time,
        flux,
        quality=quality,
        flux_err=flux_err,
        cfg=FluxDetrendConfig(gap_split_d=0.0),
    )
    split = flux_space_detrend_result(
        time,
        flux,
        quality=quality,
        flux_err=flux_err,
        cfg=FluxDetrendConfig(gap_split_d=0.5),
    )

    assert no_split.cotrend_status == "poly_spline_failed"
    assert split.n_segments == 2
    assert split.cotrend_status == "spline;spline"

    split_det = split.det_flux - np.nanmedian(split.det_flux) + 1.0
    no_split_det = no_split.det_flux - np.nanmedian(no_split.det_flux) + 1.0
    assert _binned_amp(time, split_det) < 0.35 * _binned_amp(time, no_split_det)


def test_cotrend_does_not_extrapolate_flagged_segment_edges() -> None:
    rng = np.random.default_rng(7)
    time = np.linspace(0.0, 6.0, 1200)
    trend = 400.0 * np.sin(0.8 * time) + 50.0 * time
    flux = trend + rng.normal(0.0, 30.0, time.size)
    quality = np.zeros_like(time, dtype=int)
    quality[:80] = 1
    quality[-80:] = 1

    result = flux_space_detrend_result(
        time,
        flux,
        quality=quality,
        flux_err=np.full_like(time, 30.0),
        cfg=FluxDetrendConfig(),
    )

    good = quality == 0
    assert result.cotrend_status == "spline"
    assert np.isfinite(result.cotrend).all()
    assert np.nanmax(result.cotrend[:80]) <= np.nanmax(result.cotrend[good])
    assert np.nanmin(result.cotrend[:80]) >= np.nanmin(result.cotrend[good])
    assert np.nanmax(result.cotrend[-80:]) <= np.nanmax(result.cotrend[good])
    assert np.nanmin(result.cotrend[-80:]) >= np.nanmin(result.cotrend[good])
