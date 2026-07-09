"""astropy BoxLeastSquares wrapper for the per-sector BLS first pass.

Runs BLS on a single HLSP light curve, extracts the top-N SDE peaks, and
returns a `BLSResult`. Detrending is not done here; active S56 searches use
the already-detrended ADP small/primary pair and median-normalize before BLS.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.timeseries import BoxLeastSquares

from twirl.io.hlsp import APERTURES, BJDREFI, HLSPLightCurve, quality_mask
from twirl.search.candidates import (
    BLSPeak,
    BLSResult,
    compute_sde,
    walk_peaks,
)
from twirl.search.grids import build_period_grid, duration_grid_days


@dataclass
class BLSConfig:
    apertures: tuple[str, ...] = ("DET_FLUX_ADP_SML", "DET_FLUX_ADP")
    # p_min_d raised from 2 h (v1) to 0.12 d (~2.9 h) to cut the empirical
    # alias wall pile-up (26% of best-peaks landed in [0.083, 0.10] d in v1)
    # and the cadence-sampling floor where even a central 2 R_jup transit
    # has < 4 cadences at 200 s.
    p_min_d: float = 0.12
    max_period_fraction: float = 0.45
    p_max_cap_d: float = 15.0
    # Denser duration grid: resolves the 5-13 min WD regime more finely than
    # the v1 [3, 5, 8, 12, 20] grid, and extends to 30 min so PCEBs do not
    # peg at the 21-min ceiling (an artifact we saw in v1 — most top-of-list
    # FAs had duration_min == 21).
    durations_min: tuple[float, ...] = (
        3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 13.0, 16.0, 20.0, 30.0
    )
    n_periods: int = 200_000
    # n_peaks dropped from 20 to 10: peaks 11-20 contributed no surviving
    # candidates after the heuristic vetter in v1.
    n_peaks: int = 10
    period_mask_frac: float = 0.005
    # Optional search-branch control: limit how many peaks can come from each
    # period range so a crowded systematic family does not consume the whole
    # top-N list. Defaults preserve the historical global-SDE walk.
    period_bin_edges: tuple[float, ...] = ()
    max_peaks_per_period_bin: int = 0
    min_cadences: int = 200
    sigma_clip: float = 5.0  # Reject |f - median| > sigma_clip * 1.4826 * MAD before BLS.
                             # Catches scattered-light spikes that QUALITY=0 lets through.
    orbit_edge_trim_d: float = 0.0  # If > 0, drop cadences within this many days of
                                    # each orbit's first/last cadence. Standard QLP
                                    # practice for noisy CCDs (scattered-light wings).


def _empty_result(
    lc: HLSPLightCurve, aperture: str, status: str, n_total: int, n_kept: int,
    cfg: BLSConfig,
) -> BLSResult:
    return BLSResult(
        tic=lc.tic, sector=lc.sector, cam=lc.cam, ccd=lc.ccd, tmag=lc.tmag,
        aperture=aperture, n_cad_total=n_total, n_cad_kept=n_kept,
        n_cad_quality=n_kept,
        n_cad_edge_trimmed=0,
        n_cad_sigma_clipped=0,
        quality_dropout_frac=(0.0 if n_total == 0 else (n_total - n_kept) / n_total),
        dropout_frac=(0.0 if n_total == 0 else (n_total - n_kept) / n_total),
        n_orbits=0, baseline_d=0.0, status=status,
        hlsp_path=str(lc.path), peaks=[],
    )


def _maybe_with_periodogram(
    result: BLSResult,
    return_periodogram: bool,
    spectrum: dict | None = None,
) -> BLSResult | tuple[BLSResult, dict | None]:
    if return_periodogram:
        return result, spectrum
    return result


def run_bls_on_lc(lc: HLSPLightCurve, cfg: BLSConfig | None = None,
                  aperture: str = "DET_FLUX_ADP_SML",
                  return_periodogram: bool = False) -> BLSResult | tuple[BLSResult, dict | None]:
    """Run BLS on one aperture of one light curve.

    Returns a `BLSResult`. If `return_periodogram=True`, returns a tuple
    `(result, {"period": ..., "power": ..., "sde": ..., "depth": ..., ...})`
    so the caller can stash the spectrum for diagnostics.
    """
    cfg = cfg or BLSConfig()
    if aperture not in lc.flux:
        res = _empty_result(lc, aperture, "missing_aperture", len(lc.time), 0, cfg)
        return _maybe_with_periodogram(res, return_periodogram)

    t_all = lc.time
    f_all = lc.flux[aperture]
    n_total = int(t_all.size)
    mask = quality_mask(lc, aperture)
    n_cad_quality = int(mask.sum())
    if n_cad_quality < cfg.min_cadences:
        res = _empty_result(lc, aperture, "too_few_cadences", n_total, n_cad_quality, cfg)
        return _maybe_with_periodogram(res, return_periodogram)

    t = t_all[mask].astype(np.float64)
    f = f_all[mask].astype(np.float64)
    orbitid_kept = lc.orbitid[mask]
    n_kept = n_cad_quality
    n_cad_edge_trimmed = 0
    n_cad_sigma_clipped = 0

    # Optional orbit-edge trim: drop cadences within trim_d of each orbit's
    # first/last good cadence. Targets the scattered-light wings around
    # spacecraft thermal recovery and momentum dumps.
    if cfg.orbit_edge_trim_d and cfg.orbit_edge_trim_d > 0:
        edge_keep = np.ones_like(t, dtype=bool)
        for oid in np.unique(orbitid_kept):
            sel = orbitid_kept == oid
            if sel.sum() < 2:
                continue
            t_orbit = t[sel]
            t_lo = t_orbit.min() + cfg.orbit_edge_trim_d
            t_hi = t_orbit.max() - cfg.orbit_edge_trim_d
            edge_keep[sel] = (t_orbit >= t_lo) & (t_orbit <= t_hi)
        if edge_keep.sum() >= cfg.min_cadences:
            n_cad_edge_trimmed = int((~edge_keep).sum())
            t = t[edge_keep]
            f = f[edge_keep]
            orbitid_kept = orbitid_kept[edge_keep]
            n_kept = int(edge_keep.sum())

    med = float(np.nanmedian(f))
    if not np.isfinite(med) or med == 0.0:
        res = _empty_result(lc, aperture, "all_nan", n_total, n_kept, cfg)
        return _maybe_with_periodogram(res, return_periodogram)
    y = f / med

    # Upper-tail-only sigma clip to suppress scattered-light spikes (always
    # positive excursions) without touching transit dips (always negative).
    # This asymmetry matters: a 57%-deep transit at MAD-RMS=0.10 is a 5.7σ
    # negative excursion that a symmetric clip would erase.
    if cfg.sigma_clip and cfg.sigma_clip > 0:
        mad = float(np.nanmedian(np.abs(y - 1.0)))
        if mad > 0 and np.isfinite(mad):
            sigma = 1.4826 * mad
            keep = (y - 1.0) <= cfg.sigma_clip * sigma
            if keep.sum() >= cfg.min_cadences:
                n_cad_sigma_clipped = int((~keep).sum())
                t = t[keep]
                y = y[keep]
                orbitid_kept = orbitid_kept[keep]
                n_kept = int(keep.sum())

    baseline_d = float(t.max() - t.min())
    n_orbits = int(np.unique(orbitid_kept).size) if orbitid_kept.size else 0

    durations_d = duration_grid_days(cfg.durations_min)
    periods = build_period_grid(
        baseline_d=baseline_d,
        p_min_d=cfg.p_min_d,
        max_period_fraction=cfg.max_period_fraction,
        p_max_cap_d=cfg.p_max_cap_d,
        n_periods=cfg.n_periods,
        durations_d=durations_d,
    )
    if periods.size < 2:
        res = _empty_result(lc, aperture, "degenerate_grid", n_total, n_kept, cfg)
        return _maybe_with_periodogram(res, return_periodogram)

    try:
        bls = BoxLeastSquares(t, y)
        pg = bls.power(periods, durations_d, oversample=1)
    except Exception:
        res = _empty_result(lc, aperture, "bls_fail", n_total, n_kept, cfg)
        return _maybe_with_periodogram(res, return_periodogram)

    power = np.asarray(pg.power, dtype=np.float64)
    if not np.any(np.isfinite(power)):
        res = _empty_result(lc, aperture, "nan_power", n_total, n_kept, cfg)
        return _maybe_with_periodogram(res, return_periodogram)

    sde = compute_sde(power)
    extra = {
        "log_power": power,
        "duration_d": np.asarray(pg.duration, dtype=np.float64),
        "depth": np.asarray(pg.depth, dtype=np.float64),
        "depth_snr": np.asarray(pg.depth_snr, dtype=np.float64),
        "t0_rel": np.asarray(pg.transit_time, dtype=np.float64),
    }
    raw_peaks = walk_peaks(
        period=np.asarray(pg.period, dtype=np.float64),
        sde=sde,
        extra=extra,
        n_peaks=cfg.n_peaks,
        period_mask_frac=cfg.period_mask_frac,
        period_bin_edges=cfg.period_bin_edges,
        max_peaks_per_period_bin=cfg.max_peaks_per_period_bin,
    )

    peaks = [
        BLSPeak(
            peak_rank=int(p["peak_rank"]),
            period_d=float(p["period_d"]),
            t0_bjd=float(p["t0_rel"]) + float(BJDREFI),
            duration_min=float(p["duration_d"]) * 1440.0,
            depth=float(p["depth"]),
            depth_snr=float(p["depth_snr"]),
            sde=float(p["sde"]),
            log_power=float(p["log_power"]),
        )
        for p in raw_peaks
    ]

    res = BLSResult(
        tic=lc.tic, sector=lc.sector, cam=lc.cam, ccd=lc.ccd, tmag=lc.tmag,
        aperture=aperture,
        n_cad_total=n_total, n_cad_kept=n_kept,
        n_cad_quality=n_cad_quality,
        n_cad_edge_trimmed=n_cad_edge_trimmed,
        n_cad_sigma_clipped=n_cad_sigma_clipped,
        dropout_frac=(n_total - n_kept) / n_total,
        quality_dropout_frac=(n_total - n_cad_quality) / n_total,
        n_orbits=n_orbits, baseline_d=baseline_d,
        status="ok", hlsp_path=str(lc.path), peaks=peaks,
    )

    if return_periodogram:
        spectrum = {
            "period": np.asarray(pg.period, dtype=np.float32),
            "power": power.astype(np.float32),
            "sde": sde.astype(np.float32),
            "depth": np.asarray(pg.depth, dtype=np.float32),
            "duration": np.asarray(pg.duration, dtype=np.float32),
            "t0": np.asarray(pg.transit_time, dtype=np.float64),
        }
        return res, spectrum
    return res


def save_periodogram(spectrum: dict, out_path: Path) -> None:
    """Write a per-target periodogram dict to .npz for diagnostics."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out_path, **spectrum)


__all__ = ["BLSConfig", "run_bls_on_lc", "save_periodogram", "APERTURES"]
