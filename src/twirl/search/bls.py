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

from twirl.io.hlsp import (
    APERTURES,
    BJDREFI,
    HLSPLightCurve,
    quality_mask_from_arrays,
)
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


@dataclass(frozen=True)
class BLSPreparedInputs:
    """Exact pre-search state shared by BLS and upstream eligibility checks."""

    status: str
    initial_mask: np.ndarray
    time: np.ndarray
    normalized_flux: np.ndarray
    orbitid: np.ndarray
    n_total: int
    n_cad_quality: int
    n_cad_kept: int
    n_cad_edge_trimmed: int
    n_cad_sigma_clipped: int
    flux_median: float
    baseline_d: float

    @property
    def ready(self) -> bool:
        return self.status == "ok"


def prepare_bls_inputs_from_arrays(
    *,
    time: np.ndarray,
    flux: np.ndarray,
    quality: np.ndarray,
    orbitid: np.ndarray,
    cfg: BLSConfig,
) -> BLSPreparedInputs:
    """Apply the exact BLS pre-search validity and cleaning semantics.

    This deliberately stops before constructing the full configured period
    grid or evaluating :class:`~astropy.timeseries.BoxLeastSquares`.  A small
    grid probe uses the same grid builder to determine whether the cleaned
    baseline supports at least two trial periods.
    """

    time_array = np.asarray(time)
    flux_array = np.asarray(flux)
    quality_array = np.asarray(quality)
    orbitid_array = np.asarray(orbitid)
    if orbitid_array.shape != time_array.shape:
        raise ValueError("time and orbitid arrays must have identical shapes")

    initial_mask = quality_mask_from_arrays(
        time_array,
        flux_array,
        quality_array,
    )
    n_total = int(time_array.size)
    n_cad_quality = int(initial_mask.sum())
    selected_time = time_array[initial_mask].astype(np.float64)
    selected_flux = flux_array[initial_mask].astype(np.float64)
    selected_orbitid = orbitid_array[initial_mask]

    def _result(
        status: str,
        *,
        normalized_flux: np.ndarray | None = None,
        n_cad_kept: int | None = None,
        n_cad_edge_trimmed: int = 0,
        n_cad_sigma_clipped: int = 0,
        flux_median: float = np.nan,
        baseline_d: float = 0.0,
    ) -> BLSPreparedInputs:
        normalized = (
            np.asarray(normalized_flux, dtype=np.float64)
            if normalized_flux is not None
            else np.asarray([], dtype=np.float64)
        )
        return BLSPreparedInputs(
            status=status,
            initial_mask=initial_mask,
            time=selected_time,
            normalized_flux=normalized,
            orbitid=selected_orbitid,
            n_total=n_total,
            n_cad_quality=n_cad_quality,
            n_cad_kept=(
                n_cad_quality if n_cad_kept is None else int(n_cad_kept)
            ),
            n_cad_edge_trimmed=int(n_cad_edge_trimmed),
            n_cad_sigma_clipped=int(n_cad_sigma_clipped),
            flux_median=float(flux_median),
            baseline_d=float(baseline_d),
        )

    if n_cad_quality < cfg.min_cadences:
        return _result("too_few_cadences")

    n_cad_kept = n_cad_quality
    n_cad_edge_trimmed = 0
    if cfg.orbit_edge_trim_d and cfg.orbit_edge_trim_d > 0:
        edge_keep = np.ones_like(selected_time, dtype=bool)
        for oid in np.unique(selected_orbitid):
            selected = selected_orbitid == oid
            if selected.sum() < 2:
                continue
            orbit_time = selected_time[selected]
            lower = orbit_time.min() + cfg.orbit_edge_trim_d
            upper = orbit_time.max() - cfg.orbit_edge_trim_d
            edge_keep[selected] = (orbit_time >= lower) & (orbit_time <= upper)
        if edge_keep.sum() >= cfg.min_cadences:
            n_cad_edge_trimmed = int((~edge_keep).sum())
            selected_time = selected_time[edge_keep]
            selected_flux = selected_flux[edge_keep]
            selected_orbitid = selected_orbitid[edge_keep]
            n_cad_kept = int(edge_keep.sum())

    median = float(np.nanmedian(selected_flux))
    if not np.isfinite(median) or median == 0.0:
        return _result(
            "all_nan",
            n_cad_kept=n_cad_kept,
            n_cad_edge_trimmed=n_cad_edge_trimmed,
            flux_median=median,
        )
    normalized_flux = selected_flux / median

    n_cad_sigma_clipped = 0
    if cfg.sigma_clip and cfg.sigma_clip > 0:
        mad = float(np.nanmedian(np.abs(normalized_flux - 1.0)))
        if mad > 0 and np.isfinite(mad):
            sigma = 1.4826 * mad
            keep = (
                normalized_flux - 1.0
            ) <= cfg.sigma_clip * sigma
            if keep.sum() >= cfg.min_cadences:
                n_cad_sigma_clipped = int((~keep).sum())
                selected_time = selected_time[keep]
                normalized_flux = normalized_flux[keep]
                selected_orbitid = selected_orbitid[keep]
                n_cad_kept = int(keep.sum())

    baseline_d = float(selected_time.max() - selected_time.min())
    durations_d = duration_grid_days(cfg.durations_min)
    grid_probe = build_period_grid(
        baseline_d=baseline_d,
        p_min_d=cfg.p_min_d,
        max_period_fraction=cfg.max_period_fraction,
        p_max_cap_d=cfg.p_max_cap_d,
        n_periods=2,
        durations_d=durations_d,
    )
    status = "ok" if grid_probe.size >= 2 else "degenerate_grid"
    return _result(
        status,
        normalized_flux=normalized_flux,
        n_cad_kept=n_cad_kept,
        n_cad_edge_trimmed=n_cad_edge_trimmed,
        n_cad_sigma_clipped=n_cad_sigma_clipped,
        flux_median=median,
        baseline_d=baseline_d,
    )


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

    prepared = prepare_bls_inputs_from_arrays(
        time=lc.time,
        flux=lc.flux[aperture],
        quality=lc.quality,
        orbitid=lc.orbitid,
        cfg=cfg,
    )
    n_total = prepared.n_total
    if not prepared.ready:
        res = _empty_result(
            lc,
            aperture,
            prepared.status,
            n_total,
            prepared.n_cad_kept,
            cfg,
        )
        return _maybe_with_periodogram(res, return_periodogram)

    t = prepared.time
    y = prepared.normalized_flux
    orbitid_kept = prepared.orbitid
    n_kept = prepared.n_cad_kept
    n_cad_quality = prepared.n_cad_quality
    n_cad_edge_trimmed = prepared.n_cad_edge_trimmed
    n_cad_sigma_clipped = prepared.n_cad_sigma_clipped
    baseline_d = prepared.baseline_d
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


__all__ = [
    "APERTURES",
    "BLSConfig",
    "BLSPreparedInputs",
    "prepare_bls_inputs_from_arrays",
    "run_bls_on_lc",
    "save_periodogram",
]
