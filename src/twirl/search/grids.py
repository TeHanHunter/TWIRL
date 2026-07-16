"""Period and duration grid construction for the per-sector BLS first pass.

Grid conventions (defined by :class:`twirl.search.bls.BLSConfig` and mirrored
by the runtime YAML passed to ``sector_run --config``):

- duration grid: 3, 5, 8, 12, 20 minutes (WD transits are minute-scale; the
  TESS 200 s cadence sets a ~1-cadence floor).
- period grid: log-spaced between `p_min` (1 hour by default) and
  `p_max = max_period_fraction * baseline_d` (default 0.45 × baseline so ≥2
  transits fit). Density set by `frequency_factor` à la astropy
  `BoxLeastSquares.autoperiod`.
"""
from __future__ import annotations

import numpy as np

DURATIONS_MIN: tuple[float, ...] = (3.0, 5.0, 8.0, 12.0, 20.0)


def duration_grid_days(durations_min: tuple[float, ...] = DURATIONS_MIN) -> np.ndarray:
    """Return durations as days for use with astropy BLS."""
    return np.asarray(durations_min, dtype=np.float64) / 1440.0


def build_period_grid(
    baseline_d: float,
    p_min_d: float = 2.0 / 24.0,
    max_period_fraction: float = 0.45,
    p_max_cap_d: float = 15.0,
    n_periods: int = 100_000,
    durations_d: np.ndarray | None = None,
) -> np.ndarray:
    """Frequency-uniform period grid suitable for astropy BLS.

    `n_periods` trial periods spaced uniformly in frequency between
    `1/p_max` and `1/p_min_d`. 100k is the default — at this density a 27.8 d
    sector resolves the WD 1856 1.408 d period (Δf≈2e-4 cycle/d) easily and
    BLS evaluates in ~30–60 s per LC on a single CPU core.

    Going much denser (e.g., 1M) hits the astropy "autoperiod" regime, which
    places ~10⁷ trial periods at sector baselines — minutes per LC and tens of
    hours per sector across 19k targets. Stage 2 first pass keeps the grid
    coarse on purpose; multi-sector refinement later can tighten.

    `p_max` is capped at `max_period_fraction * baseline_d` (≥2 transits) and
    at `p_max_cap_d`.
    """
    if durations_d is None:
        durations_d = duration_grid_days()
    p_max = min(max_period_fraction * float(baseline_d), float(p_max_cap_d))
    if p_max <= p_min_d:
        return np.array([p_min_d], dtype=np.float64)
    f_min = 1.0 / p_max
    f_max = 1.0 / p_min_d
    n = max(int(n_periods), 100)
    freqs = np.linspace(f_min, f_max, n)
    periods = np.sort(1.0 / freqs)
    return periods.astype(np.float64)
