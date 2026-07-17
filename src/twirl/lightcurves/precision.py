from __future__ import annotations

import math

import numpy as np


def quality_good_mask(quality: np.ndarray, drop_bits: int) -> np.ndarray:
    """Return cadences whose bitwise quality flags do not intersect ``drop_bits``."""
    return (np.asarray(quality, dtype=np.int64) & int(drop_bits)) == 0


def point_to_point_mad_precision(
    flux: np.ndarray,
    quality: np.ndarray | None = None,
    *,
    bin_minutes: float = 30.0,
    cadence_s: float = 200.0,
    drop_bits: int = 0xFFFFFFFF,
    min_points: int = 30,
    mad_scale: float = 1.4826,
) -> float:
    """Estimate binned fractional precision using point-to-point MAD.

    The per-cadence estimate is ``1.4826 * median(abs(diff(flux))) / sqrt(2)``.
    It is then scaled to the requested bin size assuming white noise.
    """
    flux = np.asarray(flux, dtype=float)
    if quality is None:
        mask = np.isfinite(flux)
    else:
        mask = quality_good_mask(quality, drop_bits) & np.isfinite(flux)
    clean = flux[mask]
    if clean.size < min_points:
        return np.nan

    diff = np.abs(np.diff(clean))
    if diff.size == 0:
        return np.nan

    sigma_per_cad = mad_scale * np.nanmedian(diff) / math.sqrt(2.0)
    n_per_bin = max(1, int(round(float(bin_minutes) * 60.0 / float(cadence_s))))
    return float(sigma_per_cad / math.sqrt(n_per_bin))
