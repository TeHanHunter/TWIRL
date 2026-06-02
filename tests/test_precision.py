from __future__ import annotations

import math

import numpy as np

from twirl.lightcurves.precision import point_to_point_mad_precision, quality_good_mask


def test_point_to_point_mad_precision_scales_to_30_min_bin() -> None:
    flux = np.array([1.0, 1.02, 1.0, 1.02, 1.0, 1.02, 1.0, 1.02, 1.0, 1.02] * 4)
    quality = np.zeros(flux.size, dtype=int)

    got = point_to_point_mad_precision(
        flux,
        quality,
        bin_minutes=30.0,
        cadence_s=200.0,
        min_points=10,
    )

    expected = 1.4826 * 0.02 / math.sqrt(2.0) / math.sqrt(9.0)
    assert np.isclose(got, expected)


def test_point_to_point_mad_precision_drops_requested_quality_bits() -> None:
    flux = np.array([1.0, 1.01, 5.0, 1.02, 1.0, 1.01, 1.0, 1.01, 1.0, 1.01, 1.0])
    quality = np.zeros(flux.size, dtype=int)
    quality[2] = 2

    keep = quality_good_mask(quality, 0xFFFFFFFF)
    assert not keep[2]

    got = point_to_point_mad_precision(
        flux,
        quality,
        bin_minutes=200.0 / 60.0,
        cadence_s=200.0,
        min_points=5,
    )
    expected = 1.4826 * 0.01 / math.sqrt(2.0)
    assert np.isclose(got, expected)
