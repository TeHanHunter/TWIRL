from __future__ import annotations

import numpy as np

from twirl.search.candidates import walk_peaks


def _extra(n: int) -> dict[str, np.ndarray]:
    return {
        "duration_d": np.full(n, 0.01),
        "depth": np.full(n, 0.1),
        "depth_snr": np.full(n, 10.0),
        "t0_rel": np.zeros(n),
        "log_power": np.arange(n, dtype=float),
    }


def test_walk_peaks_global_sde_order_without_period_quota() -> None:
    period = np.array([0.2, 0.25, 0.3, 2.0, 3.0], dtype=float)
    sde = np.array([30.0, 29.0, 28.0, 20.0, 19.0], dtype=float)

    peaks = walk_peaks(
        period=period,
        sde=sde,
        extra=_extra(len(period)),
        n_peaks=3,
        period_mask_frac=0.001,
        harmonics=(),
    )

    assert [round(peak["period_d"], 2) for peak in peaks] == [0.2, 0.25, 0.3]


def test_walk_peaks_period_quota_preserves_lower_sde_period_ranges() -> None:
    period = np.array([0.2, 0.25, 0.3, 2.0, 3.0], dtype=float)
    sde = np.array([30.0, 29.0, 28.0, 20.0, 19.0], dtype=float)

    peaks = walk_peaks(
        period=period,
        sde=sde,
        extra=_extra(len(period)),
        n_peaks=3,
        period_mask_frac=0.001,
        harmonics=(),
        period_bin_edges=(0.1, 1.0, 4.0),
        max_peaks_per_period_bin=1,
    )

    assert [round(peak["period_d"], 2) for peak in peaks] == [0.2, 2.0]
