from __future__ import annotations

import numpy as np
import pandas as pd

from twirl.vetting.recovery50_cnn import (
    CnnTrainConfig,
    _grouped_splits,
    event_view,
    folded_view,
)


def test_folded_view_uses_small_primary_and_difference_channels() -> None:
    rng = np.random.default_rng(56)
    period = 1.0
    t0 = 2459000.25
    duration_min = 30.0
    time = 2459000.0 + np.linspace(0, 4, 1600)
    phase_d = ((time - t0 + 0.5 * period) % period) - 0.5 * period
    in_transit = np.abs(phase_d) < duration_min / 1440.0 / 2
    small = 1.0 + rng.normal(0, 0.002, size=time.size)
    primary = 1.0 + rng.normal(0, 0.002, size=time.size)
    small[in_transit] -= 0.10
    primary[in_transit] -= 0.04
    quality = np.zeros(time.size, dtype=np.int32)

    values, mask, counts, x = folded_view(
        time_bjd=time,
        flux_by_aperture=[small, primary],
        quality=quality,
        period_d=period,
        t0_bjd=t0,
        duration_min=duration_min,
        n_points=129,
        window_durations=4,
    )

    assert values.shape == (3, 129)
    assert mask.shape == values.shape
    assert counts.shape == values.shape
    center = np.argmin(np.abs(x))
    assert values[0, center] < -0.05
    assert values[1, center] < -0.015
    assert values[2, center] > 0.03
    assert mask[:, center].all()


def test_event_view_keeps_event_axis_and_masks_missing_slots() -> None:
    period = 1.0
    t0 = 2459000.25
    duration_min = 20.0
    time = 2459000.0 + np.linspace(0, 3.2, 900)
    phase_d = ((time - t0 + 0.5 * period) % period) - 0.5 * period
    in_transit = np.abs(phase_d) < duration_min / 1440.0 / 2
    small = np.ones_like(time)
    primary = np.ones_like(time)
    small[in_transit] -= 0.08
    primary[in_transit] -= 0.03
    quality = np.zeros(time.size, dtype=np.int32)

    values, mask, counts, centers = event_view(
        time_bjd=time,
        flux_by_aperture=[small, primary],
        quality=quality,
        period_d=period,
        t0_bjd=t0,
        duration_min=duration_min,
        n_points=65,
        window_durations=4,
        max_events=8,
        min_event_points=1,
    )

    assert values.shape == (3, 8, 65)
    assert mask.shape == values.shape
    assert counts.shape == values.shape
    assert np.isfinite(centers).sum() >= 3
    assert not mask[:, np.isfinite(centers).sum():, :].any()


def test_grouped_splits_keep_same_tic_together() -> None:
    rows = pd.DataFrame(
        {
            "tic": [1, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            "main_teacher_target": [
                "planet_like",
                "planet_like",
                "planet_like",
                "planet_like",
                "instrumental_or_systematic",
                "instrumental_or_systematic",
                "instrumental_or_systematic",
                "instrumental_or_systematic",
                "stellar_variability",
                "stellar_variability",
            ],
        }
    )
    split = _grouped_splits(
        rows,
        classes=["planet_like", "instrumental_or_systematic", "stellar_variability"],
        cfg=CnnTrainConfig(seed=11, validation_fraction=0.25, test_fraction=0.25, require_cuda=False),
    )
    for _, idx in rows.groupby("tic").groups.items():
        assert split.iloc[list(idx)].nunique() == 1
    assert set(split.unique()).issubset({"train", "validation", "test"})
