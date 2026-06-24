from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

from twirl.search.injections import (
    batman_transit_model,
    box_transit_mask,
    choose_observed_epoch,
    inject_batman_transit,
    inject_box_transit,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_make_injections_script():
    path = REPO_ROOT / "scripts" / "stage3_injections" / "make_s56_lc_injection_training_set.py"
    spec = importlib.util.spec_from_file_location("make_s56_lc_injection_training_set", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_box_transit_mask_uses_periodic_phase() -> None:
    time = np.array([-0.01, 0.0, 0.01, 0.49, 0.50, 0.51, 1.0])
    mask = box_transit_mask(time, period_d=0.5, t0_d=0.0, duration_min=30.0)

    assert mask.tolist() == [True, True, True, True, True, True, True]

    narrow = box_transit_mask(time, period_d=0.5, t0_d=0.0, duration_min=10.0)
    assert narrow.tolist() == [False, True, False, False, True, False, True]


def test_inject_box_transit_subtracts_depth_times_baseline() -> None:
    time = np.array([-0.01, 0.0, 0.01, 0.25])
    flux = np.ones_like(time)

    injected, mask = inject_box_transit(
        time,
        flux,
        period_d=1.0,
        t0_d=0.0,
        duration_min=30.0,
        depth=0.2,
        baseline=1.5,
    )

    assert mask.tolist() == [True, True, True, False]
    assert np.allclose(injected[mask], 0.7)
    assert np.allclose(injected[~mask], 1.0)
    assert np.allclose(flux, 1.0)


def test_inject_batman_transit_applies_limb_darkened_model() -> None:
    time = np.linspace(-0.08, 0.08, 200)
    flux = np.ones_like(time)

    injected, mask, model = inject_batman_transit(
        time,
        flux,
        period_d=1.0,
        t0_d=0.0,
        duration_min=90.0,
        radius_rstar=0.4,
        a_over_rstar=20.0,
        impact_b=0.1,
        baseline=1.0,
        supersample_factor=3,
    )

    assert mask.any()
    assert np.nanmin(model) < 0.9
    assert np.allclose(injected, model)
    assert np.allclose(flux, 1.0)


def test_batman_transit_model_rejects_impossible_impact() -> None:
    try:
        batman_transit_model(
            np.linspace(-0.1, 0.1, 5),
            period_d=1.0,
            t0_d=0.0,
            radius_rstar=0.2,
            a_over_rstar=2.0,
            impact_b=2.0,
        )
    except ValueError as exc:
        assert "impact_b" in str(exc)
    else:
        raise AssertionError("impact_b >= a_over_rstar should be rejected")


def test_choose_observed_epoch_requires_good_overlap() -> None:
    rng = np.random.default_rng(7)
    time = np.linspace(0.0, 3.0, 200)
    quality = np.zeros_like(time, dtype=np.int32)
    quality[:100] = 1

    t0, mask = choose_observed_epoch(
        time,
        period_d=1.0,
        duration_min=60.0,
        rng=rng,
        quality=quality,
        min_in_transit=2,
    )

    assert np.isfinite(t0)
    assert np.count_nonzero(mask & (quality == 0)) >= 2


def test_choose_observed_epoch_rejects_empty_good_cadences() -> None:
    rng = np.random.default_rng(8)
    time = np.arange(10.0)
    quality = np.ones_like(time, dtype=np.int32)

    try:
        choose_observed_epoch(
            time,
            period_d=1.0,
            duration_min=60.0,
            rng=rng,
            quality=quality,
        )
    except ValueError as exc:
        assert "no observed good cadences" in str(exc)
    else:
        raise AssertionError("all-flagged cadences should not accept an injection")


def test_physical_wd_injection_families_store_radius_metadata() -> None:
    module = _load_make_injections_script()
    rng = np.random.default_rng(11)

    params = module._draw_params("wd_giant_or_bd", rng)

    assert params["signal_family"] == "wd_giant_or_bd"
    assert 1.2 <= params["radius_rearth"] <= module.DEFAULT_GRID_RADIUS_RANGE_REARTH[1]
    assert params["radius_rwd"] > 0
    assert params["impact_b"] >= 0
    assert params["a_over_rwd"] > params["impact_b"]
    assert params["duration_model"] == "wd_density_batman"
    assert params["injection_model"] == "batman_quadratic"
    assert 0 < params["depth"] <= 0.995


def test_depth_grid_cycles_evenly_over_period_depth_cells() -> None:
    module = _load_make_injections_script()
    rng = np.random.default_rng(12)

    rows = [
        module._draw_depth_grid_params(
            idx,
            rng=rng,
            period_range=(0.1, 10.0),
            depth_range=(0.01, 0.99),
            radius_range_rearth=(0.08, 11.2),
            period_bins=5,
            depth_bins=4,
            depth_spacing="linear",
        )
        for idx in range(40)
    ]
    counts = {}
    for row in rows:
        counts[row["grid_cell_id"]] = counts.get(row["grid_cell_id"], 0) + 1
        assert row["sampling_mode"] == "period_depth_grid"
        assert 0.01 <= row["target_depth"] <= 0.99
        assert row["a_over_rwd"] > row["impact_b"]

    assert len(counts) == 20
    assert set(counts.values()) == {2}


def test_period_radius_grid_cycles_and_draws_transiting_impact_parameters() -> None:
    module = _load_make_injections_script()
    rng = np.random.default_rng(13)

    rows = [
        module._draw_grid_params(
            idx,
            rng=rng,
            period_range=(0.1, 10.0),
            radius_range_rearth=(0.2, 16.8),
            period_bins=4,
            radius_bins=5,
        )
        for idx in range(40)
    ]
    counts = {}
    for row in rows:
        counts[row["grid_cell_id"]] = counts.get(row["grid_cell_id"], 0) + 1
        assert row["sampling_mode"] == "period_radius_grid"
        assert 0.2 <= row["radius_rearth"] <= 16.8
        assert row["period_d"] >= 0.1
        assert row["impact_b"] < 1.0 + row["radius_rwd"]
        assert row["a_over_rwd"] > row["impact_b"]
        assert 0.0 < row["inclination_deg"] <= 90.0
        assert row["injection_model"] == "batman_quadratic"

    assert len(counts) == 20
    assert set(counts.values()) == {2}
