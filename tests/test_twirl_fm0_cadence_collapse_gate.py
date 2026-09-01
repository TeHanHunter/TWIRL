from __future__ import annotations

import numpy as np
import pytest

from twirl.models.fm0.cadence_collapse_gate import (
    _evaluate_cadence_collapse_math,
    summarize_cadence_difference_geometry,
)


def _continuous_geometry(batch_size: int, length: int) -> tuple[np.ndarray, ...]:
    valid = np.ones((batch_size, length), dtype=bool)
    boundary = np.zeros((batch_size, length), dtype=bool)
    boundary[:, 0] = True
    delta_time = np.ones((batch_size, length), dtype=np.float64)
    delta_time[:, 0] = 0.0
    return valid, boundary, delta_time


def test_source_constant_cadence_codes_fail_closed() -> None:
    rng = np.random.default_rng(560067)
    source_codes = rng.normal(size=(32, 1, 32))
    tokens = np.broadcast_to(source_codes, (32, 64, 32)).copy()
    valid, boundary, delta_time = _continuous_geometry(32, 64)

    result = _evaluate_cadence_collapse_math(
        tokens,
        tokens,
        valid,
        boundary,
        delta_time,
    )

    assert result["step0"]["g_delta"] == 0.0
    assert result["step2000"]["difference_effective_rank"] == 0.0
    assert result["step2000"]["difference_energy_positive"] is False
    assert result["criteria_satisfied"] is False
    assert "passed" not in result
    assert result["authoritative_artifact_gate"] is False


def test_deterministic_position_codes_do_not_fake_data_dependent_geometry() -> None:
    rng = np.random.default_rng(560068)
    source_codes = rng.normal(size=(32, 1, 32))
    position_codes = rng.normal(size=(1, 64, 32))
    tokens = source_codes + position_codes
    valid, boundary, delta_time = _continuous_geometry(32, 64)

    summary = summarize_cadence_difference_geometry(
        tokens,
        valid,
        boundary,
        delta_time,
    )
    result = _evaluate_cadence_collapse_math(
        tokens,
        tokens,
        valid,
        boundary,
        delta_time,
    )

    assert summary["centered_endpoint_energy"] > 0.0
    assert summary["g_delta"] == 0.0
    assert summary["difference_effective_rank"] == 0.0
    assert result["criteria_satisfied"] is False


def test_data_dependent_cadence_variation_passes_frozen_gate() -> None:
    rng = np.random.default_rng(560069)
    step0 = rng.normal(size=(32, 64, 32))
    step2000 = 1.75 * step0
    valid, boundary, delta_time = _continuous_geometry(32, 64)

    result = _evaluate_cadence_collapse_math(
        step0,
        step2000,
        valid,
        boundary,
        delta_time,
    )

    assert result["step0"]["diagnostic_valid"] is True
    assert result["step2000"]["diagnostic_valid"] is True
    assert result["step2000"]["g_delta"] == pytest.approx(result["step0"]["g_delta"])
    assert result["step2000"]["difference_effective_rank"] >= 26.0
    assert result["thresholds"] == {
        "g_delta_step2000_over_step0_minimum": 0.5,
        "cadence_energy_step2000_over_step0_minimum": 0.5,
        "difference_effective_rank_minimum": 26.0,
        "maximum_contiguous_delta_cadences": 1.5,
    }
    assert result["g_delta_passed"] is True
    assert result["cadence_energy_passed"] is True
    assert result["difference_effective_rank_passed"] is True
    assert result["criteria_satisfied"] is True
    assert result["cadence_representations_averaged_or_merged"] is False


def test_uniform_near_zero_scaling_fails_energy_retention_gate() -> None:
    rng = np.random.default_rng(560071)
    step0 = rng.normal(size=(32, 64, 32))
    step2000 = 1.0e-7 * step0
    valid, boundary, delta_time = _continuous_geometry(32, 64)

    result = _evaluate_cadence_collapse_math(
        step0,
        step2000,
        valid,
        boundary,
        delta_time,
    )

    assert result["g_delta_passed"] is True
    assert result["difference_effective_rank_passed"] is True
    assert result["cadence_energy_step2000_over_step0"] == pytest.approx(1.0e-14)
    assert result["cadence_energy_passed"] is False
    assert result["criteria_satisfied"] is False


def test_invalid_boundary_and_gap_edges_are_excluded() -> None:
    rng = np.random.default_rng(560070)
    tokens = rng.normal(size=(4, 7, 8))
    valid, boundary, delta_time = _continuous_geometry(4, 7)
    valid[:, 2] = False
    delta_time[:, 4] = 4.0
    boundary[:, 5] = True

    reference = summarize_cadence_difference_geometry(
        tokens,
        valid,
        boundary,
        delta_time,
    )
    corrupted_excluded = tokens.copy()
    corrupted_excluded[:, 2:5, :] = np.nan
    changed = summarize_cadence_difference_geometry(
        corrupted_excluded,
        valid,
        boundary,
        delta_time,
    )

    assert reference["eligible_edges_by_position"] == [4, 0, 0, 0, 0, 4]
    assert reference["eligible_contiguous_edges"] == 8
    assert reference["statistical_edge_positions"] == 2
    assert reference["statistical_degrees_of_freedom"] == 6
    for name in (
        "centered_difference_energy",
        "centered_endpoint_energy",
        "g_delta",
        "difference_effective_rank",
    ):
        assert changed[name] == pytest.approx(reference[name])


def test_eligible_nonfinite_endpoint_and_bad_shapes_fail_closed() -> None:
    tokens = np.zeros((3, 4, 2), dtype=np.float64)
    valid, boundary, delta_time = _continuous_geometry(3, 4)
    tokens[0, 1, 0] = np.nan
    with pytest.raises(ValueError, match="non-finite endpoint token"):
        summarize_cadence_difference_geometry(
            tokens,
            valid,
            boundary,
            delta_time,
        )

    with pytest.raises(ValueError, match="valid must have shape"):
        summarize_cadence_difference_geometry(
            np.zeros((3, 4, 2)),
            np.ones((3, 3), dtype=bool),
            boundary,
            delta_time,
        )
