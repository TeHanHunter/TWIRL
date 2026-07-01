from __future__ import annotations

import importlib.util
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_peak_builder():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "build_injection_peak_training_table.py"
    spec = importlib.util.spec_from_file_location("build_injection_peak_training_table", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_label_requires_transit_window_overlap_for_exact_period() -> None:
    module = _load_peak_builder()

    label = module.label_peak_against_injection(
        period_d=1.0,
        t0_bjd=2459825.10,
        duration_min=8.0,
        truth_period_d=1.0,
        truth_t0_bjd=2459825.0,
        truth_duration_min=8.0,
    )

    assert label["period_rel_err"] == 0.0
    assert not label["transit_window_match"]
    assert not label["exact_ephemeris_match"]
    assert not label["is_injected_signal_peak"]
    assert label["match_kind"] == "mismatch"


def test_label_accepts_overlapping_exact_period_peak() -> None:
    module = _load_peak_builder()

    label = module.label_peak_against_injection(
        period_d=1.0,
        t0_bjd=2459825.0 + 2.0 / 1440.0,
        duration_min=8.0,
        truth_period_d=1.0,
        truth_t0_bjd=2459825.0,
        truth_duration_min=8.0,
    )

    assert label["transit_window_overlap_fraction"] >= 0.5
    assert label["exact_ephemeris_match"]
    assert label["is_injected_signal_peak"]
    assert label["match_kind"] == "exact"


def test_label_accepts_overlapping_harmonic_peak() -> None:
    module = _load_peak_builder()

    label = module.label_peak_against_injection(
        period_d=0.5,
        t0_bjd=2459825.0,
        duration_min=8.0,
        truth_period_d=1.0,
        truth_t0_bjd=2459825.0,
        truth_duration_min=8.0,
    )

    assert label["nearest_harmonic_factor"] == 0.5
    assert label["transit_window_match"]
    assert label["harmonic_ephemeris_match"]
    assert label["is_injected_signal_peak"]
    assert label["match_kind"] == "harmonic"


def test_label_rejects_harmonic_without_window_overlap() -> None:
    module = _load_peak_builder()

    label = module.label_peak_against_injection(
        period_d=0.5,
        t0_bjd=2459825.25,
        duration_min=8.0,
        truth_period_d=1.0,
        truth_t0_bjd=2459825.0,
        truth_duration_min=8.0,
    )

    assert label["nearest_harmonic_factor"] == 0.5
    assert not label["transit_window_match"]
    assert not label["harmonic_ephemeris_match"]
    assert not label["is_injected_signal_peak"]
