from __future__ import annotations

import importlib.util
from pathlib import Path

import h5py
import numpy as np
import pytest


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


def _write_injection_group(path: Path, *, contract: str) -> None:
    with h5py.File(path, "w") as h5:
        h5.attrs["contract_version"] = contract
        group = h5.create_group("injections/inj_000")
        group.attrs["contract_version"] = contract
        group.attrs["injection_id"] = "inj_000"
        group.attrs["tic"] = 123
        group.attrs["sector"] = 56
        group.attrs["camera"] = 1
        group.attrs["ccd"] = 1
        group.attrs["tessmag"] = 18.0
        group.attrs["apertures"] = '["DET_FLUX_ADP_SML"]'
        group.create_dataset("time", data=np.arange(6, dtype=float))
        group.create_dataset("orbitid", data=np.ones(6, dtype=np.int16))
        group.create_dataset("quality", data=np.array([0, 0, 4, 0, 0, 0]))
        group.create_dataset(
            "DET_FLUX_ADP_SML_injected",
            data=np.ones(6, dtype=float),
        )


def test_v2_injection_peak_input_uses_effective_quality(tmp_path: Path) -> None:
    module = _load_peak_builder()
    path = tmp_path / "v2.h5"
    contract = "s56_a2v1_fresh_injection_pair_v2"
    _write_injection_group(path, contract=contract)
    with h5py.File(path, "r+") as h5:
        group = h5["injections/inj_000"]
        group.attrs[
            "epoch_quality_policy_contract"
        ] = module.EPOCH_QUALITY_POLICY_CONTRACT
        group.create_dataset(
            "external_quality",
            data=np.array([0, 2048, 0, 0, 0, 0], dtype=np.int64),
        )
        group.create_dataset(
            "effective_quality",
            data=np.array([0, 1, 1, 0, 0, 0], dtype=np.int32),
        )

    light_curve = module._injection_lc_from_group(
        path,
        "/injections/inj_000",
        "DET_FLUX_ADP_SML",
    )
    np.testing.assert_array_equal(light_curve.quality, [0, 1, 1, 0, 0, 0])


def test_v2_injection_peak_input_rejects_stale_quality_overlay(
    tmp_path: Path,
) -> None:
    module = _load_peak_builder()
    path = tmp_path / "bad_v2.h5"
    contract = "s56_a2v1_fresh_injection_pair_v2"
    _write_injection_group(path, contract=contract)
    with h5py.File(path, "r+") as h5:
        group = h5["injections/inj_000"]
        group.attrs[
            "epoch_quality_policy_contract"
        ] = module.EPOCH_QUALITY_POLICY_CONTRACT
        group.create_dataset(
            "external_quality",
            data=np.array([0, 2048, 0, 0, 0, 0], dtype=np.int64),
        )
        group.create_dataset(
            "effective_quality",
            data=np.array([0, 0, 1, 0, 0, 0], dtype=np.int32),
        )

    with pytest.raises(ValueError, match="invalid v2 effective-quality"):
        module._injection_lc_from_group(
            path,
            "/injections/inj_000",
            "DET_FLUX_ADP_SML",
        )
