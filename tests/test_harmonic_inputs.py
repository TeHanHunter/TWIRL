from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pytest

from twirl.vetting.harmonic_inputs import (
    CHANNEL_CONTRACT,
    HARMONIC_FACTORS,
    NATIVE_DATASETS,
    RAW_PAIR_CONTRACT_VERSION,
    NativeLightCurve,
    build_harmonic_views,
    build_native_channels,
    injected_raw_uncertainty,
    pad_channel_sequences,
    read_native_light_curve,
    read_native_light_curve_from_h5,
    verify_raw_pair_contract,
)


def _native_lc(n: int = 300) -> NativeLightCurve:
    time = 2459825.0 + np.arange(n) * (200.0 / 86400.0)
    period = 0.25
    phase = ((time - 2459825.02 + 0.5 * period) % period) / period - 0.5
    transit = np.abs(phase) < 0.02
    raw_small = np.linspace(-20.0, 40.0, n)
    raw_primary = raw_small + 5.0
    det_small = np.ones(n)
    det_primary = np.ones(n)
    det_small[transit] -= 0.15
    det_primary[transit] -= 0.07
    return NativeLightCurve(
        time=time,
        cadenceno=np.arange(n),
        orbitid=np.where(np.arange(n) < n // 2, 119, 120),
        quality=np.where(np.arange(n) % 53 == 0, 1, 0),
        raw_flux_small=raw_small,
        raw_flux_err_small=np.full(n, 5.0),
        raw_flux_primary=raw_primary,
        raw_flux_err_primary=np.full(n, 7.0),
        det_flux_adp_sml=det_small,
        det_flux_adp=det_primary,
        attrs={"tic": 1},
    )


def test_native_channels_keep_every_cadence_and_negative_raw_flux() -> None:
    lc = _native_lc()
    channels = build_native_channels(lc)

    assert channels.small_values.shape == (10, len(lc.time))
    assert channels.supplemental_values.shape == (5, len(lc.time))
    assert channels.small_mask.shape == channels.small_values.shape
    assert channels.small_mask[0].all()
    assert np.isfinite(channels.small_values[0]).all()
    assert not np.allclose(channels.small_values[0], 0.0)
    assert channels.small_values[1].min() < -0.10
    assert channels.small_values[3, 0] == 0.0
    assert channels.small_values[3, -1] == 1.0
    assert channels.small_values[5, 0] == 0.0
    assert channels.small_values[5, -1] == 1.0
    assert set(np.unique(channels.small_values[7])) == {0.0, 1.0}
    assert set(np.unique(channels.small_values[9])) == {0.0, 1.0}


def test_seven_harmonic_views_are_complete_and_local_windows_are_unbinned() -> None:
    lc = _native_lc(420)
    views = build_harmonic_views(
        lc,
        period_d=0.25,
        t0_bjd=2459825.02,
        duration_min=15.0,
    )

    assert views.factors == HARMONIC_FACTORS
    assert len(views.full_values) == 7
    assert all(value.shape == (7, len(lc.time)) for value in views.full_values)
    for values in views.full_values:
        phase = values[5]
        assert np.all(np.diff(phase[np.isfinite(phase)]) >= 0)
    assert all(0 < value.shape[1] < len(lc.time) for value in views.primary_values)
    assert all(0 < value.shape[1] < len(lc.time) for value in views.secondary_values)


def test_padding_retains_exact_native_lengths() -> None:
    one = np.arange(15, dtype=np.float32).reshape(3, 5)
    two = np.arange(24, dtype=np.float32).reshape(3, 8)

    values, mask, lengths = pad_channel_sequences([one, two])

    assert values.shape == (2, 3, 8)
    assert mask.shape == values.shape
    assert lengths.tolist() == [5, 8]
    assert np.array_equal(values[0, :, :5], one)
    assert not mask[0, :, 5:].any()


def test_injected_uncertainty_preserves_floor_and_scales_source_poisson() -> None:
    error = np.full(4, 10.0)
    model = np.asarray([1.0, 0.5, 0.1, 0.0])
    out = injected_raw_uncertainty(
        error,
        model,
        source_flux_rate=100.0,
        cadence_s=2.0,
    )

    assert np.isclose(out[0], 10.0)
    assert np.all(np.diff(out) < 0)
    assert np.isclose(out[-1], np.sqrt(50.0))
    with pytest.raises(ValueError, match="source_flux_rate"):
        injected_raw_uncertainty(error, model, source_flux_rate=np.nan, cadence_s=2.0)


def test_hdf5_contract_reader_and_verifier(tmp_path: Path) -> None:
    lc = _native_lc(32)
    path = tmp_path / "native.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["contract_version"] = RAW_PAIR_CONTRACT_VERSION
        h5.attrs["time_system"] = "BJD"
        for name, channels in CHANNEL_CONTRACT.items():
            h5.attrs[name] = json.dumps(channels)
        group = h5.create_group("targets/0000000000000001")
        payload = {
            "time": lc.time,
            "cadenceno": lc.cadenceno,
            "orbitid": lc.orbitid,
            "quality": lc.quality,
            "raw_flux_small": lc.raw_flux_small,
            "raw_flux_err_small": lc.raw_flux_err_small,
            "raw_flux_primary": lc.raw_flux_primary,
            "raw_flux_err_primary": lc.raw_flux_err_primary,
            "det_flux_adp_sml": lc.det_flux_adp_sml,
            "det_flux_adp": lc.det_flux_adp,
        }
        for name in NATIVE_DATASETS:
            group.create_dataset(name, data=payload[name])

    verification = verify_raw_pair_contract(path)
    loaded = read_native_light_curve(path, group_path="targets/0000000000000001")
    with h5py.File(path, "r") as h5:
        loaded_from_handle = read_native_light_curve_from_h5(
            h5,
            group_path="targets/0000000000000001",
        )

    assert verification["passed"], verification["failures"]
    assert verification["counts"] == {"targets": 1, "injections": 0}
    assert np.array_equal(loaded.cadenceno, lc.cadenceno)
    assert np.array_equal(loaded.raw_flux_small, lc.raw_flux_small)
    assert np.array_equal(loaded_from_handle.raw_flux_small, lc.raw_flux_small)
