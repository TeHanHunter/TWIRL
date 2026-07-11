from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest

from twirl.vetting import harmonic_export
from twirl.vetting.harmonic_export import (
    RAW_SOURCE_CONTRACT_VERSION,
    align_raw_by_cadence,
    build_raw_pair_export,
    compact_adp_detector_lookup,
    discover_tglc_paths,
    merge_raw_pair_shards,
    shared_bls_periodogram,
    shared_bls_spectrum,
)
from twirl.vetting.harmonic_inputs import verify_raw_pair_contract
from twirl.vetting.harmonic_inputs import CHANNEL_CONTRACT


def test_raw_alignment_is_by_cadence_and_preserves_negative_flux() -> None:
    raw = {
        "time": np.asarray([2459825.2, 2459825.0, 2459825.1]),
        "cadenceno": np.asarray([12, 10, 11]),
        "raw_flux_small": np.asarray([2.0, -5.0, 1.0]),
    }

    aligned = align_raw_by_cadence(
        raw,
        cadenceno=np.asarray([10, 11, 12]),
        time=np.asarray([2459825.0, 2459825.1, 2459825.2]),
    )

    assert aligned["raw_flux_small"].tolist() == [-5.0, 1.0, 2.0]


def test_raw_alignment_accepts_observed_one_point_five_second_orbit_offset() -> None:
    time = np.asarray([2459825.0, 2459825.1])
    raw = {
        "time": time + 1.5 / 86400.0,
        "cadenceno": np.asarray([10, 11]),
        "raw_flux_small": np.asarray([1.0, 2.0]),
    }

    aligned = align_raw_by_cadence(raw, cadenceno=np.asarray([10, 11]), time=time)

    assert np.allclose(aligned["_time_delta_s"], 1.5, atol=5.0e-5)
    with pytest.raises(ValueError, match="timestamps differ"):
        align_raw_by_cadence(
            raw,
            cadenceno=np.asarray([10, 11]),
            time=time,
            max_time_error_s=1.0,
        )


def test_raw_path_discovery_can_recover_missing_detector_metadata(tmp_path: Path) -> None:
    expected = tmp_path / "orbit-119/ffi/cam3/ccd2/LC/42.h5"
    expected.parent.mkdir(parents=True)
    expected.touch()

    found = discover_tglc_paths(
        raw_root=tmp_path,
        tic=42,
        camera=None,
        ccd=None,
        orbits=(119, 120),
    )

    assert found == [expected]


def test_compact_adp_detector_lookup_uses_model_facing_product(tmp_path: Path) -> None:
    path = tmp_path / "compact_adp.h5"
    with h5py.File(path, "w") as h5:
        group = h5.create_group("targets/0000000000000042")
        group.attrs["camera"] = 3
        group.attrs["ccd"] = 2

    lookup = compact_adp_detector_lookup(path, tics=[42, 99])

    assert lookup == {42: (3, 2)}


def test_raw_source_export_uses_compact_detector_and_clears_stale_failure(
    tmp_path: Path, monkeypatch
) -> None:
    tic = 42
    training = tmp_path / "training.csv"
    pd.DataFrame(
        {
            "tic": [tic],
            "cam": [np.nan],
            "ccd": [np.nan],
            "morphology_include_v1": [True],
            "preserve_include_v1": [True],
            "harmonic_include_v1": [False],
        }
    ).to_csv(training, index=False)
    compact = tmp_path / "compact.h5"
    with h5py.File(compact, "w") as h5:
        group = h5.create_group(f"targets/{tic:016d}")
        group.attrs["camera"] = 3
        group.attrs["ccd"] = 2
    raw_root = tmp_path / "raw"
    paths = []
    for orbit in (119, 120):
        path = raw_root / f"orbit-{orbit}/ffi/cam3/ccd2/LC/{tic}.h5"
        path.parent.mkdir(parents=True)
        path.touch()
        paths.append(path)
    n = 8
    payload = {
        "time": 2459825.0 + np.arange(n) * 200.0 / 86400.0,
        "cadenceno": np.arange(n),
        "orbitid": np.repeat([119, 120], n // 2),
        "quality": np.zeros(n),
        "raw_flux_small": np.ones(n),
        "raw_flux_err_small": np.ones(n),
        "raw_flux_primary": np.ones(n),
        "raw_flux_err_primary": np.ones(n),
    }
    observed_paths: list[Path] = []

    def fake_merge(found: list[Path]) -> dict[str, np.ndarray]:
        observed_paths.extend(found)
        return payload

    monkeypatch.setattr(harmonic_export, "merge_tglc_raw_paths", fake_merge)
    output = tmp_path / "raw_sources.h5"
    output.with_suffix(".failures.csv").write_text("tic,error\n42,stale\n")

    summary = harmonic_export.export_tglc_raw_sources(
        training_table=training,
        raw_root=raw_root,
        out_h5=output,
        compact_adp_h5=compact,
    )

    assert summary["n_compact_adp_detectors"] == 1
    assert observed_paths == paths
    assert not output.with_suffix(".failures.csv").exists()
    with h5py.File(output, "r") as h5:
        group = h5[f"targets/{tic:016d}"]
        assert (group.attrs["camera"], group.attrs["ccd"]) == (3, 2)
        assert group.attrs["detector_source"] == "compact_adp_attrs"


def test_shared_bls_periodogram_is_finite_for_transit_signal() -> None:
    time = 2459825.0 + np.arange(800) * 200.0 / 86400.0
    period = 0.5
    phase = ((time - 2459825.1 + 0.5 * period) % period) / period - 0.5
    flux = np.ones(len(time))
    flux[np.abs(phase) < 0.015] -= 0.15
    grid = np.geomspace(0.12, 1.5, 256)

    sde = shared_bls_periodogram(
        time=time,
        quality=np.zeros(len(time), dtype=np.int32),
        flux=flux,
        periods=grid,
    )

    assert sde.shape == grid.shape
    assert np.isfinite(sde).all()
    assert grid[int(np.nanargmax(sde))] < 1.05
    spectrum = shared_bls_spectrum(
        time=time,
        quality=np.zeros(len(time), dtype=np.int32),
        flux=flux,
        periods=grid,
    )
    assert set(spectrum) == {"power", "sde"}
    assert np.isfinite(spectrum["power"]).all()


def test_real_only_raw_pair_export_writes_full_contract(tmp_path: Path) -> None:
    tic = 123
    n = 240
    cadence = np.arange(1000, 1000 + n)
    time = 2459825.0 + np.arange(n) * 200.0 / 86400.0
    raw_source = tmp_path / "raw.h5"
    with h5py.File(raw_source, "w") as h5:
        h5.attrs["contract_version"] = RAW_SOURCE_CONTRACT_VERSION
        group = h5.create_group(f"targets/{tic:016d}")
        payload = {
            "time": time,
            "cadenceno": cadence,
            "orbitid": np.full(n, 119),
            "quality": np.zeros(n),
            "raw_flux_small": np.linspace(-2.0, 5.0, n),
            "raw_flux_err_small": np.full(n, 1.0),
            "raw_flux_primary": np.linspace(0.0, 8.0, n),
            "raw_flux_err_primary": np.full(n, 2.0),
        }
        for name, values in payload.items():
            group.create_dataset(name, data=values)

    adp_source = tmp_path / "adp.h5"
    with h5py.File(adp_source, "w") as h5:
        group = h5.create_group(f"targets/{tic:016d}")
        group.attrs["tic"] = tic
        group.create_dataset("time", data=time - 2457000.0)
        group.create_dataset("cadenceno", data=cadence)
        group.create_dataset("orbitid", data=np.full(n, 119))
        group.create_dataset("quality", data=np.zeros(n))
        group.create_dataset("DET_FLUX_ADP_SML", data=np.ones(n))
        group.create_dataset("DET_FLUX_ADP", data=np.ones(n))

    table = tmp_path / "training.csv"
    pd.DataFrame(
        [
            {
                "review_id": "real:123",
                "tic": tic,
                "source_kind": "real_candidate",
                "morphology_include_v1": True,
                "preserve_include_v1": True,
                "harmonic_include_v1": False,
            }
        ]
    ).to_csv(table, index=False)
    output = tmp_path / "native.h5"

    summary = build_raw_pair_export(
        training_table=table,
        raw_source_h5=raw_source,
        compact_adp_h5=adp_source,
        injection_pair_h5=None,
        out_h5=output,
        repo_root=tmp_path,
        n_periods=64,
    )
    verification = verify_raw_pair_contract(
        output,
        require_errors=True,
        require_periodograms=True,
    )

    assert summary["n_real_targets"] == 1
    assert verification["passed"], verification["failures"]
    with h5py.File(output, "r") as h5:
        group = h5[f"targets/{tic:016d}"]
        assert np.nanmedian(group["time"][:]) > 1.0e5
        assert group["raw_flux_small"][:].min() < 0
        assert len(group["bls_power_small"]) == 64
        assert len(group["bls_sde_small"]) == 64


def test_shard_merge_rejects_no_data_loss(tmp_path: Path) -> None:
    shards = []
    for shard_index, tic in enumerate((1, 2)):
        path = tmp_path / f"shard_{shard_index}.h5"
        with h5py.File(path, "w") as h5:
            h5.attrs["contract_version"] = "s56_adp_raw_pair_v1"
            h5.attrs["time_system"] = "BJD"
            for name, channels in CHANNEL_CONTRACT.items():
                h5.attrs[name] = json.dumps(channels)
            h5.create_group("injections")
            group = h5.create_group(f"targets/{tic:016d}")
            for name in (
                "time",
                "cadenceno",
                "orbitid",
                "quality",
                "raw_flux_small",
                "raw_flux_err_small",
                "raw_flux_primary",
                "raw_flux_err_primary",
                "det_flux_adp_sml",
                "det_flux_adp",
            ):
                values = 2459825.0 + np.arange(4) if name == "time" else np.ones(4)
                group.create_dataset(name, data=values)
            for name in (
                "bls_log_period_grid",
                "bls_power_small",
                "bls_sde_small",
                "bls_power_primary",
                "bls_sde_primary",
            ):
                group.create_dataset(name, data=np.ones(3))
        shards.append(path)

    output = tmp_path / "merged.h5"
    summary = merge_raw_pair_shards(shard_paths=shards, out_h5=output)

    assert summary["counts"] == {"targets": 2, "injections": 0}
    assert verify_raw_pair_contract(
        output, require_errors=True, require_periodograms=True
    )["passed"]
