from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest

from twirl.lightcurves.a2v1_cadence_reference import (
    AUTHORITY_EXCLUSION_EXTERNAL_BIT,
    AUTHORITY_EXCLUSION_POLICY,
    AUTHORITY_EXCLUSION_POLICY_CONTRACT,
    CADENCE_REFERENCE_COLUMNS,
    authority_exclusions_sha256,
    file_sha256,
)
from twirl.lightcurves.external_quality import (
    EFFECTIVE_QUALITY_POLICY,
    EXTERNAL_QUALITY_POLICY_CONTRACT,
)
from twirl.vetting import harmonic_export
from twirl.vetting.harmonic_export import (
    RAW_SOURCE_CONTRACT_VERSION,
    align_raw_by_cadence,
    build_raw_pair_export,
    compact_adp_detector_lookup,
    discover_tglc_paths,
    merge_raw_pair_shards,
    read_candidate_table,
    shared_bls_periodogram,
    shared_bls_spectrum,
)
from twirl.vetting.harmonic_inputs import (
    CANDIDATE_PROVENANCE_CONTRACT_VERSION,
    CHANNEL_CONTRACT,
    verify_raw_pair_contract,
)


def _write_table(frame: pd.DataFrame, path: Path) -> None:
    if path.suffix == ".parquet":
        frame.to_parquet(path, compression="zstd", index=False)
    else:
        frame.to_csv(path, index=False)


def _write_s56_reference(
    root: Path,
    *,
    camera1_orbit119_cadences: np.ndarray,
    external_bad_cadences: set[int] | None = None,
) -> tuple[Path, Path]:
    external_bad_cadences = external_bad_cadences or set()
    rows: list[dict[str, int]] = []
    camera_orbit_cadences: dict[tuple[int, int], list[int]] = {}
    for camera in range(1, 5):
        for orbit in (119, 120):
            if camera == 1 and orbit == 119:
                cadences = [int(value) for value in camera1_orbit119_cadences]
            else:
                cadences = [9_000_000 + orbit * 100 + camera]
            camera_orbit_cadences[(camera, orbit)] = cadences
            for ccd in range(1, 5):
                for cadence in cadences:
                    spoc = (
                        8
                        if camera == 1
                        and ccd == 1
                        and cadence in external_bad_cadences
                        else 0
                    )
                    rows.append(
                        {
                            "sector": 56,
                            "orbitid": orbit,
                            "camera": camera,
                            "ccd": ccd,
                            "cadenceno": cadence,
                            "spoc_quality": spoc,
                            "qlp_quality": 0,
                            "external_quality": spoc,
                        }
                    )
    frame = pd.DataFrame(rows, columns=CADENCE_REFERENCE_COLUMNS).sort_values(
        ["orbitid", "camera", "ccd", "cadenceno"], kind="stable"
    )
    table = root / "cadence_reference.csv"
    frame.to_csv(table, index=False)

    sources: list[dict[str, object]] = []
    source_index = 1

    def add_source(role: str, **metadata: int) -> None:
        nonlocal source_index
        path = (root / f"authority_{source_index:03d}.txt").resolve()
        digest = f"{source_index:064x}"
        sources.append(
            {"role": role, "path": str(path), "sha256": digest, **metadata}
        )
        source_index += 1

    for orbit in (119, 120):
        for camera in range(1, 5):
            add_source("qlp_cam_quat", orbitid=orbit, camera=camera)
            for ccd in range(1, 5):
                add_source(
                    "qlp_detector_qflag",
                    orbitid=orbit,
                    camera=camera,
                    ccd=ccd,
                    n_rows=len(camera_orbit_cadences[(camera, orbit)]),
                )
    for camera in range(1, 5):
        for ccd in range(1, 5):
            add_source("spoc_flag_file", camera=camera, ccd=ccd)
    add_source("spoc_quality_table")
    add_source("spoc_quality_provenance")
    detector_names = [
        f"cam{camera}_ccd{ccd}"
        for camera in range(1, 5)
        for ccd in range(1, 5)
    ]
    authority_exclusions = {
        "contract_version": AUTHORITY_EXCLUSION_POLICY_CONTRACT,
        "policy": AUTHORITY_EXCLUSION_POLICY,
        "external_bit": AUTHORITY_EXCLUSION_EXTERNAL_BIT,
        "n_rows": 0,
        "by_detector": {
            detector: {"n_rows": 0, "rows": []}
            for detector in detector_names
        },
    }
    manifest = {
        "contract_version": "s56_a2v1_cadence_reference_v1",
        "builder_version": "a2v1_cadence_reference_builder_v3",
        "sector": 56,
        "cadence_authority": "qlp_cam_quat",
        "quality_authority": "spoc_and_qlp_quality_flags",
        "quality_composition": {
            "external_quality": "spoc_quality | (qlp_quality << 30)",
            "qlp_quality_raw_values": [0, 1],
            "qlp_quality_external_bit": 30,
        },
        "table_sha256": file_sha256(table),
        "table_columns": list(CADENCE_REFERENCE_COLUMNS),
        "n_rows": len(frame),
        "detectors": detector_names,
        "orbits": [119, 120],
        "n_rows_by_detector": {
            f"cam{int(camera)}_ccd{int(ccd)}": int(len(group))
            for (camera, ccd), group in frame.groupby(["camera", "ccd"])
        },
        "n_nonzero_spoc_quality": int(np.count_nonzero(frame["spoc_quality"])),
        "n_nonzero_qlp_quality": 0,
        "n_nonzero_external_quality": int(
            np.count_nonzero(frame["external_quality"])
        ),
        "n_spoc_authority_files_verified": 16,
        "n_qlp_qflag_files_verified": 32,
        "n_spoc_rows_excluded_by_quat": 0,
        "authority_exclusions": authority_exclusions,
        "authority_exclusions_sha256": authority_exclusions_sha256(
            authority_exclusions
        ),
        "source_file_sha256": {
            str(source["path"]): str(source["sha256"]) for source in sources
        },
        "sources": sources,
    }
    manifest_path = root / "cadence_reference.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return table, manifest_path


def test_candidate_table_reader_rejects_unknown_format(tmp_path: Path) -> None:
    path = tmp_path / "training.json"
    path.write_text("[]\n")

    with pytest.raises(ValueError, match="unsupported table format"):
        read_candidate_table(path)


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


@pytest.mark.parametrize("table_suffix", [".csv", ".parquet"])
def test_raw_source_export_uses_compact_detector_and_clears_stale_failure(
    tmp_path: Path, monkeypatch, table_suffix: str
) -> None:
    tic = 42
    training = tmp_path / f"training{table_suffix}"
    training_rows = pd.DataFrame(
        {
            "tic": [tic],
            "cam": [np.nan],
            "ccd": [np.nan],
            "morphology_include_v1": [True],
            "preserve_include_v1": [True],
            "harmonic_include_v1": [False],
        }
    )
    _write_table(training_rows, training)
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


@pytest.mark.parametrize("table_suffix", [".csv", ".parquet"])
def test_real_only_raw_pair_export_writes_full_contract(
    tmp_path: Path, table_suffix: str
) -> None:
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
        group.attrs["sector"] = 56
        group.attrs["camera"] = 1
        group.attrs["ccd"] = 1
        group.create_dataset("time", data=time - 2457000.0)
        group.create_dataset("cadenceno", data=cadence)
        group.create_dataset("orbitid", data=np.full(n, 119))
        internal_quality = np.zeros(n, dtype=np.int32)
        internal_quality[3] = 4
        group.create_dataset("quality", data=internal_quality)
        group.create_dataset("DET_FLUX_ADP_SML", data=np.ones(n))
        group.create_dataset("DET_FLUX_ADP", data=np.ones(n))

    table = tmp_path / f"training{table_suffix}"
    training_rows = pd.DataFrame(
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
    )
    _write_table(training_rows, table)
    output = tmp_path / "native.h5"
    cadence_table, cadence_manifest = _write_s56_reference(
        tmp_path,
        camera1_orbit119_cadences=cadence,
        external_bad_cadences={int(cadence[5])},
    )
    training_summary = tmp_path / "training.summary.json"
    training_summary.write_text(
        json.dumps(
            {
                "provenance_contract_version": (
                    CANDIDATE_PROVENANCE_CONTRACT_VERSION
                ),
                "candidate_table_sha256": file_sha256(table),
                "n_candidate_rows": 1,
                "tier1_target_eligibility_sha256": "4" * 64,
                "tier1_gate_json_sha256": "5" * 64,
                "adp_peaks_sha256": "6" * 64,
                "adp_peaks_summary_sha256": "7" * 64,
                "compact_lc_sha256": file_sha256(adp_source),
                "cadence_reference_sha256": file_sha256(cadence_table),
                "cadence_reference_manifest_sha256": file_sha256(cadence_manifest),
                "bls_search_contract_version": "s56_a2v1_teacher_bls_search_v1",
                "bls_config_sha256": "9" * 64,
                "tier1_gate": {"enrichment_ready": True},
                "bls_evidence": {"status": "pass"},
            }
        )
        + "\n"
    )

    summary = build_raw_pair_export(
        training_table=table,
        training_summary=training_summary,
        raw_source_h5=raw_source,
        compact_adp_h5=adp_source,
        injection_pair_h5=None,
        cadence_reference_table=cadence_table,
        cadence_reference_manifest=cadence_manifest,
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
        assert np.flatnonzero(group["quality"][:]).tolist() == [3, 5]
        assert group.attrs["n_cad_internal_bad"] == 1
        assert group.attrs["n_cad_external_bad"] == 1
        assert group.attrs["n_cad_authority_excluded"] == 0
        assert h5.attrs["external_quality_policy_contract"] == (
            EXTERNAL_QUALITY_POLICY_CONTRACT
        )
        assert h5.attrs["effective_quality_policy"] == EFFECTIVE_QUALITY_POLICY
        assert h5.attrs["cadence_reference_table_sha256"] == file_sha256(
            cadence_table
        )
        assert h5.attrs["authority_exclusion_policy_contract"] == (
            AUTHORITY_EXCLUSION_POLICY_CONTRACT
        )
        assert h5.attrs["authority_exclusion_external_bit"] == (
            AUTHORITY_EXCLUSION_EXTERNAL_BIT
        )
        assert len(h5.attrs["authority_exclusions_sha256"]) == 64
        assert h5.attrs["n_authority_exclusions"] == 0
        assert h5.attrs["compact_adp_h5_sha256"] == file_sha256(adp_source)
        assert h5.attrs["compact_lc_sha256"] == file_sha256(adp_source)
        assert h5.attrs["quality_overlay_n_cad_authority_excluded"] == 0
    assert summary["external_quality"]["counts"]["n_cad_effective_bad"] == 2
    assert summary["external_quality"]["counts"]["n_cad_authority_excluded"] == 0


def test_raw_pair_export_rejects_compact_not_bound_by_candidate_summary(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    training = tmp_path / "training.csv"
    training.write_text("review_id,tic,native_input_include\nreal:1,1,true\n")
    training_summary = tmp_path / "training.summary.json"
    training_summary.write_text("{}\n")
    raw_source = tmp_path / "raw.h5"
    compact = tmp_path / "compact.h5"
    raw_source.write_bytes(b"raw")
    compact.write_bytes(b"compact")

    monkeypatch.setattr(
        harmonic_export,
        "candidate_provenance_from_summary",
        lambda **_: {
            "training_table_sha256": file_sha256(training),
            "training_summary_sha256": file_sha256(training_summary),
            "compact_lc_sha256": "0" * 64,
        },
    )
    with pytest.raises(ValueError, match="compact ADP HDF5 SHA256"):
        build_raw_pair_export(
            training_table=training,
            training_summary=training_summary,
            raw_source_h5=raw_source,
            compact_adp_h5=compact,
            injection_pair_h5=None,
            cadence_reference_table=tmp_path / "cadence.csv",
            cadence_reference_manifest=tmp_path / "cadence.json",
            out_h5=tmp_path / "native.h5",
            repo_root=tmp_path,
        )


def test_raw_pair_export_accepts_parquet_training_table(tmp_path: Path) -> None:
    tic = 456
    n = 220
    cadence = np.arange(2000, 2000 + n)
    time = 2459825.0 + np.arange(n) * 200.0 / 86400.0
    raw_source = tmp_path / "raw.h5"
    with h5py.File(raw_source, "w") as h5:
        h5.attrs["contract_version"] = RAW_SOURCE_CONTRACT_VERSION
        group = h5.create_group(f"targets/{tic:016d}")
        for name, values in {
            "time": time,
            "cadenceno": cadence,
            "orbitid": np.full(n, 119),
            "quality": np.zeros(n),
            "raw_flux_small": np.linspace(-1.0, 5.0, n),
            "raw_flux_err_small": np.ones(n),
            "raw_flux_primary": np.linspace(0.0, 7.0, n),
            "raw_flux_err_primary": np.full(n, 2.0),
        }.items():
            group.create_dataset(name, data=values)
    adp_source = tmp_path / "adp.h5"
    with h5py.File(adp_source, "w") as h5:
        group = h5.create_group(f"targets/{tic:016d}")
        group.attrs["tic"] = tic
        group.attrs["sector"] = 56
        group.attrs["camera"] = 1
        group.attrs["ccd"] = 1
        for name, values in {
            "time": time - 2457000.0,
            "cadenceno": cadence,
            "orbitid": np.full(n, 119),
            "quality": np.zeros(n),
            "DET_FLUX_ADP_SML": np.ones(n),
            "DET_FLUX_ADP": np.ones(n),
        }.items():
            group.create_dataset(name, data=values)
    table = tmp_path / "training.parquet"
    pd.DataFrame(
        {
            "review_id": ["real:456"],
            "tic": [tic],
            "source_kind": ["real_candidate"],
            "native_input_include": [True],
        }
    ).to_parquet(table, index=False)
    output = tmp_path / "native.h5"
    cadence_table, cadence_manifest = _write_s56_reference(
        tmp_path, camera1_orbit119_cadences=cadence
    )

    summary = build_raw_pair_export(
        training_table=table,
        raw_source_h5=raw_source,
        compact_adp_h5=adp_source,
        injection_pair_h5=None,
        cadence_reference_table=cadence_table,
        cadence_reference_manifest=cadence_manifest,
        out_h5=output,
        repo_root=tmp_path,
        n_periods=32,
    )

    assert summary["n_real_targets"] == 1
    assert verify_raw_pair_contract(
        output, require_errors=True, require_periodograms=True
    )["passed"]

def test_injection_raw_pair_uses_same_authoritative_quality_overlay(
    tmp_path: Path,
) -> None:
    tic = 789
    injection_id = "inj-0001"
    n = 220
    cadence = np.arange(3000, 3000 + n)
    time_relative = 2825.0 + np.arange(n) * 200.0 / 86400.0
    time_absolute = time_relative + 2457000.0
    internal_quality = np.zeros(n, dtype=np.int32)
    internal_quality[4] = 16

    raw_source = tmp_path / "raw.h5"
    with h5py.File(raw_source, "w") as h5:
        h5.attrs["contract_version"] = RAW_SOURCE_CONTRACT_VERSION
        group = h5.create_group(f"targets/{tic:016d}")
        for name, values in {
            "time": time_absolute,
            "cadenceno": cadence,
            "orbitid": np.full(n, 119),
            "quality": internal_quality,
            "raw_flux_small": np.full(n, 100.0),
            "raw_flux_err_small": np.full(n, 2.0),
            "raw_flux_primary": np.full(n, 200.0),
            "raw_flux_err_primary": np.full(n, 3.0),
        }.items():
            group.create_dataset(name, data=values)
    adp_source = tmp_path / "adp.h5"
    with h5py.File(adp_source, "w") as h5:
        h5.create_group("targets")

    canonical_path = tmp_path / "canonical.h5"
    with h5py.File(canonical_path, "w") as h5:
        group = h5.create_group(f"injections/{injection_id}")
        for name, values in {
            "RAW_FLUX_Small_injected": np.full(n, 99.0),
            "RAW_FLUX_Small_original": np.full(n, 100.0),
            "RAW_FLUX_Primary_injected": np.full(n, 198.0),
            "RAW_FLUX_Primary_original": np.full(n, 200.0),
        }.items():
            group.create_dataset(name, data=values)
    pair_path = tmp_path / "pair.h5"
    det_flux = np.ones(n)
    det_flux[np.arange(n) % 50 < 4] -= 0.1
    with h5py.File(pair_path, "w") as h5:
        h5.attrs["source_injection_h5"] = str(canonical_path)
        group = h5.create_group(f"injections/{injection_id}")
        group.attrs["tic"] = tic
        group.attrs["sector"] = 56
        group.attrs["camera"] = 1
        group.attrs["ccd"] = 1
        group.attrs["cadence_s"] = 200.0
        group.attrs["injection_baseline_Small"] = 100.0
        group.attrs["injection_baseline_Primary"] = 200.0
        group.attrs["source_injection_h5"] = str(canonical_path)
        for name, values in {
            "time": time_relative,
            "cadenceno": cadence,
            "orbitid": np.full(n, 119),
            "quality": internal_quality,
            "transit_model": np.ones(n),
            "DET_FLUX_ADP_SML_injected": det_flux,
            "DET_FLUX_ADP_injected": det_flux,
            "DET_FLUX_ADP_SML_original": det_flux,
            "DET_FLUX_ADP_original": det_flux,
        }.items():
            group.create_dataset(name, data=values)

    training = tmp_path / "training.csv"
    pd.DataFrame(
        {
            "review_id": [f"injection:{injection_id}"],
            "tic": [tic],
            "source_kind": ["injection_recovery"],
            "is_injected_row": [True],
            "injection_id": [injection_id],
            "native_input_include": [True],
        }
    ).to_csv(training, index=False)
    cadence_table, cadence_manifest = _write_s56_reference(
        tmp_path,
        camera1_orbit119_cadences=cadence,
        external_bad_cadences={int(cadence[7])},
    )
    output = tmp_path / "native_injection.h5"

    summary = build_raw_pair_export(
        training_table=training,
        raw_source_h5=raw_source,
        compact_adp_h5=adp_source,
        injection_pair_h5=pair_path,
        cadence_reference_table=cadence_table,
        cadence_reference_manifest=cadence_manifest,
        out_h5=output,
        repo_root=tmp_path,
        n_periods=32,
    )

    assert summary["n_injections"] == 1
    with h5py.File(output, "r") as h5:
        group = h5[f"injections/{injection_id}"]
        assert np.flatnonzero(group["quality"][:]).tolist() == [4, 7]
        assert group.attrs["n_cad_external_only_bad"] == 1
        assert group.attrs["n_cad_effective_bad"] == 2
    verification = verify_raw_pair_contract(
        output, require_errors=True, require_periodograms=True
    )
    assert verification["passed"], verification["failures"]


def test_shard_merge_rejects_no_data_loss(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    shards = []
    for shard_index, tic in enumerate((1, 2)):
        path = tmp_path / f"shard_{shard_index}.h5"
        with h5py.File(path, "w") as h5:
            h5.attrs["contract_version"] = harmonic_export.RAW_PAIR_CONTRACT_VERSION
            h5.attrs["time_system"] = "BJD"
            h5.attrs["external_quality_policy_contract"] = (
                EXTERNAL_QUALITY_POLICY_CONTRACT
            )
            h5.attrs["effective_quality_policy"] = EFFECTIVE_QUALITY_POLICY
            h5.attrs["cadence_reference_contract_version"] = (
                "s56_a2v1_cadence_reference_v1"
            )
            h5.attrs["cadence_reference_cadence_authority"] = "qlp_cam_quat"
            h5.attrs["cadence_reference_quality_authority"] = (
                "spoc_and_qlp_quality_flags"
            )
            h5.attrs["cadence_reference_table"] = "/authority/reference.csv"
            h5.attrs["cadence_reference_manifest"] = "/authority/reference.json"
            h5.attrs["cadence_reference_table_sha256"] = "1" * 64
            h5.attrs["cadence_reference_manifest_sha256"] = "2" * 64
            h5.attrs["cadence_reference_source_declaration_sha256"] = "3" * 64
            h5.attrs["authority_exclusion_policy_contract"] = (
                AUTHORITY_EXCLUSION_POLICY_CONTRACT
            )
            h5.attrs["authority_exclusion_external_bit"] = (
                AUTHORITY_EXCLUSION_EXTERNAL_BIT
            )
            h5.attrs["authority_exclusions_sha256"] = "4" * 64
            h5.attrs["n_authority_exclusions"] = 0
            h5.attrs["training_table_sha256"] = "a" * 64
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
            group.attrs["quality_policy_contract"] = (
                EXTERNAL_QUALITY_POLICY_CONTRACT
            )
            quality_counts = {
                "n_cad_total": 4,
                "n_cad_internal_bad": 4,
                "n_cad_external_bad": 0,
                "n_cad_external_only_bad": 0,
                "n_cad_authority_excluded": 0,
                "n_cad_effective_bad": 4,
            }
            for name, value in quality_counts.items():
                group.attrs[name] = value
                h5.attrs[f"quality_overlay_{name}"] = value
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
    assert summary["external_quality_counts"]["n_cad_total"] == 8
    assert verify_raw_pair_contract(
        output, require_errors=True, require_periodograms=True
    )["passed"]

    real_sha256 = harmonic_export._file_sha256
    calls = 0

    def changing_sha256(path: Path) -> str:
        nonlocal calls
        calls += 1
        if calls == 3:
            return "f" * 64
        return real_sha256(path)

    monkeypatch.setattr(harmonic_export, "_file_sha256", changing_sha256)
    with pytest.raises(RuntimeError, match="shards changed during merge"):
        merge_raw_pair_shards(
            shard_paths=shards,
            out_h5=tmp_path / "changed_during_merge.h5",
        )
    monkeypatch.setattr(harmonic_export, "_file_sha256", real_sha256)

    merged_training = tmp_path / "merged_candidates.csv"
    merged_training.write_text("review_id,tic\nreal:1,1\nreal:2,2\n")
    for shard_index, shard in enumerate(shards):
        local_training = tmp_path / f"candidates_{shard_index:02d}.csv"
        local_training.write_text(
            f"review_id,tic\nreal:{shard_index + 1},{shard_index + 1}\n"
        )
        local_pair = tmp_path / f"injections_{shard_index:02d}.h5"
        local_pair.write_bytes(f"pair-{shard_index}".encode())
        with h5py.File(shard, "r+") as h5:
            h5.attrs["training_table"] = str(local_training.resolve())
            h5.attrs["training_table_sha256"] = file_sha256(local_training)
            h5.attrs["injection_pair_h5"] = str(local_pair.resolve())
            h5.attrs["injection_pair_h5_sha256"] = file_sha256(local_pair)

    aggregate_output = tmp_path / "merged_aggregate.h5"
    aggregate = merge_raw_pair_shards(
        shard_paths=shards,
        out_h5=aggregate_output,
        merged_training_table=merged_training,
    )
    assert len(aggregate["shard_local_provenance"]) == 2
    with h5py.File(aggregate_output, "r") as h5:
        assert h5.attrs["training_table"] == str(merged_training.resolve())
        assert h5.attrs["training_table_sha256"] == file_sha256(merged_training)
        assert h5.attrs["injection_pair_h5"] == ""

    with pytest.raises(ValueError, match="native-input provenance mismatch"):
        merge_raw_pair_shards(
            shard_paths=shards,
            out_h5=tmp_path / "strict_local_mismatch.h5",
        )

    with h5py.File(shards[1], "r+") as h5:
        h5.attrs["cadence_reference_table_sha256"] = "f" * 64
    with pytest.raises(ValueError, match="external-quality provenance mismatch"):
        merge_raw_pair_shards(
            shard_paths=shards,
            out_h5=tmp_path / "mismatched.h5",
            merged_training_table=merged_training,
        )
