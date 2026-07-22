from __future__ import annotations

import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pytest

from twirl.lightcurves.external_quality import (
    EFFECTIVE_QUALITY_POLICY,
    EXTERNAL_QUALITY_POLICY_CONTRACT,
)
from twirl.vetting.harmonic_inputs import (
    CHANNEL_CONTRACT,
    CANDIDATE_PROVENANCE_CONTRACT_VERSION,
    HARMONIC_FACTORS,
    NATIVE_DATASETS,
    RAW_PAIR_CONTRACT_VERSION,
    RAW_PAIR_CANDIDATE_PROVENANCE_ATTRS,
    NativeLightCurve,
    build_harmonic_views,
    build_native_channels,
    candidate_provenance_from_summary,
    injected_raw_uncertainty,
    pad_channel_sequences,
    read_native_light_curve,
    read_native_light_curve_from_h5,
    verify_native_candidate_binding,
    verify_raw_pair_contract,
)


def _set_quality_contract(h5: h5py.File, group: h5py.Group, quality: np.ndarray) -> None:
    n_total = len(quality)
    n_bad = int(np.count_nonzero(quality))
    h5.attrs["external_quality_policy_contract"] = EXTERNAL_QUALITY_POLICY_CONTRACT
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
    for name, digest in (
        ("cadence_reference_table_sha256", "1" * 64),
        ("cadence_reference_manifest_sha256", "2" * 64),
        ("cadence_reference_source_declaration_sha256", "3" * 64),
    ):
        h5.attrs[name] = digest
    counts = {
        "n_cad_total": n_total,
        "n_cad_internal_bad": n_bad,
        "n_cad_external_bad": 0,
        "n_cad_external_only_bad": 0,
        "n_cad_effective_bad": n_bad,
    }
    group.attrs["quality_policy_contract"] = EXTERNAL_QUALITY_POLICY_CONTRACT
    for name, value in counts.items():
        group.attrs[name] = value
        h5.attrs[f"quality_overlay_{name}"] = value


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
        _set_quality_contract(h5, group, lc.quality)

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

    with h5py.File(path, "r+") as h5:
        group = h5["targets/0000000000000001"]
        original = int(group.attrs["n_cad_effective_bad"])
        group.attrs["n_cad_effective_bad"] = original + 1
    corrupted = verify_raw_pair_contract(path)
    assert not corrupted["passed"]
    assert any("effective-bad arithmetic" in value for value in corrupted["failures"])
    with h5py.File(path, "r+") as h5:
        group = h5["targets/0000000000000001"]
        group.attrs["n_cad_effective_bad"] = original
        h5.attrs["quality_overlay_n_cad_effective_bad"] = original + 1
    corrupted = verify_raw_pair_contract(path)
    assert not corrupted["passed"]
    assert any(
        "root effective-bad arithmetic" in value
        or "disagrees with group total" in value
        for value in corrupted["failures"]
    )


def test_native_candidate_binding_requires_exact_groups_and_provenance(
    tmp_path: Path,
) -> None:
    lc = _native_lc(32)
    candidate_table = tmp_path / "candidates.csv"
    candidate_table.write_text(
        "review_id,tic,source_kind,native_input_include,period_d,t0_bjd,duration_min\n"
        "real:1,1,real_candidate,true,1.0,2459825.0,10.0\n"
    )
    table_digest = hashlib.sha256(candidate_table.read_bytes()).hexdigest()
    summary_path = tmp_path / "summary.json"
    summary_path.write_text(
        json.dumps(
            {
                "provenance_contract_version": (
                    CANDIDATE_PROVENANCE_CONTRACT_VERSION
                ),
                "candidate_table_sha256": table_digest,
                "n_candidate_rows": 1,
                "tier1_target_eligibility_sha256": "4" * 64,
                "tier1_gate_json_sha256": "5" * 64,
                "adp_peaks_sha256": "6" * 64,
                "adp_peaks_summary_sha256": "7" * 64,
                "compact_lc_sha256": "8" * 64,
                "cadence_reference_sha256": "1" * 64,
                "cadence_reference_manifest_sha256": "2" * 64,
                "bls_search_contract_version": "s56_a2v1_teacher_bls_search_v1",
                "bls_config_sha256": "9" * 64,
                "tier1_gate": {"enrichment_ready": True},
                "bls_evidence": {"status": "pass"},
            }
        )
        + "\n"
    )
    provenance = candidate_provenance_from_summary(
        candidate_table=candidate_table,
        candidate_summary=summary_path,
    )
    path = tmp_path / "native_bound.h5"
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
        _set_quality_contract(h5, group, lc.quality)
        h5.attrs["training_table_sha256"] = provenance["training_table_sha256"]
        for name in RAW_PAIR_CANDIDATE_PROVENANCE_ATTRS:
            h5.attrs[name] = provenance[name]

    verification = verify_native_candidate_binding(
        path,
        candidate_table=candidate_table,
        candidate_summary=summary_path,
        require_periodograms=False,
    )
    assert verification["passed"], verification["failures"]

    with h5py.File(path, "a") as h5:
        h5.create_group("targets/0000000000000002")
    verification = verify_native_candidate_binding(
        path,
        candidate_table=candidate_table,
        candidate_summary=summary_path,
        require_periodograms=False,
    )
    assert not verification["passed"]
    assert any("unexpected candidate groups" in value for value in verification["failures"])
