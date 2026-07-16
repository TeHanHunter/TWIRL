from __future__ import annotations

import h5py
import numpy as np
import pandas as pd

from twirl.lightcurves.a2v1_qa import (
    WD1856_PERIOD_D,
    WD1856_T0_BJD,
    WD1856_TIC,
    audit_compact_adp,
    audit_bls_coverage,
    audit_wd1856_bls,
    compare_raw_extractions,
    evaluate_photometric_gates,
    stratified_target_sample,
)
from twirl.vetting.adp_only import ADP_ONLY_CONTRACT_VERSION


def _write_raw_h5(path, *, offset: float = 0.0) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cadence = np.arange(100, dtype=np.int64)
    signal = 1000.0 + 20.0 * np.sin(np.linspace(0, 6 * np.pi, len(cadence))) + offset
    with h5py.File(path, "w") as h5:
        h5.attrs["TIC ID"] = 123
        h5.attrs["Sector"] = 56
        h5.attrs["Orbit"] = 119
        h5.attrs["Camera"] = 1
        h5.attrs["CCD"] = 2
        h5.attrs["TessMag"] = 17.0
        h5.attrs["BJDoffset"] = 2457000
        lc = h5.create_group("LightCurve")
        lc.create_dataset("BJD", data=2800.0 + cadence / 432.0)
        lc.create_dataset("Cadence", data=cadence)
        lc.create_dataset("QualityFlag", data=np.zeros(len(cadence), dtype=np.int32))
        lc.create_dataset("X", data=np.zeros(len(cadence)))
        lc.create_dataset("Y", data=np.zeros(len(cadence)))
        background = lc.create_group("Background")
        background.create_dataset("Value", data=np.zeros(len(cadence)))
        background.create_dataset("Error", data=np.ones(len(cadence)))
        photometry = lc.create_group("AperturePhotometry")
        for aperture, scale in (("Small", 1.0), ("Primary", 1.2), ("Large", 1.5)):
            group = photometry.create_group(f"{aperture}Aperture")
            flux = signal * scale
            group.create_dataset("RawFlux", data=flux)
            group.create_dataset("RawFluxError", data=np.full(len(cadence), 3.0))
            group.create_dataset("RawMagnitude", data=-2.5 * np.log10(flux))
            group.create_dataset("RawMagnitudeError", data=np.full(len(cadence), 0.01))
            group.create_dataset("X", data=np.zeros(len(cadence)))
            group.create_dataset("Y", data=np.zeros(len(cadence)))


def test_stratified_sample_is_deterministic_across_magnitude() -> None:
    catalog = pd.DataFrame(
        {
            "tic": np.arange(100),
            "sector": 56,
            "camera": 1,
            "ccd": 2,
            "tmag": np.linspace(14, 21, 100),
        }
    )
    first = stratified_target_sample(catalog, sample_size=25, seed=56)
    second = stratified_target_sample(catalog, sample_size=25, seed=56)
    assert first["tic"].tolist() == second["tic"].tolist()
    assert len(first) == 25
    assert first["tmag"].min() < 15.5
    assert first["tmag"].max() > 19.5


def test_compact_adp_audit_preserves_finite_and_negative_cadences(tmp_path) -> None:
    compact = tmp_path / "compact.h5"
    with h5py.File(compact, "w") as h5:
        targets = h5.create_group("targets")
        group = targets.create_group(f"{123:016d}")
        group.create_dataset("time", data=np.arange(100) / 432.0)
        group.create_dataset("cadenceno", data=np.arange(100))
        group.create_dataset("quality", data=np.zeros(100, dtype=np.int32))
        group.create_dataset("orbitid", data=np.full(100, 119, dtype=np.int16))
        small = 1.0 + 0.01 * np.sin(np.linspace(0, 4 * np.pi, 100))
        small[5] = -0.2
        group.create_dataset("DET_FLUX_ADP_SML", data=small)
        group.create_dataset("DET_FLUX_ADP", data=small * 1.01)
    sample = pd.DataFrame(
        {"tic": [123], "sector": [56], "camera": [1], "ccd": [2], "tmag": [17.0]}
    )
    audited = audit_compact_adp(compact, sample)
    assert audited.loc[0, "finite_q0_fraction_small"] == 1.0
    assert audited.loc[0, "negative_q0_fraction_small"] == 0.01
    assert audited.loc[0, "aperture_correlation"] > 0.99


def test_raw_extraction_comparison_matches_native_flux(tmp_path) -> None:
    a2v1_root = tmp_path / "a2v1"
    reference_root = tmp_path / "reference"
    relative = "orbit-119/ffi/cam1/ccd2/LC/123.h5"
    _write_raw_h5(a2v1_root / relative, offset=0.0)
    _write_raw_h5(reference_root / relative, offset=2.0)
    sample = pd.DataFrame({"tic": [123], "tmag": [17.0], "camera": [1], "ccd": [2]})
    compared = compare_raw_extractions(
        sample,
        a2v1_root=a2v1_root,
        reference_root=reference_root,
        orbits=[119],
        workers=2,
    )
    assert len(compared) == 2
    assert compared["status"].eq("ok").all()
    assert compared["separate_file_identity"].all()
    assert compared["normalized_flux_correlation"].min() > 0.999


def test_wd1856_and_sector_qa_gate_use_current_adp_bls() -> None:
    peaks = pd.DataFrame(
        {
            "tic": [WD1856_TIC, WD1856_TIC],
            "aperture": ["DET_FLUX_ADP_SML", "DET_FLUX_ADP"],
            "peak_rank": [1, 1],
            "period_d": [WD1856_PERIOD_D * 1.0001, WD1856_PERIOD_D * 0.9999],
            "t0_bjd": [WD1856_T0_BJD, WD1856_T0_BJD + 5 * WD1856_PERIOD_D],
            "sde": [60.0, 45.0],
        }
    )
    wd1856 = audit_wd1856_bls(peaks)
    assert wd1856["passed"]

    n = 100
    targets = pd.DataFrame(
        {
            "finite_q0_fraction_small": np.full(n, 0.99),
            "finite_q0_fraction_primary": np.full(n, 0.98),
            "quality0_fraction": np.full(n, 0.90),
            "mad_ppm_small": np.full(n, 20_000.0),
            "mad_ppm_primary": np.full(n, 25_000.0),
            "aperture_correlation": np.full(n, 0.80),
            "primary_to_small_mad_ratio": np.full(n, 1.25),
            "aperture_difference_mad_ppm": np.full(n, 10_000.0),
        }
    )
    raw = pd.DataFrame(
        {
            "status": ["ok"] * 100,
            "separate_file_identity": [True] * 100,
            "normalized_flux_correlation": np.full(100, 0.95),
            "normalized_difference_mad_ppm": np.full(100, 5_000.0),
        }
    )
    schema = {
        "sector": 56,
        "ok": True,
        "ok_h5": True,
        "ok_fits": True,
        "h5": {"n_present_h5": 200, "n_missing_h5_non_edge": 0, "n_zero_byte_h5": 0},
        "fits": {
            "n_found_unique_tics": 100,
            "n_missing_fits_non_edge_tics": 0,
            "n_bad_checked_fits": 0,
        },
    }
    evaluated = evaluate_photometric_gates(
        sector=56,
        schema=schema,
        n_compact_targets=100,
        bls_coverage={"passed": True, "n_expected_pairs": 200, "n_observed_pairs": 200},
        target_metrics=targets,
        raw_comparison=raw,
        wd1856=wd1856,
        reference_benchmark=None,
    )
    assert evaluated["passed"]
    assert set(evaluated["gates"]) == {
        "schema_and_completeness",
        "bls_target_aperture_coverage",
        "sampled_adp_photometry",
        "aperture_consistency",
        "prior_extraction_tree_comparison",
        "benchmark",
    }


def test_bls_coverage_rejects_truncated_table_and_accepts_accounted_failures() -> None:
    catalog = pd.DataFrame({"tic": [1, 2], "sector": [56, 56]})
    common = {
        "source_product_tag": "A2v1",
        "adp_only_contract_version": ADP_ONLY_CONTRACT_VERSION,
        "bls_n_periods": 50_000,
        "bls_n_peaks": 10,
        "bls_p_min_d": 0.12,
        "bls_p_max_cap_d": 15.0,
        "bls_max_period_fraction": 0.45,
        "bls_sigma_clip": 5.0,
        "bls_orbit_edge_trim_d": 0.0,
    }
    truncated = pd.DataFrame(
        [
            {"tic": 1, "aperture": aperture, "status": "ok", **common}
            for aperture in ("DET_FLUX_ADP_SML", "DET_FLUX_ADP")
        ]
    )
    failed = audit_bls_coverage(truncated, catalog)
    assert failed["passed"] is False
    assert failed["n_missing_pairs"] == 2

    complete = pd.concat(
        [
            truncated,
            pd.DataFrame(
                [
                    {
                        "tic": 2,
                        "aperture": aperture,
                        "status": "too_few_cadences",
                        **common,
                    }
                    for aperture in ("DET_FLUX_ADP_SML", "DET_FLUX_ADP")
                ]
            ),
        ],
        ignore_index=True,
    )
    passed = audit_bls_coverage(complete, catalog)
    assert passed["passed"] is True
    assert passed["n_missing_pairs"] == 0
    assert passed["n_failed_pairs"] == 2
    assert passed["failed_pair_fraction"] == 0.5
