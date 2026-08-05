from __future__ import annotations

from pathlib import Path

import numpy as np

from twirl.io.hlsp import HLSPLightCurve
from twirl.search.a2v1_bls_contract import (
    approved_a2v1_teacher_bls_config,
    approved_a2v1_teacher_bls_runtime_config,
)
from twirl.search.bls import (
    BLSConfig,
    prepare_bls_inputs_from_arrays,
    run_bls_on_lc,
)
from twirl.search.candidates import result_to_rows


def test_bls_defaults_use_active_adp_pair() -> None:
    assert BLSConfig().apertures == ("DET_FLUX_ADP_SML", "DET_FLUX_ADP")


def test_locked_json_and_runtime_bls_configs_match_exactly() -> None:
    payload = approved_a2v1_teacher_bls_config()
    runtime = approved_a2v1_teacher_bls_runtime_config()

    for key, value in vars(runtime).items():
        observed = list(value) if isinstance(value, tuple) else value
        assert payload[key] == observed
    assert set(payload) == {*vars(runtime), "source_product_tag"}


def _synthetic_lc() -> HLSPLightCurve:
    n = 300
    time = np.linspace(0.0, 2.0, n)
    flux = 1.0 + 0.001 * np.sin(2.0 * np.pi * time)
    flux[80] = 1.25
    flux[220] = 1.30
    quality = np.zeros(n, dtype=np.int32)
    quality[10:20] = 1
    orbitid = np.where(time < 1.0, 119, 120).astype(np.int16)
    return HLSPLightCurve(
        tic=123,
        tmag=17.5,
        sector=56,
        cam=1,
        ccd=1,
        ra=np.nan,
        dec=np.nan,
        time=time,
        cadenceno=np.arange(n),
        orbitid=orbitid,
        quality=quality,
        flux={"DET_FLUX": flux},
        path=Path("synthetic.fits"),
    )


def test_bls_reports_cadence_cleaning_breakdown() -> None:
    lc = _synthetic_lc()
    cfg = BLSConfig(
        apertures=("DET_FLUX",),
        p_min_d=0.2,
        p_max_cap_d=1.0,
        max_period_fraction=0.45,
        durations_min=(8.0,),
        n_periods=200,
        n_peaks=3,
        min_cadences=50,
        sigma_clip=3.0,
        orbit_edge_trim_d=0.05,
    )

    result = run_bls_on_lc(lc, cfg, aperture="DET_FLUX")

    assert result.status == "ok"
    assert result.n_cad_total == 300
    assert result.n_cad_quality == 290
    assert result.n_cad_edge_trimmed > 0
    assert result.n_cad_sigma_clipped >= 1
    assert result.n_cad_kept == 290 - result.n_cad_edge_trimmed - result.n_cad_sigma_clipped
    assert np.isclose(result.quality_dropout_frac, 10 / 300)
    assert result.dropout_frac > result.quality_dropout_frac

    rows = result_to_rows(
        result,
        "synthetic",
        run_config={
            "search_branch": "quota_smoke",
            "n_periods": 200,
            "n_peaks": 3,
            "period_bin_edges": (0.2, 0.5, 1.0),
            "max_peaks_per_period_bin": 1,
            "durations_min": (8.0,),
        },
    )
    assert rows[0]["n_cad_quality"] == 290
    assert rows[0]["n_cad_edge_trimmed"] == result.n_cad_edge_trimmed
    assert rows[0]["n_cad_sigma_clipped"] == result.n_cad_sigma_clipped
    assert rows[0]["quality_dropout_frac"] == result.quality_dropout_frac
    assert rows[0]["bls_search_branch"] == "quota_smoke"
    assert rows[0]["bls_n_periods"] == 200
    assert rows[0]["bls_max_peaks_per_period_bin"] == 1
    assert rows[0]["bls_period_bin_edges"] == "[0.2, 0.5, 1.0]"


def test_candidate_rows_have_stable_bls_config_defaults() -> None:
    lc = _synthetic_lc()
    cfg = BLSConfig(
        apertures=("DET_FLUX",),
        p_min_d=0.2,
        p_max_cap_d=1.0,
        max_period_fraction=0.45,
        durations_min=(8.0,),
        n_periods=200,
        n_peaks=3,
        min_cadences=50,
    )
    result = run_bls_on_lc(lc, cfg, aperture="DET_FLUX")

    row = result_to_rows(result, "legacy-call-without-config")[0]

    assert row["bls_search_branch"] == ""
    assert row["bls_n_periods"] == 0
    assert row["bls_n_peaks"] == 0
    assert row["bls_max_peaks_per_period_bin"] == 0
    assert row["bls_period_bin_edges"] == ""
    assert np.isnan(row["bls_p_min_d"])


def _preparation_config() -> BLSConfig:
    return BLSConfig(
        apertures=("DET_FLUX",),
        p_min_d=0.12,
        p_max_cap_d=15.0,
        max_period_fraction=0.45,
        durations_min=(3.0,),
        n_periods=100,
        n_peaks=1,
        min_cadences=200,
        sigma_clip=5.0,
        orbit_edge_trim_d=0.0,
    )


def _lc_from_arrays(
    time: np.ndarray,
    flux: np.ndarray,
    quality: np.ndarray,
) -> HLSPLightCurve:
    return HLSPLightCurve(
        tic=123,
        tmag=17.5,
        sector=56,
        cam=1,
        ccd=1,
        ra=np.nan,
        dec=np.nan,
        time=np.asarray(time, dtype=float),
        cadenceno=np.arange(len(time)),
        orbitid=np.ones(len(time), dtype=np.int16),
        quality=np.asarray(quality, dtype=np.int32),
        flux={"DET_FLUX": np.asarray(flux, dtype=float)},
        path=Path("synthetic.fits"),
    )


def _prepare(lc: HLSPLightCurve, cfg: BLSConfig):
    return prepare_bls_inputs_from_arrays(
        time=lc.time,
        flux=lc.flux["DET_FLUX"],
        quality=lc.quality,
        orbitid=lc.orbitid,
        cfg=cfg,
    )


def test_bls_preparation_ignores_flagged_flux_and_requires_finite_time() -> None:
    cfg = _preparation_config()
    time = np.linspace(0.0, 2.0, 500)
    quality = np.r_[
        np.ones(300, dtype=np.int32),
        np.zeros(200, dtype=np.int32),
    ]
    flux = np.r_[
        np.zeros(300),
        1.0 + 1.0e-3 * np.sin(np.linspace(0.0, 8.0 * np.pi, 200)),
    ]
    lc = _lc_from_arrays(time, flux, quality)

    prepared = _prepare(lc, cfg)
    result = run_bls_on_lc(lc, cfg, aperture="DET_FLUX")

    assert prepared.status == "ok"
    assert prepared.n_cad_quality == 200
    assert np.isclose(prepared.flux_median, 1.0)
    assert result.status == prepared.status

    lc.time[-1] = np.nan
    prepared = _prepare(lc, cfg)
    result = run_bls_on_lc(lc, cfg, aperture="DET_FLUX")
    assert prepared.status == "too_few_cadences"
    assert prepared.n_cad_quality == 199
    assert result.status == prepared.status


def test_bls_preparation_rejects_zero_good_flux_median() -> None:
    cfg = _preparation_config()
    time = np.linspace(0.0, 2.0, 250)
    flux = np.zeros(250)
    quality = np.zeros(250, dtype=np.int32)
    lc = _lc_from_arrays(time, flux, quality)

    prepared = _prepare(lc, cfg)
    result = run_bls_on_lc(lc, cfg, aperture="DET_FLUX")

    assert prepared.status == "all_nan"
    assert prepared.n_cad_quality == 250
    assert prepared.flux_median == 0.0
    assert result.status == prepared.status


def test_bls_preparation_uses_post_clip_baseline_for_grid_validity() -> None:
    cfg = _preparation_config()
    time = np.r_[np.linspace(0.0, 0.2, 200), 2.0]
    flux = np.r_[
        1.0 + 1.0e-3 * np.sin(np.linspace(0.0, 8.0 * np.pi, 200)),
        2.0,
    ]
    quality = np.zeros(201, dtype=np.int32)
    lc = _lc_from_arrays(time, flux, quality)

    prepared = _prepare(lc, cfg)
    result = run_bls_on_lc(lc, cfg, aperture="DET_FLUX")

    assert prepared.status == "degenerate_grid"
    assert prepared.n_cad_quality == 201
    assert prepared.n_cad_sigma_clipped == 1
    assert prepared.n_cad_kept == 200
    assert np.isclose(prepared.baseline_d, 0.2)
    assert result.status == prepared.status


def test_bls_preparation_declines_clip_that_crosses_cadence_floor() -> None:
    cfg = _preparation_config()
    time = np.r_[np.linspace(0.0, 0.2, 199), 2.0]
    flux = np.r_[
        1.0 + 1.0e-3 * np.sin(np.linspace(0.0, 8.0 * np.pi, 199)),
        2.0,
    ]
    quality = np.zeros(200, dtype=np.int32)
    lc = _lc_from_arrays(time, flux, quality)

    prepared = _prepare(lc, cfg)
    result = run_bls_on_lc(lc, cfg, aperture="DET_FLUX")

    assert prepared.status == "ok"
    assert prepared.n_cad_sigma_clipped == 0
    assert prepared.n_cad_kept == 200
    assert np.isclose(prepared.baseline_d, 2.0)
    assert result.status == prepared.status
