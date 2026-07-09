from __future__ import annotations

from pathlib import Path

import numpy as np

from twirl.io.hlsp import HLSPLightCurve
from twirl.search.bls import BLSConfig, run_bls_on_lc
from twirl.search.candidates import result_to_rows


def test_bls_defaults_use_active_adp_pair() -> None:
    assert BLSConfig().apertures == ("DET_FLUX_ADP_SML", "DET_FLUX_ADP")


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
