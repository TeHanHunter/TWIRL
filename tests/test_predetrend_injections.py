from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

from twirl.lightcurves.tglc_h5_reader import APERTURE_KEYS, TGLCAperture, TGLCLightCurve


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_predetrend_script():
    path = REPO_ROOT / "scripts" / "stage3_injections" / "make_s56_predetrend_injection_set.py"
    spec = importlib.util.spec_from_file_location("make_s56_predetrend_injection_set", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _synthetic_tglc() -> TGLCLightCurve:
    rng = np.random.default_rng(123)
    n = 720
    time = np.linspace(0.0, 4.0, n)
    quality = np.zeros(n, dtype=np.int64)
    trend = 1200.0 + 60.0 * np.sin(2.0 * np.pi * time / 4.0)
    apertures = {}
    for idx, ap_key in enumerate(APERTURE_KEYS, start=1):
        flux = trend * idx + rng.normal(0.0, 10.0 * idx, n)
        err = np.full(n, 10.0 * idx)
        apertures[ap_key] = TGLCAperture(
            name=ap_key,
            raw_flux=flux,
            raw_flux_err=err,
            raw_magnitude=np.full(n, np.nan),
            raw_magnitude_err=np.full(n, np.nan),
            centroid_x=np.full(n, 10.0 + idx),
            centroid_y=np.full(n, 20.0 + idx),
            flux_was_synthesized=False,
        )
    lc = TGLCLightCurve(
        tic=123456789,
        sector=56,
        orbit=-1,
        cam=4,
        ccd=1,
        tmag=17.0,
        ra=1.0,
        dec=2.0,
        bjd_offset=2457000,
        time=time,
        cadence=np.arange(n),
        quality=quality,
        centroid_x=np.full(n, 10.0),
        centroid_y=np.full(n, 20.0),
        background=np.zeros(n),
        background_err=np.ones(n),
        apertures=apertures,
        path=None,
    )
    lc.orbitid = np.ones(n, dtype=np.int32) * 119
    return lc


def test_clone_with_raw_flux_preserves_original_and_orbitid() -> None:
    module = _load_predetrend_script()
    lc = _synthetic_tglc()
    new_primary = lc.apertures["Primary"].raw_flux - 50.0

    cloned = module._clone_with_raw_flux(lc, {"Primary": new_primary})

    assert np.allclose(cloned.apertures["Primary"].raw_flux, new_primary)
    assert not np.allclose(lc.apertures["Primary"].raw_flux, new_primary)
    assert np.array_equal(cloned.orbitid, lc.orbitid)
    assert np.array_equal(cloned.apertures["Small"].raw_flux, lc.apertures["Small"].raw_flux)


def test_detrended_output_columns_include_canonical_and_adp() -> None:
    module = _load_predetrend_script()
    lc = _synthetic_tglc()

    columns, diagnostics = module._detrended_output_columns(
        lc,
        adaptive_bkspace=0.3,
        adaptive_gap_split=0.2,
    )

    for name in (
        "DET_FLUX",
        "DET_FLUX_SML",
        "DET_FLUX_LAG",
        "DET_FLUX_ADP",
        "DET_FLUX_ADP_SML",
        "DET_FLUX_ADP_LAG",
    ):
        assert name in columns
        assert columns[name].shape == lc.time.shape
        assert np.isfinite(columns[name][lc.quality == 0]).any()
        assert name in diagnostics


def test_positive_baseline_requires_positive_good_median() -> None:
    module = _load_predetrend_script()
    quality = np.array([0, 0, 1, 0])

    assert module._positive_baseline(np.array([10.0, 11.0, 100.0, 12.0]), quality) == 11.0
    assert not np.isfinite(module._positive_baseline(np.array([-2.0, -1.0, 5.0, -3.0]), quality))
