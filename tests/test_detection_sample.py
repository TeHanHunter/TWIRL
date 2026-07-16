"""Deterministic, local end-to-end smoke for periodic candidate generation."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from twirl.io.hlsp import HLSPLightCurve
from twirl.search.bls import BLSConfig, run_bls_on_lc
from twirl.search.candidates import result_to_rows, write_parquet
from twirl.search.consolidate import ConsolidateConfig, consolidate_candidates


def test_synthetic_periodic_search_and_consolidation(tmp_path: Path) -> None:
    period_d = 1.5
    duration_d = 20.0 / 1440.0
    time = np.arange(0.0, 9.0, 5.0 / 1440.0)
    phase = ((time - 0.2 + 0.5 * period_d) % period_d) - 0.5 * period_d
    flux = np.ones(time.size, dtype=float)
    flux[np.abs(phase) <= 0.5 * duration_d] -= 0.2
    flux += 2.0e-4 * np.sin(2.0 * np.pi * time / 0.37)

    lc = HLSPLightCurve(
        tic=1856,
        tmag=16.0,
        sector=56,
        cam=4,
        ccd=1,
        ra=np.nan,
        dec=np.nan,
        time=time,
        cadenceno=np.arange(time.size),
        orbitid=np.where(time < 4.5, 119, 120).astype(np.int16),
        quality=np.zeros(time.size, dtype=np.int32),
        flux={"DET_FLUX_SML": flux, "DET_FLUX": flux.copy()},
        path=Path("synthetic.fits"),
    )
    config = BLSConfig(
        apertures=("DET_FLUX_SML", "DET_FLUX"),
        p_min_d=0.5,
        p_max_cap_d=3.0,
        max_period_fraction=0.45,
        durations_min=(16.0, 20.0, 30.0),
        n_periods=2_000,
        n_peaks=3,
        min_cadences=200,
        sigma_clip=5.0,
    )

    rows: list[dict[str, object]] = []
    for aperture in config.apertures:
        result = run_bls_on_lc(lc, config, aperture=aperture)
        assert result.status == "ok"
        assert result.peaks
        assert np.isclose(result.peaks[0].period_d, period_d, rtol=0.01)
        rows.extend(result_to_rows(result, "synthetic-detection-smoke"))

    candidates = tmp_path / "candidates.parquet"
    write_parquet(rows, candidates)
    consolidated = consolidate_candidates(
        candidates,
        ConsolidateConfig(sde_min=5.0),
    ).to_pandas()

    assert len(consolidated) >= 1
    best = consolidated.iloc[0]
    assert best["n_apertures_agree"] == 2
    assert set(best["apertures_agree"].split(",")) == {
        "DET_FLUX_SML",
        "DET_FLUX",
    }
    assert np.isclose(best["period_d"], period_d, rtol=0.01)
