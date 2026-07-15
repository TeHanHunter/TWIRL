from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from twirl.io.hlsp import HLSPLightCurve  # noqa: E402
from twirl.vetting.two_aperture import (  # noqa: E402
    _plot_even_odd_bins,
    _plot_folded_bins,
    _plot_full_lc,
    _plot_full_phase_fold,
    _plot_harmonic_fold,
    _plot_periodogram,
    _shared_light_curve_ylim,
)


def _two_aperture_lc() -> HLSPLightCurve:
    time = np.arange(2825.0, 2845.0, 1.0 / 1440.0)
    phase = ((time - 2825.2 + 0.5) % 1.0) - 0.5
    secondary_phase = np.abs(np.abs(phase) - 0.5)
    small = np.ones_like(time)
    primary = np.ones_like(time)
    small[np.abs(phase) < 0.02] = 0.72
    primary[np.abs(phase) < 0.02] = 0.64
    primary[secondary_phase < 0.01] = 0.82
    return HLSPLightCurve(
        tic=1,
        tmag=17.0,
        sector=56,
        cam=1,
        ccd=1,
        ra=0.0,
        dec=0.0,
        time=time,
        cadenceno=np.arange(len(time)),
        orbitid=np.ones(len(time), dtype=int),
        quality=np.zeros(len(time), dtype=int),
        flux={"DET_FLUX_ADP_SML": small, "DET_FLUX_ADP": primary},
        path=Path("toy.fits"),
    )


def test_periodogram_distinguishes_review_period_from_aperture_bls_peak() -> None:
    fig, ax = plt.subplots()
    spec = {
        "period": np.geomspace(0.1, 12.0, 100),
        "sde": np.linspace(0.0, 8.0, 100),
    }

    _plot_periodogram(
        ax,
        spec,
        bls_peak={"period_d": 6.1796258},
        review_peak={"period_d": 0.2786347},
    )

    labeled_lines = {
        line.get_label(): float(line.get_xdata()[0])
        for line in ax.lines
        if line.get_label() and not line.get_label().startswith("_")
    }
    assert labeled_lines["review P=0.2786 d"] == 0.2786347
    assert labeled_lines["aperture BLS max=6.18 d"] == 6.1796258
    assert ax.get_xlabel() == "Period (d)"
    plt.close(fig)


def test_periodogram_uses_one_marker_when_review_and_bls_period_match() -> None:
    fig, ax = plt.subplots()
    spec = {
        "period": np.geomspace(0.1, 12.0, 100),
        "sde": np.linspace(0.0, 8.0, 100),
    }

    _plot_periodogram(
        ax,
        spec,
        bls_peak={"period_d": 0.27863466},
        review_peak={"period_d": 0.27863467},
    )

    labeled_lines = [
        line for line in ax.lines if line.get_label() and not line.get_label().startswith("_")
    ]
    assert len(labeled_lines) == 1
    assert labeled_lines[0].get_label() == "review P=0.2786 d"
    plt.close(fig)


def test_review_and_odd_even_folds_use_same_phase_sampling_density() -> None:
    lc = _two_aperture_lc()
    aperture = "DET_FLUX_ADP_SML"
    fig, axes = plt.subplots(1, 2)

    _plot_folded_bins(
        axes[0],
        lc,
        aperture,
        period_d=1.0,
        t0_bjd=2459825.2,
        duration_min=60.0,
        title="Review fold",
    )
    _plot_even_odd_bins(
        axes[1],
        lc,
        aperture,
        period_d=1.0,
        t0_bjd=2459825.2,
        duration_min=60.0,
        title="Odd/even fold",
        metrics={},
        prefix="adp_sml",
    )

    review_bins = [line for line in axes[0].lines if line.get_marker() == "o"]
    parity_bins = [line for line in axes[1].lines if line.get_marker() in {"o", "s"}]
    assert len(review_bins) == 1
    assert len(parity_bins) == 2
    assert all(len(line.get_xdata()) == len(review_bins[0].get_xdata()) for line in parity_bins)
    assert len(axes[0].collections[0].get_offsets()) == len(axes[1].collections[0].get_offsets())
    assert axes[0].collections[0].get_sizes().tolist() == axes[1].collections[0].get_sizes().tolist()
    assert axes[0].collections[0].get_alpha() == axes[1].collections[0].get_alpha()
    plt.close(fig)


def test_every_light_curve_panel_accepts_one_shared_two_aperture_scale() -> None:
    lc = _two_aperture_lc()
    apertures = ("DET_FLUX_ADP_SML", "DET_FLUX_ADP")
    peak = {
        "period_d": 1.0,
        "t0_bjd": 2459825.2,
        "duration_min": 60.0,
        "depth": 0.36,
    }
    shared_ylim = _shared_light_curve_ylim(
        {aperture: lc for aperture in apertures},
        {aperture: peak for aperture in apertures},
    )
    assert shared_ylim[0] < 0.64
    assert shared_ylim[1] > 1.0

    fig, axes = plt.subplots(1, 5)
    _plot_full_lc(
        axes[0],
        lc,
        apertures[0],
        branch_name="current_adp",
        peak=peak,
        shared_ylim=shared_ylim,
    )
    _plot_full_phase_fold(
        axes[1],
        lc,
        apertures[0],
        period_d=1.0,
        t0_bjd=2459825.2,
        duration_min=60.0,
        depth=0.36,
        title="Full phase fold",
        shared_ylim=shared_ylim,
    )
    _plot_folded_bins(
        axes[2],
        lc,
        apertures[0],
        period_d=1.0,
        t0_bjd=2459825.2,
        duration_min=60.0,
        title="Review fold",
        shared_ylim=shared_ylim,
    )
    _plot_even_odd_bins(
        axes[3],
        lc,
        apertures[0],
        period_d=1.0,
        t0_bjd=2459825.2,
        duration_min=60.0,
        title="Odd/even fold",
        metrics={},
        prefix="adp_sml",
        shared_ylim=shared_ylim,
    )
    _plot_harmonic_fold(
        axes[4],
        lc,
        apertures[1],
        period_d=1.0,
        t0_bjd=2459825.2,
        duration_min=60.0,
        depth=0.36,
        title="P",
        show_ylabel=True,
        shared_ylim=shared_ylim,
    )

    assert all(np.allclose(ax.get_ylim(), shared_ylim) for ax in axes)
    plt.close(fig)
