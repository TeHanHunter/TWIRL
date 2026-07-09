from __future__ import annotations

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from twirl.vetting.two_aperture import _plot_periodogram  # noqa: E402


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
