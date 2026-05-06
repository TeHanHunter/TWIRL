"""Diagnostic plotting for the per-sector BLS first pass.

Standard 3-panel layout produced for every BLS search where we want a visual:
    Panel A — full SDE periodogram (period vs SDE), rank-1 peak + 1/2x and 2x
              harmonics marked.
    Panel B — zoom on the rank-1 peak (period vs SDE), ±1% in period.
    Panel C — phase-folded normalized flux at the rank-1 (P, T0), with the
              BLS box (depth, duration) overlaid.

Designed to be called from the sector_run orchestrator (per top candidate) or
from a standalone CLI on a saved periodogram NPZ.
"""
from __future__ import annotations

from pathlib import Path
from typing import Mapping

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from twirl.io.hlsp import BJDREFI, HLSPLightCurve, quality_mask, read_hlsp
from twirl.plotting.style import apply_twirl_style


def _phase_fold(time: np.ndarray, period: float, t0: float) -> np.ndarray:
    """Phase fold time array to [-0.5, 0.5] with t0 at phase 0."""
    return ((time - t0 + 0.5 * period) % period) / period - 0.5


def _bin_phase(phase: np.ndarray, flux: np.ndarray, n_bins: int = 80,
               window: float = 0.5) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Median + 1.4826*MAD/sqrt(n) bins over phase ∈ [-window, window]."""
    edges = np.linspace(-window, window, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    idx = np.digitize(phase, edges) - 1
    meds = np.full(n_bins, np.nan)
    errs = np.full(n_bins, np.nan)
    for i in range(n_bins):
        sel = idx == i
        n = sel.sum()
        if n == 0:
            continue
        y = flux[sel]
        meds[i] = np.median(y)
        if n > 1:
            mad = np.median(np.abs(y - meds[i]))
            errs[i] = 1.4826 * mad / np.sqrt(n)
    return centers, meds, errs


def plot_bls_diagnostic(
    spectrum: Mapping[str, np.ndarray],
    lc: HLSPLightCurve,
    aperture: str,
    peak: Mapping[str, float],
    out_path: Path,
    title: str | None = None,
) -> Path:
    """Render the standard 3-panel BLS diagnostic plot.

    Parameters
    ----------
    spectrum
        Periodogram dict produced by `run_bls_on_lc(..., return_periodogram=True)`
        (or loaded from the NPZ written by `save_periodogram`). Must have
        keys 'period', 'sde'.
    lc
        The HLSP light curve the BLS was run on.
    aperture
        Flux column name in `lc.flux` (e.g. 'DET_FLUX').
    peak
        Mapping with keys: 'period_d', 't0_bjd' (absolute BJD), 'duration_min',
        'depth', 'sde'. Typically a `BLSPeak` converted via dataclasses.asdict.
    out_path
        PNG path. Parent directories are created.
    title
        Optional figure suptitle. Default: "TIC <id> S<NN> cam/ccd Tmag — BLS rank 1".
    """
    apply_twirl_style("full_page")

    period = np.asarray(spectrum["period"], dtype=np.float64)
    sde = np.asarray(spectrum["sde"], dtype=np.float64)
    p_peak = float(peak["period_d"])
    t0_peak = float(peak["t0_bjd"])
    dur_min = float(peak["duration_min"])
    depth = float(peak["depth"])
    sde_peak = float(peak["sde"])

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(7.1, 6.6))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 1.1], hspace=0.45)
    ax_full = fig.add_subplot(gs[0])
    ax_zoom = fig.add_subplot(gs[1])
    ax_fold = fig.add_subplot(gs[2])

    # Panel A: full SDE periodogram on log-period axis.
    ax_full.semilogx(period, sde, lw=0.5, color="0.25", rasterized=True)
    ax_full.axhline(0, color="0.7", lw=0.5)
    ax_full.axvline(p_peak, color="#c4452f", lw=1.0, ls="--",
                    label=f"rank-1 P={p_peak:.5f} d")
    for mult, lbl in ((0.5, "1/2x"), (2.0, "2x")):
        ax_full.axvline(p_peak * mult, color="#c4452f", lw=0.5, ls=":", alpha=0.6,
                        label=lbl)
    ax_full.set_xlabel("period (d)")
    ax_full.set_ylabel("SDE")
    ax_full.set_xlim(period.min(), period.max())
    ax_full.legend(loc="upper right", fontsize=6)
    ax_full.set_title(f"BLS periodogram (n={period.size:,} trial periods)",
                      fontsize=8, loc="left")

    # Panel B: ±1% zoom on rank-1 peak.
    pad = 0.01 * p_peak
    sel = (period >= p_peak - pad) & (period <= p_peak + pad)
    ax_zoom.plot(period[sel], sde[sel], lw=0.7, color="0.2")
    ax_zoom.axvline(p_peak, color="#c4452f", lw=1.0, ls="--")
    ax_zoom.axhline(sde_peak, color="#c4452f", lw=0.5, ls=":")
    ax_zoom.set_xlabel("period (d)")
    ax_zoom.set_ylabel("SDE")
    ax_zoom.set_title(
        f"rank-1 peak: P={p_peak:.6f} d, SDE={sde_peak:.2f}", fontsize=8, loc="left"
    )

    # Panel C: phase-folded LC with the BLS box.
    mask = quality_mask(lc, aperture)
    t_rel = lc.time[mask]                       # BJD - 2457000
    flux = lc.flux[aperture][mask]
    med = np.nanmedian(flux)
    flux = flux / med if med > 0 else flux
    t_abs = t_rel + BJDREFI
    phase = _phase_fold(t_abs, p_peak, t0_peak)

    inside = np.abs(phase * p_peak * 1440.0) < 60.0  # ±60 min window
    ax_fold.scatter(phase[inside] * p_peak * 1440.0, flux[inside], s=1.7,
                    color="0.55", alpha=0.4, rasterized=True, zorder=1)
    centers, meds, errs = _bin_phase(phase[inside], flux[inside], n_bins=60,
                                     window=60.0 / (p_peak * 1440.0))
    ok = np.isfinite(meds)
    ax_fold.errorbar(centers[ok] * p_peak * 1440.0, meds[ok], yerr=errs[ok],
                     fmt="o", ms=2.8, mfc="#1f3b6b", mec="#1f3b6b",
                     ecolor="#1f3b6b", elinewidth=0.7, capsize=0, zorder=3,
                     label="binned")
    # BLS box overlay: a step-function with depth `depth` and width `dur_min`.
    ax_fold.axhline(1.0, color="0.65", lw=0.6, zorder=0)
    ax_fold.axhspan(1.0 - depth, 1.0, xmin=0, xmax=0, alpha=0)  # placeholder
    half = dur_min / 2.0
    box_x = [-60, -half, -half, half, half, 60]
    box_y = [1.0, 1.0, 1.0 - depth, 1.0 - depth, 1.0, 1.0]
    ax_fold.plot(box_x, box_y, color="#c4452f", lw=1.2, ls="-",
                 label=f"BLS box: depth={depth:.3f}, dur={dur_min:.1f} min")
    ax_fold.set_xlim(-60, 60)
    y_lo = float(np.nanpercentile(flux[inside], 1))
    y_hi = float(np.nanpercentile(flux[inside], 99.5))
    ax_fold.set_ylim(min(y_lo, 1.0 - depth) - 0.05, y_hi + 0.05)
    ax_fold.set_xlabel("phase from t0 (min)")
    ax_fold.set_ylabel("normalized flux")
    ax_fold.set_title(
        f"phase fold: P={p_peak:.6f} d  T0={t0_peak:.4f} BJD  ({aperture})",
        fontsize=8, loc="left",
    )
    ax_fold.legend(loc="lower right", fontsize=6)

    if title is None:
        title = (f"TIC {lc.tic}  S{lc.sector} cam{lc.cam}/ccd{lc.ccd}  "
                 f"T={lc.tmag:.2f}  —  BLS rank-1")
    fig.suptitle(title, fontsize=10, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return out_path


def render_from_periodogram_npz(
    npz_path: Path,
    hlsp_path: Path,
    out_path: Path,
    aperture: str = "DET_FLUX",
    peak_overrides: Mapping[str, float] | None = None,
) -> Path:
    """Convenience: load a saved periodogram + HLSP and render the diagnostic.

    The rank-1 peak is the argmax of the SDE array unless overridden via
    `peak_overrides` (useful for plotting a non-rank-1 peak from the
    candidates.parquet table).
    """
    with np.load(npz_path) as z:
        spectrum = {k: np.asarray(z[k]) for k in z.files}
    lc = read_hlsp(Path(hlsp_path))
    if lc is None:
        raise RuntimeError(f"could not read HLSP {hlsp_path}")
    sde = spectrum["sde"]
    period = spectrum["period"]
    idx = int(np.nanargmax(sde))
    peak = {
        "period_d": float(period[idx]),
        "t0_bjd": float(spectrum["t0"][idx]) + float(BJDREFI)
            if "t0" in spectrum else float("nan"),
        "duration_min": float(spectrum["duration"][idx]) * 1440.0
            if "duration" in spectrum else float("nan"),
        "depth": float(spectrum["depth"][idx]) if "depth" in spectrum else 0.0,
        "sde": float(sde[idx]),
    }
    if peak_overrides:
        peak.update({k: float(v) for k, v in peak_overrides.items()})
    return plot_bls_diagnostic(spectrum, lc, aperture, peak, out_path)
