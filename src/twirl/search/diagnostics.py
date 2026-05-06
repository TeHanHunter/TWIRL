"""Diagnostic plotting for the per-sector BLS first pass.

Two output products:

1. `plot_bls_diagnostic` — single-aperture 3-panel summary (full periodogram,
   peak zoom, phase fold). Used for quick checks on one (target, aperture).

2. `plot_vet_sheet` — TOI-style per-target vet page that stacks all three
   apertures plus an odd-vs-even transit overlay. This is the main per-target
   product saved alongside the candidates table; reviewers can scan it
   without rerunning BLS.

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


def plot_vet_sheet(
    lc: HLSPLightCurve,
    spectra: dict[str, dict],
    peaks: dict[str, dict],
    best_aperture: str,
    out_path: Path,
    cluster_summary: dict | None = None,
) -> Path:
    """TOI-style per-target vet sheet.

    Parameters
    ----------
    lc
        The HLSP light curve (with all 3 aperture columns loaded).
    spectra
        Dict {aperture -> spectrum dict (period/sde/...)} for each searched
        aperture. Apertures missing here are skipped in the plot.
    peaks
        Dict {aperture -> peak dict with period_d/t0_bjd/duration_min/depth/sde}.
        Comes from the candidates.parquet rank-1 row per aperture.
    best_aperture
        Which aperture's (P, T0) defines the canonical phase fold and the
        odd/even split. Usually the audit's `passing_aperture`, or the
        highest-SDE aperture when no audit was run.
    out_path
        PNG path.
    cluster_summary
        Optional dict from consolidated.parquet — adds n_apertures_agree and
        per-aperture SDE/rank to the bottom summary panel.
    """
    apply_twirl_style("full_page")

    aps = [a for a in ("DET_FLUX_SML", "DET_FLUX", "DET_FLUX_LAG") if a in spectra]
    if not aps:
        raise RuntimeError("plot_vet_sheet: no apertures provided")
    if best_aperture not in peaks:
        best_aperture = aps[0]
    p_best = float(peaks[best_aperture]["period_d"])
    t0_best = float(peaks[best_aperture]["t0_bjd"])
    dur_best = float(peaks[best_aperture]["duration_min"])
    depth_best = float(peaks[best_aperture]["depth"])

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Layout: 5 rows.
    #   row 0: full LC (single panel spanning all 3 columns)
    #   row 1: 3 periodograms
    #   row 2: 3 phase folds
    #   row 3: odd vs even fold (best aperture, single panel)
    #   row 4: summary text
    fig = plt.figure(figsize=(11.0, 12.0))
    gs = fig.add_gridspec(
        5, 3,
        height_ratios=[1.0, 1.1, 1.1, 1.0, 0.6],
        hspace=0.55, wspace=0.30,
    )

    # --- Row 0: full LC time series (best aperture only, with QUALITY mask)
    ax_lc = fig.add_subplot(gs[0, :])
    mask = quality_mask(lc, best_aperture)
    t_rel = lc.time[mask]
    flux = lc.flux[best_aperture][mask]
    med = np.nanmedian(flux)
    flux_n = flux / med if med > 0 else flux
    ax_lc.scatter(t_rel, flux_n, s=1.5, c="0.3", alpha=0.5, rasterized=True)
    ax_lc.axhline(1.0, color="0.6", lw=0.5)
    # Mark predicted transit times for the best (P, T0).
    n_lo = int(np.ceil((t_rel.min() + BJDREFI - t0_best) / p_best))
    n_hi = int(np.floor((t_rel.max() + BJDREFI - t0_best) / p_best))
    for n in range(n_lo, n_hi + 1):
        tt = (t0_best + n * p_best) - BJDREFI
        ax_lc.axvline(tt, color="#c4452f", lw=0.4, alpha=0.5, zorder=0)
    y_lo = float(np.nanpercentile(flux_n, 1))
    y_hi = float(np.nanpercentile(flux_n, 99.5))
    pad = 0.05 * (y_hi - y_lo) if y_hi > y_lo else 0.05
    ax_lc.set_ylim(y_lo - pad, y_hi + pad)
    ax_lc.set_xlim(t_rel.min(), t_rel.max())
    ax_lc.set_xlabel("BJD - 2457000 (d)")
    ax_lc.set_ylabel(f"{best_aperture}\nnorm. flux")
    ax_lc.set_title(
        f"Full sector light curve  ({best_aperture}, n_kept={mask.sum()})  "
        f"[red lines = predicted transits at P={p_best:.5f} d]",
        fontsize=8, loc="left",
    )

    # --- Row 1: periodograms per aperture
    for i, ap in enumerate(aps):
        ax = fig.add_subplot(gs[1, i])
        spec = spectra[ap]
        period = np.asarray(spec["period"], dtype=np.float64)
        sde = np.asarray(spec["sde"], dtype=np.float64)
        ax.semilogx(period, sde, lw=0.5, color="0.25", rasterized=True)
        ax.axhline(0, color="0.7", lw=0.4)
        if ap in peaks:
            p = float(peaks[ap]["period_d"])
            ax.axvline(p, color="#c4452f", lw=0.8, ls="--",
                       label=f"P={p:.5f} d")
            ax.axvline(p_best, color="#1f3b6b", lw=0.5, ls=":", alpha=0.7,
                       label=f"best ap P={p_best:.5f}")
            ax.legend(fontsize=6, loc="upper right")
        ax.set_xlim(period.min(), period.max())
        ax.set_xlabel("period (d)")
        ax.set_ylabel(f"{ap}\nSDE")
        title = ap
        if ap in peaks:
            title += f"  (rank-1 SDE={float(peaks[ap]['sde']):.1f})"
        ax.set_title(title, fontsize=8, loc="left")

    # --- Row 2: phase folds per aperture, all using best (P, T0)
    for i, ap in enumerate(aps):
        ax = fig.add_subplot(gs[2, i])
        m = quality_mask(lc, ap)
        t_abs = lc.time[m] + BJDREFI
        f = lc.flux[ap][m] / np.nanmedian(lc.flux[ap][m])
        phase = _phase_fold(t_abs, p_best, t0_best)
        in_win = np.abs(phase * p_best * 1440.0) < 60.0
        ax.scatter(phase[in_win] * p_best * 1440.0, f[in_win],
                   s=1.7, color="0.55", alpha=0.4, rasterized=True, zorder=1)
        c, m_, e_ = _bin_phase(phase[in_win], f[in_win], n_bins=60,
                                window=60.0 / (p_best * 1440.0))
        ok = np.isfinite(m_)
        ax.errorbar(c[ok] * p_best * 1440.0, m_[ok], yerr=e_[ok],
                    fmt="o", ms=2.6, mfc="#1f3b6b", mec="#1f3b6b",
                    ecolor="#1f3b6b", elinewidth=0.6, capsize=0, zorder=3)
        # BLS box from best aperture.
        half = dur_best / 2.0
        ax.plot([-60, -half, -half, half, half, 60],
                [1.0, 1.0, 1.0 - depth_best, 1.0 - depth_best, 1.0, 1.0],
                color="#c4452f", lw=1.0, zorder=2)
        ax.axhline(1.0, color="0.65", lw=0.4)
        ax.set_xlim(-60, 60)
        if in_win.any():
            ymin = float(min(np.nanpercentile(f[in_win], 1), 1.0 - depth_best))
            ymax = float(np.nanpercentile(f[in_win], 99.5))
            ax.set_ylim(ymin - 0.05, ymax + 0.05)
        ax.set_xlabel("phase from t0 (min)")
        ax.set_ylabel(f"{ap}\nflux")
        ax.set_title(f"phase fold @ best (P,T0) — {ap}", fontsize=8, loc="left")

    # --- Row 3: odd vs even transit fold (best aperture)
    ax_oe = fig.add_subplot(gs[3, :])
    m = quality_mask(lc, best_aperture)
    t_abs = lc.time[m] + BJDREFI
    f = lc.flux[best_aperture][m] / np.nanmedian(lc.flux[best_aperture][m])
    transit_index = np.floor((t_abs - t0_best) / p_best + 0.5).astype(int)
    is_odd = (transit_index % 2) != 0
    phase_min = ((t_abs - t0_best + 0.5 * p_best) % p_best - 0.5 * p_best) * 1440.0
    in_win = np.abs(phase_min) < 60.0

    odd_sel = in_win & is_odd
    even_sel = in_win & ~is_odd
    ax_oe.scatter(phase_min[odd_sel], f[odd_sel], s=2, color="#c4452f",
                  alpha=0.35, label="odd transits", rasterized=True, zorder=1)
    ax_oe.scatter(phase_min[even_sel], f[even_sel], s=2, color="#1f3b6b",
                  alpha=0.35, label="even transits", rasterized=True, zorder=1)

    def _bin(phase, flux, label, color):
        c, mvals, evals = _bin_phase(phase, flux, n_bins=40,
                                      window=60.0 / (p_best * 1440.0))
        ok = np.isfinite(mvals)
        ax_oe.errorbar(c[ok] * p_best * 1440.0, mvals[ok], yerr=evals[ok],
                       fmt="o", ms=4, mfc=color, mec=color, ecolor=color,
                       elinewidth=0.8, capsize=0, zorder=3, label=f"{label} binned")
        # Median in-transit depth for this subset.
        in_tr = np.abs(c * p_best * 1440.0) < dur_best / 2
        if (in_tr & ok).any():
            return float(np.nanmedian(mvals[in_tr & ok]))
        return float("nan")

    odd_depth = _bin(_phase_fold(t_abs[odd_sel], p_best, t0_best),
                      f[odd_sel], "odd", "#c4452f")
    even_depth = _bin(_phase_fold(t_abs[even_sel], p_best, t0_best),
                       f[even_sel], "even", "#1f3b6b")
    half = dur_best / 2.0
    ax_oe.plot([-60, -half, -half, half, half, 60],
               [1.0, 1.0, 1.0 - depth_best, 1.0 - depth_best, 1.0, 1.0],
               color="0.4", lw=0.8, ls="--", zorder=2)
    ax_oe.axhline(1.0, color="0.65", lw=0.4)
    ax_oe.set_xlim(-60, 60)
    ax_oe.set_xlabel("phase from t0 (min)")
    ax_oe.set_ylabel(f"{best_aperture}\nodd / even")
    ax_oe.legend(fontsize=7, loc="lower right", ncol=2)
    n_odd = int(odd_sel.sum())
    n_even = int(even_sel.sum())
    odd_minus_even = (odd_depth - even_depth) if (
        np.isfinite(odd_depth) and np.isfinite(even_depth)
    ) else float("nan")
    ax_oe.set_title(
        f"odd/even transit fold (EB check)  "
        f"odd_med={odd_depth:.3f} (n={n_odd})  "
        f"even_med={even_depth:.3f} (n={n_even})  "
        f"Δ={odd_minus_even:.3f}",
        fontsize=8, loc="left",
    )

    # --- Row 4: summary text panel
    ax_tx = fig.add_subplot(gs[4, :])
    ax_tx.axis("off")
    lines = [
        f"TIC {lc.tic}   Sector {lc.sector}   cam{lc.cam}/ccd{lc.ccd}   "
        f"Tmag={lc.tmag:.2f}",
        f"best aperture: {best_aperture}",
        f"  P     = {p_best:.6f} d",
        f"  T0    = {t0_best:.6f} BJD",
        f"  dur   = {dur_best:.2f} min",
        f"  depth = {depth_best:.3f}",
        f"  SDE   = {float(peaks[best_aperture]['sde']):.2f}",
    ]
    if cluster_summary is not None:
        lines.append(
            f"cross-aperture: n_apertures_agree="
            f"{cluster_summary.get('n_apertures_agree', '?')}  "
            f"({cluster_summary.get('apertures_agree', '')})"
        )
        for ap in ("DET_FLUX_SML", "DET_FLUX", "DET_FLUX_LAG"):
            sde_ap = cluster_summary.get(f"sde_{ap}")
            if sde_ap is not None and np.isfinite(sde_ap):
                lines.append(
                    f"  {ap:<14} SDE={float(sde_ap):.2f}  "
                    f"rank={int(cluster_summary.get(f'rank_{ap}', 0))}"
                )
    if np.isfinite(odd_minus_even):
        verdict = ("ok" if abs(odd_minus_even) < 0.05
                   else "*** SUSPICIOUS odd/even Δ — possible EB ***")
        lines.append(f"odd-even Δdepth = {odd_minus_even:+.3f}   [{verdict}]")
    ax_tx.text(0.0, 1.0, "\n".join(lines), transform=ax_tx.transAxes,
               va="top", ha="left", fontsize=8, family="monospace")

    fig.suptitle(
        f"TIC {lc.tic}  S{lc.sector} cam{lc.cam}/ccd{lc.ccd}  T={lc.tmag:.2f}  —  TWIRL BLS vet sheet",
        fontsize=11, y=0.997,
    )
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    # Also save the same figure as a PDF alongside the PNG. PDFs are lighter
    # for archival vetting (vector text + scatter rendered as raster), and
    # downstream collaborators tend to prefer them for paper figures.
    pdf_path = out_path.with_suffix(".pdf")
    try:
        fig.savefig(pdf_path, bbox_inches="tight")
    except Exception:
        pass
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
