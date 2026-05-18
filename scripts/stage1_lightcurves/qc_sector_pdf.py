#!/usr/bin/env python3
"""Per-sector HLSP quality-check PDF.

Pages:
  1. Sector summary: counts by Tmag bin, Tmag histogram, sky map.
  2. Photometric precision (1-hr CDPP-like RMS) vs Tmag, with Sullivan 2015
     model overlay and binned medians.
  3..N. Example detrended LCs grouped by Tmag bin. 200-s cadences after
     QUALITY mask in gray; 30-min binned median overlay in color.

Usage:
  qc_sector_pdf.py \\
      --hlsp-root /pdo/users/tehan/tglc-gpu-production/hlsp_s0056 \\
      --sector 56 \\
      --output reports/stage1_lightcurves/qc_pdf/s56_qc.pdf \\
      --n-precision 5000 --n-per-bin 8

Sullivan+2015 1-hr precision model:
  sigma_1hr^2 = sigma_phot^2 + sigma_sky^2 + sigma_read^2 + sigma_sys^2
  defaults from their Eq. 1, 1-hr equivalent (10 effective pixels).
"""

from __future__ import annotations

import argparse
import math
import random
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.io.hlsp import (  # noqa: E402
    iter_hlsp_fits,
    read_hlsp,
    tglc_mad_error,
)


# ----------------------------- Sullivan 2015 model -----------------------------

def sullivan_1hr_ppm(tmag: np.ndarray,
                     n_pix: float = 10.0,
                     read_noise: float = 10.0,  # e- per read
                     subexp_s: float = 2.0,     # 2-s subexposures
                     integration_s: float = 3600.0,
                     star_e_per_s_at_T10: float = 1.5e6,
                     sky_e_per_pix_per_s: float = 1.0e2,
                     sys_floor_ppm: float = 60.0) -> np.ndarray:
    """Sullivan+2015 Eq. 1 1-hour precision in ppm.

    Defaults: 10-pixel effective aperture, read noise 10 e-, 2-s subexposures,
    100 e-/pix/s sky, 1.5e6 e-/s star at Tmag=10, 60 ppm systematic floor.
    """
    tmag = np.asarray(tmag, dtype=float)
    f_star = star_e_per_s_at_T10 * 10.0 ** (-0.4 * (tmag - 10.0)) * integration_s
    f_sky = sky_e_per_pix_per_s * n_pix * integration_s
    n_reads = integration_s / subexp_s
    sigma_phot = np.sqrt(f_star) / f_star
    sigma_sky = np.sqrt(f_sky) / f_star
    sigma_read = read_noise * np.sqrt(n_reads * n_pix) / f_star
    sigma_sys = sys_floor_ppm * 1e-6
    sigma = np.sqrt(sigma_phot ** 2 + sigma_sky ** 2 + sigma_read ** 2 + sigma_sys ** 2)
    return sigma * 1e6


# ----------------------------- HLSP I/O -----------------------------
# Reader is centralized in twirl.io.hlsp; this script adapts it into the dict
# shape the rest of the QC PDF code expects, with a bitmask-overridable
# QUALITY filter so we can run the legacy "drop any nonzero" or a custom mask.


def good_mask(quality: np.ndarray, drop_bits: int = 0xFFFFFFFF) -> np.ndarray:
    """Default: any nonzero QUALITY is bad. Override with --quality-drop-bits."""
    return (np.asarray(quality, dtype=np.int64) & drop_bits) == 0


def read_lc(path: Path, drop_bits: int):
    """Adapter: load via twirl.io.hlsp.read_hlsp, return the dict the PDF
    code expects. ``ferr`` is the MAD-based per-cadence error from
    tglc_mad_error (HLSP's DET_FLUX_ERR is unreliable; many sectors ship it
    all-NaN or as a single constant value).

    ``HLSPLightCurve`` does not expose ``ra``/``dec`` on the dataclass, so
    we pull ``RA_OBJ``/``DEC_OBJ`` directly from the primary FITS header
    here. NaN if either is missing (preserved downstream by the precision
    panel; only used for the sky-coverage page, which falls back gracefully).
    """
    lc = read_hlsp(path)
    if lc is None:
        return None
    flux = lc.flux.get("DET_FLUX")
    if flux is None:
        flux = lc.flux.get("SAP_FLUX", np.zeros(len(lc.time)))
    # Keep negative-flux cadences. They are legitimate Poisson + background-
    # subtraction excursions at faint magnitudes; dropping them here exactly
    # reproduces the upstream TGLC `flux[flux <= 0] = nan` bug we patched out
    # in `twirl/preserve-negative-flux`. The 5-sigma MAD clip in
    # `sigma_clip_mask` catches genuine outliers; we don't need a hard `>0`
    # floor at the mask level. (cf. progress log 2026-05-13 Plan A entry.)
    mask = good_mask(lc.quality, drop_bits) & np.isfinite(flux)
    ferr = tglc_mad_error(lc, "DET_FLUX") if "DET_FLUX" in lc.flux else np.full(len(lc.time), np.nan)
    try:
        hdr = fits.getheader(path, ext=0)
        # QLP uses RA_OBJ / DEC_OBJ; TWIRL v3 HLSPs prior to commit XXXX used
        # RA / DEC (the writer now emits both, but old files may have only RA).
        ra_v = hdr.get("RA_OBJ", hdr.get("RA", np.nan))
        dec_v = hdr.get("DEC_OBJ", hdr.get("DEC", np.nan))
        ra = float(ra_v) if ra_v is not None else np.nan
        dec = float(dec_v) if dec_v is not None else np.nan
    except Exception:
        ra, dec = np.nan, np.nan
    return dict(path=lc.path, tic=lc.tic, tmag=lc.tmag, ra=ra, dec=dec,
                time=lc.time, flux=flux, ferr=ferr, quality=lc.quality, mask=mask)


def per_lc_precision_dmad(lc: dict, bin_minutes: float = 30.0,
                          cadence_s: float = 200.0) -> float:
    """Estimated photometric precision via the differential-MAD recipe used
    in the TGLC paper:

        sigma_per_cad = 1.48 * median(|diff(flux)|) / sqrt(2)
        sigma_bin     = sigma_per_cad / sqrt(N_per_bin)

    Insensitive to long-period stellar variability; isolates per-cadence
    photon+systematic noise. DET_FLUX is already median-normalized, so the
    output is a fractional precision (no division by mag scaling needed).
    """
    f = lc["flux"][lc["mask"]]
    if len(f) < 30:
        return np.nan
    df = np.abs(np.diff(f))
    if len(df) == 0:
        return np.nan
    sigma_per_cad = 1.48 * np.nanmedian(df) / math.sqrt(2.0)
    n_per_bin = max(1, int(round(bin_minutes * 60.0 / cadence_s)))
    return sigma_per_cad / math.sqrt(n_per_bin)


def per_lc_rms_ppm(lc: dict, cadence_s: float = 200.0) -> float:
    """1-hour-equivalent fractional RMS in ppm — used in example-LC titles.

    sigma_1hr = 1.4826 * MAD(flux/median) / sqrt(N_per_hr).
    """
    f = lc["flux"][lc["mask"]]
    if len(f) < 30:
        return np.nan
    med = np.median(f)
    if med <= 0:
        return np.nan
    mad = np.median(np.abs(f - med))
    sigma_per_cad = 1.4826 * mad / med
    n_per_hr = 3600.0 / cadence_s
    return sigma_per_cad / math.sqrt(n_per_hr) * 1e6


def bin_lc(time_: np.ndarray, flux: np.ndarray, mask: np.ndarray,
           bin_minutes: float = 30.0):
    """Median-bin a masked LC into bin_minutes-wide bins. Returns (t_b, f_b)."""
    t = time_[mask]
    f = flux[mask]
    if len(t) == 0:
        return np.array([]), np.array([])
    t0 = t.min()
    bin_d = bin_minutes / 60.0 / 24.0
    idx = ((t - t0) / bin_d).astype(int)
    edges = np.unique(idx)
    t_b = np.array([np.median(t[idx == e]) for e in edges])
    f_b = np.array([np.median(f[idx == e]) for e in edges])
    return t_b, f_b


# ----------------------------- Pages -----------------------------

TMAG_BINS = [(8, 12), (12, 14), (14, 16), (16, 18), (18, 20)]


def _aitoff_xy(ra_deg: np.ndarray, dec_deg: np.ndarray):
    coords = SkyCoord(ra=np.asarray(ra_deg) * u.deg,
                      dec=np.asarray(dec_deg) * u.deg, frame="icrs").galactic
    lon = -coords.l.wrap_at(180 * u.deg).radian
    lat = coords.b.radian
    return lon, lat


def _decorate_aitoff(ax, *, label_size=8, tick_size=7, grid_lw=0.55):
    xtk = np.array([-180, -120, -60, 0, 60, 120], dtype=float)
    ytk = np.array([-45, -30, -15, 0, 15, 30, 45, 60, 75], dtype=float)
    labels = ["", "120°", "60°", "0°", "300°", "240°"]
    ax.set_xticks(np.radians(xtk)); ax.set_xticklabels([])
    ax.set_yticks(np.radians(ytk))
    ax.tick_params(axis="both", labelsize=tick_size, colors="0.28", pad=4)
    ax.grid(color="0.85", linewidth=grid_lw)
    lats = np.linspace(-np.pi / 2, np.pi / 2, 361)
    for x in np.radians(xtk):
        ax.plot(np.full_like(lats, x), lats, color="black", lw=0.55,
                ls=":", alpha=0.75, zorder=5.2)
    ax.axhline(0.0, color="black", lw=0.75, ls="--", alpha=0.95, zorder=5.0)
    for x, lab in zip(np.radians(xtk), labels):
        if not lab:
            continue
        t = ax.text(x, 0.03, lab, ha="center", va="center", color="black",
                    fontsize=tick_size - 1, zorder=5.5)
        t.set_path_effects([pe.withStroke(linewidth=1.3,
                                          foreground=(1.0, 1.0, 1.0, 0.8))])


def _sector_hull_xy(ra_deg: np.ndarray, dec_deg: np.ndarray):
    """Convex-hull boundary of this-sector points in galactic-Aitoff space.
    Aitoff is non-Euclidean so we hull in (l_wrap_deg, b_deg) Cartesian then
    convert. Robust enough as a visual border, not a precise FOV polygon.
    """
    if len(ra_deg) < 4:
        return None
    coords = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs").galactic
    l = coords.l.wrap_at(180 * u.deg).deg
    b = coords.b.deg
    pts = np.column_stack([l, b])
    try:
        from scipy.spatial import ConvexHull
        hull = ConvexHull(pts)
        l_h = pts[hull.vertices, 0]
        b_h = pts[hull.vertices, 1]
    except Exception:
        return None
    l_h = np.append(l_h, l_h[0]); b_h = np.append(b_h, b_h[0])
    return -np.radians(l_h), np.radians(b_h)


def page_summary(pdf, sector: int, headers: list[dict], n_total: int,
                 catalog_path: Path | None = None,
                 observations_path: Path | None = None):
    """Page 1: Tmag distribution, table of counts. Page 2: Aitoff sky map
    showing all WDs faint, other-sector WDs as a lighter background, and the
    current sector highlighted with a convex-hull border."""
    tmags = np.array([h["tmag"] for h in headers if np.isfinite(h["tmag"])])
    ras_obs = np.array([h["ra"] for h in headers if np.isfinite(h["ra"])])
    decs_obs = np.array([h["dec"] for h in headers if np.isfinite(h["dec"])])

    # ---- Page 1: histogram + counts table
    fig = plt.figure(figsize=(8.5, 6.0))
    fig.suptitle(f"TWIRL Sector {sector} HLSP QC summary", fontsize=14)
    ax1 = fig.add_subplot(1, 1, 1)
    bins = np.arange(8, 21, 0.25)
    ax1.hist(tmags, bins=bins, color="#3b6ea5", alpha=0.85,
             edgecolor="white", linewidth=0.4)
    ax1.set_xlabel("TESS magnitude")
    ax1.set_ylabel("# of WD light curves")
    ax1.set_title(f"Tmag distribution  (N = {n_total})")
    ax1.set_xlim(8, 20.5)
    counts_lines = ["Tmag bin counts:"]
    for lo, hi in TMAG_BINS:
        n = int(((tmags >= lo) & (tmags < hi)).sum())
        counts_lines.append(f"  T∈[{lo:>2}, {hi:>2}): {n:>6}")
    ax1.text(0.02, 0.98, "\n".join(counts_lines), transform=ax1.transAxes,
             ha="left", va="top", fontsize=8, family="monospace",
             bbox=dict(facecolor="white", edgecolor="0.7", alpha=0.9))
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    pdf.savefig(fig); plt.close(fig)

    # ---- Page 2: Aitoff sky map
    fig = plt.figure(figsize=(8.5, 5.5))
    ax = fig.add_subplot(1, 1, 1, projection="aitoff")

    # Layer 1: full master catalog (lightest)
    if catalog_path is not None and catalog_path.is_file():
        try:
            from astropy.table import Table
            cat = Table.read(catalog_path, hdu=1, unit_parse_strict="silent")
            ra_all = np.asarray(cat["ra"], dtype=float)
            dec_all = np.asarray(cat["dec"], dtype=float)
            ok = np.isfinite(ra_all) & np.isfinite(dec_all)
            x_all, y_all = _aitoff_xy(ra_all[ok], dec_all[ok])
            ax.scatter(x_all, y_all, c="0.88", s=0.5, alpha=0.20,
                       linewidths=0, rasterized=True, zorder=1,
                       label="all WD catalog")
        except Exception as e:
            print(f"[qc-pdf] catalog load fail ({catalog_path}): {e}")

    # Layer 2: WDs in other sectors (sector != current)
    other_ra, other_dec = None, None
    if observations_path is not None and observations_path.is_file():
        try:
            from astropy.table import Table
            obs = Table.read(observations_path, hdu=1, unit_parse_strict="silent")
            sec_arr = np.asarray(obs["sector"], dtype=int)
            ra_o = np.asarray(obs["ra"], dtype=float)
            dec_o = np.asarray(obs["dec"], dtype=float)
            mask_other = (sec_arr != sector) & np.isfinite(ra_o) & np.isfinite(dec_o)
            other_ra = ra_o[mask_other]; other_dec = dec_o[mask_other]
            x_o, y_o = _aitoff_xy(other_ra, other_dec)
            ax.scatter(x_o, y_o, c="0.62", s=0.6, alpha=0.30,
                       linewidths=0, rasterized=True, zorder=2,
                       label="other sectors")
        except Exception as e:
            print(f"[qc-pdf] observations load fail: {e}")

    # Layer 3: this sector's WDs (from HLSP headers we already read)
    x_s, y_s = _aitoff_xy(ras_obs, decs_obs)
    ax.scatter(x_s, y_s, c="#c44d10", s=2.0, alpha=0.85,
               linewidths=0, rasterized=True, zorder=3,
               label=f"S{sector}  (N={len(ras_obs)})")

    # Layer 4: convex-hull border outline of the highlighted sector
    hull_drawn = False
    hull = _sector_hull_xy(ras_obs, decs_obs)
    if hull is not None:
        # Drop hull edges that wrap across the antimeridian (l = ±180°),
        # otherwise the polygon draws a long ribbon across the projection.
        x_h, y_h = hull
        for i in range(len(x_h) - 1):
            if abs(x_h[i + 1] - x_h[i]) < np.pi:
                ax.plot([x_h[i], x_h[i + 1]],
                        [y_h[i], y_h[i + 1]],
                        color="#7a2c08", lw=1.2, zorder=4)
                hull_drawn = True

    legend_handles = []
    from matplotlib.lines import Line2D
    legend_handles.append(Line2D([0], [0], marker="o", color="w",
                                  markerfacecolor="0.88", markersize=5,
                                  label="all WD catalog"))
    if other_ra is not None and len(other_ra) > 0:
        legend_handles.append(Line2D([0], [0], marker="o", color="w",
                                      markerfacecolor="0.62", markersize=5,
                                      label="other sectors"))
    legend_handles.append(Line2D([0], [0], marker="o", color="w",
                                  markerfacecolor="#c44d10", markersize=6,
                                  label=f"S{sector}  (N={len(ras_obs)})"))
    if hull_drawn:
        legend_handles.append(Line2D([0], [0], color="#7a2c08", lw=1.2,
                                      label="sector hull"))

    _decorate_aitoff(ax)
    ax.set_title(f"S{sector} sky coverage  (galactic Aitoff)", pad=14)
    ax.set_xlabel("Galactic longitude")
    ax.set_ylabel("Galactic latitude")
    ax.legend(handles=legend_handles, loc="lower left", fontsize=8,
              bbox_to_anchor=(0.0, -0.10))
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


def load_noisemodel(path: Path | None):
    """Sullivan+2015 σ_base(T) reference curve (Tmag, fractional precision)."""
    if path is None or not path.is_file():
        return None
    try:
        data = np.loadtxt(path)
        return data[:, 0], data[:, 1]
    except Exception as e:
        print(f"[qc-pdf] noisemodel load fail: {e}")
        return None


def page_precision(pdf, sector: int, sample_lcs: list[dict],
                   noisemodel_path: Path | None = None,
                   bin_minutes: float = 30.0):
    """30-min DMAD photometric precision vs Tmag, with σ_base baseline + ratio.

    Top panel: per-LC precision scatter (very faint), running median, σ_base.
    Bottom panel: precision / σ_base(T) ratio + running median.
    """
    from scipy.ndimage import median_filter

    prec = np.array([per_lc_precision_dmad(lc, bin_minutes=bin_minutes)
                     for lc in sample_lcs])
    tmag = np.array([lc["tmag"] for lc in sample_lcs])
    ok = np.isfinite(prec) & np.isfinite(tmag) & (prec > 0)
    prec = prec[ok]; tmag = tmag[ok]

    # Sort by Tmag for a smooth running median.
    order = np.argsort(tmag)
    tmag_s = tmag[order]; prec_s = prec[order]

    nm = load_noisemodel(noisemodel_path)
    sigma_base_at_T = None
    if nm is not None:
        sigma_base_at_T = np.interp(tmag_s, nm[0], nm[1])

    fig, axes = plt.subplots(2, 1, figsize=(8.5, 8.5), sharex=True,
                              gridspec_kw={"height_ratios": [3, 2], "hspace": 0.08})
    ax_t, ax_b = axes

    # Top — scatter + running median + σ_base
    ax_t.scatter(tmag_s, prec_s, s=2, alpha=0.05, color="#d2691e",
                 rasterized=True, label=f"TGLC HLSP (N={len(prec_s)})")
    win = max(50, len(prec_s) // 60)
    rmed = median_filter(prec_s, size=win, mode="nearest")
    ax_t.plot(tmag_s, rmed, "-", color="#c44d10", lw=2.0,
              path_effects=[pe.Stroke(linewidth=3.5, foreground="k"), pe.Normal()],
              label="running median")
    if nm is not None:
        ax_t.plot(nm[0], nm[1], "k-", lw=1.6, label=r"$\sigma_{\rm base}(T)$")
    ax_t.set_yscale("log")
    ax_t.set_ylabel("Estimated Photometric Precision (frac, "
                    f"{int(bin_minutes)}-min bin)")
    ax_t.set_ylim(1e-4, 1.0)
    ax_t.set_xlim(7, 20.5)
    ax_t.set_title(f"S{sector}  {int(bin_minutes)}-min bin")
    ax_t.grid(True, which="both", alpha=0.3)
    ax_t.legend(loc="upper left", fontsize=9)

    # Bottom — ratio
    if sigma_base_at_T is not None:
        ratio = prec_s / sigma_base_at_T
        ax_b.scatter(tmag_s, ratio, s=2, alpha=0.05, color="#d2691e",
                     rasterized=True)
        rmed_r = median_filter(ratio, size=win, mode="nearest")
        ax_b.plot(tmag_s, rmed_r, "-", color="#c44d10", lw=2.0,
                  path_effects=[pe.Stroke(linewidth=3.5, foreground="k"),
                                pe.Normal()],
                  label="running median")
        ax_b.axhline(1.0, color="k", lw=1.0)
        ax_b.set_ylim(0.5, 5.0)
        ax_b.set_yscale("log")
        ax_b.set_yticks([0.5, 1, 2, 5])
        ax_b.set_yticklabels(["0.5", "1", "2", "5"])
        ax_b.legend(loc="upper left", fontsize=9)
    ax_b.set_xlabel("TESS magnitude")
    ax_b.set_ylabel(r"Precision / $\sigma_{\rm base}(T)$")
    ax_b.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def sigma_clip_mask(flux: np.ndarray, base_mask: np.ndarray,
                    n_sigma: float = 5.0) -> np.ndarray:
    """Add a robust 5-sigma clip on top of the QUALITY mask using MAD."""
    out = base_mask.copy()
    if out.sum() < 10:
        return out
    f = flux[out]
    med = np.median(f)
    mad = np.median(np.abs(f - med))
    sig = 1.4826 * mad
    if sig <= 0:
        return out
    keep = np.abs(flux - med) < n_sigma * sig
    return out & keep


def page_lc_examples(pdf, sector: int, lcs: list[dict], bin_label: str):
    n = len(lcs)
    if n == 0:
        return
    cols = 2
    rows = int(math.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(8.5, max(2.5 * rows, 4)),
                              squeeze=False)
    fig.suptitle(f"S{sector} example LCs — {bin_label}  "
                 f"(QUALITY≠0 + 5σ clip removed)",
                 fontsize=12)
    for ax, lc in zip(axes.flatten(), lcs):
        m_q = lc["mask"]
        t = lc["time"]
        f = lc["flux"]
        m_clean = sigma_clip_mask(f, m_q, n_sigma=5.0)
        med = np.median(f[m_clean]) if m_clean.sum() else (np.median(f[m_q]) if m_q.sum() else 1.0)
        # Three layers:
        #   red    — QUALITY-flagged cadences (removed by flag)
        #   orange — sigma-clipped outliers (removed by outlier filter)
        #   gray   — kept 200 s cadences (after both filters)
        flagged = ~m_q & np.isfinite(f) & (f > 0)
        outlier = m_q & ~m_clean & np.isfinite(f) & (f > 0)
        if flagged.sum():
            ax.plot(t[flagged], f[flagged] / med, ".", ms=1.2, color="#c33",
                    alpha=0.6, label="QUALITY flag")
        if outlier.sum():
            ax.plot(t[outlier], f[outlier] / med, "x", ms=2.5, color="#e08e1e",
                    alpha=0.8, label="5σ outlier")
        ax.plot(t[m_clean], f[m_clean] / med, ".", ms=1.0, color="0.55",
                alpha=0.7, label="kept 200 s")
        # 30-min binned overlay (using cleaned mask)
        tb, fb = bin_lc(t, f, m_clean, bin_minutes=30.0)
        if len(tb):
            ax.plot(tb, fb / med, "-", color="#1f77b4", lw=1.3,
                    label="30-min median")
        rms = per_lc_rms_ppm(lc)
        title = (f"TIC {lc['tic']}  T={lc['tmag']:.2f}"
                 + (f"  1-hr RMS={rms:.0f} ppm" if np.isfinite(rms) else ""))
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("BJD - 2457000 (d)", fontsize=8)
        ax.set_ylabel("relative flux", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)
        # Scale y-axis from clean data only (flagged/outliers may extend off).
        if m_clean.sum() >= 5:
            f_clean_norm = f[m_clean] / med
            lo = np.percentile(f_clean_norm, 0.5)
            hi = np.percentile(f_clean_norm, 99.5)
            pad = max(0.02, 0.15 * (hi - lo))
            ax.set_ylim(lo - pad, hi + pad)
    for ax in axes.flatten()[len(lcs):]:
        ax.set_visible(False)
    # Single legend in figure
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=8,
                   bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    pdf.savefig(fig)
    plt.close(fig)


# ----------------------------- main -----------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hlsp-root", type=Path, required=True)
    ap.add_argument("--sector", type=int, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--n-precision", type=int, default=5000,
                    help="Random LCs sampled for precision-vs-Tmag panel.")
    ap.add_argument("--n-per-bin", type=int, default=8,
                    help="Example LCs per Tmag bin page.")
    ap.add_argument("--quality-drop-bits", type=lambda s: int(s, 0),
                    default=0xFFFFFFFF,
                    help="Bitmask of QUALITY bits considered 'bad' "
                         "(default: any nonzero). Pass 0 to keep all.")
    ap.add_argument("--noisemodel",
                    type=Path,
                    default=Path("/pdo/users/tehan/TWIRL/data_local/refs/noisemodel.dat"),
                    help="Sullivan+2015 σ_base(T) reference table (Tmag, frac).")
    ap.add_argument("--bin-minutes", type=float, default=30.0,
                    help="Bin width for the photometric-precision panel.")
    ap.add_argument("--catalog",
                    type=Path,
                    default=Path("/pdo/users/tehan/TWIRL/data_local/catalogs/twirl_master_catalog/twirl_wd_master_catalog_v0_tesscoverage.fits"),
                    help="Master WD catalog FITS for the sky-coverage background layer.")
    ap.add_argument("--observations",
                    type=Path,
                    default=Path("/pdo/users/tehan/TWIRL/data_local/catalogs/twirl_master_catalog/twirl_wd_tess_observations_v0.fits"),
                    help="WD TESS observations FITS for other-sector layer.")
    ap.add_argument("--data-output",
                    type=Path,
                    default=None,
                    help="Optional .npz path to save the per-LC precision arrays "
                         "alongside the PDF (default: <output>.npz).")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    random.seed(args.seed)

    paths = list(iter_hlsp_fits(args.hlsp_root))
    print(f"[qc-pdf] {len(paths)} HLSP FITS in {args.hlsp_root}")
    if not paths:
        print("[qc-pdf] no HLSP found; aborting")
        return 1

    # Pass 1: read all headers (cheap) for the summary page.
    headers = []
    for i, p in enumerate(paths):
        try:
            hdr = fits.getheader(p, ext=0)
            headers.append(dict(
                tmag=float(hdr.get("TESSMAG", np.nan)),
                ra=float(hdr.get("RA_OBJ", hdr.get("RA", np.nan))),
                dec=float(hdr.get("DEC_OBJ", hdr.get("DEC", np.nan))),
                path=p))
        except Exception:
            continue
        if (i + 1) % 2000 == 0:
            print(f"[qc-pdf]   header pass {i+1}/{len(paths)}")
    print(f"[qc-pdf] read {len(headers)} headers")

    # Sample for precision panel.
    n_prec = min(args.n_precision, len(headers))
    prec_sample_paths = [h["path"] for h in random.sample(headers, n_prec)]
    print(f"[qc-pdf] reading {n_prec} LCs for precision panel...")
    sample_lcs = []
    for i, p in enumerate(prec_sample_paths):
        lc = read_lc(p, args.quality_drop_bits)
        if lc is not None:
            sample_lcs.append(lc)
        if (i + 1) % 500 == 0:
            print(f"[qc-pdf]   precision read {i+1}/{n_prec}")
    print(f"[qc-pdf] precision sample ready: {len(sample_lcs)} LCs")

    # Per-bin example LCs.
    by_bin = {b: [] for b in TMAG_BINS}
    for h in headers:
        for lo, hi in TMAG_BINS:
            if np.isfinite(h["tmag"]) and lo <= h["tmag"] < hi:
                by_bin[(lo, hi)].append(h["path"])
                break
    bin_examples: dict[tuple, list[dict]] = {}
    for bin_, plist in by_bin.items():
        sample = random.sample(plist, min(args.n_per_bin, len(plist)))
        bin_examples[bin_] = [r for p in sample if (r := read_lc(p, args.quality_drop_bits)) is not None]
        print(f"[qc-pdf]   bin T∈[{bin_[0]},{bin_[1]}): {len(plist)} total, "
              f"sampled {len(bin_examples[bin_])}")

    # Render PDF.
    with PdfPages(args.output) as pdf:
        page_summary(pdf, args.sector, headers, len(paths),
                     catalog_path=args.catalog,
                     observations_path=args.observations)
        page_precision(pdf, args.sector, sample_lcs,
                       noisemodel_path=args.noisemodel,
                       bin_minutes=args.bin_minutes)
        for (lo, hi), lcs in bin_examples.items():
            page_lc_examples(pdf, args.sector, lcs, f"T∈[{lo}, {hi})")
    print(f"[qc-pdf] wrote {args.output}")

    # Companion NPZ — raw inputs to the precision plot, so the PDF can be
    # regenerated locally without re-reading 19k FITS.
    npz_path = args.data_output or args.output.with_suffix(".npz")
    prec_arr = np.array([per_lc_precision_dmad(lc, bin_minutes=args.bin_minutes)
                         for lc in sample_lcs])
    tmag_prec = np.array([lc["tmag"] for lc in sample_lcs])
    nm_data = load_noisemodel(args.noisemodel)
    np.savez_compressed(
        npz_path,
        sector=args.sector,
        bin_minutes=args.bin_minutes,
        n_total_hlsp=len(paths),
        tmag_all=np.array([h["tmag"] for h in headers]),
        ra_all=np.array([h["ra"] for h in headers]),
        dec_all=np.array([h["dec"] for h in headers]),
        prec_30min=prec_arr,
        tmag_prec_sample=tmag_prec,
        sigma_base_tmag=(nm_data[0] if nm_data is not None else np.array([])),
        sigma_base_value=(nm_data[1] if nm_data is not None else np.array([])),
    )
    print(f"[qc-pdf] wrote {npz_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
