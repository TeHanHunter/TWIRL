import os
from glob import glob
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

# ---------------- Helpers ----------------

def fold_phase(t, t0, P):
    # Phase in [-0.5, 0.5)
    return ((t - t0) / P + 0.5) % 1.0 - 0.5

def phase_bin_equal(phase, y, yerr, binsize_phase=None, nbins=100):
    """
    Bin AFTER phase fold on the SAME domain you plot: [-0.5, 0.5).
    Each binned point sits at the inverse-variance weighted mean phase in that bin.
    """
    phase = np.asarray(phase); y = np.asarray(y); yerr = np.asarray(yerr)
    m = np.isfinite(phase) & np.isfinite(y) & np.isfinite(yerr) & (yerr > 0)
    if not np.any(m): return np.array([]), np.array([]), np.array([])
    ph = phase[m]; yy = y[m]; ye = yerr[m]

    if binsize_phase is not None:
        nbins = max(2, int(np.ceil(1.0 / binsize_phase)))
    edges = np.linspace(-0.5, 0.5, nbins + 1)

    idx = np.digitize(ph, edges, right=False) - 1
    idx = np.clip(idx, 0, nbins - 1)

    w = 1.0 / (ye**2)
    sumw   = np.bincount(idx, weights=w,    minlength=nbins)
    sumwy  = np.bincount(idx, weights=w*yy, minlength=nbins)
    sumwph = np.bincount(idx, weights=w*ph, minlength=nbins)

    good = sumw > 0
    x  = np.full(nbins, np.nan); yb = np.full(nbins, np.nan); eb = np.full(nbins, np.nan)
    x[good]  = sumwph[good] / sumw[good]
    yb[good] = sumwy[good]  / sumw[good]
    eb[good] = 1.0 / np.sqrt(sumw[good])

    x, yb, eb = x[good], yb[good], eb[good]
    order = np.argsort(x)
    return x[order], yb[order], eb[order]

def snr_from_phase_window(phase, y, yerr, half_width_phase, ootr_buffer_phase):
    """
    Compute SNR in the style of snr = sqrt(N_points) * depth / noise.
    - Counts in-transit points, weights them by instantaneous depth
      relative to mid-transit depth.
    - Uses the mid-transit depth as the 'signal' amplitude.
    - Assumes roughly constant error bars (take median yerr).
    Returns: depth_center, ferr, snr, N_eff
    """
    phase = np.asarray(phase); y = np.asarray(y); yerr = np.asarray(yerr)
    m = np.isfinite(phase) & np.isfinite(y) & np.isfinite(yerr) & (yerr > 0)
    if not np.any(m):
        return np.nan, np.nan, np.nan, np.nan
    ph, yy, ye = phase[m], y[m], yerr[m]

    # Define windows
    in_tr = np.abs(ph) <= half_width_phase
    oot   = np.abs(ph) >= (half_width_phase + ootr_buffer_phase)
    if (np.sum(in_tr) < 1) or (np.sum(oot) < 1):
        return np.nan, np.nan, np.nan, np.nan

    # Out-of-transit baseline
    mu_oot = np.mean(yy[oot])

    # Mid-transit depth (use points nearest phase=0)
    i_mid = np.argmin(np.abs(ph))
    depth_center = mu_oot - yy[i_mid]

    # Weight each in-transit point by its instantaneous relative depth
    inst_depths = mu_oot - yy[in_tr]
    rel_depths = np.clip(inst_depths / depth_center, 0.0, 1.0)
    N_eff = np.sum(rel_depths)  # effective number of points

    # Per-point error (assume constant across points)
    ferr = np.median(ye[in_tr])

    # SNR in the style of theory: sqrt(N_eff) * depth_center / ferr
    snr = np.sqrt(N_eff) * depth_center / ferr
    return depth_center, ferr, snr, N_eff

def extract_sector(hdul):
    try:
        return int(hdul[0].header["SECTOR"])
    except Exception:
        return None

# ---------------- Core plotting ----------------

def plot_pf_lc(local_directory, period, mid_transit_tbjd, kind="cal_aper_flux",
               fits_file=None, outname=None, ylimits=(0.3, 1.5),
               binsize_seconds=100.0,
               transit_duration_min=7.998,  # full T14 (min)
               transit_duration_min_err=0.023,
               zoom_n_halfwidths=3.0):
    files = [fits_file] if fits_file is not None else glob(os.path.join(local_directory, "*.fits"))
    plots_dir = os.path.join(local_directory, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    for j, fpath in enumerate(files):
        with fits.open(fpath, mode="denywrite") as hdul:
            sec = extract_sector(hdul) or "unknown"
            data = hdul[1].data
            if len(data[kind]) != len(data["time"]):
                print(f"Skip {fpath}: {kind} length != time length"); continue

            q = (data["TESS_flags"] == 0) & (data["TGLC_flags"] == 0)
            t = np.asarray(data["time"][q], float)
            f = np.asarray(data[kind][q], float)

            ferr_val = float(hdul[1].header["CAPE_ERR"])
            ferr = np.full_like(f, ferr_val, dtype=float)

            # Phase fold (for SNR windowing)
            phi = fold_phase(t, mid_transit_tbjd, period)  # [-0.5, 0.5)

            # Display binning (in phase)
            binsize_phase = (binsize_seconds/86400.0) / period
            ph_b, f_b, fe_b = phase_bin_equal(phi, f, ferr, binsize_phase=binsize_phase)

            # Window in phase for SNR; convert to minutes only for plotting
            half_width_phase = ( (transit_duration_min/2.0) / 1440.0 ) / period
            ootr_buffer_phase = 3 * half_width_phase

            depth, sigma_depth, snr, mu_oot = snr_from_phase_window(
                phi, f, ferr,
                half_width_phase=half_width_phase,
                ootr_buffer_phase=ootr_buffer_phase
            )

            # ---- X axis in minutes from mid-transit ----
            day2min = 1440.0
            minutes = phi * period * day2min
            minutes_b = ph_b * period * day2min if len(ph_b) else np.array([])

            half_width_min = transit_duration_min / 2.0
            xlim = zoom_n_halfwidths * half_width_min

            fig = plt.figure(figsize=(8,5))
            ax = plt.gca()

            # raw
            ax.errorbar(minutes, f, ferr, c="silver", ls="", elinewidth=0.1, marker=".", ms=3,
                        zorder=2, label="200-s Raw")
            # binned
            if len(minutes_b):
                ax.errorbar(minutes_b, f_b, fe_b, c=f"C{j%10}", ls="", elinewidth=1.2, marker=".", ms=8,
                            zorder=3, label=f"{int(binsize_seconds)}-s Bin")

            # Guides in minutes
            ax.axvline(0.0, ymin=0, ymax=1, ls="dotted", color="grey", zorder=1)
            ax.axvspan(-half_width_min, +half_width_min, color="0.9", zorder=1, alpha=0.5)

            ax.set_xlim(-60, 60)
            if ylimits: ax.set_ylim(*ylimits)
            ax.set_xlabel("Time from mid-transit [min]")
            ax.set_ylabel("Normalized flux")
            if len(minutes_b): ax.legend(loc="best")

            if np.isfinite(snr):
                txt = (r"SNR$_\mathrm{tr}$ = " + f"{snr:.1f}   "
                       r"$\Delta$ = " + f"{depth:.4f} ± {sigma_depth:.4f}")
                ax.text(0.98, 0.03, txt, ha="right", va="bottom", transform=ax.transAxes)

            ax.set_title(f"WD 1856 b, TGLC Sector {sec}")

            save_as = outname or os.path.join(plots_dir, f"WD_1856_s{sec}.pdf")
            plt.savefig(save_as, dpi=300)
            plt.close(fig)

# ---------------- Driver ----------------

def main(local_directory, period, mid_transit_tbjd, kind="cal_aper_flux", outdir="plots"):
    os.makedirs(outdir, exist_ok=True)
    fits_files = sorted(glob(os.path.join(local_directory, "lc", "*.fits")))
    if not fits_files:
        fits_files = sorted(glob(os.path.join(local_directory, "*.fits")))

    for f in fits_files:
        print(f"Processing {f}")
        try:
            plot_pf_lc(
                local_directory=os.path.dirname(f),
                period=period,
                mid_transit_tbjd=mid_transit_tbjd,
                kind=kind,
                fits_file=f,
                outname=None,
                binsize_seconds=100.0,       # display bin width (does not affect SNR)
                transit_duration_min=7.998,  # <-- half measured T14
                transit_duration_min_err=0.023,
                zoom_n_halfwidths=3.0
            )
        except Exception as e:
            print(f"Failed on {f}: {e}")

if __name__ == "__main__":
    local_directory = "/Users/tehan/Documents/TWIRL/WD_1856/"
    period = 1.40794050                # days
    mid_transit_tbjd = 1779.375082800 - 0.001*period      # TBJD/BTJD consistent with your FITS
    # [BLS] Refined period: 1.407933459 d (from 1.407940500); t0: 1779.381033931 (from 1779.375082800)
    kind = "cal_aper_flux"
    main(local_directory, period, mid_transit_tbjd, kind=kind)