from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from docutils.nodes import legend
from matplotlib.colors import Normalize
import matplotlib.patheffects as pe


def load_wd_data(path, t_limit=19.0, poe_min=None, clip_pct=99.5):
    t = Table.read(path, hdu="Joined", unit_parse_strict="silent")
    t.convert_bytestring_to_unicode()
    ra = np.asarray(t["ra"])
    dec = np.asarray(t["dec"])
    plx = np.asarray(t["parallax"])
    g = np.asarray(t["phot_g_mean_mag"])
    bp = np.asarray(t["phot_bp_mean_mag"]) if "phot_bp_mean_mag" in t.colnames else np.full_like(g, np.nan, dtype=float)
    rp = np.asarray(t["phot_rp_mean_mag"]) if "phot_rp_mean_mag" in t.colnames else np.full_like(g, np.nan, dtype=float)
    poe = np.asarray(t["parallax_over_error"]) if "parallax_over_error" in t.colnames else None
    designation = np.array([
        int(val.split()[-1]) if val.split()[-1].isdigit() else -1
        for val in t["designation"]
    ], dtype=int)
    # distance (needed for mask later)
    plx_ok = np.isfinite(plx) & (plx > 0)

    # Gaia -> TESS (fallback if NaN)
    dif = bp - rp
    tess = g - 0.00522555 * dif ** 3 + 0.0891337 * dif ** 2 - 0.633923 * dif + 0.0324473
    bad = ~np.isfinite(tess)
    if np.any(bad):
        tess[bad] = g[bad] - 0.430

    # mask
    m = np.isfinite(ra) & np.isfinite(dec) & plx_ok & np.isfinite(g) & (tess < t_limit)
    if poe is not None and poe_min is not None:
        m &= np.isfinite(poe) & (poe > poe_min)

    # apply mask
    ra, dec, plx, g, bp, rp, tess, designation = ra[m], dec[m], plx[m], g[m], bp[m], rp[m], tess[m], designation[m]

    # distance in pc (apply clip after mask)
    d_pc = 1000.0 / plx
    if clip_pct is not None:
        d_pc = np.minimum(d_pc, np.percentile(d_pc, clip_pct))
    print(f'Total number of WDs: {len(ra)}')
    return ra, dec, plx, d_pc, tess, designation


def plot_2d_polar(dec_deg, dist_pc, tess_mag,
                  rmax=160, r_pct=99, max_points=250_000,
                  s_min=2, s_max=4, alpha=0.6, tick_count=5,
                  cmap="viridis", wd_bp=None, wd_rp=None, wd_g=16.958):
    """Declination as angle, distance as radius. Seaborn theme; colorbar in TESS mag."""
    sns.set_theme(style="whitegrid", context="talk",
                  rc={'font.family': 'serif', 'font.serif': ['DejaVu Serif'], 'font.size': 8,
                      'axes.edgecolor': '0.2', 'axes.labelcolor': '0.', 'xtick.color': '0.', 'ytick.color': '0.',
                      'axes.facecolor': '1', 'grid.color': '0.9'})

    n = dist_pc.size
    if n > max_points:
        idx = np.random.default_rng(42).choice(n, size=max_points, replace=False)
        dec_deg, dist_pc, tess_mag = dec_deg[idx], dist_pc[idx], tess_mag[idx]

    theta = np.deg2rad(dec_deg)
    r = dist_pc

    if rmax is None:
        rmax = np.percentile(r, r_pct)
    r = np.clip(r, 0, rmax)

    # color/size from TESS magnitude
    t_clip = np.clip(tess_mag, tess_mag.min(), tess_mag.max())
    sizes = np.interp(t_clip, (t_clip.min(), t_clip.max()), (s_max, s_min))
    colors = t_clip

    fig, ax = plt.subplots(figsize=(11, 6), subplot_kw={'projection': 'polar'})
    ax.set_thetalim(-np.pi / 2, np.pi / 2)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(1)

    sc = ax.scatter(theta, r, c=colors, cmap=cmap,
                    s=sizes, alpha=alpha, linewidths=0, rasterized=True)

    # --- WD 1856+534 using same cmap & norm (compute T if bp/rp given; else fallback G-0.430) ---
    wd_dec = 53 + 30 / 60 + 33.3 / 3600  # ≈ 53.509°
    wd_dist = 1000 / 40.3983  # ≈ 24.77 pc
    if (wd_bp is not None) and (wd_rp is not None):
        wd_dif = float(wd_bp - wd_rp)
        wd_t = wd_g - 0.00522555 * wd_dif ** 3 + 0.0891337 * wd_dif ** 2 - 0.633923 * wd_dif + 0.0324473
        if not np.isfinite(wd_t): wd_t = wd_g - 0.430
    else:
        wd_t = wd_g - 0.430
    wd_t = float(np.clip(wd_t, t_clip.min(), t_clip.max()))
    wd_theta = np.deg2rad(wd_dec)

    ax.scatter(wd_theta, wd_dist,
               c=[wd_t], cmap=sc.cmap, norm=sc.norm,
               s=60, edgecolor="k", linewidth=0.8,
               marker="o", alpha=0.95, zorder=5, label="WD 1856")
    ax.text(wd_theta + 0.15, wd_dist + 5, "WD 1856", ha="right", va="center", fontsize=10)

    # ticks
    dec_ticks = [90, 60, 30, 0, -30, -60, -90]
    ax.set_thetagrids(dec_ticks, labels=[f"{d}°" for d in dec_ticks])
    ax.set_rlim(0, rmax)
    rticks = np.linspace(0, rmax, tick_count)
    ax.set_rticks(rticks)
    ax.set_rlabel_position(0)
    ax.tick_params(labelsize=10)

    ax.grid(color="0.7", linewidth=0.8)
    ax.spines['polar'].set_color("0.4")
    ax.text(np.deg2rad(102), 0.5 * rmax, "Distance (pc)", ha="center", va="center", fontsize=10)
    ax.text(np.deg2rad(80), 1.15 * rmax, "Dec (degree)", ha="center", va="bottom", rotation=90, fontsize=10)

    # colorbar (TESS mag)
    # --- horizontal colorbar aligned with main plot ---
    pos = ax.get_position()
    cbar_ax = fig.add_axes([pos.x0, pos.y0 + 0.08, pos.width, 0.03])
    cb = fig.colorbar(sc, cax=cbar_ax, orientation="horizontal")

    # Force ticks including T=18
    cb.set_ticks([10, 12, 14, 16, 18])
    cb.set_label("TESS magnitude (T)", fontsize=10)
    cb.ax.tick_params(labelsize=10)

    plt.savefig("plots/wd_dist_dec_Tmag.pdf", dpi=300, bbox_inches="tight")
    return fig, ax


def _to_numpy(x):
    # unwrap masked arrays -> regular ndarray
    return np.asanyarray(getattr(x, "data", x))


def _finite_pos(x):
    x = _to_numpy(x).astype(float)
    return np.isfinite(x) & (x > 0)


def _finite_any(x):
    x = _to_numpy(x).astype(float)
    return np.isfinite(x)


def plot_2d_polar_two(
        dec_deg_L, dist_pc_L, tess_mag_L,
        dec_deg_R, dist_pc_R, tess_mag_R, new_Tmag,
        rmax=None, r_pct=99, max_points=250_000,
        s_min=0.1, s_max=1, alpha=0.9, tick_count=5,
        cmap="viridis", wd_bp=None, wd_rp=None, wd_g=16.958,
        savepath="plots/wd_dist_dec_Tmag_dual.pdf", titles=("", "")
):
    sns.set_theme(style="whitegrid", context="talk",
                  rc={'font.family': 'serif', 'font.serif': ['DejaVu Serif'], 'font.size': 8,
                      'axes.edgecolor': '0.2', 'axes.labelcolor': '0.', 'xtick.color': '0.', 'ytick.color': '0.',
                      'axes.facecolor': '1', 'grid.color': '0.9'})

    # Convert/clean arrays
    dec_deg_L = _to_numpy(dec_deg_L);
    dist_pc_L = _to_numpy(dist_pc_L);
    tess_mag_L = _to_numpy(tess_mag_L)
    dec_deg_R = _to_numpy(dec_deg_R);
    dist_pc_R = _to_numpy(dist_pc_R);
    tess_mag_R = _to_numpy(tess_mag_R)

    # Per-panel finite masks (require finite dec, dist, T)
    mL = _finite_any(dec_deg_L) & _finite_pos(dist_pc_L) & _finite_any(tess_mag_L)
    mR = _finite_any(dec_deg_R) & _finite_pos(dist_pc_R) & _finite_any(tess_mag_R)

    dec_deg_L, dist_pc_L, tess_mag_L = dec_deg_L[mL], dist_pc_L[mL], tess_mag_L[mL]
    dec_deg_R, dist_pc_R, tess_mag_R = dec_deg_R[mR], dist_pc_R[mR], tess_mag_R[mR]

    # Shared rmax from finite, positive distances across both panels
    if rmax is None:
        all_dist = np.concatenate([dist_pc_L, dist_pc_R])
        if all_dist.size == 0:
            raise ValueError("No finite, positive distances found for either panel.")
        rmax = np.nanpercentile(all_dist, r_pct)
        if not np.isfinite(rmax) or rmax <= 0:
            rmax = float(np.nanmax(all_dist))  # fallback
        if not np.isfinite(rmax) or rmax <= 0:
            raise ValueError("Could not determine a valid rmax (all distances invalid).")

    # Clip radii
    rL = np.clip(dist_pc_L, 0, rmax)
    rR = np.clip(dist_pc_R, 0, rmax)

    # Angles
    thetaL = np.deg2rad(dec_deg_L)
    thetaR = np.deg2rad(dec_deg_R)

    # Shared color normalization using finite T only
    all_T = np.concatenate([tess_mag_L, tess_mag_R])
    all_T = all_T[np.isfinite(all_T)]
    if all_T.size == 0:
        Tmin, Tmax = 10.0, 20.0  # sane default
    else:
        Tmin, Tmax = float(np.nanmin(all_T)), float(np.nanmax(all_T))
        if not np.isfinite(Tmin) or not np.isfinite(Tmax) or Tmin == Tmax:
            Tmin, Tmax = 10.0, 20.0
    norm = Normalize(vmin=Tmin, vmax=Tmax)

    def sizes_from_T(T):
        T = np.clip(T, Tmin, Tmax)
        return np.interp(T, (Tmin, Tmax), (s_max, s_min))

    sL = sizes_from_T(tess_mag_L)
    sR = sizes_from_T(tess_mag_R)

    fig = plt.figure(figsize=(12, 5))

    ax1 = plt.subplot2grid((1, 12), (0, 0), colspan=4, projection="polar")
    ax2 = plt.subplot2grid((1, 12), (0, 3), colspan=4, projection="polar")
    ax3 = plt.subplot2grid((1, 12), (0, 8), colspan=4)  # normal axis

    plt.subplots_adjust(wspace=-0.3)

    def setup_ax(ax, rotate="none"):
        ax.set_thetalim(-np.pi / 2, np.pi / 2)
        if rotate == "left":  # rotate counterclockwise 90
            ax.set_theta_zero_location("W")  # put 0° at east
            ax.set_theta_direction(-1)  # CCW
        elif rotate == "right":  # rotate clockwise 90
            ax.set_theta_zero_location("E")  # put 0° at west
            ax.set_theta_direction(1)  # CW
        else:  # default (north up, CCW)
            ax.set_theta_zero_location("N")
            ax.set_theta_direction(1)

        dec_ticks = [90, 60, 30, 0, -30, -60, -90]
        ax.set_thetagrids(dec_ticks, labels=[f"{d}°" for d in dec_ticks])
        ax.set_rlim(0, float(rmax))
        rticks = np.linspace(0, float(rmax), tick_count)
        ax.set_rticks(rticks)
        ax.set_rlabel_position(0)
        ax.tick_params(labelsize=10)
        ax.grid(color="0.7", linewidth=0.8)
        ax.spines['polar'].set_color("0.4")
        # ax.text(np.deg2rad(102), 0.5 * float(rmax), "Distance (pc)",
        #         ha="center", va="center", fontsize=10)
        ax.text(np.deg2rad(0), 1.3 * float(rmax), "Declination",
                ha="center", va="center", rotation=90, fontsize=10)

    # LEFT scatter
    setup_ax(ax1, rotate="left")
    ax1.text(np.deg2rad(102), 0.5 * float(rmax), "Distance (pc)",
             ha="center", va="center", rotation=90, fontsize=10)
    # Rotate radial tick labels
    for label in ax1.get_yticklabels():
        label.set_rotation(90)  # or any angle
    orderL = np.argsort(tess_mag_L)[::-1]  # reverse: dimmer (large mag) first
    scL = ax1.scatter(thetaL[orderL], rL[orderL], c=tess_mag_L[orderL],
                      cmap=cmap, norm=norm, s=sL[orderL], alpha=alpha,
                      linewidths=0, rasterized=True)

    # RIGHT scatter
    setup_ax(ax2, rotate="right")
    ax2.set_yticklabels([])
    orderR = np.argsort(tess_mag_R)[::-1]  # same for right panel
    scR = ax2.scatter(thetaR[orderR], rR[orderR], c=tess_mag_R[orderR],
                      cmap=cmap, norm=norm, s=sR[orderR], alpha=alpha,
                      linewidths=0, rasterized=True)

    # WD 1856+534 on both panels
    wd_dec = 53 + 30 / 60 + 33.3 / 3600
    wd_dist = 1000 / 40.3983
    if (wd_bp is not None) and (wd_rp is not None):
        wd_dif = float(wd_bp - wd_rp)
        wd_t = wd_g - 0.00522555 * wd_dif ** 3 + 0.0891337 * wd_dif ** 2 - 0.633923 * wd_dif + 0.0324473
        if not np.isfinite(wd_t): wd_t = wd_g - 0.430
    else:
        wd_t = wd_g - 0.430
    wd_t = float(np.clip(wd_t, Tmin, Tmax))
    wd_theta = np.deg2rad(wd_dec)

    for ax in (ax1, ax2):
        ax.scatter(wd_theta, wd_dist, c=[wd_t], cmap=scR.cmap, norm=norm,
                   s=30, edgecolor="w", linewidth=0.8, marker="o", alpha=0.95, zorder=5)
        ax.text(
            0.15, 220, "WD 1856",
            ha="center", va="center", fontsize=10,
            color="black",
            path_effects=[pe.withStroke(linewidth=2, foreground="white")]
        )

    # Titles
    if titles and len(titles) == 2:
        if titles[0]: ax1.set_title(titles[0], fontsize=12, pad=10)
        if titles[1]: ax2.set_title(titles[1], fontsize=12, pad=10)

    # Shared horizontal colorbar
    posL = ax1.get_position();
    posR = ax2.get_position()
    cbar_ax = fig.add_axes([posL.x0 + 0.075, posL.y0 - 0.1, posR.x0 + posR.width - posL.x0 - 0.15, 0.02])
    cb = fig.colorbar(scR, cax=cbar_ax, orientation="horizontal")
    # keep your tick preference but clamp inside [Tmin, Tmax]
    desired = np.array([10, 12, 14, 16, 18, 20], float)
    ticks = desired[(desired >= Tmin) & (desired <= Tmax)]
    if ticks.size < 2:
        ticks = np.linspace(Tmin, Tmax, 5)
    cb.set_ticks(ticks)
    cb.set_label("TESS magnitude (T)", fontsize=10)
    cb.ax.tick_params(labelsize=10)
    # Clamp colorbar to [14, 20]
    scR.set_clim(14, Tmax)
    cb = fig.colorbar(scR, cax=cbar_ax, orientation="horizontal")

    # Define ticks including the left boundary
    ticks = np.arange(14, 21, 2)
    cb.set_ticks(ticks)

    # Replace the first tick label with "<14"
    tick_labels = [f"{t:.0f}" for t in ticks]
    tick_labels[0] = "<14"
    cb.set_ticklabels(tick_labels)

    cb.set_label("TESS magnitude (T)", fontsize=10)
    cb.ax.tick_params(labelsize=10)

    # histogram
    print(len(tess_mag_R))
    print(len(tess_mag_L))
    print(len(new_Tmag))
    ax3.hist(tess_mag_R, bins=np.linspace(10, 22, 13), histtype='step', color="0.9", edgecolor="C0", linewidth=2, alpha=0.9,
             label='Gaia WD', zorder=1)
    ax3.hist(tess_mag_L, bins=np.linspace(10, 22, 13), histtype='step', color="0.9", edgecolor="C1", linewidth=2, alpha=0.9,
             label='TOI WD', zorder=2)
    ax3.hist(new_Tmag, bins=np.linspace(10, 22, 13), histtype='step', color="0.9", edgecolor="k", linewidth=2, alpha=0.9,
             label='Unexplored WD', zorder=3)
    ax3.legend(loc="upper left", fontsize=10)
    ax3.tick_params(labelsize=10)
    ax3.set_xlabel("TESS magnitude (T)", fontsize=10)
    ax3.set_ylabel("Number of Stars", fontsize=10)
    ax3.set_yscale("log")
    ax3.set_xlim(12,22)
    plt.savefig(savepath, dpi=300, bbox_inches="tight")
    return fig, (ax1, ax2, ax3)


if __name__ == "__main__":
    # --- Right panel: Gaia WD dataset ---
    path = "/Users/tehan/Downloads/GaiaEDR3_WD_main.fits"
    ra, dec, plx, d_pc, Tmag, designation = load_wd_data(path, t_limit=100)
    print(len(Tmag))
    # --- Left panel: TIC dataset ---
    tic_path = "/Users/tehan/PycharmProjects/TWIRL/Vanderburg_TOI/vanderburg_TOIs.csv"
    tic = Table.read(tic_path)

    # Extract arrays
    ra_L = tic["ra"].data
    dec_L = tic["dec"].data
    d_pc_L = tic["distance_pc"].data
    designation_L = tic["dr3_source_id"].data
    phot_g_mean_mag = tic["phot_g_mean_mag"].data
    phot_bp_mean_mag = tic["phot_bp_mean_mag"].data
    phot_rp_mean_mag = tic["phot_rp_mean_mag"].data
    Tmag_L = np.zeros(np.shape(tic))
    for i in range(len(tic)):
        if (phot_bp_mean_mag[i] is not None) and (phot_rp_mean_mag[i] is not None):
            wd_dif = float(phot_bp_mean_mag[i] - phot_rp_mean_mag[i])
            Tmag_L[i] = phot_g_mean_mag[
                            i] - 0.00522555 * wd_dif ** 3 + 0.0891337 * wd_dif ** 2 - 0.633923 * wd_dif + 0.0324473
            if not np.isfinite(Tmag_L[i]): Tmag_L[i] = phot_g_mean_mag[i] - 0.430
        else:
            Tmag_L[i] = phot_g_mean_mag[i] - 0.430
    # new_designation
    mask = ~np.isin(designation, designation_L)
    idx = np.where(mask)[0]  # indices of those elements

    # Single-panel check (optional)
    # plot_2d_polar(dec, d_pc, Tmag, wd_bp=17.5032, wd_rp=16.2780)

    # Dual-panel plot
    fig, axes = plot_2d_polar_two(
        dec_L, d_pc_L, Tmag_L,  # Left panel
        dec, d_pc, Tmag,  # Right panel
        Tmag[idx],
        rmax=800,
        titles=("2-min TOI WD", "200-s FFI WD"),
        wd_bp=17.5032, wd_rp=16.2780
    )
