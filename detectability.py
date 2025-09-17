import os
import glob
import numpy as np
from astropy import units as u
from astropy.constants import G, M_sun, R_earth
from astropy.constants import M_sun, R_earth as R_earth_const
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.ticker import PercentFormatter

# --- Styling: apply BEFORE creating any figures ---
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update({
    # Fonts
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "STIX Two Text", "Times New Roman"],
    "mathtext.fontset": "dejavuserif",

    # Global font sizes
    "font.size": 12,          # base size
    "axes.titlesize": 12,     # panel titles
    "axes.labelsize": 12,     # x/y labels
    "xtick.labelsize": 12,    # tick labels
    "ytick.labelsize": 12,
    "legend.fontsize": 11,
    "figure.titlesize": 14,   # overall fig title if you add one

    # Colors / grid
    "axes.edgecolor": "0.2",
    "axes.labelcolor": "0.",
    "xtick.color": "0.",
    "ytick.color": "0.",
    "axes.facecolor": "1",
    "grid.color": "0.5",

    # Tick directions
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.bottom": True,
    "ytick.left": True,

})
# White dwarf: 0.6 Msun, 1.0 R_earth
M_star = 0.6 * M_sun
R_star = 1.0 * R_earth

def _to_scalar(x):
    """Return a Python float for scalar-like inputs (Quantity | 0-D ndarray | float)."""
    if hasattr(x, "value"):
        x = x.value
    a = np.asarray(x)
    return a.item() if a.shape == () else a  # returns float for scalar, ndarray otherwise

def _to_scalar_float(x):
    """Strict float scalar (or np.nan) — use when you KNOW x should be scalar."""
    s = _to_scalar(x)
    if isinstance(s, np.ndarray):
        # fall back conservatively; you can also raise if you prefer
        return s.ravel()[0]
    return float(s)

def _qty(x, assumed_unit):
    q = u.Quantity(x)
    return q if q.unit != u.dimensionless_unscaled else (q.value * assumed_unit)

def _uniform_depth(k, b):
    """
    Fractional flux drop for a uniform stellar disk (R*=1), planet radius k, separation b.
    Returns values in [0, 1]. Vectorized.
    """
    k = np.asarray(k, dtype=float)
    d = np.asarray(b, dtype=float)

    depth = np.zeros_like(d, dtype=float)

    # Full overlap of the smaller disk by the larger:
    full = d <= np.abs(1.0 - k)
    # If k <= 1: planet fully inside stellar disk → area = π k^2 → depth = k^2
    # If k > 1: planet fully covers the star → depth = 1
    depth[full] = np.where(k <= 1.0, k**2, 1.0)

    # Partial (grazing) overlap
    part = (d < 1.0 + k) & (~full)
    if np.any(part):
        dp = d[part]
        c1 = np.clip((dp**2 + 1.0 - k**2) / (2.0 * dp), -1.0, 1.0)
        c2 = np.clip((dp**2 + k**2 - 1.0) / (2.0 * dp * k), -1.0, 1.0)
        term1 = np.arccos(c1)                    # r1 = 1
        term2 = np.arccos(c2)                    # r2 = k
        rad = (-dp + 1.0 + k) * (dp + 1.0 - k) * (dp - 1.0 + k) * (dp + 1.0 + k)
        lens = term1 + k**2 * term2 - 0.5 * np.sqrt(np.clip(rad, 0.0, None))
        depth[part] = lens / np.pi               # normalize by π R_*^2 (R_*=1)

    # No-overlap stays 0. Extra safety clip:
    return np.clip(depth, 0.0, 1.0)

def n_pts_in_transit(
    P, k, inc, W=27*u.day, cadence=200*u.s,
    M_star=None, R_star=None, G_const=G, return_depth=True,
    weighted=True, nsamp=2048
):
    if M_star is None:
        from detectability import M_star as M_star_global
        M_star = M_star_global
    if R_star is None:
        from detectability import R_star as R_star_global
        R_star = R_star_global

    P = _qty(P, u.day).to(u.s)
    i = _qty(inc, u.deg).to(u.rad)
    W = _qty(W, u.day)
    cadence = _qty(cadence, u.s)

    # Broadcast
    P_b, k_b, i_b = np.broadcast_arrays(np.atleast_1d(P.value),
                                        np.atleast_1d(k),
                                        np.atleast_1d(i.value))
    out_n   = np.zeros(P_b.shape, float)
    out_dep = np.zeros(P_b.shape, float)

    for idx in np.ndindex(P_b.shape):
        P_s = (P_b[idx]) * u.s
        k_s = float(k_b[idx])
        i_s = (i_b[idx]) * u.rad

        a = (G_const * M_star * P_s**2 / (4*np.pi**2))**(1/3)
        b = (a / R_star) * np.cos(i_s.value)
        sin_i = max(np.sin(i_s.value), 1e-12)
        inside = (1.0 + k_s)**2 - b**2

        if (inside <= 0.0) or (b > (1.0 + k_s)):
            out_n[idx] = 0.0; out_dep[idx] = 0.0
            continue

        chord = np.sqrt(max(0.0, inside))
        T14 = (P_s/np.pi) * (R_star / a) * (chord / sin_i)
        Ntr = (W.to(u.s) / P_s).value
        depth_max = _uniform_depth(k_s, b.value)
        out_dep[idx] = depth_max

        if not weighted:
            T_tot = (W.to(u.s) / np.pi) * (R_star / a) * (chord / sin_i)
            out_n[idx] = (T_tot / cadence).to(u.dimensionless_unscaled).value
            continue

        if depth_max <= 0.0:
            out_n[idx] = 0.0
            continue

        # Depth-weighted effective duration per transit
        ns = int(max(32, nsamp))
        t = np.linspace(-0.5, 0.5, ns) * T14.value
        v_perp = (2.0 * chord) / T14.value
        x = v_perp * t
        d = np.sqrt(b.value**2 + x**2)
        depth_t = _uniform_depth(k_s, d)
        w_t = np.clip(depth_t / depth_max, 0.0, 1.0)
        w_avg = np.trapz(w_t, t) / T14.value

        T_eff_total = Ntr * (w_avg * T14.value)
        out_n[idx] = T_eff_total / cadence.value

    n_points = out_n if out_n.shape != () else float(out_n)
    depth    = out_dep if out_dep.shape != () else float(out_dep)
    if not return_depth:
        return n_points
    return n_points, depth


def calculate_snr(Tmag, n_points, depth, noisemodel="/Users/tehan/Documents/TGLC/noisemodel.dat"):
    noise_model = np.loadtxt(noisemodel).T
    mags, noise_hr = noise_model[0], noise_model[1]
    noise_interp = np.interp(Tmag, mags, noise_hr)

    # If noise_hr is 1-hr precision, the 200 s scaling is sqrt(3600/200) = sqrt(18)
    noise_200s = noise_interp * np.sqrt(18)  # change to sqrt(9) only if your file is 30 min

    n_val = _to_scalar(n_points)  # supports Quantity/0-D/float
    d_val = _to_scalar(depth)
    return np.sqrt(n_val) * d_val / noise_200s

# --- assumes n_pts_in_transit(P, k, inc) and calculate_snr(Tmag, n_pts, depth) exist ---

def detection_limit_curve(
    Tmag,
    periods=np.linspace(0.2, 20.0, 400, endpoint=True),
    incs=np.linspace(89.5, 90.0, 150, endpoint=True),
    snr_threshold=7.0,
    k=1.0,
    right_censor="nan"
):
    snr_grid  = np.zeros((len(periods), len(incs)), dtype=float)
    dept_grid = np.zeros_like(snr_grid, dtype=float)

    for i, P in enumerate(periods):
        for j, inc in enumerate(incs):
            n_pts, depth = n_pts_in_transit(P, k, inc)  # expect numerics
            n_pts = _to_scalar_float(n_pts)
            depth = _to_scalar_float(depth)
            snr_grid[i, j]  = _to_scalar_float(calculate_snr(Tmag, n_pts, depth))
            dept_grid[i, j] = depth

    P_list, D_list = [], []
    for j, inc in enumerate(incs):
        s = snr_grid[:, j]
        i_cross = None
        for i in range(len(periods)-2, -1, -1):
            if (s[i] >= snr_threshold) and (s[i+1] < snr_threshold):
                i_cross = i
                break

        if i_cross is None:
            if np.all(s < snr_threshold):
                P_list.append(np.nan); D_list.append(np.nan)
            else:
                if right_censor == "nan":
                    P_list.append(np.nan); D_list.append(np.nan)
                else:
                    P0, P1, S0, S1 = periods[-2], periods[-1], s[-2], s[-1]
                    if np.isclose(S0, S1):
                        P_star = np.nan; D_star = np.nan
                    else:
                        P_star = P0 + (snr_threshold - S0)/(S1 - S0)*(P1 - P0)
                        _, D_star = n_pts_in_transit(P_star, k, inc)
                        D_star = _to_scalar_float(D_star)
                    P_list.append(_to_scalar_float(P_star))
                    D_list.append(D_star)
        else:
            P0, P1, S0, S1 = periods[i_cross], periods[i_cross+1], s[i_cross], s[i_cross+1]
            if np.isclose(S0, S1):
                P_star = P0
            else:
                P_star = P0 + (snr_threshold - S0)/(S1 - S0)*(P1 - P0)
            _, D_star = n_pts_in_transit(P_star, k, inc)
            P_list.append(_to_scalar_float(P_star))
            D_list.append(_to_scalar_float(D_star))

    P_curve = np.asarray(P_list, dtype=float)
    D_curve = np.asarray(D_list, dtype=float)

    good = np.isfinite(P_curve) & np.isfinite(D_curve)
    o = np.argsort(P_curve[good])
    return {
        "Tmag": float(Tmag),
        "P_curve": P_curve[good][o],
        "D_curve": D_curve[good][o],
        "k": float(k),
        "snr_threshold": float(snr_threshold),
        "periods_grid": np.asarray(periods, dtype=float),
        "incs_grid": np.asarray(incs, dtype=float),
        "right_censor": right_censor,
    }

def save_curve(curve, outdir="/Users/tehan/PycharmProjects/TWIRL/plots"):
    os.makedirs(outdir, exist_ok=True)
    k = curve["k"]; snr = curve["snr_threshold"]; Tmag = curve["Tmag"]
    fname = f"detcurve_Tmag{int(Tmag)}_SNR{snr:g}_k{k:.3f}.npz"
    path = os.path.join(outdir, fname)
    np.savez(
        path,
        **{k:v for k,v in curve.items()}
    )
    return path

def batch_compute_save(
    Tmags=range(15,21),
    periods=np.linspace(0.2, 20.0, 400, endpoint=True),
    incs=np.linspace(89.5, 90.0, 150, endpoint=True),
    snr_threshold=7.0,
    k=1.0,
    outdir="/Users/tehan/PycharmProjects/TWIRL/plots",
    right_censor="nan"
):
    paths=[]
    for Tmag in Tmags:
        curve = detection_limit_curve(
            Tmag=Tmag, periods=periods, incs=incs,
            snr_threshold=snr_threshold, k=k, right_censor=right_censor
        )
        p = save_curve(curve, outdir=outdir); paths.append(p)
    return paths

def _smooth_curve(P, D, enable=True, target_width_days=0.3, min_win=5):
    """
    Length-preserving moving-average smoother.
    - Keeps length == len(P) always.
    - Skips smoothing if series is too short or sampling already coarse.
    """
    P = np.asarray(P).ravel()
    D = np.asarray(D).ravel()
    n = len(P)
    if (not enable) or n < min_win + 2:
        return P, D

    # Require strictly increasing P for a meaningful dP; otherwise sort
    if not np.all(np.diff(P) > 0):
        o = np.argsort(P)
        P, D = P[o], D[o]

    dP = np.median(np.diff(P))
    if not np.isfinite(dP) or dP <= 0:
        return P, D

    # If already coarse compared to target width, skip smoothing
    if dP >= target_width_days / 2.0:
        return P, D

    # Choose an odd window size; cap by data length
    win = int(round(target_width_days / max(dP, 1e-12)))
    win = max(win, min_win)
    if win % 2 == 0:
        win += 1
    if win >= n:
        # too wide; skip smoothing
        return P, D

    # Moving average with symmetric padding to preserve length
    pad = win // 2
    Dpad = np.pad(D, pad_width=pad, mode="edge")
    kernel = np.ones(win, dtype=float) / win
    Dsmooth = np.convolve(Dpad, kernel, mode="valid")  # length n
    return P, Dsmooth

def plot_from_saved(
    ax,
    pattern,
    title="Detection-limit curves vs TESS magnitude",
    smooth=True,
    smooth_width_days=0.5,
    fill_alpha=0.12,
    remove_flat=True,
    eps=1e-3,
    roche=1,
):
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files match: {pattern}")

    footer_parts = ''

    for f in files:
        data = np.load(f, allow_pickle=True)
        P = np.asarray(data["P_curve"]).ravel()
        D = np.asarray(data["D_curve"]).ravel()
        Tmag = float(data["Tmag"])
        k = float(data["k"])
        snr = float(data["snr_threshold"])
        incs = np.asarray(data["incs_grid"]).ravel()
        periods = np.asarray(data["periods_grid"]).ravel()
        censor = str(data["right_censor"])

        footer_parts = rf"$R_p/R_*={k:.0f}$, Required transit SNR={snr:g}, Sector=27 days, inc=[{incs.min():.0f},{incs.max():.0f}]°"

        if len(P) == 0:
            continue

        good = np.isfinite(P) & np.isfinite(D)
        P, D = P[good], D[good]
        if len(P) < 2:
            continue

        Pp, Dp = _smooth_curve(P, D, enable=smooth, target_width_days=smooth_width_days)

        n = min(len(Pp), len(Dp))
        Pp, Dp = Pp[:n], Dp[:n]

        cap = min(k**2, 1.0)

        if remove_flat:
            mask = Dp < (cap - eps)
            if len(mask) != len(Pp):
                m = min(len(mask), len(Pp))
                mask = mask[:m]; Pp = Pp[:m]; Dp = Dp[:m]
            if not np.any(mask):
                continue
            Pp, Dp = Pp[mask], Dp[mask]
        else:
            Dp = np.minimum(Dp, cap - eps)

        if len(Pp) < 2:
            continue

        line, = ax.plot(Pp, Dp, lw=2, label=f"Tmag={int(Tmag)}")
        color = line.get_color()
        ax.fill_between(Pp, Dp, cap, alpha=fill_alpha, linewidth=0, color=color)
    ax.vlines(roche, 0, 1, linestyles="dotted", color="k")
    ax.text(roche*1.1, 0.5, 'Roche-limit', fontsize=9, rotation=90,va='center', )
    ax.set_xlabel("Orbital period [days]")
    ax.set_ylabel("Transit depth")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=True, ncol=2, loc="upper right", fontsize=10)
    ax.set_xlim(0.2, 15.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xscale("log")
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1.0, decimals=0))

    # plot WD 1856 b
    if title[:3] == 'Jup':
        mag_1856 = 16.338
        p_1856 = 1.4079405
        d_1856 = 0.5665
        k_1856 = 7.28
        i_1856 = 88.778
        # mag_1856 = 17
        # p_1856 = 2
        # d_1856 = 0.2
        # k_1856 = 10
        # i_1856 = 89.02
        snr_1856, param = estimate_snr_single(mag_1856, p_1856, k_1856, i_1856, M_star=0.518*M_sun, R_star=1.429*R_earth)
        print(snr_1856, param)
        ax.scatter(p_1856, param['depth'], marker="o", s=50, c="g", label="WD 1856 b")
    footer_txt = footer_parts
    return footer_txt

def estimate_snr_single(
        Tmag, P_days, k, inc_deg,
        W_days=27.0, cadence_s=200.0,
        M_star=M_star, R_star=R_star,
        noisemodel="/Users/tehan/Documents/TGLC/noisemodel.dat"
):
    """
    Compute transit SNR using your n_pts_in_transit() and calculate_snr().

    Parameters
    ----------
    Tmag : float
        TESS magnitude
    P_days : float
        Orbital period in days
    k : float
        Rp/R*
    inc_deg : float
        Inclination in degrees
    """
    n_pts, depth = n_pts_in_transit(
        P_days, k, inc_deg,
        W=W_days * u.day, cadence=cadence_s * u.s,
        M_star=M_star, R_star=R_star,
        return_depth=True
    )
    snr = calculate_snr(Tmag, n_pts, depth, noisemodel=noisemodel)
    return _to_scalar_float(snr), {
        "n_pts": _to_scalar_float(n_pts),
        "depth": _to_scalar_float(depth),
        "k": float(k),
        "inc": float(inc_deg),
    }


if __name__ == "__main__":
    mag_1856 = 16.338
    p_1856 = 1.4079405
    d_1856 = 0.5665
    k_1856 = 7.28
    i_1856 = 88.735
    snr_1856, param = estimate_snr_single(mag_1856, p_1856, k_1856, i_1856, M_star=0.518 * M_sun,
                                          R_star=1.429 * R_earth)
    print(snr_1856, param)
    # 1) Compute & save curves for Tmag=15..20
    k = 1.0
    paths = batch_compute_save(
        Tmags=range(14, 20),
        periods=np.linspace(0.2, 20.0, 400),
        incs=np.linspace(88, 90.0, 300),
        snr_threshold=7.0,
        k=k,
        outdir="/Users/tehan/PycharmProjects/TWIRL/plots/curves/",
        right_censor="nan"
    )

    k = 10.0
    paths = batch_compute_save(
        Tmags=range(14, 20),
        periods=np.linspace(0.2, 20.0, 400),
        incs=np.linspace(86, 90.0, 300),
        snr_threshold=7.0,
        k=k,
        outdir="/Users/tehan/PycharmProjects/TWIRL/plots/curves/",
        right_censor="nan"
    )

    # 2) Plot from the saved files (auto-discovers via glob pattern)
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    # Left panel: k=1, smooth=False
    k = 1.0
    footer1 = plot_from_saved(
        ax=axes[0],
        pattern=f"/Users/tehan/PycharmProjects/TWIRL/plots/curves/detcurve_Tmag*_SNR*_k{k:.3f}.npz",
        title="Earth-like WD planet detectability map",
        smooth=False,
        roche=0.218,
    )

    # Right panel: k=10, smooth=True
    k = 10.0
    footer2 = plot_from_saved(
        ax=axes[1],
        pattern=f"/Users/tehan/PycharmProjects/TWIRL/plots/curves/detcurve_Tmag*_SNR*_k{k:.3f}.npz",
        title="Jupiter-like WD planet detectability map",
        smooth=False,
        roche=0.44,
    )

    # Put a shared footer text
    fig.text(
        0.05, 0.01,
        footer1,
        ha="left", va="bottom", fontsize=10
    )
    fig.text(
        0.535, 0.01,
        footer2,
        ha="left", va="bottom", fontsize=10
    )
    outpdf = "/Users/tehan/PycharmProjects/TWIRL/plots/depth_vs_period_multiTmag_combined.pdf"
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(outpdf, dpi=300)
    plt.close(fig)
    print(f"Saved: {outpdf}")