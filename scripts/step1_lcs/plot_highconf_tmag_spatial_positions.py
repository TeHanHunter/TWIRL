#!/usr/bin/env python

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.plotting.style import apply_twirl_style


DEFAULT_CATALOG = Path("data_local/catalogs/GaiaEDR3_WD_main.fits")
DEFAULT_OUTDIR = Path("reports/step1_lcs")


class QuantileNormalize(Normalize):
    """Map values through an empirical CDF so dense ranges use more of the colormap."""

    def __init__(self, values: np.ndarray, n_knots: int = 1024):
        finite_values = np.asarray(values, dtype=float)
        finite_values = finite_values[np.isfinite(finite_values)]
        if finite_values.size == 0:
            raise ValueError("QuantileNormalize requires at least one finite value.")

        probs = np.linspace(0.0, 1.0, n_knots)
        self._value_knots = np.quantile(finite_values, probs)
        self._cdf_knots = probs
        super().__init__(
            vmin=float(self._value_knots[0]),
            vmax=float(self._value_knots[-1]),
            clip=True,
        )

    def __call__(self, value, clip=None):
        result, is_scalar = self.process_value(value)
        data = np.asarray(result.data, dtype=float)
        mapped = np.interp(data, self._value_knots, self._cdf_knots)
        mapped = np.ma.masked_array(mapped, mask=result.mask, copy=False)
        if is_scalar:
            mapped = mapped[0]
        return mapped

    def inverse(self, value):
        return np.interp(value, self._cdf_knots, self._value_knots)


class PwdLogNormalize(Normalize):
    """Log-stretch the contamination fraction so values near P_wd ~= 1 remain visible."""

    def __init__(self, vmin: float = 0.75, vmax: float = 1.0, epsilon: float = 1e-4):
        self.epsilon = epsilon
        self.max_contamination = max(1.0 - vmin, epsilon)
        self.log_span = np.log10(self.max_contamination) - np.log10(self.epsilon)
        super().__init__(vmin=vmin, vmax=vmax, clip=True)

    def __call__(self, value, clip=None):
        result, is_scalar = self.process_value(value)
        data = np.asarray(result.data, dtype=float)
        data = np.clip(data, self.vmin, self.vmax)
        contamination = np.clip(1.0 - data, self.epsilon, self.max_contamination)
        mapped = (np.log10(self.max_contamination) - np.log10(contamination)) / self.log_span
        mapped = np.ma.masked_array(mapped, mask=result.mask, copy=False)
        if is_scalar:
            mapped = mapped[0]
        return mapped

    def inverse(self, value):
        contamination = 10 ** (np.log10(self.max_contamination) - value * self.log_span)
        return 1.0 - contamination


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot the high-confidence white dwarf sample (P_wd > 0.75) with a "
            "TESS-magnitude color scale in a polar distance view and an Aitoff sky map."
        )
    )
    parser.add_argument(
        "--catalog",
        type=Path,
        default=DEFAULT_CATALOG,
        help="Path to the local Gaia EDR3 WD FITS catalog.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help="Directory for generated figure products.",
    )
    parser.add_argument(
        "--distance-clip-pct",
        type=float,
        default=99.5,
        help="Upper percentile used to clip distances in the polar panel.",
    )
    parser.add_argument(
        "--template",
        choices=("column", "full_page"),
        default="full_page",
        help="Named TWIRL plotting template from doc/plotting_style.md.",
    )
    parser.add_argument(
        "--dense-marker-scale",
        type=float,
        default=0.55,
        help="Scale factor applied to the default dense scatter marker size.",
    )
    parser.add_argument(
        "--output-stem",
        type=str,
        default="wd_highconf_tmag_spatial_positions",
        help="Base filename for output figure products.",
    )
    parser.add_argument(
        "--color-scale",
        choices=("quantile", "linear"),
        default="quantile",
        help="Color normalization for TESS magnitude. 'quantile' spreads dense magnitude ranges across the colormap.",
    )
    return parser.parse_args()


def load_highconf_sample(path: Path) -> dict[str, np.ndarray]:
    table = Table.read(path, hdu="Joined", unit_parse_strict="silent")
    table.convert_bytestring_to_unicode()

    ra = np.asarray(table["ra"], dtype=float)
    dec = np.asarray(table["dec"], dtype=float)
    parallax = np.asarray(table["parallax"], dtype=float)
    g_mag = np.asarray(table["phot_g_mean_mag"], dtype=float)
    bp_mag = np.asarray(table["phot_bp_mean_mag"], dtype=float)
    rp_mag = np.asarray(table["phot_rp_mean_mag"], dtype=float)
    pwd = np.asarray(table["Pwd"], dtype=float)

    # Han & Brandt (2023), Eq. (1), adopting the Stassun et al. (2019) Gaia->TESS conversion.
    color = bp_mag - rp_mag
    tmag = (
        g_mag
        - 0.00522555 * color**3
        + 0.0891337 * color**2
        - 0.633923 * color
        + 0.0324473
    )
    bad_tmag = ~np.isfinite(tmag)
    tmag[bad_tmag] = g_mag[bad_tmag] - 0.430

    mask = (
        np.isfinite(ra)
        & np.isfinite(dec)
        & np.isfinite(parallax)
        & (parallax > 0)
        & np.isfinite(g_mag)
        & np.isfinite(tmag)
        & np.isfinite(pwd)
        & (pwd > 0.75)
    )

    distance_pc = 1000.0 / parallax[mask]
    sample = {
        "ra": ra[mask],
        "dec": dec[mask],
        "distance_pc": distance_pc,
        "pwd": pwd[mask],
        "tmag": tmag[mask],
    }

    # Plot faint objects first so the brighter subset is not buried.
    order = np.argsort(sample["tmag"])[::-1]
    for key in sample:
        sample[key] = sample[key][order]
    return sample


def galactic_aitoff_coords(ra_deg: np.ndarray, dec_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    coords = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs").galactic
    lon_wrap = coords.l.wrap_at(180 * u.deg).radian
    return -lon_wrap, coords.b.radian


def wd1856_tmag() -> float:
    wd_g = 16.958
    wd_bp = 17.5032
    wd_rp = 16.2780
    color = wd_bp - wd_rp
    return (
        wd_g
        - 0.00522555 * color**3
        + 0.0891337 * color**2
        - 0.633923 * color
        + 0.0324473
    )


def plot_sample(
    sample: dict[str, np.ndarray],
    outdir: Path,
    distance_clip_pct: float,
    template_name: str,
    dense_marker_scale: float,
    output_stem: str,
    color_scale: str,
) -> tuple[Path, Path]:
    template = apply_twirl_style(template_name)
    marker_size = template["dense_marker_size"] * dense_marker_scale

    distance_pc = sample["distance_pc"]
    dec_deg = sample["dec"]
    ra_deg = sample["ra"]
    pwd = sample["pwd"]
    tmag = sample["tmag"]

    clipped_rmax = float(np.nanpercentile(distance_pc, distance_clip_pct))
    clipped_distance = np.clip(distance_pc, 0, clipped_rmax)
    theta = np.deg2rad(dec_deg)

    # Use a robust display range while preserving the true TESS magnitudes in the data arrays.
    tmag_vmin = max(17.0, float(np.nanpercentile(tmag, 1)))
    tmag_vmax = float(np.nanpercentile(tmag, 99))
    if color_scale == "quantile":
        color_values = np.clip(tmag, tmag_vmin, tmag_vmax)
        tmag_norm = QuantileNormalize(color_values)
    else:
        color_values = np.clip(tmag, tmag_vmin, tmag_vmax)
        tmag_norm = Normalize(vmin=tmag_vmin, vmax=tmag_vmax)
    pwd_norm = PwdLogNormalize(vmin=0.75, vmax=1.0, epsilon=1e-4)
    pwd_cmap = plt.get_cmap("viridis")
    tmag_cmap = plt.get_cmap("viridis_r")

    fig = plt.figure(figsize=template["figsize"])
    grid = fig.add_gridspec(
        1,
        2,
        width_ratios=[1.0, 1.45],
        wspace=max(template["panel_wspace"], 0.12),
    )
    ax = fig.add_subplot(grid[0, 0], projection="polar")
    ax_sky = fig.add_subplot(grid[0, 1], projection="aitoff")

    ax.set_thetalim(-np.pi / 2, np.pi / 2)
    ax.set_theta_zero_location("W")
    ax.set_theta_direction(-1)

    sc = ax.scatter(
        theta,
        clipped_distance,
        c=pwd,
        cmap=pwd_cmap,
        norm=pwd_norm,
        s=marker_size,
        alpha=0.10,
        linewidths=0,
        rasterized=True,
    )

    sky_x, sky_y = galactic_aitoff_coords(ra_deg, dec_deg)
    ax_sky.scatter(
        sky_x,
        sky_y,
        c=color_values,
        cmap=tmag_cmap,
        norm=tmag_norm,
        s=max(0.45, marker_size - 0.15),
        alpha=0.20,
        linewidths=0,
        rasterized=True,
    )

    wd_ra = 18 * 15 + 57 / 60 * 15 + 39.0 / 3600 * 15
    wd_dec = 53 + 30 / 60 + 33.3 / 3600
    wd_dist = 1000.0 / 40.3983
    wd_theta = np.deg2rad(wd_dec)
    wd_tmag = float(np.clip(wd1856_tmag(), tmag_vmin, tmag_vmax))
    wd_pwd = 0.999
    ax.scatter(
        [wd_theta],
        [wd_dist],
        c=[wd_pwd],
        cmap=pwd_cmap,
        norm=pwd_norm,
        s=42,
        marker="o",
        edgecolor="k",
        linewidth=0.9,
        alpha=0.92,
        zorder=5,
    )
    ax.text(
        wd_theta,
        wd_dist + 200,
        "WD 1856",
        ha="right",
        va="center",
        fontsize=template["annotation_size"],
    )

    wd_x, wd_y = galactic_aitoff_coords(np.array([wd_ra]), np.array([wd_dec]))
    ax_sky.scatter(
        wd_x,
        wd_y,
        c=[wd_tmag],
        cmap=tmag_cmap,
        norm=tmag_norm,
        s=42,
        marker="o",
        edgecolor="k",
        linewidth=0.9,
        alpha=0.92,
        zorder=6,
    )
    ax_sky.text(
        wd_x[0] + np.radians(6),
        wd_y[0] + np.radians(2),
        "WD 1856",
        ha="left",
        va="center",
        fontsize=template["annotation_size"],
    )

    dec_ticks = [90, 60, 30, 0, -30, -60, -90]
    ax.set_thetagrids(dec_ticks, labels=[f"{deg}°" for deg in dec_ticks])
    ax.set_rlim(0, clipped_rmax)
    radial_ticks = np.array([500.0, 1000.0, 1500.0])
    radial_ticks = radial_ticks[radial_ticks < clipped_rmax]
    ax.set_rticks(radial_ticks)
    ax.set_yticklabels([f"{tick / 1000.0:g}" for tick in radial_ticks])
    ax.set_rlabel_position(135)
    ax.tick_params(axis="x", labelsize=template["tick_size"], pad=-3)
    ax.tick_params(axis="y", labelsize=template["tick_size"], pad=2)
    ax.grid(color="0.72", linewidth=template["grid_linewidth"])
    ax.spines["polar"].set_color("0.4")
    ax.text(
        0.80,
        0.82,
        "Distance (kpc)",
        ha="center",
        va="center",
        rotation=90,
        fontsize=template["label_size"],
        transform=ax.transAxes,
    )
    ax.text(
        0.10,
        0.52,
        "Dec (degree)",
        ha="center",
        va="center",
        rotation=90,
        fontsize=template["label_size"],
        transform=ax.transAxes,
    )

    ax_sky.set_xlabel("Galactic longitude (degree)", fontsize=template["label_size"], labelpad=10)
    ax_sky.set_ylabel("Galactic latitude (degree)", fontsize=template["label_size"], labelpad=6)
    ax_sky.yaxis.set_label_coords(-0.08, 0.5)
    xticks_deg = np.array([-180, -120, -60, 0, 60, 120])
    ax_sky.set_xticks(np.radians(xticks_deg))
    ax_sky.set_xticklabels([])
    ax_sky.tick_params(axis="both", labelsize=template["tick_size"])
    ax_sky.grid(color="0.85", linewidth=template["grid_linewidth"])
    latitudes = np.linspace(-np.pi / 2, np.pi / 2, 361)
    for x_tick in np.radians(xticks_deg):
        ax_sky.plot(
            np.full_like(latitudes, x_tick),
            latitudes,
            color="black",
            linewidth=0.55,
            linestyle=":",
            alpha=0.75,
            zorder=5.2,
        )
    ax_sky.axhline(0.0, color="black", linewidth=0.75, linestyle="--", alpha=0.95, zorder=5.0)
    longitude_labels = ["", "120°", "60°", "0°", "300°", "240°"]
    for x_tick, label in zip(np.radians(xticks_deg), longitude_labels):
        ax_sky.text(
            x_tick,
            0.0,
            label,
            ha="center",
            va="bottom",
            fontsize=template["tick_size"],
            color="0.25",
            zorder=6,
        )

    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.17, top=0.96)
    left_pos = ax.get_position()
    right_pos = ax_sky.get_position()
    cbar_height = 0.026
    cbar_bottom = 0.138

    cbar_pwd_ax = fig.add_axes([left_pos.x0, cbar_bottom, left_pos.width, cbar_height])
    cbar_pwd = fig.colorbar(
        ScalarMappable(norm=pwd_norm, cmap=pwd_cmap),
        cax=cbar_pwd_ax,
        orientation="horizontal",
    )
    cbar_pwd.set_ticks([0.75, 0.95, 0.99, 1.00])
    cbar_pwd.set_label(r"$P_{\mathrm{wd}}$", fontsize=template["label_size"])
    cbar_pwd.ax.tick_params(labelsize=template["tick_size"])

    cbar_t_ax = fig.add_axes([right_pos.x0, cbar_bottom, right_pos.width, cbar_height])
    cbar = fig.colorbar(
        ScalarMappable(norm=tmag_norm, cmap=tmag_cmap),
        cax=cbar_t_ax,
        orientation="horizontal",
    )
    tick_start = int(np.floor(tmag_vmin))
    tick_stop = int(np.ceil(tmag_vmax))
    cbar.set_ticks(np.arange(tick_start, tick_stop + 1))
    cbar.set_label("TESS magnitude (T)", fontsize=template["label_size"])
    cbar.ax.tick_params(labelsize=template["tick_size"])

    pdf_path = outdir / f"{output_stem}.pdf"
    png_path = outdir / f"{output_stem}.png"
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return pdf_path, png_path


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    sample = load_highconf_sample(args.catalog)
    pdf_path, png_path = plot_sample(
        sample=sample,
        outdir=args.outdir,
        distance_clip_pct=args.distance_clip_pct,
        template_name=args.template,
        dense_marker_scale=args.dense_marker_scale,
        output_stem=args.output_stem,
        color_scale=args.color_scale,
    )

    print(f"Catalog path: {args.catalog}")
    print(f"High-confidence sample size: {sample['ra'].size}")
    print(f"Figure PDF: {pdf_path}")
    print(f"Figure PNG: {png_path}")


if __name__ == "__main__":
    main()
