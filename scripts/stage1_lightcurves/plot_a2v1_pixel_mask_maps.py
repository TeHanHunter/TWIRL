#!/usr/bin/env python3
"""Render representative A2v1 saturated-pixel mask diagnostics by sector.

Each map uses the median of the first requested staged FFIs for one
orbit/camera/CCD. Red contours show the corresponding TGLC bad-pixel threshold
proxy: saturated, low-flux, and nonfinite pixels plus their four neighbors.
The companion grid shows the cadence-median central ePSF coefficient for each
cutout, normalized by the CCD median. Red cell outlines identify an on-disk
A2v1 ePSF rather than a legacy empty-mask ePSF symlink.
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass, asdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LogNorm
from matplotlib.patches import Rectangle
import numpy as np
from astropy.io import fits

from twirl.plotting.style import apply_twirl_style


DEFAULT_A2V1_ROOT = Path("/pdo/users/tehan/tglc-gpu-production-A2v1")
DEFAULT_OUTPUT_DIR = Path("reports/stage1_lightcurves/a2v1_pixel_mask_maps")
DEFAULT_SECTOR_ORBITS = {
    56: 119,
    57: 121,
    58: 123,
    59: 125,
    60: 127,
    61: 129,
    62: 131,
    63: 133,
}
SOURCE_RE = re.compile(r"^source_(?P<x>\d+)_(?P<y>\d+)\.pkl$")
BACKGROUND_PARAMETER_COUNT = 6


@dataclass(frozen=True)
class MapSummary:
    sector: int
    orbit: int
    camera: int
    ccd: int
    n_cadences: int
    n_source_tiles: int
    n_local_epsf_tiles: int
    n_proxy_bad_pixels: int
    epsf_center_index: int
    epsf_relative_scale_median: float
    png: str
    pdf: str


def parse_sector_orbit(value: str) -> tuple[int, int]:
    try:
        sector_text, orbit_text = value.split(":", maxsplit=1)
        sector = int(sector_text)
        orbit = int(orbit_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "sector orbit must have format '<sector>:<first-orbit>'"
        ) from exc
    if sector <= 0 or orbit <= 0:
        raise argparse.ArgumentTypeError("sector and orbit must be positive")
    return sector, orbit


def parse_source_grid(path: Path) -> tuple[int, int]:
    match = SOURCE_RE.match(path.name)
    if match is None:
        raise ValueError(f"Unexpected source-pickle filename: {path}")
    return int(match.group("x")), int(match.group("y"))


def central_epsf_index(n_parameters: int) -> int:
    """Return the flattened center index of the ePSF component block."""
    n_epsf_parameters = n_parameters - BACKGROUND_PARAMETER_COUNT
    side_length = int(np.sqrt(n_epsf_parameters))
    if side_length * side_length != n_epsf_parameters or side_length % 2 == 0:
        raise ValueError(
            "ePSF parameters must contain an odd square block followed by "
            f"{BACKGROUND_PARAMETER_COUNT} background parameters; got {n_parameters}"
        )
    return (side_length // 2) * side_length + side_length // 2


def percentile_limits(image: np.ndarray) -> tuple[float, float]:
    positive = image[np.isfinite(image) & (image > 0)]
    if positive.size == 0:
        raise ValueError("Stitched image has no finite positive pixels")
    vmin, vmax = np.nanpercentile(positive, (5.0, 99.8))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin <= 0 or vmax <= vmin:
        raise ValueError("Could not derive positive finite image display limits")
    return float(vmin), float(vmax)


def science_pixel_slices(scipixs: str) -> tuple[slice, slice, tuple[float, float, float, float]]:
    """Return science-image slices and zero-indexed CCD display bounds."""
    try:
        x_range, y_range = scipixs.strip("[]").split(",")
        min_x, max_x = (int(value) for value in x_range.split(":"))
        min_y, max_y = (int(value) for value in y_range.split(":"))
    except ValueError as exc:
        raise ValueError(f"Invalid SCIPIXS header value: {scipixs!r}") from exc
    x0, y0 = min_x - 1, min_y - 1
    return (
        slice(y0, max_y),
        slice(x0, max_x),
        (float(x0), float(max_x), float(y0), float(max_y)),
    )


def bad_pixel_proxy(median_flux: np.ndarray) -> np.ndarray:
    """Apply TGLC's bad-pixel thresholds to a first-cadence median image."""
    bad = ~np.isfinite(median_flux)
    finite = median_flux[np.isfinite(median_flux)]
    if finite.size == 0:
        return bad
    bad |= median_flux > 0.8 * np.nanmax(finite)
    bad |= median_flux < 0.2 * np.nanmedian(finite)
    expanded = bad.copy()
    expanded[1:, :] |= bad[:-1, :]
    expanded[:-1, :] |= bad[1:, :]
    expanded[:, 1:] |= bad[:, :-1]
    expanded[:, :-1] |= bad[:, 1:]
    return expanded


def mask_edges(mask: np.ndarray) -> np.ndarray:
    """Return one-pixel outlines for a binary mask without contourpy."""
    padded = np.pad(mask, 1, constant_values=False)
    interior = (
        mask
        & padded[:-2, 1:-1]
        & padded[2:, 1:-1]
        & padded[1:-1, :-2]
        & padded[1:-1, 2:]
    )
    return mask & ~interior


def load_sector_map(
    *,
    a2v1_root: Path,
    sector: int,
    orbit: int,
    camera: int,
    ccd: int,
    n_cadences: int,
) -> tuple[
    np.ndarray,
    np.ndarray,
    tuple[float, float, float, float],
    np.ndarray,
    np.ndarray,
    int,
    int,
]:
    source_dir = a2v1_root / f"orbit-{orbit}" / "ffi" / f"cam{camera}" / f"ccd{ccd}" / "source"
    epsf_dir = source_dir.parent / "epsf"
    ffi_dir = source_dir.parent / "ffi"
    source_paths = sorted(source_dir.glob("source_*.pkl"), key=parse_source_grid)
    if not source_paths:
        raise FileNotFoundError(f"No source pickles found: {source_dir}")
    ffi_paths = sorted(path for path in ffi_dir.glob("*.fits") if path.is_file())
    if len(ffi_paths) < n_cadences:
        raise ValueError(f"{ffi_dir} has {len(ffi_paths)} FFIs; requested {n_cadences}")

    images: list[np.ndarray] = []
    science_extent: tuple[float, float, float, float] | None = None
    for ffi_path in ffi_paths[:n_cadences]:
        with fits.open(ffi_path, memmap=True) as hdus:
            y_slice, x_slice, extent = science_pixel_slices(hdus[0].header["SCIPIXS"])
            if science_extent is None:
                science_extent = extent
            elif science_extent != extent:
                raise ValueError(f"Inconsistent SCIPIXS value in {ffi_path}")
            images.append(np.asarray(hdus[0].data[y_slice, x_slice], dtype=np.float32))
    median_flux = np.nanmedian(np.stack(images, axis=0), axis=0)
    proxy_mask = bad_pixel_proxy(median_flux)

    epsf_scale: dict[tuple[int, int], float] = {}
    local_epsf_by_grid: dict[tuple[int, int], bool] = {}
    n_source_tiles = 0
    epsf_center = -1

    for source_path in source_paths:
        grid_x, grid_y = parse_source_grid(source_path)
        n_source_tiles += 1

        epsf_path = epsf_dir / f"epsf_{grid_x}_{grid_y}.npy"
        epsf = np.load(epsf_path, mmap_mode="r")
        if epsf.ndim != 2 or epsf.shape[0] < n_cadences:
            raise ValueError(f"Unexpected ePSF shape {epsf.shape} in {epsf_path}")
        center = central_epsf_index(epsf.shape[1])
        if epsf_center not in {-1, center}:
            raise ValueError(f"Inconsistent ePSF parameterization in {epsf_path}")
        epsf_center = center
        epsf_scale[(grid_x, grid_y)] = float(np.nanmedian(epsf[:n_cadences, center]))
        local_epsf_by_grid[(grid_x, grid_y)] = not epsf_path.is_symlink()

    assert science_extent is not None
    max_x = max(grid_x for grid_x, _ in epsf_scale) + 1
    max_y = max(grid_y for _, grid_y in epsf_scale) + 1
    scale_grid = np.full((max_y, max_x), np.nan, dtype=float)
    local_epsf = np.zeros((max_y, max_x), dtype=bool)
    for (grid_x, grid_y), scale in epsf_scale.items():
        scale_grid[grid_y, grid_x] = scale
        local_epsf[grid_y, grid_x] = local_epsf_by_grid[(grid_x, grid_y)]
    scale_median = float(np.nanmedian(scale_grid))
    if not np.isfinite(scale_median) or scale_median == 0:
        raise ValueError("Cannot normalize ePSF center coefficients")
    return (
        median_flux,
        proxy_mask,
        science_extent,
        scale_grid / scale_median,
        local_epsf,
        n_source_tiles,
        epsf_center,
    )


def render_sector_map(
    *,
    a2v1_root: Path,
    output_dir: Path,
    sector: int,
    orbit: int,
    camera: int,
    ccd: int,
    n_cadences: int,
) -> MapSummary:
    (
        median_flux,
        proxy_mask,
        science_extent,
        relative_scale,
        local_epsf,
        n_tiles,
        epsf_center,
    ) = load_sector_map(
        a2v1_root=a2v1_root,
        sector=sector,
        orbit=orbit,
        camera=camera,
        ccd=ccd,
        n_cadences=n_cadences,
    )
    apply_twirl_style("full_page")
    figure, (image_axis, scale_axis) = plt.subplots(
        1,
        2,
        figsize=(10.8, 4.75),
        gridspec_kw={"width_ratios": (1.13, 1.0), "wspace": 0.24},
    )
    vmin, vmax = percentile_limits(median_flux)
    image_axis.imshow(
        median_flux,
        origin="lower",
        cmap="gray",
        norm=LogNorm(vmin=vmin, vmax=vmax),
        interpolation="nearest",
        extent=science_extent,
    )
    if np.any(proxy_mask):
        image_axis.imshow(
            np.ma.masked_where(~mask_edges(proxy_mask), proxy_mask),
            origin="lower",
            cmap=ListedColormap(["#d1495b"]),
            interpolation="nearest",
            extent=science_extent,
        )
    image_axis.set(
        xlabel="CCD x [pixel]",
        ylabel="CCD y [pixel]",
        title=(
            f"S{sector:02d} orbit {orbit}, cam{camera}/ccd{ccd}\n"
            f"median of first {n_cadences} cadences; TGLC bad-pixel proxy"
        ),
    )
    image_axis.grid(False)

    finite_scale = relative_scale[np.isfinite(relative_scale)]
    deviation = max(0.02, float(np.nanmax(np.abs(finite_scale - 1.0))))
    scale_image = scale_axis.imshow(
        relative_scale,
        origin="lower",
        cmap="coolwarm",
        vmin=1.0 - deviation,
        vmax=1.0 + deviation,
        interpolation="nearest",
    )
    for grid_y, grid_x in np.ndindex(relative_scale.shape):
        value = relative_scale[grid_y, grid_x]
        if np.isfinite(value):
            scale_axis.text(
                grid_x,
                grid_y,
                f"{value:.2f}",
                ha="center",
                va="center",
                fontsize=4.7,
                color="#202020",
            )
        if local_epsf[grid_y, grid_x]:
            scale_axis.add_patch(
                Rectangle(
                    (grid_x - 0.5, grid_y - 0.5),
                    1,
                    1,
                    fill=False,
                    edgecolor="#d1495b",
                    linewidth=0.7,
                )
            )
    scale_axis.set(
        xlabel="cutout x",
        ylabel="cutout y",
        title="central ePSF coefficient\nrelative to CCD median; red = local A2v1 fit",
        xticks=np.arange(relative_scale.shape[1]),
        yticks=np.arange(relative_scale.shape[0]),
    )
    scale_axis.grid(False)
    colorbar = figure.colorbar(scale_image, ax=scale_axis, fraction=0.046, pad=0.045)
    colorbar.set_label("relative ePSF coefficient")
    figure.subplots_adjust(left=0.07, right=0.96, top=0.91, bottom=0.12)

    output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"s{sector:04d}_o{orbit}_cam{camera}_ccd{ccd}_pixel_mask_map"
    png = output_dir / f"{stem}.png"
    pdf = output_dir / f"{stem}.pdf"
    figure.savefig(png, dpi=240)
    figure.savefig(pdf)
    plt.close(figure)
    return MapSummary(
        sector=sector,
        orbit=orbit,
        camera=camera,
        ccd=ccd,
        n_cadences=n_cadences,
        n_source_tiles=n_tiles,
        n_local_epsf_tiles=int(local_epsf.sum()),
        n_proxy_bad_pixels=int(proxy_mask.sum()),
        epsf_center_index=epsf_center,
        epsf_relative_scale_median=float(np.nanmedian(relative_scale)),
        png=str(png),
        pdf=str(pdf),
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a2v1-root", type=Path, default=DEFAULT_A2V1_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--sector-orbit", type=parse_sector_orbit, action="append")
    parser.add_argument("--camera", type=int, default=1, choices=range(1, 5))
    parser.add_argument("--ccd", type=int, default=1, choices=range(1, 5))
    parser.add_argument("--n-cadences", type=int, default=20)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.n_cadences < 1:
        raise ValueError("--n-cadences must be positive")
    sector_orbits = args.sector_orbit or list(DEFAULT_SECTOR_ORBITS.items())
    summaries = [
        render_sector_map(
            a2v1_root=args.a2v1_root,
            output_dir=args.output_dir,
            sector=sector,
            orbit=orbit,
            camera=args.camera,
            ccd=args.ccd,
            n_cadences=args.n_cadences,
        )
        for sector, orbit in sector_orbits
    ]
    summary_path = args.output_dir / "summary.json"
    summary_path.write_text(json.dumps([asdict(item) for item in summaries], indent=2) + "\n")
    for summary in summaries:
        print(
            "[a2v1-mask-map] "
            f"S{summary.sector:02d} orbit={summary.orbit} tiles={summary.n_source_tiles} "
            f"local_epsf_tiles={summary.n_local_epsf_tiles} "
            f"proxy_bad_pixels={summary.n_proxy_bad_pixels} "
            f"png={summary.png}"
        )


if __name__ == "__main__":
    main()
