#!/usr/bin/env python3
"""Render representative A2v1 saturated-pixel mask diagnostics by sector.

Each map stitches the median of the first requested cadences from the 14x14
TGLC source-cutout grid for one orbit/camera/CCD. Red contours show the static
``Source.mask.mask`` pixels passed to the A2v1 ePSF fitter. The companion grid
shows the cadence-median central ePSF coefficient for each cutout, normalized
by the CCD median. It is a relative fit-scale diagnostic, not photometry.
"""

from __future__ import annotations

import argparse
import json
import pickle
import re
from dataclasses import dataclass, asdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

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
    n_masked_tiles: int
    n_masked_pixels: int
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


def load_sector_map(
    *,
    a2v1_root: Path,
    sector: int,
    orbit: int,
    camera: int,
    ccd: int,
    n_cadences: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, int, int]:
    source_dir = a2v1_root / f"orbit-{orbit}" / "ffi" / f"cam{camera}" / f"ccd{ccd}" / "source"
    epsf_dir = source_dir.parent / "epsf"
    source_paths = sorted(source_dir.glob("source_*.pkl"), key=parse_source_grid)
    if not source_paths:
        raise FileNotFoundError(f"No source pickles found: {source_dir}")

    mosaic_sum: np.ndarray | None = None
    mosaic_count: np.ndarray | None = None
    masked_count: np.ndarray | None = None
    epsf_scale: dict[tuple[int, int], float] = {}
    n_masked_tiles = 0
    n_source_tiles = 0
    epsf_center = -1

    for source_path in source_paths:
        grid_x, grid_y = parse_source_grid(source_path)
        with source_path.open("rb") as source_handle:
            source = pickle.load(source_handle)
        flux = np.asarray(source.flux[:n_cadences], dtype=float)
        if flux.shape[0] < n_cadences:
            raise ValueError(
                f"{source_path} has only {flux.shape[0]} cadences; requested {n_cadences}"
            )
        tile = np.nanmedian(flux, axis=0)
        source_mask = np.ma.getmaskarray(source.mask).astype(bool, copy=False)
        if tile.shape != source_mask.shape:
            raise ValueError(f"Flux/mask shape mismatch in {source_path}")

        x0 = int(source.ccd_x)
        y0 = int(source.ccd_y)
        y1 = y0 + tile.shape[0]
        x1 = x0 + tile.shape[1]
        if mosaic_sum is None:
            mosaic_sum = np.zeros((y1, x1), dtype=float)
            mosaic_count = np.zeros((y1, x1), dtype=np.uint16)
            masked_count = np.zeros((y1, x1), dtype=np.uint8)
        elif y1 > mosaic_sum.shape[0] or x1 > mosaic_sum.shape[1]:
            new_shape = (max(y1, mosaic_sum.shape[0]), max(x1, mosaic_sum.shape[1]))
            expanded_sum = np.zeros(new_shape, dtype=float)
            expanded_count = np.zeros(new_shape, dtype=np.uint16)
            expanded_mask = np.zeros(new_shape, dtype=np.uint8)
            expanded_sum[: mosaic_sum.shape[0], : mosaic_sum.shape[1]] = mosaic_sum
            expanded_count[: mosaic_count.shape[0], : mosaic_count.shape[1]] = mosaic_count
            expanded_mask[: masked_count.shape[0], : masked_count.shape[1]] = masked_count
            mosaic_sum, mosaic_count, masked_count = expanded_sum, expanded_count, expanded_mask

        finite = np.isfinite(tile)
        mosaic_sum[y0:y1, x0:x1] += np.where(finite, tile, 0.0)
        mosaic_count[y0:y1, x0:x1] += finite.astype(np.uint16)
        masked_count[y0:y1, x0:x1] += source_mask.astype(np.uint8)
        n_masked_tiles += int(source_mask.any())
        n_source_tiles += 1

        epsf_path = epsf_dir / f"epsf_{grid_x}_{grid_y}.npy"
        epsf = np.load(epsf_path)
        if epsf.ndim != 2 or epsf.shape[0] < n_cadences:
            raise ValueError(f"Unexpected ePSF shape {epsf.shape} in {epsf_path}")
        center = central_epsf_index(epsf.shape[1])
        if epsf_center not in {-1, center}:
            raise ValueError(f"Inconsistent ePSF parameterization in {epsf_path}")
        epsf_center = center
        epsf_scale[(grid_x, grid_y)] = float(np.nanmedian(epsf[:n_cadences, center]))

    assert mosaic_sum is not None
    assert mosaic_count is not None
    assert masked_count is not None
    mosaic = np.divide(
        mosaic_sum,
        mosaic_count,
        out=np.full(mosaic_sum.shape, np.nan, dtype=float),
        where=mosaic_count > 0,
    )
    mask = masked_count > 0
    max_x = max(grid_x for grid_x, _ in epsf_scale) + 1
    max_y = max(grid_y for _, grid_y in epsf_scale) + 1
    scale_grid = np.full((max_y, max_x), np.nan, dtype=float)
    for (grid_x, grid_y), scale in epsf_scale.items():
        scale_grid[grid_y, grid_x] = scale
    scale_median = float(np.nanmedian(scale_grid))
    if not np.isfinite(scale_median) or scale_median == 0:
        raise ValueError("Cannot normalize ePSF center coefficients")
    return mosaic, mask, scale_grid / scale_median, n_source_tiles, n_masked_tiles, epsf_center


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
    mosaic, mask, relative_scale, n_tiles, n_masked_tiles, epsf_center = load_sector_map(
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
    vmin, vmax = percentile_limits(mosaic)
    image_axis.imshow(
        mosaic,
        origin="lower",
        cmap="gray",
        norm=LogNorm(vmin=vmin, vmax=vmax),
        interpolation="nearest",
    )
    if np.any(mask):
        image_axis.contour(
            mask.astype(float),
            levels=(0.5,),
            origin="lower",
            colors=("#d1495b",),
            linewidths=0.45,
        )
    image_axis.set(
        xlabel="CCD x [pixel]",
        ylabel="CCD y [pixel]",
        title=(
            f"S{sector:02d} orbit {orbit}, cam{camera}/ccd{ccd}\n"
            f"median of first {n_cadences} cadences; mask outlines"
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
    scale_axis.set(
        xlabel="cutout x",
        ylabel="cutout y",
        title="central ePSF coefficient\nrelative to CCD median",
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
        n_masked_tiles=n_masked_tiles,
        n_masked_pixels=int(mask.sum()),
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
            f"masked_tiles={summary.n_masked_tiles} masked_pixels={summary.n_masked_pixels} "
            f"png={summary.png}"
        )


if __name__ == "__main__":
    main()
