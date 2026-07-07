#!/usr/bin/env python3
"""Build lightweight TIC-match overlays for existing TGLC source pickles.

The TGLC ``Source`` pickle stores a large FFI cutout cube plus a small
``source.tic`` table.  For TWIRL A2v1 production, the expensive pixel payload can
be reused, but the TIC table must be regenerated from no-cap TIC catalogs.  This
script writes small per-source sidecars:

    source_tic/source_X_Y.ecsv

A tiny TGLC light-curve-stage patch can then replace ``source.tic`` with the
sidecar table after unpickling, without rewriting the large source pickle.
For requested TWIRL WD targets, sidecars can be built directly from the
observation table detector coordinates with ``--overlay-from-observations``.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.lightcurves.tglc_faint_allstar import _clean_int_column, _find_column


DEFAULT_TGLC_DATA_DIR = Path("/pdo/users/tehan/tglc-gpu-production")
DEFAULT_OBSERVATIONS_FITS = Path(
    "/pdo/users/tehan/TWIRL/data_local/catalogs/twirl_master_catalog/"
    "twirl_wd_tess_observations_v0.fits"
)
PDO_USER_ROOT = Path("/pdo/users/tehan")
SOURCE_RE = re.compile(r"^source_(?P<x>\d+)_(?P<y>\d+)\.pkl$")


@dataclass(frozen=True)
class Ccd:
    camera: int
    ccd: int

    @property
    def label(self) -> str:
        return f"cam{self.camera}_ccd{self.ccd}"


@dataclass
class CcdOverlaySummary:
    orbit: int
    camera: int
    ccd: int
    overlay_source: str
    source_dir: str
    output_source_tic_dir: str
    tic_catalog_file: str | None
    observations_fits: str | None
    ffi_file: str | None
    n_source_files: int
    n_sidecars_written: int
    n_sidecars_existing: int
    n_source_links_written: int
    n_source_links_existing: int
    n_total_overlay_rows: int
    n_nonempty_sidecars: int
    n_allowed_tics: int | None


def parse_ccd(value: str) -> Ccd:
    try:
        camera_text, ccd_text = value.split(",", maxsplit=1)
        camera = int(camera_text)
        ccd = int(ccd_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"CCD must have format camera,ccd; got {value!r}"
        ) from exc
    if camera not in {1, 2, 3, 4} or ccd not in {1, 2, 3, 4}:
        raise argparse.ArgumentTypeError(f"Camera/CCD must be in 1..4; got {value!r}")
    return Ccd(camera=camera, ccd=ccd)


def default_ccds() -> list[Ccd]:
    return [Ccd(camera=camera, ccd=ccd) for camera in range(1, 5) for ccd in range(1, 5)]


def enforce_pdo_write_path(path: Path) -> None:
    """Reject writes into shared PDO trees while allowing local dry-run tests."""
    resolved = Path(path).expanduser().resolve(strict=False)
    if not str(resolved).startswith("/pdo/"):
        return
    allowed = PDO_USER_ROOT.resolve(strict=False)
    if resolved != allowed and allowed not in resolved.parents:
        raise ValueError(f"Refusing to write outside {allowed}: {resolved}")


def parse_source_grid(path: Path) -> tuple[int, int]:
    match = SOURCE_RE.match(path.name)
    if match is None:
        raise ValueError(f"Not a TGLC source pickle name: {path.name}")
    return int(match.group("x")), int(match.group("y"))


def source_bounds(
    grid_x: int,
    grid_y: int,
    *,
    cutout_size: int = 150,
    cutout_overlap: int = 2,
    ccd_x_offset: int = 44,
) -> tuple[float, float, float, float]:
    """Return inclusive CCD pixel bounds matching ``tglc.ffi.Source``."""
    stride = cutout_size - (2 * cutout_overlap)
    left = grid_x * stride + ccd_x_offset
    bottom = grid_y * stride
    return float(left), float(left + cutout_size), float(bottom), float(bottom + cutout_size)


def normalize_tic_catalog_with_positions(
    tic_catalog: Table,
    *,
    tmag_max: float | None = None,
    allowed_tics: set[int] | None = None,
) -> Table:
    """Normalize a TGLC TIC ECSV to TIC, Gaia DR3 source ID, Tmag, RA, and Dec."""
    tic_col = _find_column(tic_catalog, ("TIC", "ID", "tic", "id"))
    gaia_col = _find_column(tic_catalog, ("gaia3", "dr3_source_id", "GAIA3", "GAIA"))
    tmag_col = _find_column(tic_catalog, ("Tmag", "tmag", "TESSMAG", "tessmag"))
    ra_col = _find_column(tic_catalog, ("ra", "RA"))
    dec_col = _find_column(tic_catalog, ("dec", "DEC"))

    tic, tic_good = _clean_int_column(tic_catalog[tic_col])
    gaia3, gaia_good = _clean_int_column(tic_catalog[gaia_col])
    tmag = np.asarray(tic_catalog[tmag_col], dtype=float)
    ra = np.asarray(tic_catalog[ra_col], dtype=float)
    dec = np.asarray(tic_catalog[dec_col], dtype=float)

    good = tic_good & gaia_good & np.isfinite(tmag) & np.isfinite(ra) & np.isfinite(dec)
    if tmag_max is not None:
        good &= tmag <= float(tmag_max)
    if allowed_tics is not None:
        allowed = np.fromiter(allowed_tics, dtype=np.int64, count=len(allowed_tics))
        good &= np.isin(tic, allowed)

    out = Table()
    out["TIC"] = tic[good]
    out["gaia3"] = gaia3[good]
    out["tmag"] = tmag[good]
    out["ra"] = ra[good]
    out["dec"] = dec[good]
    if len(out) == 0:
        return out

    order = np.lexsort((np.asarray(out["gaia3"]), np.asarray(out["TIC"])))
    out = out[order]
    keep = np.ones(len(out), dtype=bool)
    keep[1:] = np.asarray(out["TIC"][1:]) != np.asarray(out["TIC"][:-1])
    return out[keep]


def project_tic_catalog(tic_catalog: Table, wcs: WCS) -> Table:
    """Attach CCD x/y coordinates to a normalized TIC catalog."""
    if len(tic_catalog) == 0:
        projected = tic_catalog.copy()
        projected["ccd_x"] = np.asarray([], dtype=float)
        projected["ccd_y"] = np.asarray([], dtype=float)
        return projected
    coords = SkyCoord(tic_catalog["ra"], tic_catalog["dec"], unit="deg")
    x, y = wcs.world_to_pixel(coords)
    projected = tic_catalog.copy()
    projected["ccd_x"] = np.asarray(x, dtype=float)
    projected["ccd_y"] = np.asarray(y, dtype=float)
    return projected


def overlay_for_source_bounds(
    projected_tic_catalog: Table,
    bounds: tuple[float, float, float, float],
) -> Table:
    """Return the two-column ``source.tic`` replacement for one source cutout."""
    left, right, bottom, top = bounds
    if len(projected_tic_catalog) == 0:
        return Table({"TIC": np.array([], dtype=np.int64), "gaia3": np.array([], dtype=np.int64)})

    x = np.asarray(projected_tic_catalog["ccd_x"], dtype=float)
    y = np.asarray(projected_tic_catalog["ccd_y"], dtype=float)
    inside = (left <= x) & (x <= right) & (bottom <= y) & (y <= top)
    matched = projected_tic_catalog[inside]

    out = Table()
    out["TIC"] = np.asarray(matched["TIC"], dtype=np.int64)
    out["gaia3"] = np.asarray(matched["gaia3"], dtype=np.int64)
    return out


def load_observation_targets(
    observations_fits: Path,
    *,
    orbit: int,
    camera: int,
    ccd: int,
    tmag_max: float | None = None,
) -> Table:
    """Load requested TWIRL targets with detector coordinates for one orbit/CCD."""
    with fits.open(observations_fits, memmap=True) as hdul:
        rows = hdul[1].data
        tic = np.asarray(rows["tic_id"])
        gaia3 = np.asarray(rows["source_id"])
        x = np.asarray(rows["colpix"], dtype=float)
        y = np.asarray(rows["rowpix"], dtype=float)
        mask = (
            (np.asarray(rows["orbit"]) == orbit)
            & (np.asarray(rows["camera"]) == camera)
            & (np.asarray(rows["ccd"]) == ccd)
            & (tic > 0)
            & np.isfinite(x)
            & np.isfinite(y)
        )
        if tmag_max is not None:
            mask &= np.asarray(rows["tmag"], dtype=float) <= float(tmag_max)

        out = Table()
        out["TIC"] = tic[mask].astype(np.int64)
        out["gaia3"] = gaia3[mask].astype(np.int64)
        out["ccd_x"] = x[mask]
        out["ccd_y"] = y[mask]

    if len(out) == 0:
        return out

    order = np.lexsort((np.asarray(out["gaia3"]), np.asarray(out["TIC"])))
    out = out[order]
    keep = np.ones(len(out), dtype=bool)
    keep[1:] = np.asarray(out["TIC"][1:]) != np.asarray(out["TIC"][:-1])
    return out[keep]


def load_allowed_tics(
    observations_fits: Path,
    *,
    orbit: int,
    camera: int,
    ccd: int,
) -> set[int]:
    with fits.open(observations_fits, memmap=True) as hdul:
        rows = hdul[1].data
        mask = (
            (np.asarray(rows["orbit"]) == orbit)
            & (np.asarray(rows["camera"]) == camera)
            & (np.asarray(rows["ccd"]) == ccd)
            & (np.asarray(rows["tic_id"]) > 0)
        )
        return set(np.unique(np.asarray(rows["tic_id"])[mask]).astype(np.int64).tolist())


def first_ffi_file(ffi_dir: Path) -> Path:
    files = sorted(ffi_dir.glob("*.fits"))
    if not files:
        raise FileNotFoundError(f"No FITS files found in {ffi_dir}")
    return files[0]


def write_table_if_requested(
    table: Table,
    path: Path,
    *,
    apply: bool,
    overwrite: bool,
) -> str:
    if path.exists() and not overwrite:
        return "existing"
    if not apply:
        return "dry_run"
    enforce_pdo_write_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    table.write(path, format="ascii.ecsv", overwrite=True)
    return "written"


def link_source_if_requested(
    source_file: Path,
    output_source_file: Path,
    *,
    apply: bool,
    overwrite: bool,
) -> str:
    if source_file.resolve(strict=False) == output_source_file.resolve(strict=False):
        return "same_root"
    if output_source_file.exists() or output_source_file.is_symlink():
        if output_source_file.is_symlink() and output_source_file.resolve(strict=False) == source_file.resolve(strict=False):
            return "existing"
        if not overwrite:
            return "existing"
        if apply:
            output_source_file.unlink()
    if not apply:
        return "dry_run"
    enforce_pdo_write_path(output_source_file)
    output_source_file.parent.mkdir(parents=True, exist_ok=True)
    output_source_file.symlink_to(source_file.resolve())
    return "written"


def build_one_ccd_overlays(
    *,
    source_tglc_data_dir: Path,
    catalog_tglc_data_dir: Path,
    output_tglc_data_dir: Path,
    orbit: int,
    ccd: Ccd,
    catalog_orbit: int | None,
    observations_fits: Path | None,
    overlay_from_observations: bool,
    restrict_to_observation_tics: bool,
    tmag_max: float | None,
    apply: bool,
    overwrite: bool,
    link_sources: bool,
    limit_source_files: int | None,
    cutout_size: int,
    cutout_overlap: int,
    ccd_x_offset: int,
) -> CcdOverlaySummary:
    input_ccd_root = source_tglc_data_dir / f"orbit-{orbit}" / "ffi" / f"cam{ccd.camera}" / f"ccd{ccd.ccd}"
    output_ccd_root = output_tglc_data_dir / f"orbit-{orbit}" / "ffi" / f"cam{ccd.camera}" / f"ccd{ccd.ccd}"
    source_dir = input_ccd_root / "source"
    output_source_dir = output_ccd_root / "source"
    output_source_tic_dir = output_ccd_root / "source_tic"
    tic_catalog_file = (
        catalog_tglc_data_dir
        / f"orbit-{catalog_orbit if catalog_orbit is not None else orbit}"
        / "ffi"
        / "catalogs"
        / f"TIC_cam{ccd.camera}_ccd{ccd.ccd}.ecsv"
    )

    if not source_dir.exists():
        raise FileNotFoundError(f"Missing source directory: {source_dir}")

    allowed_tics: set[int] | None = None
    if restrict_to_observation_tics and not overlay_from_observations:
        if observations_fits is None:
            raise ValueError("--restrict-to-observation-tics requires --observations-fits")
        allowed_tics = load_allowed_tics(
            observations_fits,
            orbit=orbit,
            camera=ccd.camera,
            ccd=ccd.ccd,
        )

    if overlay_from_observations:
        if observations_fits is None:
            raise ValueError("--overlay-from-observations requires --observations-fits")
        projected = load_observation_targets(
            observations_fits,
            orbit=orbit,
            camera=ccd.camera,
            ccd=ccd.ccd,
            tmag_max=tmag_max,
        )
        tic_catalog_path: Path | None = None
        ffi_path: Path | None = None
    else:
        if not tic_catalog_file.exists():
            raise FileNotFoundError(f"Missing TIC catalog: {tic_catalog_file}")
        tic_catalog = Table.read(tic_catalog_file)
        normalized = normalize_tic_catalog_with_positions(
            tic_catalog,
            tmag_max=tmag_max,
            allowed_tics=allowed_tics,
        )
        ffi_path = first_ffi_file(input_ccd_root / "ffi")
        with fits.open(ffi_path, memmap=False) as hdul:
            wcs = WCS(hdul[0].header)
        projected = project_tic_catalog(normalized, wcs)
        tic_catalog_path = tic_catalog_file

    source_files = sorted(source_dir.glob("source_*.pkl"))
    if limit_source_files is not None:
        source_files = source_files[: int(limit_source_files)]

    n_sidecars_written = 0
    n_sidecars_existing = 0
    n_source_links_written = 0
    n_source_links_existing = 0
    n_total_overlay_rows = 0
    n_nonempty_sidecars = 0

    for source_file in source_files:
        grid_x, grid_y = parse_source_grid(source_file)
        bounds = source_bounds(
            grid_x,
            grid_y,
            cutout_size=cutout_size,
            cutout_overlap=cutout_overlap,
            ccd_x_offset=ccd_x_offset,
        )
        overlay = overlay_for_source_bounds(projected, bounds)
        n_total_overlay_rows += len(overlay)
        n_nonempty_sidecars += int(len(overlay) > 0)

        sidecar_path = output_source_tic_dir / source_file.with_suffix(".ecsv").name
        status = write_table_if_requested(
            overlay,
            sidecar_path,
            apply=apply,
            overwrite=overwrite,
        )
        n_sidecars_written += int(status == "written")
        n_sidecars_existing += int(status == "existing")

        if link_sources:
            link_status = link_source_if_requested(
                source_file,
                output_source_dir / source_file.name,
                apply=apply,
                overwrite=overwrite,
            )
            n_source_links_written += int(link_status == "written")
            n_source_links_existing += int(link_status == "existing")

    return CcdOverlaySummary(
        orbit=orbit,
        camera=ccd.camera,
        ccd=ccd.ccd,
        overlay_source="observations" if overlay_from_observations else "tic_catalog",
        source_dir=str(source_dir),
        output_source_tic_dir=str(output_source_tic_dir),
        tic_catalog_file=None if tic_catalog_path is None else str(tic_catalog_path),
        observations_fits=str(observations_fits) if overlay_from_observations else None,
        ffi_file=None if ffi_path is None else str(ffi_path),
        n_source_files=len(source_files),
        n_sidecars_written=n_sidecars_written,
        n_sidecars_existing=n_sidecars_existing,
        n_source_links_written=n_source_links_written,
        n_source_links_existing=n_source_links_existing,
        n_total_overlay_rows=n_total_overlay_rows,
        n_nonempty_sidecars=n_nonempty_sidecars,
        n_allowed_tics=len(projected) if overlay_from_observations else (
            None if allowed_tics is None else len(allowed_tics)
        ),
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tglc-data-dir", type=Path, default=DEFAULT_TGLC_DATA_DIR)
    parser.add_argument(
        "--source-tglc-data-dir",
        type=Path,
        default=None,
        help="Root containing existing source pickles and FFI symlinks. Defaults to --tglc-data-dir.",
    )
    parser.add_argument(
        "--catalog-tglc-data-dir",
        type=Path,
        default=None,
        help="Root containing max-99 TIC catalogs. Defaults to --tglc-data-dir.",
    )
    parser.add_argument(
        "--output-tglc-data-dir",
        type=Path,
        default=None,
        help="Root for sidecars/source symlinks. Defaults to --tglc-data-dir.",
    )
    parser.add_argument("--orbit", type=int, required=True)
    parser.add_argument(
        "--catalog-orbit",
        type=int,
        default=None,
        help="Orbit whose TIC_cam*_ccd*.ecsv catalogs should be used.",
    )
    parser.add_argument("--ccd", type=parse_ccd, action="append", default=None)
    parser.add_argument("--observations-fits", type=Path, default=DEFAULT_OBSERVATIONS_FITS)
    parser.add_argument(
        "--overlay-from-observations",
        action="store_true",
        help=(
            "Build source_tic sidecars directly from TWIRL observation-table detector "
            "positions, bypassing no-cap TIC catalog generation."
        ),
    )
    parser.add_argument(
        "--restrict-to-observation-tics",
        action="store_true",
        help="Only include TIC IDs requested by the TWIRL observation table for this orbit/CCD.",
    )
    parser.add_argument("--tmag-max", type=float, default=None)
    parser.add_argument("--cutout-size", type=int, default=150)
    parser.add_argument("--cutout-overlap", type=int, default=2)
    parser.add_argument("--ccd-x-offset", type=int, default=44)
    parser.add_argument("--limit-source-files", type=int, default=None)
    parser.add_argument("--link-sources", action="store_true")
    parser.add_argument("--apply", action="store_true", help="Actually write sidecars/symlinks.")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--summary-json", type=Path, default=None)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    ccds: list[Ccd] = args.ccd if args.ccd is not None else default_ccds()
    source_root = args.source_tglc_data_dir or args.tglc_data_dir
    catalog_root = args.catalog_tglc_data_dir or args.tglc_data_dir
    output_root = args.output_tglc_data_dir or args.tglc_data_dir

    summaries: list[CcdOverlaySummary] = []
    for ccd in ccds:
        summary = build_one_ccd_overlays(
            source_tglc_data_dir=source_root,
            catalog_tglc_data_dir=catalog_root,
            output_tglc_data_dir=output_root,
            orbit=args.orbit,
            ccd=ccd,
            catalog_orbit=args.catalog_orbit,
            observations_fits=args.observations_fits,
            overlay_from_observations=args.overlay_from_observations,
            restrict_to_observation_tics=args.restrict_to_observation_tics,
            tmag_max=args.tmag_max,
            apply=args.apply,
            overwrite=args.overwrite,
            link_sources=args.link_sources,
            limit_source_files=args.limit_source_files,
            cutout_size=args.cutout_size,
            cutout_overlap=args.cutout_overlap,
            ccd_x_offset=args.ccd_x_offset,
        )
        summaries.append(summary)
        print(
            f"[source-tic] orbit={summary.orbit} cam{summary.camera}/ccd{summary.ccd} "
            f"sources={summary.n_source_files} nonempty={summary.n_nonempty_sidecars} "
            f"rows={summary.n_total_overlay_rows} written={summary.n_sidecars_written}",
            flush=True,
        )

    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "apply": bool(args.apply),
        "tglc_data_dir": str(args.tglc_data_dir),
        "source_tglc_data_dir": str(source_root),
        "catalog_tglc_data_dir": str(catalog_root),
        "output_tglc_data_dir": str(output_root),
        "orbit": args.orbit,
        "catalog_orbit": args.catalog_orbit,
        "overlay_from_observations": bool(args.overlay_from_observations),
        "restrict_to_observation_tics": bool(args.restrict_to_observation_tics),
        "tmag_max": args.tmag_max,
        "link_sources": bool(args.link_sources),
        "summaries": [asdict(summary) for summary in summaries],
    }
    if args.summary_json is not None:
        if args.apply:
            enforce_pdo_write_path(args.summary_json)
        args.summary_json.parent.mkdir(parents=True, exist_ok=True)
        args.summary_json.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    else:
        print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
