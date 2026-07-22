"""Build an independent WD 1856 S56 light curve from MAST TESSCut pixels.

The current A2v1 compact product is used only as a timestamp-to-cadence
authority.  All photometry in the emitted reference product comes directly
from the official TESSCut ``FLUX`` pixel cube.  This separation is deliberate:
the product is suitable as the non-TGLC/TWIRL reference consumed by the
Tier-1 independent-extraction comparison.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import tempfile
from typing import Any, Callable, Mapping, Sequence
from urllib.parse import urlparse
import uuid

import h5py
import numpy as np
import pandas as pd

from twirl.lightcurves.a2v1_qa import WD1856_TIC, file_sha256
from twirl.lightcurves.external_quality import load_external_quality_reference


SECTOR = 56
CAMERA = 4
CCD = 1
TIME_OFFSET_BJD = 2_457_000.0
# TIC/TOI position used for the MAST TESSCut request.  This is intentionally
# the TESS catalog position rather than the older, unpropagated J2000 position.
WD1856_TESSCUT_RA_DEG = 284.415675
WD1856_TESSCUT_DEC_DEG = 53.509024
DEFAULT_POSITION_TOLERANCE_ARCSEC = 1.0
DEFAULT_MAX_TIME_DELTA_SECONDS = 60.0
DEFAULT_ROLLING_WINDOW_CADENCES = 433
BUILDER_VERSION = "wd1856_s56_tesscut_independent_extraction_v2"
MASKED_TREND_FILL_POLICY = (
    "linear interpolation by row index between supported rolling medians, "
    "with nearest-supported edge carry; permitted only where effective quality != 0"
)
MASKED_TREND_FILL_HEADER = "linear row interpolation; quality-masked only"
CADENCE_MAP_HASH_ALGORITHM = (
    "sha256(contract UTF-8 + row-count <u8 + TESSCut TIME <f8 + "
    "compact TIME <f8 + CADENCENO <i8)"
)


@dataclass(frozen=True)
class _SourceState:
    path: Path
    sha256: str
    size_bytes: int
    mtime_ns: int


def _capture_source_state(path: Path) -> _SourceState:
    resolved = Path(path).resolve()
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    stat = resolved.stat()
    return _SourceState(
        path=resolved,
        sha256=file_sha256(resolved),
        size_bytes=int(stat.st_size),
        mtime_ns=int(stat.st_mtime_ns),
    )


def _assert_sources_unchanged(states: Sequence[_SourceState]) -> None:
    for initial in states:
        if not initial.path.is_file():
            raise RuntimeError(f"source disappeared during extraction: {initial.path}")
        stat = initial.path.stat()
        observed = _SourceState(
            path=initial.path,
            sha256=file_sha256(initial.path),
            size_bytes=int(stat.st_size),
            mtime_ns=int(stat.st_mtime_ns),
        )
        if observed != initial:
            raise RuntimeError(f"source mutated during extraction: {initial.path}")


def _validate_source_url(source_url: str) -> str:
    value = str(source_url).strip()
    parsed = urlparse(value)
    hostname = (parsed.hostname or "").lower()
    official_host = hostname == "stsci.edu" or hostname.endswith(".stsci.edu")
    if parsed.scheme.lower() != "https" or not official_host:
        raise ValueError("source_url must be an HTTPS STScI/MAST URL")
    if not parsed.path:
        raise ValueError("source_url must include a nonempty path")
    return value


def _header_values(headers: Sequence[Any], keyword: str) -> list[Any]:
    return [header[keyword] for header in headers if keyword in header]


def _required_integer_header(
    headers: Sequence[Any], keyword: str, expected: int
) -> int:
    raw = _header_values(headers, keyword)
    if not raw:
        raise ValueError(f"TESSCut FITS is missing required {keyword} metadata")
    values: list[int] = []
    for value in raw:
        number = float(value)
        if not np.isfinite(number) or number != np.floor(number):
            raise ValueError(f"TESSCut {keyword} metadata is not an integer")
        values.append(int(number))
    if set(values) != {expected}:
        raise ValueError(
            f"TESSCut {keyword} must be {expected}; observed {sorted(set(values))}"
        )
    return expected


def _position_separation_arcsec(
    ra_deg: float, dec_deg: float, target_ra_deg: float, target_dec_deg: float
) -> float:
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    left = SkyCoord(ra_deg * u.deg, dec_deg * u.deg, frame="icrs")
    right = SkyCoord(target_ra_deg * u.deg, target_dec_deg * u.deg, frame="icrs")
    return float(left.separation(right).arcsec)


def _validate_header_position(
    headers: Sequence[Any], *, tolerance_arcsec: float
) -> list[float]:
    coordinate_pairs: list[tuple[Any, Any]] = []
    for header in headers:
        has_ra = "RA_OBJ" in header
        has_dec = "DEC_OBJ" in header
        if has_ra != has_dec:
            raise ValueError("TESSCut RA_OBJ/DEC_OBJ metadata are not paired")
        if has_ra:
            coordinate_pairs.append((header["RA_OBJ"], header["DEC_OBJ"]))
    if not coordinate_pairs:
        raise ValueError("TESSCut FITS must declare RA_OBJ and DEC_OBJ")
    separations: list[float] = []
    for ra_raw, dec_raw in coordinate_pairs:
        ra = float(ra_raw)
        dec = float(dec_raw)
        if not np.isfinite([ra, dec]).all() or not (0 <= ra < 360 and -90 <= dec <= 90):
            raise ValueError("TESSCut RA_OBJ/DEC_OBJ metadata are invalid")
        separation = _position_separation_arcsec(
            ra,
            dec,
            WD1856_TESSCUT_RA_DEG,
            WD1856_TESSCUT_DEC_DEG,
        )
        if separation > tolerance_arcsec:
            raise ValueError(
                "TESSCut target position is not WD 1856: "
                f"separation={separation:.6g} arcsec exceeds "
                f"{tolerance_arcsec:.6g} arcsec"
            )
        separations.append(separation)
    return separations


def _validate_time_metadata(headers: Sequence[Any]) -> None:
    references = [
        float(header["BJDREFI"]) + float(header.get("BJDREFF", 0.0))
        for header in headers
        if "BJDREFI" in header
    ]
    if not references:
        raise ValueError("TESSCut FITS must declare BJDREFI")
    if not np.isfinite(references).all() or any(
        abs(value - TIME_OFFSET_BJD) > 1.0e-9 for value in references
    ):
        raise ValueError("TESSCut TIME must use BJD - 2457000")
    units = {str(value).strip().lower() for value in _header_values(headers, "TIMEUNIT")}
    if units and not units <= {"d", "day", "days"}:
        raise ValueError(f"TESSCut TIMEUNIT must be days; observed {sorted(units)}")
    systems = {str(value).strip().upper() for value in _header_values(headers, "TIMESYS")}
    if systems and systems != {"TDB"}:
        raise ValueError(f"TESSCut TIMESYS must be TDB; observed {sorted(systems)}")


def _extract_flux_wcs(header: Any, flux_column_number: int) -> Any:
    from astropy.wcs import WCS

    try:
        wcs = WCS(
            header,
            keysel=["binary"],
            colsel=[int(flux_column_number)],
            relax=True,
        )
    except Exception as exc:
        raise ValueError("TESSCut FLUX column lacks a valid celestial WCS") from exc
    if int(wcs.naxis) != 2 or not bool(wcs.has_celestial):
        raise ValueError("TESSCut FLUX column lacks a two-dimensional celestial WCS")
    return wcs


def _read_tesscut(path: Path, *, tolerance_arcsec: float) -> dict[str, Any]:
    from astropy.io import fits

    with fits.open(path, mode="readonly", memmap=False, checksum=True) as hdul:
        matches = []
        for hdu in hdul:
            columns = getattr(hdu, "columns", None)
            names = (
                {str(name).upper() for name in (columns.names or [])}
                if columns is not None
                else set()
            )
            if {"TIME", "QUALITY", "FLUX"} <= names:
                matches.append(hdu)
        if len(matches) != 1:
            raise ValueError(
                "TESSCut FITS must contain exactly one table with TIME, QUALITY, "
                f"and FLUX; found {len(matches)}"
            )
        pixels = matches[0]
        headers = (hdul[0].header, pixels.header)
        _required_integer_header(headers, "SECTOR", SECTOR)
        _required_integer_header(headers, "CAMERA", CAMERA)
        _required_integer_header(headers, "CCD", CCD)
        position_separations = _validate_header_position(
            headers, tolerance_arcsec=tolerance_arcsec
        )
        _validate_time_metadata(headers)

        origin = " ".join(str(value) for value in _header_values(headers, "ORIGIN"))
        creator = " ".join(str(value) for value in _header_values(headers, "CREATOR"))
        if not any(token in origin.upper() for token in ("STSCI", "MAST")):
            raise ValueError("TESSCut FITS ORIGIN does not identify STScI/MAST")
        if not any(token in creator.upper() for token in ("ASTROCUT", "TESSCUT")):
            raise ValueError("TESSCut FITS CREATOR does not identify Astrocut/TESSCut")

        column_names = [str(name).upper() for name in pixels.columns.names]
        flux_column_number = column_names.index("FLUX") + 1
        wcs = _extract_flux_wcs(pixels.header, flux_column_number)
        time = np.asarray(pixels.data["TIME"], dtype=np.float64)
        quality_raw = np.asarray(pixels.data["QUALITY"])
        quality_float = np.asarray(quality_raw, dtype=np.float64)
        flux = np.asarray(pixels.data["FLUX"], dtype=np.float64)

        if time.ndim != 1 or quality_float.ndim != 1 or flux.ndim != 3:
            raise ValueError("TESSCut TIME, QUALITY, or FLUX has an invalid rank")
        if not len(time) or len(time) != len(quality_float) or len(time) != len(flux):
            raise ValueError("TESSCut TIME, QUALITY, and FLUX lengths are inconsistent")
        if not np.isfinite(time).all() or not np.isfinite(flux).all():
            raise ValueError("TESSCut TIME and FLUX arrays must be entirely finite")
        if not np.isfinite(quality_float).all() or not np.equal(
            quality_float, np.floor(quality_float)
        ).all():
            raise ValueError("TESSCut QUALITY must contain finite integers")
        quality = quality_float.astype(np.int64)
        if (quality < 0).any():
            raise ValueError("TESSCut QUALITY values must be nonnegative")
        if len(np.unique(time)) != len(time):
            raise ValueError("TESSCut TIME contains duplicate timestamps")

        _, ny, nx = flux.shape
        if nx != ny or nx < 3 or nx % 2 != 1:
            raise ValueError(
                "TESSCut FLUX must be an odd square NxN cutout with N >= 3"
            )
        target_pixel = np.asarray(
            wcs.all_world2pix(
                [[WD1856_TESSCUT_RA_DEG, WD1856_TESSCUT_DEC_DEG]], 0
            ),
            dtype=float,
        )
        if target_pixel.shape != (1, 2) or not np.isfinite(target_pixel).all():
            raise ValueError("TESSCut WCS produced an invalid WD 1856 pixel position")
        pixel_x_float, pixel_y_float = target_pixel[0]
        pixel_x = int(np.rint(pixel_x_float))
        pixel_y = int(np.rint(pixel_y_float))
        if not (1 <= pixel_x < nx - 1 and 1 <= pixel_y < ny - 1):
            raise ValueError(
                "the WCS-nearest WD 1856 pixel cannot support a centered 3x3 aperture"
            )
        center_world = np.asarray(
            wcs.all_pix2world([[float(pixel_x), float(pixel_y)]], 0), dtype=float
        )
        if center_world.shape != (1, 2) or not np.isfinite(center_world).all():
            raise ValueError("TESSCut WCS produced an invalid center-pixel sky position")
        center_separation = _position_separation_arcsec(
            float(center_world[0, 0]),
            float(center_world[0, 1]),
            WD1856_TESSCUT_RA_DEG,
            WD1856_TESSCUT_DEC_DEG,
        )

        small_raw = flux[:, pixel_y, pixel_x]
        primary_raw = np.sum(
            flux[:, pixel_y - 1 : pixel_y + 2, pixel_x - 1 : pixel_x + 2],
            axis=(1, 2),
            dtype=np.float64,
        )
        if not np.isfinite(small_raw).all() or not np.isfinite(primary_raw).all():
            raise ValueError("TESSCut aperture sums contain nonfinite values")
        return {
            "time_btjd": time,
            "quality": quality,
            "small_raw": small_raw,
            "primary_raw": primary_raw,
            "shape": (int(ny), int(nx)),
            "pixel_x": pixel_x,
            "pixel_y": pixel_y,
            "pixel_x_float": float(pixel_x_float),
            "pixel_y_float": float(pixel_y_float),
            "center_separation_arcsec": center_separation,
            "header_position_separation_arcsec": position_separations,
            "origin": origin,
            "creator": creator,
            "procver": [
                str(value) for value in _header_values(headers, "PROCVER")
            ],
            "flux_column_number": flux_column_number,
        }


def _read_compact_cadence_authority(path: Path) -> dict[str, np.ndarray]:
    target_key = f"{WD1856_TIC:016d}"
    with h5py.File(path, "r") as h5:
        if int(h5.attrs.get("sector", -1)) != SECTOR:
            raise ValueError("compact cadence authority is not Sector 56")
        time_unit = str(h5.attrs.get("time_unit", "")).replace(" ", "").lower()
        if time_unit != "bjd-2457000":
            raise ValueError("compact cadence authority must use BJD - 2457000")
        if "targets" not in h5 or target_key not in h5["targets"]:
            raise KeyError(f"compact cadence authority lacks WD 1856 TIC {WD1856_TIC}")
        group = h5["targets"][target_key]
        identity = {
            "tic": int(group.attrs.get("tic", -1)),
            "sector": int(group.attrs.get("sector", -1)),
            "camera": int(group.attrs.get("camera", -1)),
            "ccd": int(group.attrs.get("ccd", -1)),
        }
        if identity != {
            "tic": WD1856_TIC,
            "sector": SECTOR,
            "camera": CAMERA,
            "ccd": CCD,
        }:
            raise ValueError(f"compact WD 1856 identity is invalid: {identity}")
        missing = [
            name
            for name in ("time", "cadenceno", "quality", "orbitid")
            if name not in group
        ]
        if missing:
            raise KeyError(f"compact cadence authority lacks datasets: {missing}")
        time = np.asarray(group["time"], dtype=np.float64)
        cadence_raw = np.asarray(group["cadenceno"])
        quality_raw = np.asarray(group["quality"])
        orbit_raw = np.asarray(group["orbitid"])

    cadence_float = np.asarray(cadence_raw, dtype=np.float64)
    quality_float = np.asarray(quality_raw, dtype=np.float64)
    orbit_float = np.asarray(orbit_raw, dtype=np.float64)
    lengths = {len(time), len(cadence_float), len(quality_float), len(orbit_float)}
    if (
        time.ndim != 1
        or cadence_float.ndim != 1
        or quality_float.ndim != 1
        or orbit_float.ndim != 1
        or len(lengths) != 1
    ):
        raise ValueError("compact cadence-authority arrays have invalid shapes")
    if not len(time) or not np.isfinite(time).all():
        raise ValueError("compact TIME must be a nonempty finite array")
    if not np.isfinite(cadence_float).all() or not np.equal(
        cadence_float, np.floor(cadence_float)
    ).all():
        raise ValueError("compact CADENCENO must contain finite integers")
    cadence = cadence_float.astype(np.int64)
    for values, label in (
        (quality_float, "QUALITY"),
        (orbit_float, "ORBITID"),
    ):
        if not np.isfinite(values).all() or not np.equal(
            values, np.floor(values)
        ).all():
            raise ValueError(f"compact {label} must contain finite integers")
    quality = quality_float.astype(np.int64)
    orbitid = orbit_float.astype(np.int64)
    if (quality < 0).any() or (orbitid < 0).any():
        raise ValueError("compact QUALITY and ORBITID must be nonnegative")
    if len(np.unique(time)) != len(time):
        raise ValueError("compact WD 1856 TIME contains duplicate timestamps")
    if len(np.unique(cadence)) != len(cadence):
        raise ValueError("compact WD 1856 CADENCENO contains duplicates")
    return {
        "time_btjd": time,
        "cadenceno": cadence,
        "quality": quality,
        "orbitid": orbitid,
    }


def _nearest_unique_cadence_map(
    tesscut_time: np.ndarray,
    compact_time: np.ndarray,
    compact_cadence: np.ndarray,
    *,
    max_time_delta_seconds: float,
) -> dict[str, np.ndarray]:
    order = np.argsort(compact_time, kind="stable")
    sorted_time = compact_time[order]
    insertion = np.searchsorted(sorted_time, tesscut_time, side="left")
    right = np.clip(insertion, 0, len(sorted_time) - 1)
    left = np.clip(insertion - 1, 0, len(sorted_time) - 1)
    left_delta = np.abs(tesscut_time - sorted_time[left])
    right_delta = np.abs(tesscut_time - sorted_time[right])
    chosen_sorted = np.where(left_delta <= right_delta, left, right)
    chosen = order[chosen_sorted]
    matched_time = compact_time[chosen]
    matched_cadence = compact_cadence[chosen]
    delta_seconds = np.abs(tesscut_time - matched_time) * 86_400.0
    if not np.isfinite(delta_seconds).all():
        raise ValueError("timestamp matching produced nonfinite offsets")
    if len(np.unique(chosen)) != len(chosen) or len(np.unique(matched_cadence)) != len(
        matched_cadence
    ):
        raise ValueError("TESSCut timestamps do not map one-to-one to compact cadences")
    maximum = float(np.max(delta_seconds))
    if maximum > max_time_delta_seconds:
        raise ValueError(
            f"maximum timestamp mismatch is {maximum:.6g} s, above "
            f"{max_time_delta_seconds:.6g} s"
        )
    return {
        "compact_time_btjd": matched_time,
        "cadenceno": matched_cadence,
        "delta_seconds": delta_seconds,
        "compact_index": chosen,
    }


def _cadence_map_sha256(
    tesscut_time: np.ndarray,
    compact_time: np.ndarray,
    cadence: np.ndarray,
) -> str:
    digest = hashlib.sha256()
    digest.update(BUILDER_VERSION.encode("utf-8"))
    digest.update(np.asarray([len(cadence)], dtype="<u8").tobytes())
    for values, dtype in (
        (tesscut_time, "<f8"),
        (compact_time, "<f8"),
        (cadence, "<i8"),
    ):
        digest.update(np.ascontiguousarray(values, dtype=dtype).tobytes())
    return digest.hexdigest()


def _centered_rolling_median_detrend_with_diagnostics(
    flux: np.ndarray, quality: np.ndarray, *, window_cadences: int
) -> tuple[np.ndarray, np.ndarray, int, int]:
    """Detrend and report the count of quality-masked trend fills."""

    if (
        isinstance(window_cadences, bool)
        or not isinstance(window_cadences, (int, np.integer))
        or int(window_cadences) < 3
        or int(window_cadences) % 2 != 1
    ):
        raise ValueError("rolling window must be an odd integer >= 3")
    values = np.asarray(flux, dtype=np.float64)
    flags = np.asarray(quality)
    if values.ndim != 1 or flags.ndim != 1 or len(values) != len(flags) or not len(values):
        raise ValueError("flux and quality must be nonempty, equal-length vectors")
    if not np.isfinite(values).all():
        raise ValueError("aperture flux must be entirely finite")
    minimum_periods = max(2, (int(window_cadences) + 3) // 4)
    baseline = np.where(flags == 0, values, np.nan)
    trend = (
        pd.Series(baseline, dtype="float64")
        .rolling(
            int(window_cadences),
            center=True,
            min_periods=minimum_periods,
        )
        .median()
        .to_numpy(dtype=np.float64)
    )
    unsupported = ~np.isfinite(trend)
    unsupported_on_usable = unsupported & (flags == 0)
    if unsupported_on_usable.any():
        raise ValueError(
            "quality-zero rolling median is undefined at one or more "
            "quality-zero cadences"
        )
    n_masked_trend_filled = int(np.count_nonzero(unsupported))
    if n_masked_trend_filled:
        supported_index = np.flatnonzero(~unsupported)
        if not len(supported_index):
            raise ValueError("quality-zero rolling median has no supported cadences")
        unsupported_index = np.flatnonzero(unsupported)
        trend = trend.copy()
        trend[unsupported_index] = np.interp(
            unsupported_index,
            supported_index,
            trend[supported_index],
        )
    if not np.isfinite(trend).all():
        raise RuntimeError("masked-cadence rolling-median interpolation failed")
    if (trend <= 0).any():
        raise ValueError("quality-zero rolling-median trend must remain positive")
    detrended = values / trend
    if not np.isfinite(detrended).all():
        raise ValueError("rolling-median detrending produced nonfinite flux")
    return detrended, trend, minimum_periods, n_masked_trend_filled


def centered_rolling_median_detrend(
    flux: np.ndarray, quality: np.ndarray, *, window_cadences: int
) -> tuple[np.ndarray, np.ndarray, int]:
    """Divide by a centered quality-zero rolling-median trend.

    Undefined trend samples may be interpolated only where the effective
    quality is nonzero.  The public three-value return contract is preserved;
    product builders use the internal diagnostic helper to record fill counts.
    """

    detrended, trend, minimum_periods, _ = (
        _centered_rolling_median_detrend_with_diagnostics(
            flux,
            quality,
            window_cadences=window_cadences,
        )
    )
    return detrended, trend, minimum_periods


def _float_summary(values: np.ndarray) -> dict[str, float]:
    data = np.asarray(values, dtype=float)
    return {
        "min": float(np.min(data)),
        "median": float(np.median(data)),
        "max": float(np.max(data)),
    }


def _build_output_hdul(
    *,
    time_btjd: np.ndarray,
    cadence: np.ndarray,
    quality: np.ndarray,
    small_flux: np.ndarray,
    primary_flux: np.ndarray,
    raw_sha256: str,
    compact_sha256: str,
    cadence_map_sha256: str,
    cadence_reference_table_sha256: str,
    cadence_reference_manifest_sha256: str,
    quality_policy_contract: str,
    source_url: str,
    rolling_window_cadences: int,
    rolling_min_periods: int,
    rolling_masked_fill_count: int,
    max_time_delta_seconds: float,
    position_tolerance_arcsec: float,
    pixel_x: int,
    pixel_y: int,
    created_utc: str,
) -> Any:
    from astropy.io import fits

    primary = fits.PrimaryHDU()
    identity: Mapping[str, Any] = {
        "TICID": (WD1856_TIC, "WD 1856 TIC identifier"),
        "SECTOR": (SECTOR, "TESS observing sector"),
        "CAMERA": (CAMERA, "TESS camera"),
        "CCD": (CCD, "TESS CCD"),
        "RA_OBJ": (WD1856_TESSCUT_RA_DEG, "locked target right ascension, deg"),
        "DEC_OBJ": (WD1856_TESSCUT_DEC_DEG, "locked target declination, deg"),
        "ORIGIN": ("TWIRL", "reference-product publisher"),
        "CREATOR": (BUILDER_VERSION, "builder contract"),
        "DATE": (created_utc, "UTC creation timestamp"),
        "RAWSHA": raw_sha256,
        "COMPSHA": compact_sha256,
        "MAPSHA": cadence_map_sha256,
        "QREFSHA": cadence_reference_table_sha256,
        "QMANISHA": cadence_reference_manifest_sha256,
        "QPOLICY": quality_policy_contract,
        "SRCURL": (source_url, "official MAST TESSCut source URL"),
        "ROLLWIN": (rolling_window_cadences, "centered rolling-median window"),
        "ROLLMIN": (rolling_min_periods, "minimum quality-zero baseline rows"),
        "ROLLFILL": (rolling_masked_fill_count, "quality-masked trend fills"),
        "ROLLPOL": MASKED_TREND_FILL_HEADER,
        "FILLUSE": (False, "filled trend rows are not science-eligible"),
        "MAXDTS": (max_time_delta_seconds, "allowed map mismatch, seconds"),
        "POSTOL": (position_tolerance_arcsec, "target-position tolerance, arcsec"),
        "CENTX": (pixel_x, "zero-based central aperture x pixel"),
        "CENTY": (pixel_y, "zero-based central aperture y pixel"),
        "APERSML": ("1x1 central pixel", "small reference aperture"),
        "APERBIG": ("3x3 centered sum", "primary reference aperture"),
        "FLUXSRC": ("TESSCut FLUX pixels", "photometric pixel authority"),
        "DETREND": ("centered quality-zero rolling median division", "method"),
        "CURFLUX": (False, "compact/TGLC flux was not read or used"),
    }
    for keyword, value in identity.items():
        primary.header[keyword] = value

    n_rows = len(time_btjd)
    table = fits.BinTableHDU.from_columns(
        [
            fits.Column(
                name="TICID",
                format="K",
                array=np.full(n_rows, WD1856_TIC, dtype=np.int64),
            ),
            fits.Column(
                name="SECTOR",
                format="I",
                array=np.full(n_rows, SECTOR, dtype=np.int16),
            ),
            fits.Column(name="CADENCENO", format="K", array=cadence.astype(np.int64)),
            fits.Column(name="TIME", format="D", unit="d", array=time_btjd),
            fits.Column(name="QUALITY", format="J", array=quality.astype(np.int32)),
            fits.Column(name="TESSCUT_FLUX_SML", format="D", array=small_flux),
            fits.Column(name="TESSCUT_FLUX", format="D", array=primary_flux),
        ],
        name="LIGHTCURVE",
    )
    table.header["TICID"] = WD1856_TIC
    table.header["SECTOR"] = SECTOR
    table.header["CAMERA"] = CAMERA
    table.header["CCD"] = CCD
    table.header["TIMESYS"] = "TDB"
    table.header["BJDREFI"] = 2457000
    table.header["BJDREFF"] = 0.0
    table.header["TIMEUNIT"] = "d"
    table.header["TREFPOS"] = "BARYCENTER"
    return fits.HDUList([primary, table])


def _temporary_path(destination: Path, suffix: str) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        prefix=f".{destination.name}.", suffix=suffix, dir=destination.parent, delete=False
    )
    handle.close()
    os.chmod(handle.name, 0o640)
    return Path(handle.name)


def _publish_pair(
    *,
    temporary_fits: Path,
    temporary_json: Path,
    output_fits: Path,
    manifest_json: Path,
    expected_fits_sha256: str,
    source_check: Callable[[], None],
) -> None:
    destinations = (output_fits, manifest_json)
    backups: dict[Path, Path] = {}
    published: list[Path] = []
    try:
        source_check()
        for destination in destinations:
            if destination.exists():
                backup = destination.with_name(
                    f".{destination.name}.backup-{uuid.uuid4().hex}"
                )
                os.replace(destination, backup)
                backups[destination] = backup
        for temporary, destination in (
            (temporary_fits, output_fits),
            (temporary_json, manifest_json),
        ):
            os.replace(temporary, destination)
            published.append(destination)
        source_check()
        if file_sha256(output_fits) != expected_fits_sha256:
            raise RuntimeError("published FITS hash differs from the staged product")
        json.loads(manifest_json.read_text(encoding="utf-8"))
    except Exception:
        for destination in reversed(published):
            if destination.exists():
                destination.unlink()
        for destination, backup in backups.items():
            if backup.exists():
                os.replace(backup, destination)
        raise
    else:
        for backup in backups.values():
            if backup.exists():
                backup.unlink()


def build_wd1856_tesscut_independent_extraction(
    *,
    tesscut_fits: Path,
    compact_lc: Path,
    cadence_reference_table: Path,
    cadence_reference_manifest: Path,
    output_fits: Path,
    manifest_json: Path,
    source_url: str,
    rolling_window_cadences: int = DEFAULT_ROLLING_WINDOW_CADENCES,
    position_tolerance_arcsec: float = DEFAULT_POSITION_TOLERANCE_ARCSEC,
    max_time_delta_seconds: float = DEFAULT_MAX_TIME_DELTA_SECONDS,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Build and publish the independent pair with rollback on errors."""

    source_url = _validate_source_url(source_url)
    floating = np.asarray(
        [position_tolerance_arcsec, max_time_delta_seconds], dtype=float
    )
    if not np.isfinite(floating).all() or (floating <= 0).any():
        raise ValueError("position and timestamp tolerances must be finite and positive")
    if position_tolerance_arcsec > DEFAULT_POSITION_TOLERANCE_ARCSEC:
        raise ValueError("position tolerance cannot exceed the locked 1 arcsec limit")
    if max_time_delta_seconds > DEFAULT_MAX_TIME_DELTA_SECONDS:
        raise ValueError("timestamp tolerance cannot exceed the locked 60 second limit")
    if (
        isinstance(rolling_window_cadences, bool)
        or not isinstance(rolling_window_cadences, (int, np.integer))
        or int(rolling_window_cadences) < 3
        or int(rolling_window_cadences) % 2 != 1
    ):
        raise ValueError("rolling window must be an odd integer >= 3")

    inputs = (
        Path(tesscut_fits).resolve(),
        Path(compact_lc).resolve(),
        Path(cadence_reference_table).resolve(),
        Path(cadence_reference_manifest).resolve(),
    )
    output_fits = Path(output_fits).resolve()
    manifest_json = Path(manifest_json).resolve()
    if output_fits.suffix.lower() not in {".fits", ".fit", ".fts"}:
        raise ValueError("independent reference output must be a FITS file")
    if manifest_json.suffix.lower() != ".json":
        raise ValueError("independent reference sidecar must be a JSON file")
    if output_fits == manifest_json:
        raise ValueError("FITS and JSON outputs must be distinct")
    collisions = {str(path) for path in inputs if path in {output_fits, manifest_json}}
    if collisions:
        raise ValueError(f"inputs and outputs collide: {sorted(collisions)}")
    if not overwrite:
        existing = [str(path) for path in (output_fits, manifest_json) if path.exists()]
        if existing:
            raise FileExistsError(f"outputs already exist; pass overwrite=True: {existing}")

    states = tuple(_capture_source_state(path) for path in inputs)
    (
        tesscut_state,
        compact_state,
        cadence_reference_state,
        cadence_manifest_state,
    ) = states
    tesscut = _read_tesscut(
        tesscut_state.path, tolerance_arcsec=float(position_tolerance_arcsec)
    )
    compact = _read_compact_cadence_authority(compact_state.path)
    quality_reference = load_external_quality_reference(
        table_path=cadence_reference_state.path,
        manifest_path=cadence_manifest_state.path,
        sector=SECTOR,
        expected_orbits=(119, 120),
    )
    compact_overlay = quality_reference.apply(
        sector=SECTOR,
        camera=CAMERA,
        ccd=CCD,
        cadenceno=compact["cadenceno"],
        orbitid=compact["orbitid"],
        internal_quality=compact["quality"],
        context="WD 1856 compact cadence authority",
    )
    cadence_map = _nearest_unique_cadence_map(
        tesscut["time_btjd"],
        compact["time_btjd"],
        compact["cadenceno"],
        max_time_delta_seconds=float(max_time_delta_seconds),
    )
    cadence_map_sha256 = _cadence_map_sha256(
        tesscut["time_btjd"],
        cadence_map["compact_time_btjd"],
        cadence_map["cadenceno"],
    )
    mapped_compact_index = cadence_map["compact_index"]
    tesscut_overlay = quality_reference.apply(
        sector=SECTOR,
        camera=CAMERA,
        ccd=CCD,
        cadenceno=cadence_map["cadenceno"],
        orbitid=compact["orbitid"][mapped_compact_index],
        internal_quality=tesscut["quality"],
        context="WD 1856 mapped TESSCut reference",
    )
    effective_tesscut_quality = tesscut_overlay.quality
    (
        small_flux,
        small_trend,
        rolling_min_periods,
        small_masked_trend_fill_count,
    ) = _centered_rolling_median_detrend_with_diagnostics(
        tesscut["small_raw"],
        effective_tesscut_quality,
        window_cadences=int(rolling_window_cadences),
    )
    (
        primary_flux,
        primary_trend,
        primary_min_periods,
        primary_masked_trend_fill_count,
    ) = _centered_rolling_median_detrend_with_diagnostics(
        tesscut["primary_raw"],
        effective_tesscut_quality,
        window_cadences=int(rolling_window_cadences),
    )
    if primary_min_periods != rolling_min_periods:
        raise RuntimeError("aperture detrending used inconsistent baseline minima")
    if primary_masked_trend_fill_count != small_masked_trend_fill_count:
        raise RuntimeError("aperture detrending used inconsistent masked trend fills")

    _assert_sources_unchanged(states)
    quality_reference.assert_unchanged()
    created_utc = datetime.now(timezone.utc).isoformat()
    temporary_fits = _temporary_path(output_fits, ".fits")
    temporary_json = _temporary_path(manifest_json, ".json")
    try:
        hdul = _build_output_hdul(
            time_btjd=tesscut["time_btjd"],
            cadence=cadence_map["cadenceno"],
            quality=effective_tesscut_quality,
            small_flux=small_flux,
            primary_flux=primary_flux,
            raw_sha256=tesscut_state.sha256,
            compact_sha256=compact_state.sha256,
            cadence_map_sha256=cadence_map_sha256,
            cadence_reference_table_sha256=quality_reference.table_sha256,
            cadence_reference_manifest_sha256=quality_reference.manifest_sha256,
            quality_policy_contract=quality_reference.provenance[
                "policy_contract"
            ],
            source_url=source_url,
            rolling_window_cadences=int(rolling_window_cadences),
            rolling_min_periods=rolling_min_periods,
            rolling_masked_fill_count=small_masked_trend_fill_count,
            max_time_delta_seconds=float(max_time_delta_seconds),
            position_tolerance_arcsec=float(position_tolerance_arcsec),
            pixel_x=tesscut["pixel_x"],
            pixel_y=tesscut["pixel_y"],
            created_utc=created_utc,
        )
        try:
            hdul.writeto(temporary_fits, overwrite=True, checksum=True)
        finally:
            hdul.close()
        output_sha256 = file_sha256(temporary_fits)
        delta = cadence_map["delta_seconds"]
        manifest: dict[str, Any] = {
            "builder_version": BUILDER_VERSION,
            "created_utc": created_utc,
            "target": {
                "tic": WD1856_TIC,
                "sector": SECTOR,
                "camera": CAMERA,
                "ccd": CCD,
                "ra_deg": WD1856_TESSCUT_RA_DEG,
                "dec_deg": WD1856_TESSCUT_DEC_DEG,
                "position_tolerance_arcsec": float(position_tolerance_arcsec),
            },
            "tesscut_source": {
                "path": str(tesscut_state.path),
                "source_url": source_url,
                "sha256": tesscut_state.sha256,
                "size_bytes": tesscut_state.size_bytes,
                "mtime_ns": tesscut_state.mtime_ns,
                "origin": tesscut["origin"],
                "creator": tesscut["creator"],
                "procver": tesscut["procver"],
                "flux_cube_shape": [len(tesscut["time_btjd"]), *tesscut["shape"]],
                "flux_column_number": tesscut["flux_column_number"],
                "header_position_separation_arcsec": tesscut[
                    "header_position_separation_arcsec"
                ],
            },
            "compact_cadence_authority": {
                "path": str(compact_state.path),
                "sha256": compact_state.sha256,
                "size_bytes": compact_state.size_bytes,
                "mtime_ns": compact_state.mtime_ns,
                "target_group": f"targets/{WD1856_TIC:016d}",
                "datasets_read": ["time", "cadenceno", "quality", "orbitid"],
                "flux_read_or_used": False,
                "n_available_cadences": int(len(compact["cadenceno"])),
                "cadence_map_sha256": cadence_map_sha256,
                "cadence_map_hash_algorithm": CADENCE_MAP_HASH_ALGORITHM,
            },
            "external_quality_overlay": {
                **quality_reference.provenance,
                "applied_before_reference_detrending": True,
                "compact_full_audit_counts": dict(compact_overlay.counts),
                "tesscut_mapped_audit_counts": dict(tesscut_overlay.counts),
                "tesscut_native_nonzero_quality": int(
                    np.count_nonzero(tesscut["quality"])
                ),
                "output_quality_definition": (
                    "binary (TESSCut native QUALITY != 0) OR "
                    "(authoritative external_quality != 0)"
                ),
            },
            "extraction": {
                "pixel_source": "official MAST TESSCut calibrated FLUX cube",
                "center_pixel_zero_based": {
                    "x": tesscut["pixel_x"],
                    "y": tesscut["pixel_y"],
                    "wcs_x": tesscut["pixel_x_float"],
                    "wcs_y": tesscut["pixel_y_float"],
                    "sky_separation_arcsec": tesscut[
                        "center_separation_arcsec"
                    ],
                },
                "apertures": {
                    "TESSCUT_FLUX_SML": "single WCS-nearest central pixel",
                    "TESSCUT_FLUX": "centered 3x3 sum",
                },
                "detrending": {
                    "method": "centered rolling median division",
                    "baseline_mask": (
                        "effective quality == 0 after TESSCut native OR "
                        "authoritative external overlay"
                    ),
                    "window_cadences": int(rolling_window_cadences),
                    "minimum_baseline_cadences": rolling_min_periods,
                    "edge_policy": "centered partial windows with the same minimum",
                    "masked_unsupported_trend_policy": MASKED_TREND_FILL_POLICY,
                    "n_masked_trend_filled": small_masked_trend_fill_count,
                    "filled_samples_science_eligible": False,
                },
            },
            "diagnostics": {
                "n_rows": int(len(tesscut["time_btjd"])),
                "n_quality_zero": int(
                    np.count_nonzero(effective_tesscut_quality == 0)
                ),
                "quality_zero_fraction": float(
                    np.mean(effective_tesscut_quality == 0)
                ),
                "time_delta_seconds": {
                    "median": float(np.median(delta)),
                    "p95": float(np.quantile(delta, 0.95)),
                    "max": float(np.max(delta)),
                },
                "raw_1x1": _float_summary(tesscut["small_raw"]),
                "raw_3x3": _float_summary(tesscut["primary_raw"]),
                "trend_1x1": _float_summary(small_trend),
                "trend_3x3": _float_summary(primary_trend),
                "detrended_1x1": _float_summary(small_flux),
                "detrended_3x3": _float_summary(primary_flux),
            },
            "output": {
                "path": str(output_fits),
                "sha256": output_sha256,
                "n_rows": int(len(tesscut["time_btjd"])),
                "columns": [
                    "TICID",
                    "SECTOR",
                    "CADENCENO",
                    "TIME",
                    "QUALITY",
                    "TESSCUT_FLUX_SML",
                    "TESSCUT_FLUX",
                ],
            },
        }
        temporary_json.write_text(
            json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        _publish_pair(
            temporary_fits=temporary_fits,
            temporary_json=temporary_json,
            output_fits=output_fits,
            manifest_json=manifest_json,
            expected_fits_sha256=output_sha256,
            source_check=lambda: (
                _assert_sources_unchanged(states),
                quality_reference.assert_unchanged(),
            ),
        )
    finally:
        for temporary in (temporary_fits, temporary_json):
            if temporary.exists():
                temporary.unlink()
    return manifest


__all__ = [
    "BUILDER_VERSION",
    "CADENCE_MAP_HASH_ALGORITHM",
    "DEFAULT_MAX_TIME_DELTA_SECONDS",
    "DEFAULT_POSITION_TOLERANCE_ARCSEC",
    "DEFAULT_ROLLING_WINDOW_CADENCES",
    "WD1856_TESSCUT_DEC_DEG",
    "WD1856_TESSCUT_RA_DEG",
    "build_wd1856_tesscut_independent_extraction",
    "centered_rolling_median_detrend",
]
