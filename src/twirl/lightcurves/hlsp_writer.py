"""TWIRL-FS HLSP FITS writer.

Emits a 2-HDU FITS file (PRIMARY metadata + LIGHTCURVE BinTable) in the
column schema produced by QLP ``lctools hlsp``, so that the existing
:mod:`twirl.io.hlsp` reader and all downstream code (BLS, heuristic
vetter, LEO adapter) consume TWIRL-FS HLSPs unchanged.

Difference from QLP HLSPs:

* ``SAP_FLUX`` / ``DET_FLUX`` are emitted in **linear relative-flux
  units** using the TWIRL-FS robust positive flux scale, not
  magnitudes-converted-back. Cadences where QLP would have silently
  dropped a measurement (``flux <= 0`` pre-detrend) are preserved here
  as real negative values.
* ``DET_FLUX_ERR`` is the MAD-broadcast per-cadence sigma from our
  flux-space detrender, not propagated through a magnitude space.

Filename convention:
``hlsp_twirlfs_tess_ffi_s<sector>-<tic16>_tess_v01_llc.fits``.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits

from .flux_detrend import FluxDetrendConfig


PRODUCT_PREFIX = "twirlfs"
PIPELINE_NAME = "TWIRL-FS"
METHOD_VERSION = "twirl-fs-v2"

HLSP_COLUMN_SCHEMA: tuple[tuple[str, str, str | None], ...] = (
    ("TIME",         "D", "BJD-2457000, days"),
    ("CADENCENO",    "J", None),
    ("SAP_FLUX",     "E", None),
    ("DET_FLUX",     "E", None),
    ("DET_FLUX_ERR", "E", None),
    ("QUALITY",      "J", None),
    ("ORBITID",      "J", None),
    ("SAP_X",        "E", "pixel"),
    ("SAP_Y",        "E", "pixel"),
    ("SAP_BKG",      "E", None),
    ("SAP_BKG_ERR",  "E", None),
    ("DET_FLUX_SML", "E", None),
    ("DET_FLUX_LAG", "E", None),
    ("SYS_RM_FLUX",  "E", None),
)


@dataclass
class HLSPTarget:
    """Per-target metadata for the PRIMARY HDU header."""
    tic: int
    sector: int
    cam: int
    ccd: int
    tmag: float
    ra: float
    dec: float
    version: str = "v01"


def hlsp_filename(target: HLSPTarget) -> str:
    """``hlsp_twirlfs_tess_ffi_s<sector>-<tic16>_tess_<version>_llc.fits``."""
    return (
        f"hlsp_{PRODUCT_PREFIX}_tess_ffi_s{target.sector:04d}-{target.tic:016d}"
        f"_tess_{target.version}_llc.fits"
    )


def hlsp_path(out_root: Path, target: HLSPTarget) -> Path:
    """Mirror the QLP HLSP path convention (4-deep TIC-prefix directories)."""
    s = f"{target.tic:016d}"
    return (
        Path(out_root) / s[0:4] / s[4:8] / s[8:12] / s[12:16]
        / hlsp_filename(target)
    )


def write_twirl_hlsp(
    out_path: Path,
    target: HLSPTarget,
    *,
    time_btjd: np.ndarray,
    cadenceno: np.ndarray,
    sap_flux: np.ndarray,        # raw (pre-detrend) robust-scale flux
    det_flux: np.ndarray,        # detrended robust-scale flux
    det_flux_err: np.ndarray,    # MAD-broadcast per-cadence sigma
    quality: np.ndarray,
    orbitid: np.ndarray,
    sap_x: np.ndarray,
    sap_y: np.ndarray,
    sap_bkg: np.ndarray,
    sap_bkg_err: np.ndarray,
    det_flux_sml: np.ndarray,    # detrended small-aperture flux
    det_flux_lag: np.ndarray,    # detrended large-aperture flux
    sys_rm_flux: np.ndarray | None = None,   # optional systematics-removed
    detrend_config: FluxDetrendConfig | None = None,
    detrend_diagnostics: Mapping[str, object] | None = None,
    method_version: str = METHOD_VERSION,
    extra_flux_columns: Mapping[str, np.ndarray] | None = None,
    extra_header: Mapping[str, object] | None = None,
    include_canonical_det_flux: bool = True,
    include_sys_rm_flux: bool = True,
) -> Path:
    """Write a single TWIRL-FS HLSP FITS file in QLP-compatible schema.

    Negative values in any FLUX column are preserved as real numbers.
    """
    detrend_config = detrend_config or FluxDetrendConfig()
    out_path = Path(out_path)
    arrays: dict[str, np.ndarray] = {
        "time_btjd": np.asarray(time_btjd),
        "cadenceno": np.asarray(cadenceno),
        "sap_flux": np.asarray(sap_flux),
        "det_flux": np.asarray(det_flux),
        "det_flux_err": np.asarray(det_flux_err),
        "quality": np.asarray(quality),
        "orbitid": np.asarray(orbitid),
        "sap_x": np.asarray(sap_x),
        "sap_y": np.asarray(sap_y),
        "sap_bkg": np.asarray(sap_bkg),
        "sap_bkg_err": np.asarray(sap_bkg_err),
        "det_flux_sml": np.asarray(det_flux_sml),
        "det_flux_lag": np.asarray(det_flux_lag),
    }
    time_values = arrays["time_btjd"]
    if time_values.ndim != 1:
        raise ValueError(
            f"FITS column 'time_btjd' must be one-dimensional; got shape {time_values.shape}"
        )
    n = len(time_values)
    for name, values in arrays.items():
        if values.ndim != 1:
            raise ValueError(
                f"FITS column {name!r} must be one-dimensional; got shape {values.shape}"
            )
        if len(values) != n:
            raise ValueError(
                f"FITS column {name!r} has length {len(values)}; expected {n}"
            )

    validated_extra: dict[str, np.ndarray] = {}
    if extra_flux_columns:
        for name, values in extra_flux_columns.items():
            values = np.asarray(values)
            if values.ndim != 1:
                raise ValueError(
                    f"extra FITS column {name!r} must be one-dimensional; "
                    f"got shape {values.shape}"
                )
            if len(values) != n:
                raise ValueError(
                    f"extra FITS column {name!r} has length {len(values)}; expected {n}"
                )
            validated_extra[str(name)] = values

    if include_sys_rm_flux:
        if sys_rm_flux is None:
            sys_rm_values = np.full(n, np.nan, dtype=np.float32)
        else:
            sys_rm_values = np.asarray(sys_rm_flux)
            if sys_rm_values.ndim != 1:
                raise ValueError(
                    "FITS column 'sys_rm_flux' must be one-dimensional; "
                    f"got shape {sys_rm_values.shape}"
                )
            if len(sys_rm_values) != n:
                raise ValueError(
                    f"FITS column 'sys_rm_flux' has length {len(sys_rm_values)}; expected {n}"
                )
    else:
        sys_rm_values = None

    time_btjd = arrays["time_btjd"]
    cadenceno = arrays["cadenceno"]
    sap_flux = arrays["sap_flux"]
    det_flux = arrays["det_flux"]
    det_flux_err = arrays["det_flux_err"]
    quality = arrays["quality"]
    orbitid = arrays["orbitid"]
    sap_x = arrays["sap_x"]
    sap_y = arrays["sap_y"]
    sap_bkg = arrays["sap_bkg"]
    sap_bkg_err = arrays["sap_bkg_err"]
    det_flux_sml = arrays["det_flux_sml"]
    det_flux_lag = arrays["det_flux_lag"]
    out_path.parent.mkdir(parents=True, exist_ok=True)

    cols = [
        fits.Column(name="TIME", format="D", unit="BJD-2457000, days", array=time_btjd),
        fits.Column(name="CADENCENO", format="J", array=cadenceno.astype(np.int32)),
        fits.Column(name="SAP_FLUX", format="E", array=sap_flux.astype(np.float32)),
    ]
    if include_canonical_det_flux:
        cols.extend([
            fits.Column(name="DET_FLUX", format="E", array=det_flux.astype(np.float32)),
            fits.Column(name="DET_FLUX_ERR", format="E", array=det_flux_err.astype(np.float32)),
        ])
    cols.extend([
        fits.Column(name="QUALITY", format="J", array=quality.astype(np.int32)),
        fits.Column(name="ORBITID", format="J", array=orbitid.astype(np.int32)),
        fits.Column(name="SAP_X", format="E", unit="pixel", array=sap_x.astype(np.float32)),
        fits.Column(name="SAP_Y", format="E", unit="pixel", array=sap_y.astype(np.float32)),
        fits.Column(name="SAP_BKG", format="E", array=sap_bkg.astype(np.float32)),
        fits.Column(name="SAP_BKG_ERR", format="E", array=sap_bkg_err.astype(np.float32)),
    ])
    if include_canonical_det_flux:
        cols.extend([
            fits.Column(name="DET_FLUX_SML", format="E", array=det_flux_sml.astype(np.float32)),
            fits.Column(name="DET_FLUX_LAG", format="E", array=det_flux_lag.astype(np.float32)),
        ])
    if include_sys_rm_flux:
        assert sys_rm_values is not None
        cols.append(
            fits.Column(name="SYS_RM_FLUX", format="E", array=sys_rm_values.astype(np.float32))
        )
    if validated_extra:
        for name, values in validated_extra.items():
            cols.append(
                fits.Column(
                    name=name,
                    format="E",
                    array=values.astype(np.float32),
                )
            )
    lc_hdu = fits.BinTableHDU.from_columns(cols, name="LIGHTCURVE")

    primary = fits.PrimaryHDU()
    h = primary.header
    h["TICID"]   = (int(target.tic), "TIC identifier")
    h["TESSMAG"] = (float(target.tmag), "TESS magnitude")
    h["SECTOR"]  = (int(target.sector), "TESS sector")
    h["CAMERA"]  = (int(target.cam), "TESS camera")
    h["CCD"]     = (int(target.ccd), "TESS CCD")
    # QLP convention uses RA_OBJ / DEC_OBJ for the primary-header sky position.
    # We write RA_OBJ / DEC_OBJ as the canonical keys so QC + reader code that
    # follows QLP conventions (e.g. qc_sector_pdf.py) finds them, and keep
    # RA / DEC as backward-compat aliases for earlier TWIRL-only consumers.
    h["RA_OBJ"]  = (float(target.ra), "Right ascension (deg)")
    h["DEC_OBJ"] = (float(target.dec), "Declination (deg)")
    h["RA"]      = (float(target.ra), "Right ascension (deg), alias of RA_OBJ")
    h["DEC"]     = (float(target.dec), "Declination (deg), alias of DEC_OBJ")
    h["ORIGIN"]  = ("TWIRL", "Producer pipeline")
    h["PIPELINE"] = (PIPELINE_NAME, "TWIRL Flux-Spline detrend")
    h["METHOD"] = (str(method_version), "Detrend method version")
    h["VERSION"] = (str(target.version), "HLSP version tag")
    h["BKSPACE"] = (float(detrend_config.bkspace_d), "Spline knot spacing (d)")
    h["BSPLK"] = (int(detrend_config.k), "BSpline degree")
    h["SIGCLIP"] = (float(detrend_config.sigma_clip), "Spline rejection threshold")
    h["MAXITER"] = (int(detrend_config.max_iter), "Maximum spline clip iterations")
    h["EDGEPAD"] = (float(detrend_config.edge_pad_d), "Knot edge padding (d)")
    h["GAPSPLIT"] = (float(detrend_config.gap_split_d), "Independent fit gap threshold (d)")
    h["KNOTSTR"] = (str(detrend_config.knot_strategy), "Spline knot placement strategy")
    h["OUTMODE"] = (str(detrend_config.output_mode), "Detrended output mode")
    h["SCALE"] = (str(detrend_config.scale_strategy), "Flux-scale strategy")
    h["MINSNR"] = (float(detrend_config.min_scale_snr), "Auto-scale median S/N floor")
    if detrend_diagnostics:
        h["NSEG"] = (int(detrend_diagnostics.get("n_segments", -1)), "Cotrend segment count")
        h["FITCNT"] = (int(detrend_diagnostics.get("fit_count", -1)), "Good cadences used for fit")
        h["SCALESRC"] = (
            str(detrend_diagnostics.get("scale_source", ""))[:68],
            "Subtractive scale source",
        )
        h["COTSTAT"] = (
            str(detrend_diagnostics.get("cotrend_status", ""))[:68],
            "Cotrend fit status",
        )
    if extra_header:
        for key, value in extra_header.items():
            key = str(key)
            if isinstance(value, tuple) and len(value) == 2:
                h[key] = (_fits_header_value(value[0]), str(value[1])[:68])
            else:
                h[key] = _fits_header_value(value)
    h["BJDREFI"] = (2457000, "Integer BJD epoch reference")

    fits.HDUList([primary, lc_hdu]).writeto(out_path, overwrite=True)
    return out_path


def _fits_header_value(value: object) -> object:
    """Return a FITS-header-safe scalar value."""
    if isinstance(value, np.generic):
        value = value.item()
    if value is None:
        return ""
    if isinstance(value, str):
        return value[:68]
    if isinstance(value, (bool, int, float)):
        return value
    return str(value)[:68]
