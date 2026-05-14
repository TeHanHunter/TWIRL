"""HLSP (`hlsp_qlp_tess_ffi_*_v01_llc.fits`) reader for TWIRL Stage 2.

Centralizes the load/quality/aperture pattern that previously lived in
`reports/exploratory/_scripts/plot_wd1856.py` and `scripts/stage1_lightcurves/qc_plot_hlsp_sample.py`.

HLSP files are produced by `qlp lctools hlsp` at the sector level (vstack of
all orbits in the sector). Per-target schema:
    TIME, CADENCENO, SAP_FLUX, DET_FLUX, DET_FLUX_ERR, QUALITY, ORBITID,
    SAP_X, SAP_Y, SAP_BKG, SAP_BKG_ERR, DET_FLUX_SML, DET_FLUX_LAG,
    SYS_RM_FLUX

DET_FLUX is QLP-detrended (BSpline cotrending). TIME is BJD - 2457000.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Sequence

import numpy as np
from astropy.io import fits

APERTURES: tuple[str, ...] = ("DET_FLUX_SML", "DET_FLUX", "DET_FLUX_LAG")
BJDREFI: int = 2457000


@dataclass
class HLSPLightCurve:
    tic: int
    tmag: float
    sector: int
    cam: int
    ccd: int
    ra: float
    dec: float
    time: np.ndarray
    cadenceno: np.ndarray
    orbitid: np.ndarray
    quality: np.ndarray
    flux: dict[str, np.ndarray]
    path: Path


def iter_hlsp_fits(hlsp_root: Path) -> Iterator[Path]:
    """Yield every `hlsp_{qlp,twirl}_tess_ffi_*.fits` under `hlsp_root`, sorted.

    Matches both the QLP-produced HLSPs (`hlsp_qlp_*`) and the TWIRL v3
    flux-space-detrended HLSPs (`hlsp_twirl_*`) so downstream code reads
    whichever tree it's pointed at without changes.
    """
    yield from sorted(Path(hlsp_root).rglob("hlsp_*_tess_ffi_*.fits"))


def discover_sector_targets(hlsp_root: Path, sector: int) -> list[Path]:
    """All HLSP files for `sector` under `hlsp_root`.

    Files are named `hlsp_{qlp,twirl}_tess_ffi_s{sector:04d}-{tic:016d}_tess_v01_llc.fits`;
    we filter by the `s{NNNN}` prefix to avoid loading unrelated sectors that
    might co-exist in the tree. The wildcard between `hlsp_` and `_tess_ffi_`
    matches both QLP (`qlp`) and TWIRL v3 (`twirl`) HLSP variants.
    """
    pat = f"hlsp_*_tess_ffi_s{sector:04d}-*.fits"
    return sorted(Path(hlsp_root).rglob(pat))


def read_hlsp(
    path: Path,
    columns: Sequence[str] = APERTURES,
) -> HLSPLightCurve | None:
    """Read one HLSP FITS into an `HLSPLightCurve`. Returns None on failure.

    `columns` chooses which aperture flux columns to load. SAP_FLUX is also
    loaded if requested. Returns None when the file cannot be opened or the
    table is empty.
    """
    try:
        with fits.open(path, memmap=False) as hdul:
            hdr = hdul[0].header
            tab = hdul[1].data
            n = len(tab)
            if n == 0:
                return None
            tic = int(hdr.get("TICID", -1))
            tmag = float(hdr.get("TESSMAG", np.nan))
            sector = int(hdr.get("SECTOR", -1))
            cam = int(hdr.get("CAMERA", -1))
            ccd = int(hdr.get("CCD", -1))
            ra = float(hdr.get("RA_OBJ", np.nan))
            dec = float(hdr.get("DEC_OBJ", np.nan))
            time = np.asarray(tab["TIME"], dtype=np.float64)
            cadenceno = np.asarray(
                tab["CADENCENO"] if "CADENCENO" in tab.columns.names else np.arange(n)
            )
            orbitid = np.asarray(
                tab["ORBITID"] if "ORBITID" in tab.columns.names else np.zeros(n, dtype=np.int32)
            )
            quality = np.asarray(
                tab["QUALITY"] if "QUALITY" in tab.columns.names else np.zeros(n, dtype=np.int32)
            )
            flux: dict[str, np.ndarray] = {}
            for c in columns:
                if c in tab.columns.names:
                    flux[c] = np.asarray(tab[c], dtype=np.float64)
        return HLSPLightCurve(
            tic=tic, tmag=tmag, sector=sector, cam=cam, ccd=ccd,
            ra=ra, dec=dec,
            time=time, cadenceno=cadenceno, orbitid=orbitid, quality=quality,
            flux=flux, path=Path(path),
        )
    except Exception:
        return None


def quality_mask(lc: HLSPLightCurve, aperture: str = "DET_FLUX") -> np.ndarray:
    """Boolean mask: QUALITY==0 AND finite flux for the named aperture.

    QUALITY bit semantics (from `qlp.lctools.bin.hlsp:424`): SPOC bits 0..15,
    TGLC bit 29, QLP bit 30. We require QUALITY==0 (any flagged cadence is
    excluded), matching the convention in the Stage 1 QA scripts.
    """
    if aperture not in lc.flux:
        return np.zeros(len(lc.time), dtype=bool)
    f = lc.flux[aperture]
    return (lc.quality == 0) & np.isfinite(f) & np.isfinite(lc.time)


def tglc_mad_error(
    lc: HLSPLightCurve,
    aperture: str = "DET_FLUX",
) -> np.ndarray:
    """Canonical TWIRL per-cadence flux error: 1.4826 × MAD of the
    good-quality detrended flux, broadcast across all good cadences.

    This is the TGLC-recommended estimator and should be used unconditionally
    — many QLP/HLSP FITS in our production tree (S56 included) ship with
    ``DET_FLUX_ERR`` as all-NaN due to an upstream wrapper issue, so the
    column is unreliable. The MAD method is robust to outliers (transit dips
    and flags) and matches the TGLC photometric noise model.

    Returns a same-length array as ``lc.time`` with ``sigma`` filled at every
    good-quality cadence and ``NaN`` elsewhere, so it can be used as
    per-cadence ``flux_err`` directly.
    """
    if aperture not in lc.flux:
        raise KeyError(f"aperture {aperture!r} not loaded in HLSPLightCurve")
    f = lc.flux[aperture]
    keep = quality_mask(lc, aperture)
    fg = f[keep]
    if fg.size == 0:
        return np.full_like(f, np.nan, dtype=np.float64)
    med = np.nanmedian(fg)
    mad = np.nanmedian(np.abs(fg - med))
    sigma = 1.4826 * float(mad)
    out = np.full_like(f, np.nan, dtype=np.float64)
    out[keep] = sigma
    return out
