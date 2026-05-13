"""TGLC per-orbit HDF5 reader.

Reads the HDF5 files produced by the TGLC ``lightcurves`` stage. Two
schema variants are supported:

* **Patched (post twirl/preserve-negative-flux branch)** — has both
  ``RawFlux``/``RawFluxError`` AND ``RawMagnitude``/``RawMagnitudeError``
  for each aperture. ``RawFlux`` preserves negative cadences (the ~50% of
  T>=19 cadences that QLP would silently drop).
* **Legacy (pre-patch)** — magnitude only. The reader synthesizes a flux
  array from the magnitude on the fly (``flux = 10^(-0.4 * mag)``),
  losing the negative cadences. Cadences where the magnitude is NaN
  remain NaN in the synthesized flux. Useful for sectors that haven't
  been re-extracted yet so the downstream pipeline still works.

Schema (TGLC HDF5):
    LightCurve/
      BJD                                  shape=(N,) BTJD (BJD - 2457000)
      Cadence                              shape=(N,)
      QualityFlag                          shape=(N,)
      X, Y                                 shape=(N,) global centroid
      Background/Value, Background/Error   shape=(N,)
      AperturePhotometry/{Small,Primary,Large}Aperture/
        RawMagnitude, RawMagnitudeError    shape=(N,)
        RawFlux, RawFluxError              shape=(N,)   [patched only]
        X, Y                               shape=(N,)   per-aperture centroid
    [attrs] BJDoffset, CCD, Camera, Dec, Orbit, RA, Sector, "TIC ID", TessMag
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import h5py
import numpy as np


APERTURE_KEYS: tuple[str, ...] = ("Small", "Primary", "Large")
# Aperture-name → HLSP column suffix used by QLP and downstream consumers.
HLSP_NAME = {"Small": "SML", "Primary": "MED", "Large": "LAG"}


@dataclass
class TGLCAperture:
    name: str            # "Small" / "Primary" / "Large"
    raw_flux: np.ndarray
    raw_flux_err: np.ndarray
    raw_magnitude: np.ndarray
    raw_magnitude_err: np.ndarray
    centroid_x: np.ndarray
    centroid_y: np.ndarray
    flux_was_synthesized: bool = False  # True when read from a legacy schema


@dataclass
class TGLCLightCurve:
    tic: int
    sector: int
    orbit: int
    cam: int
    ccd: int
    tmag: float
    ra: float
    dec: float
    bjd_offset: int
    time: np.ndarray         # BTJD = BJD - bjd_offset
    cadence: np.ndarray
    quality: np.ndarray
    centroid_x: np.ndarray   # global
    centroid_y: np.ndarray
    background: np.ndarray
    background_err: np.ndarray
    apertures: dict[str, TGLCAperture] = field(default_factory=dict)
    path: Path | None = None

    def time_bjd(self) -> np.ndarray:
        """Absolute BJD (adds back the bjd_offset)."""
        return self.time + float(self.bjd_offset)


def read_tglc_h5(path: Path | str) -> TGLCLightCurve:
    """Read a TGLC per-orbit HDF5 into a :class:`TGLCLightCurve`.

    Falls back to synthesized flux from magnitude if ``RawFlux`` is absent
    (legacy schema). The ``flux_was_synthesized`` attribute on each
    :class:`TGLCAperture` records which path was used.
    """
    path = Path(path)
    with h5py.File(path, "r") as f:
        lc_attrs = dict(f.attrs)
        lc = f["LightCurve"]
        bjd = np.asarray(lc["BJD"][:], dtype=np.float64)
        cadence = np.asarray(lc["Cadence"][:], dtype=np.int64)
        quality = np.asarray(lc["QualityFlag"][:], dtype=np.int64)
        cx = np.asarray(lc["X"][:], dtype=np.float64)
        cy = np.asarray(lc["Y"][:], dtype=np.float64)
        bg = np.asarray(lc["Background"]["Value"][:], dtype=np.float64)
        bg_err = np.asarray(lc["Background"]["Error"][:], dtype=np.float64)

        apertures: dict[str, TGLCAperture] = {}
        for ap in APERTURE_KEYS:
            g = lc["AperturePhotometry"][f"{ap}Aperture"]
            mag = np.asarray(g["RawMagnitude"][:], dtype=np.float64)
            mag_err = np.asarray(g["RawMagnitudeError"][:], dtype=np.float64)
            apx = np.asarray(g["X"][:], dtype=np.float64)
            apy = np.asarray(g["Y"][:], dtype=np.float64)
            if "RawFlux" in g:
                flux = np.asarray(g["RawFlux"][:], dtype=np.float64)
                flux_err = np.asarray(g["RawFluxError"][:], dtype=np.float64)
                synthesized = False
            else:
                # Legacy schema: synthesize flux from magnitude. Negative
                # cadences are unrecoverable; NaN propagates.
                flux = np.full_like(mag, np.nan)
                ok = np.isfinite(mag)
                flux[ok] = 10.0 ** (-0.4 * mag[ok])
                # Conservative error propagation: |df/dm| = 0.4 ln(10) * f
                # but mag_err is itself a MAD-broadcast constant, so the
                # synthesized flux err is also a broadcast. Use MAD on flux.
                if ok.any():
                    sigma = 1.4826 * np.nanmedian(np.abs(flux[ok] - np.nanmedian(flux[ok])))
                else:
                    sigma = np.nan
                flux_err = np.where(ok, sigma, np.nan)
                synthesized = True
            apertures[ap] = TGLCAperture(
                name=ap,
                raw_flux=flux, raw_flux_err=flux_err,
                raw_magnitude=mag, raw_magnitude_err=mag_err,
                centroid_x=apx, centroid_y=apy,
                flux_was_synthesized=synthesized,
            )

    return TGLCLightCurve(
        tic=int(lc_attrs.get("TIC ID", -1)),
        sector=int(lc_attrs.get("Sector", -1)),
        orbit=int(lc_attrs.get("Orbit", -1)),
        cam=int(lc_attrs.get("Camera", -1)),
        ccd=int(lc_attrs.get("CCD", -1)),
        tmag=float(lc_attrs.get("TessMag", np.nan)),
        ra=float(lc_attrs.get("RA", np.nan)),
        dec=float(lc_attrs.get("Dec", np.nan)),
        bjd_offset=int(lc_attrs.get("BJDoffset", 2457000)),
        time=bjd, cadence=cadence, quality=quality,
        centroid_x=cx, centroid_y=cy,
        background=bg, background_err=bg_err,
        apertures=apertures, path=path,
    )
