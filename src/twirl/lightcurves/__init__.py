"""TWIRL-FS flux-space cotrending.

Replaces QLP `lctools detrend` for the TWIRL pipeline. The reason the
replacement exists is a single line in upstream TGLC photometry:

    tess-gaia-light-curve/tglc/aperture_photometry.py:114
        flux[flux <= 0] = np.nan  # Prevent runtime warnings converting to magnitude

For faint hosts (T >= 19.5) about half the QUALITY==0 cadences have
negative post-background flux from Poisson + sky noise — physically real,
necessary for unbiased noise statistics, but unrepresentable in
magnitude space. The line above silently drops them. Downstream QLP
detrend then operates on the already-gappy magnitude series; the
resulting BSpline overshoots near the gaps and biases the recovered
relative flux upward. For a transit-survey context: short, deep dips
sit precisely where the background-subtracted flux is most likely to be
negative, so the survey loses sensitivity exactly where it should be
strongest.

The TWIRL-FS (`twirl-fs-v2`) fix is structural, in three coordinated
parts:

1. **TGLC patch** (upstream, our fork at /pdo/users/tehan/tess-gaia-light-curve-twirl/):
   - `aperture_photometry.py`: do not NaN the flux array. Compute the
     magnitude on a `np.where(flux > 0, flux, np.nan)` copy so the
     `tess_magnitude` call still works, but keep `flux` itself intact.
   - `aperture_light_curve.py`: add `RawFlux` and `RawFluxError`
     datasets to the per-aperture HDF5 group, alongside the existing
     `RawMagnitude` / `RawMagnitudeError`.

2. **This module** (`src/twirl/lightcurves/`):
   - `flux_detrend.py`: flux-space BSpline cotrend. The current
     production default is `bkspace = 0.8 d`, chosen from S56 injection
     sweeps to preserve 5 min-6 hr transit/eclipse signals better than
     QLP's historical 0.3 d spacing. It fits independent cotrends across
     large time gaps so inter-orbit gaps do not force a polynomial
     fallback. NaN inputs are interpolated through; negatives are kept
     and weighted normally.
   - `hlsp_writer.py`: emits FITS files in the QLP HLSP schema
     (TIME, SAP_FLUX, DET_FLUX, DET_FLUX_ERR, QUALITY, SAP_X, SAP_Y,
     SAP_BKG, DET_FLUX_SML, DET_FLUX_LAG, ...) so all existing
     downstream code (BLS, heuristic vetter, LEO adapter) reads
     these unchanged. Naming convention: `hlsp_twirlfs_tess_ffi_*` to
     coexist with QLP's `hlsp_qlp_tess_ffi_*`.
   - `tglc_h5_reader.py`: reads the (post-patch) TGLC HDF5 with
     `RawFlux` columns; backward-compatible path for the legacy
     magnitude-only schema converts mag -> flux for the bright regime
     where the loss is invisible.

3. **Driver** (`scripts/stage1_lightcurves/build_twirl_hlsp.py`):
   Per TIC, merge orbits in the sector, run flux_detrend on each
   aperture, write a TWIRL-FS HLSP FITS. Replaces the `qlp lctools
   detrend -> qlp lctools hlsp` step of the existing pipeline.

The TGLC photometry layer is untouched; the ePSF fit, aperture summing,
saturation handling, and local-background subtraction remain TGLC's.
"""
