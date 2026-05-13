"""TWIRL flux-space cotrending — Plan A.

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

The fix is structural, in three coordinated parts:

1. **TGLC patch** (upstream, our fork at /pdo/users/tehan/tess-gaia-light-curve-twirl/):
   - `aperture_photometry.py`: do not NaN the flux array. Compute the
     magnitude on a `np.where(flux > 0, flux, np.nan)` copy so the
     `tess_magnitude` call still works, but keep `flux` itself intact.
   - `aperture_light_curve.py`: add `RawFlux` and `RawFluxError`
     datasets to the per-aperture HDF5 group, alongside the existing
     `RawMagnitude` / `RawMagnitudeError`.

2. **This module** (`src/twirl/lightcurves/`):
   - `flux_detrend.py`: flux-space BSpline cotrend, mirrors QLP
     `lctools detrend` knot spacing (bkspace = 0.3 d) but operates on
     linear flux. NaN inputs are interpolated through; negatives are
     kept and weighted normally.
   - `hlsp_writer.py`: emits FITS files in the QLP HLSP schema
     (TIME, SAP_FLUX, DET_FLUX, DET_FLUX_ERR, QUALITY, SAP_X, SAP_Y,
     SAP_BKG, DET_FLUX_SML, DET_FLUX_LAG, ...) so all existing
     downstream code (BLS, heuristic vetter, LEO adapter) reads
     these unchanged. Naming convention: `hlsp_twirl_tess_ffi_*` to
     coexist with QLP's `hlsp_qlp_tess_ffi_*`.
   - `tglc_h5_reader.py`: reads the (post-patch) TGLC HDF5 with
     `RawFlux` columns; backward-compatible path for the legacy
     magnitude-only schema converts mag -> flux for the bright regime
     where the loss is invisible.

3. **Driver** (`scripts/stage1_lcs/build_twirl_hlsp.py`):
   Per TIC, merge orbits in the sector, run flux_detrend on each
   aperture, write a TWIRL HLSP FITS. Replaces the `qlp lctools
   detrend -> qlp lctools hlsp` step of the existing pipeline.

The TGLC photometry layer is untouched; the ePSF fit, aperture summing,
saturation handling, and local-background subtraction remain TGLC's.
"""
