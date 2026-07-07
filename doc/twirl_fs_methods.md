# TWIRL-FS Detrending Method

## Method Name

Use **TWIRL-FS** for the faint-target flux-space detrending product. The
method-version tag for the current implementation is `twirl-fs-v2`, short
for **TWIRL Flux-Spline version 2**.

## Motivation

Faint TGLC light curves can cross zero after local-background subtraction.
Those negative and near-zero linear-flux excursions are legitimate noise
samples and are important for white-dwarf transit searches. A
magnitude-space detrend or a divisive `flux / spline` normalization is
fragile in this regime because `flux <= 0` is dropped or the fitted spline
can pass through zero.

TWIRL-FS therefore detrends in linear flux and writes QLP-compatible HLSP
columns while preserving negative cadences.

## Current `twirl-fs-v2` Parameters

- Input flux: TGLC HDF5 `RawFlux` and `RawFluxError`, one sector merged
  across all available orbits.
- Fit model: cubic `LSQUnivariateSpline`.
- Gap handling: independent cotrend fits across time gaps larger than
  `gap_split_d = 0.5` days. This prevents the multi-day gap between TESS
  orbits from violating the spline knot conditions and silently degrading
  to a low-order polynomial.
- Knot spacing: `bkspace_d = 0.8` days.
- Edge padding: `edge_pad_d = 0.5` days.
- Iterative rejection: `sigma_clip = 5.0`, `max_iter = 5`.
- Fit mask: finite flux, finite positive flux error, and `QUALITY == 0`.
- Capped evaluation: flagged edge cadences outside the fitted quality-zero
  time span are evaluated at the nearest fitted boundary, avoiding
  unconstrained spline extrapolation in HLSP diagnostics.
- Output mode: subtractive relative flux,
  `DET_FLUX = 1 + (flux - spline) / scale`, recentered so the good-cadence
  median is exactly 1.
- Scale strategy: `auto`. The code uses the signed median only when its
  absolute value is at least `3 * robust_sigma(flux)` and positive enough
  to be stable; otherwise it falls back to a positive robust absolute-flux
  scale.
- Error model: `DET_FLUX_ERR` is the MAD-based scatter of the detrended
  relative-flux series broadcast to all cadences.

## Why Not The QLP-Like 0.3 Day Spacing

The local S56 sweep tested `bkspace_d` values of `0.3`, `0.5`, `0.8`,
`1.2`, and `2.0` days with injected 5, 10, 30, 60, 120, and 180 minute
events at depths of 3%, 10%, 30%, and 50%. The 0.3 day spacing removed
low-frequency structure most aggressively, but it also absorbed long
signals: its worst median depth retention was about 0.62 for 3 hour
events, and 10% events were already strongly suppressed at 1-2 hours.

The `0.8` day, `sigma_clip=5` configuration was the shortest tested
spacing that kept the 16th-84th percentile retention envelope within
0.90-1.10 for the 5 minute to 3 hour sweep, while retaining more trend
removal than 1.2-2.0 day spacings. A separate 3-6 hour stress test kept
the same configuration above the retention threshold, with worst
16th-percentile retention about 0.93.

## Product Labeling

TWIRL-FS FITS files use:

- Filename prefix: `hlsp_twirlfs_tess_ffi_*`.
- FITS `PIPELINE`: `TWIRL-FS`.
- FITS `METHOD`: `twirl-fs-v2`.
- FITS detrend metadata: `BKSPACE`, `BSPLK`, `SIGCLIP`, `MAXITER`,
  `EDGEPAD`, `GAPSPLIT`, `OUTMODE`, `SCALE`, `MINSNR`, `NSEG`,
  `FITCNT`, `SCALESRC`, and `COTSTAT`.

For the S56 production test, use an output tree named
`hlsp_s0056_twirl_fs_v2`.

## S56 Compare Columns

The Franklin/Michelle S56 handoff also has a compare tree named
`hlsp_s0056_twirl_fs_v2_compare`. In that tree, the normal columns remain
canonical `twirl-fs-v2` and should be the default search input:

- `DET_FLUX`, `DET_FLUX_ERR`, `DET_FLUX_SML`, `DET_FLUX_LAG`.

The same FITS files also include opt-in adaptive columns:

- `DET_FLUX_ADP`, `DET_FLUX_ADP_ERR`, `DET_FLUX_ADP_SML`,
  `DET_FLUX_ADP_LAG`.

These adaptive columns use method tag `twirl-fs-v2-adp03q`: the same
robust-auto subtractive residual, but with `bkspace_d = 0.3` days,
`gap_split_d = 0.2` days, and quantile-based knot placement. Quantile knots
are important here; a first uniform-knot 0.3 day attempt failed QA because
quality-zero cadence gaps caused many adaptive fits to fall back to the
polynomial path. The corrected compare product records adaptive diagnostics
in `ADPMETH`, `ADPBKSP`, `ADPGAP`, `ADPKNOT`, `ADPNSEG`, `ADPFIT`,
`ADPSCAL`, and `ADPCOTS`.

## Candidate S56 Search Branch: `twirl-fs-v2-adp015q`

The July 2026 raw-flux injection/BLS audit selected a modestly stronger
small-aperture search branch for continued vetting:

- method tag: `twirl-fs-v2-adp015q`
- validation/table branch key: `twirl_fs_v2_adp015q`
- compare columns: `DET_FLUX_ADP015`, `DET_FLUX_ADP015_ERR`,
  `DET_FLUX_ADP015_SML`, and `DET_FLUX_ADP015_LAG`
- parameters: `bkspace_d = 0.15` days, `gap_split_d = 0.2` days,
  quantile knots, robust-auto subtractive residuals

This is not yet the canonical survey FITS product. It is the current
candidate search/vetting branch: use `DET_FLUX_ADP015_SML` as the primary
BLS search aperture and `DET_FLUX_ADP015` as the two-aperture contamination
comparison. The existing `twirl-fs-v2-adp03q` columns remain the historical
S56 compare columns.

## Validation Gate Before Production Replacement

Before replacing the normal QLP detrend path, TWIRL-FS should pass:

- Full S56 build success with no new cadence-loss mode relative to the TGLC
  HDF5 inputs.
- Faint-end QA on negative and near-zero flux light curves.
- WD 1856 S56 recovery across `1x1`, `3x3`, and `5x5` apertures.
- Injection preservation for short signals, with explicit checks at about
  30 minutes and multi-hour durations.
- BLS and heuristic checks on the labeled TWIRL-FS S56 HLSP tree.
