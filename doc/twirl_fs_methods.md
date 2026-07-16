# TWIRL-FS Detrending Method

## Method family

Use **TWIRL-FS** for the faint-target flux-space detrending product. The
base method-version tag is `twirl-fs-v2`, short for **TWIRL Flux-Spline
version 2**. Current A2v1 products save two adaptive members of this family:

- ADP: `twirl-fs-v2-adp03q`, with `bkspace_d = 0.3 d`;
- ADP015: `twirl-fs-v2-adp015q`, with `bkspace_d = 0.15 d`.

Both use `gap_split_d = 0.2 d`, quantile knots, and robust-auto subtractive
relative flux. The older `0.8 d` base branch remains a methods comparison, not
an A2v1 output column.

## Motivation

Faint TGLC light curves can cross zero after local-background subtraction.
Those negative and near-zero linear-flux excursions are legitimate noise
samples and are important for white-dwarf transit searches. A
magnitude-space detrend or a divisive `flux / spline` normalization is
fragile in this regime because `flux <= 0` is dropped or the fitted spline
can pass through zero.

TWIRL-FS therefore detrends in linear flux and writes QLP-compatible HLSP
columns while preserving negative cadences.

## Base `twirl-fs-v2` parameters

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

## Product labeling

Historical TWIRL-FS comparison FITS use:

- Filename prefix: `hlsp_twirlfs_tess_ffi_*`.
- FITS `PIPELINE`: `TWIRL-FS`.
- FITS `METHOD`: `twirl-fs-v2`.
- FITS detrend metadata: `BKSPACE`, `BSPLK`, `SIGCLIP`, `MAXITER`,
  `EDGEPAD`, `GAPSPLIT`, `OUTMODE`, `SCALE`, `MINSNR`, `NSEG`,
  `FITCNT`, `SCALESRC`, and `COTSTAT`.

The accepted production-family tag is now `A2v1`. A2v1 sector trees use names
such as `hlsp_s0056_A2v1`, carry `METHOD=A2v1`, `PRODTAG=A2v1`, and `A2V1=True`,
and contain ADP/ADP015 branch columns without canonical `DET_FLUX*` or
`SYS_RM_FLUX`. See the [production protocol](a2v1_production_protocol.md).

## Historical S56 compare columns

The Franklin/Michelle S56 handoff has a historical compare tree named
`hlsp_s0056_twirl_fs_v2_compare`. In that tree, the normal columns remain
canonical `twirl-fs-v2`:

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

## ADP015 branch

The July 2026 raw-flux injection/BLS audit selected a modestly stronger
small-aperture search branch for continued vetting:

- method tag: `twirl-fs-v2-adp015q`
- validation/table branch key: `twirl_fs_v2_adp015q`
- compare columns: `DET_FLUX_ADP015`, `DET_FLUX_ADP015_ERR`,
  `DET_FLUX_ADP015_SML`, and `DET_FLUX_ADP015_LAG`
- parameters: `bkspace_d = 0.15` days, `gap_split_d = 0.2` days,
  quantile knots, robust-auto subtractive residuals

ADP015 is now part of every A2v1 product, but it is not the active harmonic
teacher input. The current S56 teacher and human-sheet contract uses
`DET_FLUX_ADP_SML` as the search/morphology channel and `DET_FLUX_ADP` as the
supplemental contamination channel. ADP015 remains available for transparent
search/injection comparisons; changing the production search branch requires a
versioned Stage 2/3 decision rather than silently changing the FITS product.

## Validation status and remaining gates

S56 and S57 A2v1 production now pass HDF5/FITS coverage and schema validation,
and S56 passes the current photometric/WD 1856 gate. The remaining survey-level
requirements are:

- Repeatable per-sector photometric QA rather than schema validation alone.
- A stable A2v1 index and compact-export schema.
- End-to-end recovery through the frozen periodic/dip, ranker, vetter, and
  candidate-merging contract.
- A representative pixel-level calibration of extraction, crowding, aperture,
  and centroid losses.
