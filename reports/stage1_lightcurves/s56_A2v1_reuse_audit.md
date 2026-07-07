# S56 `A2v1` Reuse Audit

Generated: 2026-07-07

## Product Name

Use **`A2v1`** for the regenerated production light-curve family.

Expansion: adaptive-two, version 1.

Sector-specific examples:

- FITS/tree root: `/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_A2v1`
- Compact export tag: `s56_A2v1`

Product contract:

- TGLC catalogs use `--max-magnitude 99`.
- ePSF fitting uses the saturated/overexposed-pixel mask from the ePSF normalization test work.
- Saved derived products are ADP and ADP015 only, for `1x1`, `3x3`, and `5x5` apertures:
  `DET_FLUX_ADP_SML`, `DET_FLUX_ADP`, `DET_FLUX_ADP_LAG`,
  `DET_FLUX_ADP015_SML`, `DET_FLUX_ADP015`, and `DET_FLUX_ADP015_LAG`.
- Canonical/default `DET_FLUX*` columns are not part of this production product except in explicitly labeled compatibility or comparison trees.

## Existing S56 State On `pdogpu6`

Existing root inspected: `/pdo/users/tehan/tglc-gpu-production/`

The existing S56 manifests record `max_magnitude=20.0`, so the old S56 products cannot be treated as no-cap products.

| Orbit | Catalog ECSV | Source pickles | ePSF npy | Active `LC/*.h5` | Archived `LC_pre_relight_v0/*.h5` |
|---:|---:|---:|---:|---:|---:|
| 119 | 32 | 3136 | 3136 | 19086 | 19086 |
| 120 | 0 top-level at inspection, sibling reuse pattern | 3136 | 3136 | 19085 | 19085 |

The missing-HDF5 audit found `31,590` observed S56 TICs, `19,072` unique TICs with existing HDF5 light curves, and `12,518` missing TICs. Of the missing TICs, `12,430` have TIC-side `Tmag > 20`, so this is a real catalog-gate expansion rather than a small targeted backfill.

## Reuse Decision

Can reuse without science changes:

- staged TICA FFI inputs and FFI symlinks;
- Gaia catalogs, because the current TGLC magnitude gate applies to TIC catalog creation;
- sector/orbit/camera/CCD detector tables and requested TIC lists, although the wrapper can cheaply regenerate the per-CCD text files.

Can reuse with a new validated overlay step:

- the pixel/time/quality/mask arrays inside `source_*.pkl` files. These are the expensive cutout payloads and do not depend on the TIC magnitude cap.
- The preferred implementation is a small `source_tic/source_X_Y.ecsv` sidecar built by `scripts/stage1_lightcurves/build_source_tic_overlays.py`, plus a TGLC light-curve-stage hook that replaces `source.tic` from the sidecar after unpickling.
- For TWIRL WD target emission, the preferred sidecar source is now the TWIRL observation table (`twirl_wd_tess_observations_v0.fits`) because it already has `tic_id`, Gaia `source_id`, orbit/camera/CCD, and detector `colpix`/`rowpix`. This bypasses max-99 TIC catalog generation for requested WD targets.
- Max-99 TIC catalogs are still useful if we want a broad all-TIC source match table, but they are not required to recover requested faint WD TICs from existing source pickles.
- This also works for S94+ default source trees, as long as the source files are readable and a matching TWIRL observation table exists for the requested targets.

Cannot reuse as-is for `A2v1`:

- TIC ECSV catalogs from the old `--max-magnitude 20` run;
- embedded `Source.tic` match tables inside existing source pickles;
- ePSF files, because saturated-pixel masking changes the fit;
- TGLC HDF5 light curves, because they depend on both the source TIC match table and the ePSF;
- TWIRL-FS FITS/compact exports and downstream BLS/vetting/injection products.

## Practical Production Boundary

Safe path for the first S56 remake:

1. Use a fresh or clearly versioned staging root to avoid skip-if-exists mixing with old `max_magnitude=20` products.
2. Build `source_tic/*.ecsv` overlays from the TWIRL observation table and symlink the existing `source_*.pkl` files into the new staging root.
3. Only rebuild TIC catalogs with `--max-magnitude 99` when a broad all-TIC source match table is needed; do not use TIC catalog generation as the default requested-WD recovery path.
4. Rerun ePSFs with saturated-pixel masking.
5. Rerun `tglc lightcurves`.
6. Build `A2v1` ADP/ADP015-only products and rerun the missing-HDF5, WD 1856, Julien-ten, and cadence/precision QA gates.

Required small TGLC-side hook:

- After `tglc lightcurves` unpickles a source, check for the matching sidecar at sibling path `source_tic/source_X_Y.ecsv`; if present, read it as an Astropy table and assign `source.tic` before calling `generate_light_curves`.
- The tracked patch artifact is `scripts/stage1_lightcurves/tglc_source_tic_overlay_hook.patch`.
- Validate the hook on TIC `1400899528` plus the three non-edge boundary TICs before trusting it for the full S56-S93 prepared-source backlog.

Example overlay build for one S56 CCD:

```bash
python scripts/stage1_lightcurves/build_source_tic_overlays.py \
  --source-tglc-data-dir /pdo/users/tehan/tglc-gpu-production \
  --output-tglc-data-dir /pdo/users/tehan/tglc-gpu-production-A2v1 \
  --orbit 119 \
  --ccd 4,4 \
  --overlay-from-observations \
  --link-sources \
  --apply \
  --summary-json /pdo/users/tehan/tglc-gpu-production-A2v1/twirl_logs/s56_o119_cam4_ccd4_source_tic_overlay_observations.json
```

## S56 Overlay Smoke On `pdogpu6`

Date: 2026-07-07

The observation-table overlay path was smoke-tested on S56 `cam4/ccd4`, cutout `source_4_12`, using missing faint target TIC `1400899528` (`TIC`-side `Tmag > 20` in the missing-HDF5 audit, TGLC HDF5 metadata `TessMag=19.8334`).

Remote edits and staging:

- Deployed `build_source_tic_overlays.py` to `/pdo/users/tehan/TWIRL/scripts/stage1_lightcurves/`.
- Installed the `source_tic` hook in the user-owned TGLC fork at `/pdo/users/tehan/tess-gaia-light-curve-twirl/tglc/scripts/light_curves.py`.
- Wrote sidecars under `/pdo/users/tehan/tglc-gpu-production-A2v1/orbit-{119,120}/ffi/cam4/ccd4/source_tic/`.

Results:

| Orbit | Source links | Sidecars | ePSF links used for smoke | H5 files | Target cadences |
|---:|---:|---:|---:|---:|---:|
| 119 | 196 | 196 | 1 | 1 | 5638 |
| 120 | 196 | 196 | 1 | 1 | 6137 |

Both orbit sidecars for `source_4_12` include TIC `1400899528` with Gaia DR3 source ID `1439208012920230656`. Both HDF5 files were written at:

- `/pdo/users/tehan/tglc-gpu-production-A2v1/orbit-119/ffi/cam4/ccd4/LC/1400899528.h5`
- `/pdo/users/tehan/tglc-gpu-production-A2v1/orbit-120/ffi/cam4/ccd4/LC/1400899528.h5`

This validates the source-pickle reuse and target-emission hook. The smoke intentionally reused the old `epsf_4_12.npy` through symlinks to isolate the overlay path, so it is **not** a saturated-mask ePSF validation.
