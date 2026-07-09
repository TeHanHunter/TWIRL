# S56 `A2v1` Reuse Audit

Generated: 2026-07-07

## Product Name

Use **`A2v1`** for the regenerated production light-curve family.

Expansion: adaptive-two, version 1.

Sector-specific examples:

- production tree root: `/pdo/users/tehan/tglc-gpu-production-A2v1/`
- S56 FITS root: `/pdo/users/tehan/tglc-gpu-production-A2v1/hlsp_s0056_A2v1`
- Compact export tag: `s56_A2v1`

Product contract:

- TGLC catalogs use `--max-magnitude 99`.
- ePSF fitting uses the saturated/overexposed-pixel mask from the ePSF normalization test work.
- Saved derived products are ADP and ADP015 only, for `1x1`, `3x3`, and `5x5` apertures:
  `DET_FLUX_ADP_SML`, `DET_FLUX_ADP`, `DET_FLUX_ADP_LAG`,
  `DET_FLUX_ADP015_SML`, `DET_FLUX_ADP015`, and `DET_FLUX_ADP015_LAG`.
- Canonical/default `DET_FLUX*` columns are not part of this production product except in explicitly labeled compatibility or comparison trees.

Clean PDO layout contract:

- `orbit-<orbit>/ffi/cam<camera>/ccd<ccd>/source/`: symlinked prepared source
  pickles from the existing production tree.
- `orbit-<orbit>/ffi/cam<camera>/ccd<ccd>/source_tic/`: A2v1 sidecar target
  tables built from the TWIRL observation table.
- `orbit-<orbit>/ffi/cam<camera>/ccd<ccd>/epsf/`: reused empty-mask ePSFs plus
  newly computed saturated-mask A2v1 ePSFs.
- `orbit-<orbit>/ffi/cam<camera>/ccd<ccd>/LC/`: regenerated A2v1 HDF5 light
  curves.
- `hlsp_s<sector>_A2v1/`: ADP/ADP015-only FITS products.
- `twirl_logs/`: reusable run, prefill, handoff, and validation logs.

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

## Full S56 HDF5 Reproduction Launch

Date: 2026-07-07

The full S56 A2v1 HDF5 reproduction was launched on `pdogpu6` in tmux session
`twirl-s56-a2v1-r2` with:

- runner: `/pdo/users/tehan/TWIRL/scripts/stage1_lightcurves/run_s56_a2v1_reproduction_pdo.sh`
- output root: `/pdo/users/tehan/tglc-gpu-production-A2v1/`
- log root: `/pdo/users/tehan/tglc-gpu-production-A2v1/twirl_logs/`
- GPU config for orbit 119: `--max-parallel-ccd-jobs 4 --gpu-list 0,1,2,3 --epsfs-nprocs 4`

The runner rebuilt all observation-table overlays for both S56 orbits:

| Orbit | Sidecars | Source links |
|---:|---:|---:|
| 119 | 3136 | 3136 |
| 120 | 3136 | 3136 |

The remote TGLC fork was checked before launch:

- `tglc/scripts/light_curves.py` contains the `source_tic` overlay hook.
- `tglc/scripts/epsfs.py` contains the saturated-pixel mask path using
  `source.mask.mask`.
- The fork is on branch `twirl/preserve-negative-flux` at commit `8235695`
  with the local `light_curves.py` hook modification.

First HDF5 evidence:

- A completed new orbit-119 cam1/ccd1 ePSF cutout (`source_0_10`) was used for
  a one-target light-curve smoke while the full CCD ePSF stage continued.
- TIC `2053615925` wrote
  `/pdo/users/tehan/tglc-gpu-production-A2v1/orbit-119/ffi/cam1/ccd1/LC/2053615925.h5`.
- The file opens with `HDF5_USE_FILE_LOCKING=FALSE`, has `5638` cadences, and
  contains the expected nested raw-flux schema:
  `LightCurve/AperturePhotometry/{SmallAperture,PrimaryAperture,LargeAperture}/RawFlux`.

Last status snapshot before local monitoring stopped:

| Time (EDT) | Orbit | ePSF files | HDF5 files |
|---|---:|---:|---:|
| 2026-07-07 17:59 | 119 | 108 | 1 |
| 2026-07-07 17:59 | 120 | 0 | 0 |

The full S56 A2v1 HDF5 reproduction is still in progress. Completion still
requires both orbits to finish all ePSFs and light curves, followed by the
ADP/ADP015-only derived product and QA gates.

## Empty-Mask ePSF Reuse

The saturated-pixel ePSF change affects only cutouts where
`source.mask.mask` has at least one masked pixel. For cutouts with an empty
mask, the old ePSF fit and the A2v1 fit are equivalent. The clean A2v1 tree
therefore uses:

- `source/`: symlinks to the existing prepared source pickles.
- `source_tic/`: A2v1 overlay sidecars built from the TWIRL observation table.
- `epsf/`: symlinks to old ePSFs for empty-mask cutouts, plus newly computed
  A2v1 ePSFs for masked cutouts.
- `LC/`: regenerated A2v1 HDF5 light curves.

This keeps the A2v1 folder structure explicit while avoiding unnecessary ePSF
refits. A future archival copy can materialize symlinked files if a fully
self-contained tree is needed.

S56 manual prefill results:

| Orbit | Existing A2v1 ePSFs kept | Empty-mask ePSFs linked | Masked cutouts left for recompute | Errors |
|---:|---:|---:|---:|---:|
| 119 | 244 | 2090 | 802 | 0 |
| 120 | 0 | 2320 | 816 | 0 |

After this prefill, production was relaunched in tmux `twirl-s56-a2v1-r3` with
`TWIRL_A2V1_PREFILL_EMPTY_MASK_EPSFS=0`. The orbit-119 cam1 ePSF logs confirm
that linked ePSFs are skipped quickly and only masked cutouts enter the slower
saturated-mask fit.

## Reusable A2v1 Launcher

The clean A2v1 path is no longer S56-specific. The generic launcher is:

```bash
bash scripts/stage1_lightcurves/run_a2v1_reproduction_pdo.sh 56 119:o1 120:o2
```

For a later S94-style run with the same source-pickle reuse policy:

```bash
bash scripts/stage1_lightcurves/run_a2v1_reproduction_pdo.sh 94 195:o1 196:o2
```

The S56 convenience wrapper delegates to this generic launcher. GPU selection
and worker counts are controlled by environment variables such as
`TWIRL_A2V1_GPU_MAX_PARALLEL`, `TWIRL_A2V1_EPSFS_NPROCS`, and
`TWIRL_A2V1_PREFILL_NPROCS`, so future sectors can use the same A2v1 folder
structure while adapting to current PDO load.

## GPU-Lane Resume Fix

The first post-prefill resume (`twirl-s56-a2v1-r3`) exposed a scheduler bug in
the orbit driver: GPU IDs were assigned by CCD list index before executor
scheduling, so when a fast CCD finished, the next queued CCD could start on a
GPU still occupied by a slower masked-cutout CCD. Orbit-119 `cam2/ccd2` hit a
CuPy out-of-memory error this way while another lane was still using the same
physical GPU.

The driver now uses one queue-consuming worker lane per GPU ID. Each lane keeps
its assigned `CUDA_VISIBLE_DEVICES` value and only pulls another CCD after the
current CCD finishes. The stale `r3` run was stopped without cleaning outputs,
and `twirl-s56-a2v1-r4` resumed from the existing A2v1 tree. The patched run
auto-skipped completed `cam1/ccd1` and `cam1/ccd3`, then resumed the incomplete
orbit-119 CCDs on distinct GPU lanes.

The generic launcher also accepts `TWIRL_A2V1_GPU_LIST` for future explicit-GPU
resumes when the free-memory picker would choose a busy device. Live checkpoint:
at `2026-07-07 21:31 EDT`, `twirl-s56-a2v1-r4` was still running orbit `119`
with `2438` ePSFs and `3143` HDF5 files under the A2v1 root.

## A2v1-Only FITS Export Mode

The ADP/ADP015-only export path is implemented as:

```bash
bash scripts/stage1_lightcurves/run_a2v1_hlsp_pdo.sh 56 119 120
```

For S56, the default output root is:

```text
/pdo/users/tehan/tglc-gpu-production-A2v1/hlsp_s0056_A2v1
```

The wrapper calls `build_twirl_hlsp.py --a2v1-only`, which writes only the
A2v1 detrended branch products:

- `DET_FLUX_ADP_SML`, `DET_FLUX_ADP`, `DET_FLUX_ADP_LAG`
- `DET_FLUX_ADP015_SML`, `DET_FLUX_ADP015`, `DET_FLUX_ADP015_LAG`

Primary-branch error columns `DET_FLUX_ADP_ERR` and `DET_FLUX_ADP015_ERR` are
also present. Canonical `DET_FLUX*` and `SYS_RM_FLUX` are omitted in this mode.
A one-target smoke on `pdogpu6` from a real orbit-119 A2v1 HDF5 passed with
`METHOD=A2v1`, `PRODTAG=A2v1`, required ADP/ADP015 columns present, and
canonical detrended columns absent.

The handoff monitor `twirl-s56-a2v1-hlsp-monitor` is running on `pdogpu6`. It
waits for the HDF5 reproduction completion marker before launching this FITS
build, so the A2v1 FITS tree is not generated from partial orbit products.

## A2v1 Product Validator

The reusable validator is:

```bash
python scripts/stage1_lightcurves/validate_a2v1_product.py \
  --a2v1-root /pdo/users/tehan/tglc-gpu-production-A2v1 \
  --sector 56 \
  --orbits 119 120 \
  --summary-json reports/stage1_lightcurves/s56_A2v1_validation_summary.json
```

For S56, it derives the no-cap expected target set from
`twirl_wd_tess_observations_v0.fits`: `63,238` orbit-level HDF5 files and
`31,590` unique TIC-level FITS files. It verifies HDF5 coverage by
orbit/camera/CCD and opens the A2v1 FITS products to enforce:

- `METHOD=A2v1`, `PRODTAG=A2v1`, and `A2V1=True`;
- required branch columns `DET_FLUX_ADP*` and `DET_FLUX_ADP015*`;
- forbidden canonical columns absent: `DET_FLUX`, `DET_FLUX_ERR`,
  `DET_FLUX_SML`, `DET_FLUX_LAG`, and `SYS_RM_FLUX`.

Local synthetic tests exercise both the accepted ADP/ADP015-only schema and a
rejected file that still contains canonical detrended columns.

## Current Remote Checkpoint

Last confirmed `pdogpu6` state before the SSH gateway began requiring
keyboard-interactive authentication:

- timestamp: `2026-07-07 21:55:14 EDT`;
- tmux: `twirl-s56-a2v1-r4` and `twirl-s56-a2v1-hlsp-monitor` active;
- orbit `119`: `2484` ePSFs and `4677` HDF5 files;
- orbit `120`: `2320` reused ePSFs and `0` HDF5 files;
- A2v1 FITS: `0` files under
  `/pdo/users/tehan/tglc-gpu-production-A2v1/hlsp_s0056_A2v1`.

After that checkpoint, non-interactive SSH through `hostess3` failed with
`Permission denied (keyboard-interactive)`. Per TWIRL policy, Codex did not
initiate interactive auth. Production may still be running under tmux, but the
next authoritative status check requires the user-opened PDO control session to
be available again.

## R5 Resume Checkpoint

PDO access was restored on `2026-07-08`. The `r4` run had exited after orbit
`119` failed on two remaining masked ePSF CCDs:

- `cam3/ccd1`: CuPy OOM during `tglc epsfs --nprocs 4`;
- `cam3/ccd4`: CuPy OOM during `tglc epsfs --nprocs 4`.

The A2v1 tree was not cleaned. A resume run was launched in tmux session
`twirl-s56-a2v1-r5` with:

```bash
TWIRL_A2V1_PREFILL_EMPTY_MASK_EPSFS=0
TWIRL_A2V1_GPU_LIST=0,1,2,3,6,7
TWIRL_A2V1_GPU_MAX_PARALLEL=4
TWIRL_A2V1_EPSFS_NPROCS=2
HDF5_USE_FILE_LOCKING=FALSE
bash scripts/stage1_lightcurves/run_a2v1_reproduction_pdo.sh 56 119:o1 120:o2
```

The resume auto-skipped the `14` completed orbit-119 CCDs and reran only
`cam3/ccd1` plus `cam3/ccd4` with the lower ePSF worker count. At
`2026-07-08 10:17:24 EDT`, both jobs were still healthy and past the previous
failure point:

- orbit `119`: `3030` ePSFs and `27149` HDF5 files;
- orbit `120`: `2320` reused ePSFs and `0` HDF5 files;
- `cam3/ccd1`: `65 / 196` ePSFs in progress;
- `cam3/ccd4`: `73 / 196` ePSFs in progress.

A one-off tmux monitor `twirl-s56-a2v1-r5-hlsp` watches the r5-specific driver
log and will run `run_a2v1_hlsp_pdo.sh 56 119 120` after the r5 HDF5
completion marker appears. It ignores the older r4 abort because it watches
`s56_a2v1_r5_driver.log`, not the cumulative reproduction log.

At `2026-07-08 11:05:20 EDT`, the reduced-parallelism r5 ePSF rerun was still
healthy:

- orbit `119`: `3086` ePSFs and `27149` HDF5 files;
- orbit `120`: `2320` reused ePSFs and `0` HDF5 files;
- A2v1 FITS: `0` files under `hlsp_s0056_A2v1`.
