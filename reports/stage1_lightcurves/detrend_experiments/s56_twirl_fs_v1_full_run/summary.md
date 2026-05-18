# S56 TWIRL-FS v1 Full-Run Summary

Run date: `2026-05-18`

Output tree:
`/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl_fs_v1`

## Build

- Host: `pdogpu6`
- Runtime: `/pdo/users/tehan/twirl-gpu-venv/bin/python` with
  `LD_LIBRARY_PATH=/sw/python-versions/python-3.11.9/lib`
- Inputs:
  - `/pdo/users/tehan/tglc-gpu-production/orbit-119/ffi`
  - `/pdo/users/tehan/tglc-gpu-production/orbit-120/ffi`
- Script: `scripts/stage1_lcs/build_twirl_hlsp.py`
- Workers: `32`
- Result: `19,072 ok / 0 fail`
- Wall time from script summary: `17.4 min`

The target union was `19,072` TICs. Orbit 119 contributed `19,072` HDF5
files and orbit 120 contributed `19,071` HDF5 files.

## Method Metadata

Smoke-test and sample-audit FITS headers confirm:

- Filename prefix: `hlsp_twirlfs_tess_ffi_*`
- `PIPELINE = TWIRL-FS`
- `METHOD = twirl-fs-v1`
- `BKSPACE = 0.8`
- `SIGCLIP = 5.0`
- `OUTMODE = subtractive`
- `SCALE = auto`

## Sample Audit

Audited WD 1856 plus seven representative faint S56 targets:

| TIC | cadences | q0 cadences | finite DET q0 | negative SAP q0 | finite DET for negative SAP q0 | DET q0 MAD |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 267574918 | 11775 | 11136 | 11136 | 164 | 164 | 0.1282 |
| 1884514391 | 11775 | 11775 | 11775 | 7079 | 7079 | 0.1961 |
| 1981551646 | 11775 | 11775 | 11775 | 6631 | 6631 | 0.1631 |
| 1981662747 | 11775 | 11735 | 11735 | 7411 | 7411 | 0.1527 |
| 2014550155 | 11775 | 11775 | 11775 | 8112 | 8112 | 0.1304 |
| 2041357363 | 11775 | 11775 | 11775 | 7818 | 7818 | 0.1512 |
| 2053290277 | 11775 | 11550 | 11550 | 7473 | 7473 | 0.2656 |
| 2053333283 | 11775 | 11419 | 11419 | 7404 | 7404 | 0.2865 |

Aggregate over the eight audited targets:

- `92,940 / 92,940` quality-zero cadences have finite `DET_FLUX`.
- `52,092 / 52,092` negative quality-zero `SAP_FLUX` cadences retain finite
  `DET_FLUX`.
- Good-cadence `DET_FLUX` medians are recentered to `1.0` for every audited
  target.

## Interpretation

The full S56 TWIRL-FS v1 build completed without new per-target failures and
the representative faint-target audit preserved every negative quality-zero
cadence as a finite detrended measurement. This supports using
`twirl-fs-v1` as the labeled S56 product for the next BLS/vetting comparison.
