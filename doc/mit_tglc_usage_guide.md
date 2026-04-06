# Using The MIT-Adapted TGLC For TWIRL

This note compares your original `TESS_Gaia_Light_Curve` repo with the MIT fork at `mit-kavli-institute/tess-gaia-light-curve` and summarizes how to use the new version effectively for the WD-planet search.

## Executive Summary

Your original repo is centered on single-target extraction through Python functions such as `tglc.quick_lc.tglc_lc()` and `target_lightcurve.epsf()`. The MIT fork is a production pipeline for orbit/camera/CCD batches, exposed through a `tglc` command-line interface and organized around pre-staged FFIs and local catalog databases.

For TWIRL, this means:

- do not treat the MIT fork as a drop-in replacement for `quick_lc.py`,
- think in terms of orbit-camera-CCD jobs,
- use the MIT fork to mass-produce HDF5 light curves,
- build TWIRL-specific indexing, QC, detection, and validation around those products.

## Additional PDO Operational References

Beyond the `tess-gaia-light-curve` code itself, the MIT Kavli `qlp-ops` wiki contains the current operator-facing PDO workflows around TGLC and downstream QLP processing:

- Photometry roadmap:
  `https://github.com/mit-kavli-institute/qlp-ops/wiki/Photometry-Roadmap`
- Planet-search roadmap:
  `https://github.com/mit-kavli-institute/qlp-ops/wiki/Planet-Search-Roadmap`

For TWIRL, the photometry roadmap is the more directly relevant reference because it documents the PDO operational context around orbit setup, `tglc all`, QLP ingestion, detrending, and CCD-level quality-flag handling on `/pdo/qlp-data`. The planet-search roadmap is less directly tied to TWIRL Stage 1, but it is useful background on how MIT's broader PDO pipeline continues from finished light curves into BLS, astronet triage, difference imaging, and report generation.

These wiki pages should be treated as operational context rather than TWIRL requirements. TWIRL still differs from standard QLP operation in several important ways:

- the scientific sample remains Gaia-first and WD-defined, not TIC-defined
- the target regime is much fainter than the bright-star defaults used in standard QLP/TGLC operations
- TWIRL Stage 1 is about producing and indexing WD light curves, not automatically adopting the downstream QLP planet-search stack

## What Changed Relative To The Original Repo

| Area | Original repo | MIT fork | Practical consequence for TWIRL |
| --- | --- | --- | --- |
| Packaging | `setup.py`, `requirements.txt`, `README.rst` | `pyproject.toml`, `README.md`, modern test suite | The MIT fork is easier to deploy reproducibly on shared infrastructure. |
| Entry point | Python functions and ad hoc scripts | `tglc` CLI with `catalogs`, `cutouts`, `epsfs`, `lightcurves`, plus `all` | Production should be scripted from the CLI, not from `quick_lc.py`. |
| Processing model | target-by-target | orbit/camera/CCD batch processing | Build job tables by orbit and detector footprint. |
| Catalog access | online lookups and local helper code | local `pyticdb` access to `tic_82` and `gaia3` | MIT PDO machines must have working local catalog DB access. |
| FFI access | `ffi_cut` / MAST-driven workflows | user-staged TICA FFIs on disk | FFIs must already exist in the expected directory tree. |
| Data layout | custom `logs/`, `lc/`, `epsf/`, `source/`, `plots/` under a local directory | strict `tglc-data/orbit-<n>/ffi/...` layout managed by `Manifest` | TWIRL should adopt the MIT data layout rather than invent a parallel one. |
| Output product | FITS light curves with PSF/aperture columns, optional saved apertures | HDF5 light curves with background plus 1x1, 3x3, 5x5 aperture photometry | Downstream TWIRL code needs new readers and should not expect `cal_psf_flux`. |
| PSF usage | user-visible PSF and prior tuning in the single-target workflow | ePSF fit is internal, output is decontaminated aperture photometry | The new workflow is less interactive but better for mass production. |
| Extra apertures | optional 5x5 handling | native small 1x1, primary 3x3, large 5x5 apertures | Useful for WD transit sensitivity experiments and vetting. |
| Tests | minimal | explicit CLI, utility, and end-to-end tests | The MIT fork is much safer for long-running survey production. |

## Important Behavioral Differences

### 1. There is no direct `quick_lc.py` equivalent

Old usage:

```python
from tglc.quick_lc import tglc_lc
tglc_lc(target="TIC 264468702", local_directory="...", sector=56, save_aper=True)
```

New usage:

- stage FFIs for an orbit/camera/CCD,
- generate catalogs,
- make cutouts,
- fit ePSFs,
- extract light curves.

That is a production-pipeline mindset, not an interactive target-first one.

### 2. The MIT fork expects local infrastructure

The code assumes:

- TICA FFIs are already on disk,
- `pyticdb` databases named `tic_82` and `gaia3` are available,
- you have a working `tglc-data` tree.

The code does not fetch FFIs for you.

### 3. Output targets are TIC-driven

The MIT fork uses Gaia for the field-star model, but the light-curve extraction stage produces files for TIC targets. For TWIRL this creates one key requirement:

- verify that the WD sample is sufficiently complete after Gaia-to-TIC matching.

If important WDs do not have usable TIC IDs, the MIT fork will need a small extension so the light-curve stage can emit Gaia-selected targets directly.

### 4. The default magnitude limits are too bright for TWIRL

The MIT CLI defaults are tailored to QLP-like bright-star work:

- `--max-magnitude 13.5`
- `--mdwarf-magnitude 15.0`

For WD work, you should explicitly raise the target magnitude limit, likely to around `20`, during catalog creation or `all` runs.

### 5. The old `prior` knob is gone from the public interface

Your original workflow exposed a `prior` argument that controlled floating field stars in the fit. The MIT CLI does not currently expose an equivalent parameter.

Implication:

- if TWIRL needs that flexibility for crowded WD fields, it will require a code change instead of a runtime option.

## The MIT Fork Data Layout

The `Manifest` class hard-codes the expected file structure:

```text
tglc-data/
  orbit-185/
    ffi/
      catalogs/
        TIC_cam1_ccd1.ecsv
        Gaia_cam1_ccd1.ecsv
      cam1/
        ccd1/
          ffi/
            hlsp_tica_tess_ffi_s0089-o1-........-cam1-ccd1_tess_v01_img.fits
          source/
            source_0_0.pkl
          epsf/
            epsf_0_0.npy
          LC/
            123456789.h5
```

If you run from inside a directory named `tglc-data`, the CLI auto-detects it. Otherwise pass `--tglc-data-dir`.

## Recommended TWIRL Usage Pattern On MIT PDO

### 1. Install the MIT fork explicitly

Use the MIT repo, not the released package from your original repo:

```bash
pip install git+https://github.com/mit-kavli-institute/tess-gaia-light-curve.git
```

For development on PDO or in a clone:

```bash
pip install -e ".[dev]"
```

### 2. Stage the FFIs first

For each orbit/camera/CCD you want to process, place the TICA calibrated FFI files in:

```text
tglc-data/orbit-<orbit>/ffi/cam<camera>/ccd<ccd>/ffi/
```

The MIT code looks for file names matching the TICA pattern:

```text
hlsp_tica_tess_ffi_s<sector>-o<orbit-within-sector>-<cadence>-cam<camera>-ccd<ccd>_tess_v01_img.fits
```

### 3. Confirm database access

Before large runs, verify that the PDO environment can resolve:

- `Databases["tic_82"]`
- `Databases["gaia3"]`

The MIT code uses those exact database names.

### 4. Run one orbit-camera-CCD job at a time

A typical production sequence is:

```bash
tglc catalogs --orbit 185 --ccd 1,1 --max-magnitude 20 --nprocs 16 --tglc-data-dir /path/to/tglc-data
tglc cutouts --orbit 185 --ccd 1,1 --nprocs 16 --tglc-data-dir /path/to/tglc-data
tglc epsfs --orbit 185 --ccd 1,1 --nprocs 16 --tglc-data-dir /path/to/tglc-data
tglc lightcurves --orbit 185 --ccd 1,1 --nprocs 16 --tglc-data-dir /path/to/tglc-data
```

In TWIRL, the first safe wrapper around this pattern is
[`scripts/stage1_lcs/run_tglc_catalogs.py`](/Users/tehan/PycharmProjects/TWIRL/scripts/stage1_lcs/run_tglc_catalogs.py).
It reads the orbit-aware TWIRL detector summary and, by default, prints or writes the exact
`tglc catalogs` commands for the selected orbit/camera/CCD jobs before any execution happens.

For the current Sector 56 full-orbit benchmark, use
[`scripts/stage1_lcs/run_tglc_orbit_pipeline.py`](/Users/tehan/PycharmProjects/TWIRL/scripts/stage1_lcs/run_tglc_orbit_pipeline.py).
That driver runs `catalogs -> cutouts -> epsfs -> lightcurves` sequentially for each CCD, launches
multiple CCD jobs concurrently at the outer level, and filters `lightcurves` to WD TIC IDs from the
TWIRL observation table.

One extra `catalogs` caveat found in the Sector 56 orbit `119` benchmark: in dense CCDs, a single
TIC-cone query can return more than `65,535` Gaia DR2 source IDs, and the stock DR2->DR3 bridge
lookup then fails with `psycopg.OperationalError` because SQLAlchemy binds each ID as one query
parameter in an `IN (...)` clause. The TWIRL user-owned TGLC fork now batches that bridge query in
chunks of `50,000` IDs to stay below the Postgres parameter limit.

Or, once the setup is trusted:

```bash
tglc all --orbit 185 --ccd 1,1 --max-magnitude 20 --nprocs 16 --tglc-data-dir /path/to/tglc-data
```

### 5. Use `--cutout` for smoke tests

Before a full CCD run, test one cutout:

```bash
tglc all --orbit 185 --ccd 1,1 --cutout 0,0 --max-magnitude 20 --tglc-data-dir /path/to/tglc-data
```

This is the fastest way to validate:

- FFI staging
- catalog DB access
- ePSF fitting
- output serialization

### 6. Use GPUs only for the ePSF stage

The MIT code uses optional `cupy` support in `epsfs`. If GPU support is missing or unstable, run:

```bash
tglc epsfs ... --no-gpu
```

The other stages are CPU-oriented.

### 7. Understand where `epsfs` parallelizes

`tglc epsfs` parallelizes over cutout pickle files **within one CCD** through `--nprocs`.
If multiple CCDs are passed in one command, the current implementation still loops over
`camera,ccd` serially and only parallelizes the `source_*.pkl` files for the active CCD.

Operational consequence for TWIRL:

- one `orbit/camera/ccd` job with `--nprocs 16` launches one parent process and sixteen worker
  processes, plus one Python multiprocessing resource-tracker helper in GPU mode
- scaling to a whole orbit through one serial `tglc epsfs` command is likely to be slow because all
  `16` CCDs are processed one after another
- before mass production, prefer a job launcher that parallelizes across independent
  `orbit/camera/ccd` units, with a modest per-job `--nprocs`, rather than one monolithic orbit job
- the first full-orbit benchmark setting is orbit `119` on `pdogpu1` with `4` CCD jobs in parallel,
  while keeping `catalogs/cutouts/epsfs/lightcurves` worker counts at `16/16/16/16` inside each CCD
  job; orbit `120` should use a different outer CCD concurrency for comparison while leaving
  `epsfs --nprocs 16` fixed
- expect uneven CCD runtimes because crowded/heavy CCDs are slower than sparse CCDs
- the current default for TWIRL one-CCD `epsfs` jobs is `--nprocs 16`: in the Sector 56 benchmark,
  orbit `120` on `pdogpu1` with `--nprocs 16` finished in `4:01:31`, while orbit `119` on `pdogpu6`
  with `--nprocs 64` finished in `4:03:04`, so `64` workers did not show a clear advantage

## Known-Good PDO Environment And Benchmark Pattern

The current PDO benchmark work established a concrete, working launch pattern for the MIT fork.

### PDO environment notes

- The orbit-local virtual environment at `/pdo/users/tehan/tglc-deep-catalogs/orbit-195/ffi/catalogs/.venv/` is not sufficient by itself on `pdogpu6`; calling its Python directly failed with `libpython3.11.so.1.0` missing.
- The usable shared runtime is the QLP environment at `/sw/qlp-environment/.venv/`.
- On PDO, the MIT fork is importable from that shared environment when the Python 3.11 library path is added through `LD_LIBRARY_PATH`.

The minimal validated import test was:

```bash
deactivate 2>/dev/null
env LD_LIBRARY_PATH=/pdo/app/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH} \
  /sw/qlp-environment/.venv/bin/python \
  -c "import sys, tglc; print(sys.executable); print(tglc.__file__)"
```

That resolved to:

- `/sw/qlp-environment/.venv/bin/python`
- `/sw/qlp-environment/.venv/lib/python3.11/site-packages/tglc/__init__.py`

For current PDO work, prefer this shared environment pattern over trying to revive an orbit-local `.venv`.

Two benchmark-specific runtime caveats found on `2026-04-02`:

- on `pdogpu6`, the shared QLP environment currently reports `HAS_CUPY=False`, so `tglc epsfs`
  is running CPU-only unless a `cupy`-enabled environment is prepared
- on `pdogpu1`, importing `numba` from `/sw/qlp-environment/.venv/bin/python` failed until
  `/pdo/app/anaconda/anaconda2-4.4.0/lib` was added ahead of the Python 3.11 library path so
  `libffi.so.6` can be found

For the current CPU-only benchmark, prefer `pdogpu1` with `--nprocs 16` for `tglc epsfs` unless a
new timing test changes that default. Do not assume `pdogpu6` is faster for `epsfs` until a
`cupy`-enabled TGLC environment is available there.

The validated `pdogpu1` import pattern for the user-owned fork is:

```bash
deactivate 2>/dev/null
env \
  LD_LIBRARY_PATH=/pdo/app/anaconda/anaconda2-4.4.0/lib:/pdo/app/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH} \
  PYTHONPATH=/pdo/users/tehan/tess-gaia-light-curve-twirl:${PYTHONPATH} \
  /sw/qlp-environment/.venv/bin/python \
  -c "import tglc; from tglc.utils._optional_deps import HAS_CUPY; print(tglc.__file__); print(HAS_CUPY)"
```

### Existing PDO template tree

The existing deep-catalog production tree that TWIRL should mirror is:

```text
/pdo/users/tehan/tglc-deep-catalogs/
```

The checked template orbit is:

```text
/pdo/users/tehan/tglc-deep-catalogs/orbit-195/ffi/
```

What was confirmed there:

- `catalogs` is complete, with `32` files in `orbit-195/ffi/catalogs/`
- `cutouts` is complete, with `3136` `source_*.pkl` files under `orbit-195/ffi/`
- the tree shape matches the MIT `Manifest` layout, including `catalogs/`, `camX/ccdY/`, and `run/`

This is the safest operational template for new TWIRL benchmark runs because it writes into the same user-owned area without touching shared PDO staging more than necessary.

The concrete user-owned tree observed for the Sector `56` benchmark currently looks like:

```text
/pdo/users/tehan/tglc-deep-catalogs/
  orbit-119/
    ffi/
      catalogs/
        Gaia_cam4_ccd1.ecsv
        TIC_cam4_ccd1.ecsv
      cam4/
        ccd1/
          ffi/
            hlsp_tica_tess_ffi_s0056-o1-......-cam4-ccd1_tess_v01_img.fits -> /pdo/qlp-data/tica-delivery/s0056/cam4-ccd1/...
          source/
            source_*.pkl
  orbit-120/
    ffi/
      catalogs/
        Gaia_cam4_ccd1.ecsv
        TIC_cam4_ccd1.ecsv
      cam4/
        ccd1/
          ffi/
            hlsp_tica_tess_ffi_s0056-o2-......-cam4-ccd1_tess_v01_img.fits -> /pdo/qlp-data/tica-delivery/s0056/cam4-ccd1/...
```

### PDO write boundary

For TWIRL PDO work, treat `/pdo/users/tehan/` as the only writable area.

- `/pdo/qlp-data/...` and other shared PDO trees are read-only inputs for TWIRL purposes.
- Do not create, edit, overwrite, move, or delete anything outside `/pdo/users/tehan/`.
- If TICA FFIs or other shared inputs are needed for a benchmark run, expose them under `/pdo/users/tehan/...` with user-owned staging or symlinks rather than changing the shared source tree.

### Sector 56 benchmark unit for WD 1856

For the fixed Sector `56` benchmark:

- WD 1856 is on `cam4/ccd1`
- Sector `56` maps to orbit `119` and orbit `120`

So the first benchmark unit is not the whole sector. It is exactly these two orbit-camera-CCD jobs:

- orbit `119`, `ccd 4,1`
- orbit `120`, `ccd 4,1`

### Known-good `catalogs` command pattern on PDO

Use the shared QLP Python explicitly and keep the benchmark outputs under `/pdo/users/tehan/tglc-deep-catalogs/`.

Orbit `119`:

```bash
mkdir -p /pdo/users/tehan/tglc-deep-catalogs/orbit-119/ffi/catalogs
deactivate 2>/dev/null
env LD_LIBRARY_PATH=/pdo/app/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH} \
  /sw/qlp-environment/.venv/bin/python \
  -m tglc catalogs \
  --orbit 119 \
  --ccd 4,1 \
  --max-magnitude 20 \
  --nprocs 16 \
  --tglc-data-dir /pdo/users/tehan/tglc-deep-catalogs
```

Orbit `120`:

```bash
mkdir -p /pdo/users/tehan/tglc-deep-catalogs/orbit-120/ffi/catalogs
deactivate 2>/dev/null
env LD_LIBRARY_PATH=/pdo/app/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH} \
  /sw/qlp-environment/.venv/bin/python \
  -m tglc catalogs \
  --orbit 120 \
  --ccd 4,1 \
  --max-magnitude 20 \
  --nprocs 16 \
  --tglc-data-dir /pdo/users/tehan/tglc-deep-catalogs
```

The expected output files for each orbit are:

```text
/pdo/users/tehan/tglc-deep-catalogs/orbit-<orbit>/ffi/catalogs/Gaia_cam4_ccd1.ecsv
/pdo/users/tehan/tglc-deep-catalogs/orbit-<orbit>/ffi/catalogs/TIC_cam4_ccd1.ecsv
```

### Operational distinction between PDO trees

Keep these two roles separate:

- `/pdo/qlp-data/tica-delivery/s0056/cam4-ccd1/...` is the raw TICA FFI delivery tree used as the Sector `56` benchmark input source
- `/pdo/qlp-data/sector-56/ffi/cam4/ccd1/...` is a downstream QLP products tree (`BLS`, `DI`, `GPU`, `LC`, `MCMC`, `REPORTS`), not the raw TGLC `cutouts` input tree
- `/pdo/users/tehan/tglc-deep-catalogs/...` is the user-owned TGLC production tree currently serving as the safest TWIRL benchmark root

For the WD 1856 benchmark, first mirror the existing `tglc-deep-catalogs` arrangement and only then proceed to `cutouts`, `epsfs`, and `lightcurves`. In particular, keep all new writes under `/pdo/users/tehan/...` and treat `/pdo/qlp-data/...` as read-only.

### Known-good `cutouts` staging pattern on PDO

For Sector `56`, the raw TICA FFIs for WD 1856's benchmark CCD were found at:

```text
/pdo/qlp-data/tica-delivery/s0056/cam4-ccd1/
```

The orbit-local `ffi/` directories under the user-owned benchmark tree can be populated with symlinks:

```bash
mkdir -p /pdo/users/tehan/tglc-deep-catalogs/orbit-119/ffi/cam4/ccd1/ffi
mkdir -p /pdo/users/tehan/tglc-deep-catalogs/orbit-120/ffi/cam4/ccd1/ffi
```

```bash
for f in /pdo/qlp-data/tica-delivery/s0056/cam4-ccd1/*-o1-*-cam4-ccd1_tess_v01_img.fits; do
  ln -sfn "$f" /pdo/users/tehan/tglc-deep-catalogs/orbit-119/ffi/cam4/ccd1/ffi/
done
```

```bash
for f in /pdo/qlp-data/tica-delivery/s0056/cam4-ccd1/*-o2-*-cam4-ccd1_tess_v01_img.fits; do
  ln -sfn "$f" /pdo/users/tehan/tglc-deep-catalogs/orbit-120/ffi/cam4/ccd1/ffi/
done
```

The orbit `119` `o1` staging count was `5638` symlinks in
`/pdo/users/tehan/tglc-deep-catalogs/orbit-119/ffi/cam4/ccd1/ffi/`.

The corresponding `cutouts` command pattern is:

```bash
deactivate 2>/dev/null
env LD_LIBRARY_PATH=/pdo/app/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH} \
  /sw/qlp-environment/.venv/bin/python \
  -m tglc cutouts \
  --orbit 119 \
  --ccd 4,1 \
  --nprocs 16 \
  --tglc-data-dir /pdo/users/tehan/tglc-deep-catalogs
```

Known orbit `119` `cutouts` warnings from the Sector 56 benchmark:

- one TICA file was skipped as invalid because the FITS header lacked `COARSE`:
  `/pdo/qlp-data/tica-delivery/s0056/cam4-ccd1/hlsp_tica_tess_ffi_s0056-o1-00694257-cam4-ccd1_tess_v01_img.fits`
- `5637 cadence gaps != 1 detected` was reported after sorting by timestamp

Despite those warnings, orbit `119` `cutouts` completed and wrote `196/196` source pickle files to
`/pdo/users/tehan/tglc-deep-catalogs/orbit-119/ffi/cam4/ccd1/source/`.

These warnings should be revisited during benchmark QA, but they do not block continuing to `epsfs` after both orbit `119` and orbit `120` finish `cutouts`.

### User-owned fork for pre-Sector-67 TICA quality headers

Sector `56` TICA FFIs predate the `FINE`, `COARSE`, `RW_DESAT`, and `STRAYLT1-4` headers that were
added starting in Sector `67`, so the stock MIT reader treated missing quality headers as invalid
FFIs.

The current TWIRL benchmark therefore uses the user-owned fork at:

```text
/pdo/users/tehan/tess-gaia-light-curve-twirl/
```

That fork carries a pre-Sector-67 fallback patch in `tglc/ffi.py` so missing TICA quality headers
default to zero only for `sector < 67`. On PDO, expose that fork ahead of the installed package with
`PYTHONPATH`:

```bash
deactivate 2>/dev/null
env \
  LD_LIBRARY_PATH=/pdo/app/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH} \
  PYTHONPATH=/pdo/users/tehan/tess-gaia-light-curve-twirl:${PYTHONPATH} \
  /sw/qlp-environment/.venv/bin/python \
  -m tglc cutouts \
  --orbit 120 \
  --ccd 4,1 \
  --nprocs 16 \
  --replace \
  --tglc-data-dir /pdo/users/tehan/tglc-deep-catalogs
```

### Live PDO process hygiene

When rerunning a benchmark stage after code changes, check for stale stopped jobs before trusting the current outputs.

In this benchmark, a pre-patch orbit `120` `cutouts` process remained stopped in `pts/10`:

```text
PID 71335  STAT Tl  /sw/qlp-environment/.venv/bin/python -m tglc cutouts --orbit 120 --ccd 4,1 --nprocs 16 --tglc-data-dir /pdo/users/tehan/tglc-deep-catalogs
```

while the patched rerun was active in a user-owned tmux session:

```text
PID 91844  STAT Rl+  ... PYTHONPATH=/pdo/users/tehan/tess-gaia-light-curve-twirl:${PYTHONPATH} ... -m tglc cutouts --orbit 120 --ccd 4,1 --nprocs 16 --replace --tglc-data-dir /pdo/users/tehan/tglc-deep-catalogs
```

Do not resume stale pre-patch jobs such as the stopped `pts/10` process, because they can overwrite user-owned benchmark products with old code.

For the active benchmark, `pdogpu6` orbit `120` `epsf/` and `LC/` were temporarily made read-only in
the main benchmark tree so the orbit `119` driver cannot overwrite orbit `120` products while a
separate `pdogpu1` CPU-only `epsfs` benchmark is running in:

```text
/pdo/users/tehan/tglc-deep-catalogs-pdogpu1-epsf-test/
```

## TWIRL-Specific Operating Advice

### Use a WD control table to drive production

Do not launch jobs directly from the raw WD catalog. First build a control table that maps each WD to:

- TIC ID
- Gaia DR3 ID
- orbit
- camera/CCD
- cutout

Then collapse that table to the unique orbit-camera-CCD jobs needed for production.

In the current TWIRL repo, the relevant planning products are:

- `data_local/catalogs/twirl_master_catalog/twirl_wd_tess_detector_summary_v0.csv`
- `data_local/catalogs/twirl_master_catalog/tess_detector_target_tables_v0/`

The summary table is the right input for job selection; the per-orbit/camera/CCD target tables are
the TWIRL-side reference for what each planned job is expected to contain.

### Keep the MIT layout, add TWIRL metadata on top

Do not change the MIT file structure unless necessary. Instead:

- keep the raw TGLC products where the MIT code expects them,
- add TWIRL index tables and metadata files alongside them or in a separate repo-controlled location.

### Filter final outputs to the WD target list when possible

The `lightcurves` subcommand accepts `--tic` IDs. That is useful when:

- you have already identified the WD TIC IDs in a given CCD,
- you want to avoid writing light curves for unrelated stars.

This only helps at the final stage. Cataloging, cutouts, and ePSF fitting still need the full field model.

### Start with WD 1856+534 and a small pilot sample

Before launching full production, require the new workflow to:

1. regenerate a strong single-sector detection of WD 1856+534 b,
2. process a small sample of quiet WDs,
3. produce a first QA summary versus magnitude and crowding.

## Old-To-New Task Mapping

| If you used to do this | Use this now |
| --- | --- |
| `tglc.quick_lc.tglc_lc()` | `tglc all` or the four CLI stages |
| `ffi_cut()` | `tglc cutouts` |
| `target_lightcurve.epsf()` | `tglc epsfs` then `tglc lightcurves` |
| `run.py` batch logic | orbit-camera-CCD job wrappers around the CLI |
| inspect FITS `cal_aper_flux` / `cal_psf_flux` | read HDF5 aperture photometry outputs |
| `save_aper=True` for extra aperture information | use the built-in small, primary, and large apertures in the HDF5 product |

## Known Constraints To Plan Around

- The MIT fork is designed for TICA FFIs, not arbitrary remote download flows.
- The default target magnitude limits are not WD-ready.
- The product format changed from FITS to HDF5.
- The public CLI no longer exposes the old `prior` option.
- Light-curve emission is TIC-based, so Gaia-only WD targets may require a small extension.

## Recommendation For TWIRL

Use the MIT fork as the Stage 1 production engine, but keep TWIRL-specific logic outside it whenever possible:

- WD catalog assembly
- orbit/camera/CCD scheduling
- HDF5 indexing and consolidation
- QA
- ML detection
- injection-recovery
- candidate vetting

Only modify the MIT fork if TWIRL truly needs one of these missing capabilities:

- Gaia-first target emission without TIC IDs
- deeper or custom target selection logic
- reintroduction of a user-controlled floating-star prior
- extra output fields needed by the downstream detection stack
