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
tglc epsfs --orbit 185 --ccd 1,1 --nprocs 4 --tglc-data-dir /path/to/tglc-data
tglc lightcurves --orbit 185 --ccd 1,1 --nprocs 16 --tglc-data-dir /path/to/tglc-data
```

In TWIRL, the first safe wrapper around this pattern is
[`scripts/stage1_lcs/run_tglc_catalogs.py`](/Users/tehan/PycharmProjects/TWIRL/scripts/stage1_lcs/run_tglc_catalogs.py).
It reads the orbit-aware TWIRL detector summary and, by default, prints or writes the exact
`tglc catalogs` commands for the selected orbit/camera/CCD jobs before any execution happens.

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
