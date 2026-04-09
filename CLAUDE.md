# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Purpose

TWIRL (TESS White dwarf Investigation of Remnant pLanets) is a scientific survey pipeline to search for transiting/occulting objects around white dwarfs using 200 s TESS FFIs (Sector ≥ 56), establish the first statistical occurrence-rate constraints on close-in WD planets (~1% precision), and produce validated candidates for post-main-sequence planetary survival studies. The benchmark target is **WD 1856+534** in Sector 56, orbits 119/120, cam4/ccd1.

## Read First

Before working on any task, consult these documents in `doc/`:
- `twirl_plan.md` — forward-looking executable plan; each stage has a status marker and the `## Current Status Snapshot` section shows where things stand
- `twirl_progress_log.md` — dated execution history; each active subsection ends with `**Next:**` showing the immediate action
- `ideas.md` — open science questions and referee risks; starts with `## Critical Before Next Stage`
- `mit_tglc_usage_guide.md` — PDO-specific operational details
- `local_data.md` — local data path conventions (data not in git)

## Environment Setup

```bash
python -m pip install -r requirements.txt
# On MIT PDO only:
python -m pip install pyticdb
```

There is no formal test suite yet. Validation is manual via script outputs and checkpoint files.

## Architecture

### Pipeline Stages
1. **Stage 1** (in progress): Generate WD light curves on MIT PDO with TGLC
2. **Stage 2** (planned): Build search stack (periodic & dip searches)
3. **Stage 3** (planned): Injection-recovery tests for completeness
4. **Stage 4** (planned): Full-sample search
5. **Stage 5** (planned): Validation and follow-up coordination

### Code Layout

- **`src/twirl/`** — importable production package
  - `catalogs/master_catalog.py` — builds the WD master catalog from the Gaia EDR3 seed FITS
  - `catalogs/tess_coverage.py` — maps TESS sector/orbit/camera/CCD coverage per target
  - `plotting/style.py` — publication-quality matplotlib RC params and figure templates
- **`scripts/stage1_lcs/`** — CLI wrappers for Stage 1 production steps; each script is a thin argparse driver over `src/twirl/` functions
- **`data_local/`** — large files not tracked in git (see `doc/local_data.md`)
- **`doc/`** — project documentation and planning; `twirl_progress_log.md` is append-only

### Data Flow (Stage 1)

```
GaiaEDR3_WD_main.fits (external seed, data_local/)
  → catalogs.master_catalog  → twirl_wd_master_catalog_v0.fits + JSON sidecar
  → catalogs.tess_coverage   → enriched catalog with tess_observations_json column
  → TIC merge (PDO-only)     → tic_id / tic_match_status columns
  → detector tables          → per-orbit job lists for MIT TGLC
  → MIT TGLC (on PDO)        → HDF5 light curves
```

### Import Pattern in Scripts

All scripts add `src/` to `sys.path` at the top rather than installing the package editably:

```python
REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))
```

### Full Stage 1 Pipeline Per orbit/camera/CCD

```
tglc catalogs → tglc cutouts → tglc epsfs → tglc lightcurves (HDF5)
  → qlp lctools detrend → qlp lctools hlsp (FITS)
```

Two separate PDO environments are used:
- **TGLC env** (`/sw/qlp-environment/.venv/bin/python` + TGLC fork via PYTHONPATH): for all `tglc` stages
- **QLP env** (`/pdo/app/qlp-environment/.venv/bin/python`, `qlp==0.13.2`): for `detrend` and `hlsp`

Activate the QLP env with `source scripts/activate_qlp_env.sh`, then use `qlp_run lctools detrend/hlsp`. The system PYTHONPATH includes FFITools `qlp==0.1` which shadows the correct version — the activation script fixes this by setting PYTHONPATH to the TGLC fork only.

The detrend step must run before hlsp — it writes `bestdmagkey` into HDF5. For Sector < 67 (including Sector 56), hlsp uses `--flag-type spoc --flag-source fits`. SPOC flags are at `/pdo/qlp-data/spocflags/` and confirmed present for all Sector 56 CCDs.

**Standard recovery rate check** after each production run: requested TIC IDs (detector table) → h5 count → FITS count, binned by Tmag. Gaia-only targets (no TIC bridge) are always 0% recovered until the fork is extended.

### Key Constants

- `FIRST_200S_SECTOR = 56` — only 200 s FFIs from this sector onward
- Default high-confidence WD cut: `Pwd > 0.75` (~359,073 targets)
- PDO defaults: `--max-magnitude 20`, `--nprocs 16` for epsfs stage

## Data Policy

- Large inputs live under `data_local/` — never commit them
- Do not create duplicate products with different filenames; use versioned names (`_v0`, `_v1`)
- Every output must carry a JSON sidecar with: seed-file path, size, mod-time, build version tag, SHA256 hash, and key assumptions

## Scientific Scope

- **Core survey**: 200 s TESS FFIs only (Sector ≥ 56)
- **Seed catalog**: Gentile Fusillo et al. (2021) Gaia EDR3 WD catalog
- **Authoritative target identifier**: Gaia DR3 source_id (TIC is secondary)
- **First-year regime**: large, deep, short-duration events (WD 1856-like)
- Do not broaden science claims or merge content from TWIRL_proposal without asking

## Working Conventions

- Plan first, then edit. Prefer minimal, local changes.
- Do not refactor files unrelated to the current task.
- Production code goes in `src/twirl/` and `scripts/`, not notebooks.
- Use `astropy.table.Table` for structured catalog data; FITS for catalogs, HDF5 for light curves, JSON for metadata.
- Plotting: use `src/twirl/plotting/style.py` templates (`column` 3.4×2.65 in, `full_page` 7.1×4.1 in); no figure titles by default.
- Update `doc/twirl_progress_log.md` after completing a milestone.

## PDO-Specific Notes

- MIT-adapted TGLC fork lives on MIT PDO machines (not in this repo)
- TGLC CLI stages: `catalogs → cutouts → epsfs → lightcurves`
- Prefer `pdogpu1` for CPU-only benchmark work; keep outputs under `/pdo/users/tehan/`
- PDO fork requires a patch for pre-Sector-67 TICA headers (see `doc/mit_tglc_usage_guide.md`)

## End-of-Day Wrap

When the user says "wrap for the day", do all of the following automatically:
1. Update `doc/twirl_plan.md` `## Current Status Snapshot` to reflect current milestone state.
2. Append dated execution notes to `doc/twirl_progress_log.md` for anything completed this session.
3. Commit all tracked changes with a concise message.
4. Push to remote.
5. Leave the next concrete step clear in the progress log (`**Next:**` line).

## Stop and Ask Before

- Selecting a new benchmark sector
- Broadening science claims beyond the WD 1856-like regime
- Merging content from TWIRL_proposal repo
- Committing to a follow-up instrument
- Changing repo scope
