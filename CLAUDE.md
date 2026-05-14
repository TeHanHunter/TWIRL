# CLAUDE.md

## Project

TWIRL (TESS White dwarf Investigation of Remnant pLanets) searches for transiting / occulting objects around white dwarfs using 200 s TESS FFIs (Sector ≥ 56). Benchmark target: **WD 1856+534**, Sector 56 orbits 119/120, cam4/ccd1.

## Read First

Before any task, consult `doc/`:

- `twirl_plan.md` — forward plan; `## Current Status Snapshot`, the active sprint section, and `## External Data Policy` for local data paths
- `twirl_progress_log.md` — dated history; active subsections end with `**Next:**`
- `ideas.md` — open questions / referee risks; starts with `## Critical Before Next Stage`
- `mit_tglc_usage_guide.md` — PDO operational details (envs, per-orbit pipeline, recovery audits, Sector < 67 flag handling)

No formal test suite. Validation is manual via script outputs and checkpoint files.

## Key Constants

- `FIRST_200S_SECTOR = 56`
- High-confidence WD cut: `Pwd > 0.75` (~359,073 targets)
- Authoritative target identifier: **Gaia DR3 `source_id`** (TIC is secondary)

## Layout

### Code Layout

- **`src/twirl/`** — importable production package
  - `catalogs/` — WD master catalog + TESS sector/orbit/camera/CCD coverage mapping
  - `io/hlsp.py` — HLSP FITS reader (both `hlsp_qlp_*` and `hlsp_twirl_*` schemas)
  - `lightcurves/` — TGLC HDF5 reader, flux-space cotrend, HLSP writer
  - `plotting/style.py` — publication-quality matplotlib templates
  - `search/bls.py` etc. — BLS-based search engine + candidate consolidation
  - `vetting/` — heuristic + LEO-Vetter adapters + centroid on-target test
- **`scripts/stage{1_lightcurves,2_search,5_validation}/`** — thin argparse CLI drivers over `src/twirl/`; `stage3_injections/` and `stage4_search/` planned but not built
- **`data_local/`** — large external + staged files, gitignored (layout in `doc/twirl_plan.md` § External Data Policy)
- **`doc/`** — `twirl_plan.md`, `twirl_progress_log.md` (append-only), `mit_tglc_usage_guide.md`, `ideas.md`

### Data Flow (Stage 1)

```
GaiaEDR3_WD_main.fits (external seed, data_local/)
  → catalogs.master_catalog  → twirl_wd_master_catalog_v0.fits + JSON sidecar
  → catalogs.tess_coverage   → enriched catalog with tess_observations_json column
  → TIC merge (PDO-only)     → tic_id / tic_match_status columns
  → detector tables          → per-orbit job lists for MIT TGLC
  → MIT TGLC (on PDO)        → HDF5 light curves
  → twirl flux_space_detrend → TWIRL HLSP FITS (hlsp_twirl_*)
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
  → twirl build_twirl_hlsp (flux-space cotrend, FITS)
```

Two separate PDO environments are used:
- **TGLC env** (`/sw/qlp-environment/.venv/bin/python` + TGLC fork via PYTHONPATH): for all `tglc` stages; has astropy 7.1 (`coordinates.representation` is a package)
- **QLP env** (`/pdo/app/qlp-environment/.venv/bin/python`, `qlp==0.13.2`): for legacy `lctools detrend` and `lctools hlsp`

For Sector < 67 (including Sector 56), hlsp uses `--flag-type spoc --flag-source fits`. SPOC flags at `/pdo/qlp-data/spocflags/`.

**HLSP operates at the sector level**, not per-orbit. `qlp/lctools/bin/hlsp.py` stacks all orbits of a sector into a single FITS per TIC. All orbits in a sector must be detrended before HLSP runs.

**Catalog reuse across orbits in the same sector.** TGLC's `catalogs` stage keys queries on `(sector, camera, ccd)`, so sibling orbits produce identical `Gaia_camC_ccdD.ecsv` / `TIC_camC_ccdD.ecsv`. The stage has native skip-if-exists; symlinking a completed sibling's catalogs lets dense-field CCDs skip multi-hour Gaia queries. `run_tglc_orbit_pipeline.py --reuse-catalogs-from-orbit N` does this automatically per CCD. Only pass sibling orbits in the same sector.

**Standard recovery rate check** after each production run: requested TIC IDs → h5 count → FITS count, binned by Tmag. Gaia-only targets (no TIC bridge) are always 0% recovered until the fork is extended.

**Mandatory cadence-alignment check before legacy QLP detrend.** Before running `qlp lctools detrend` on any new set of TGLC h5s, verify that the TGLC `LightCurve/Cadence` array is a subset of the camera's `camC_quat.txt` cadence column for that orbit. Any cadence present in TGLC but missing from the quat file will crash `quatspline.fit_trend` with a shape mismatch (first seen Sector 56 cam3/orbit-120: TGLC kept cadence `699957` that the quat file drops, breaking every cam3 LC). Fix with `scripts/stage1_lightcurves/fix_tglc_quat_cadence_mismatch.py` — it trims the TGLC-only cadence(s) from every equal-length dataset in each h5. The quat file is authoritative for spacecraft pointing; never widen it to match TGLC. (The new TWIRL flux-space-detrend path in `build_twirl_hlsp.py` doesn't use the quat file.)

## Data Policy

- Large inputs live under `data_local/` — never commit them.
- No duplicate products with different filenames; use versioned names (`_v0`, `_v1`).
- Every output carries a JSON sidecar with seed-file path, size, mod-time, build version, SHA256, and key assumptions.

## Scientific Scope

- **Core survey**: 200 s TESS FFIs only (Sector ≥ 56)
- **Seed catalog**: Gentile Fusillo et al. (2021) Gaia EDR3 WD catalog
- **Authoritative target identifier**: Gaia DR3 source_id (TIC is secondary)
- **First-year regime**: large, deep, short-duration events (WD 1856-like)
- **Secondary goal (2026-05-01)**: WD 1145+017-style disintegrating-planet candidates — variable-depth, asymmetric, possibly aperiodic dust-tail dips. Search stack must include a depth-allowed-to-vary template detector in addition to BLS.
- Do not broaden science claims or merge content from TWIRL_proposal without asking.

## Working Conventions

- Plan first, then edit. Prefer minimal, local changes; don't refactor files unrelated to the task.
- Production code goes in `src/twirl/` and `scripts/`, not notebooks.
- `astropy.table.Table` for structured catalog data; FITS for catalogs, HDF5 for light curves, JSON for metadata.
- Plotting (style authority is [`src/twirl/plotting/style.py`](src/twirl/plotting/style.py); update the code first, then any docs that reference it):
  - Use `apply_twirl_style(template)` to pick up the templates (`column` 3.4×2.65 in, `full_page` 7.1×4.1 in). Don't hardcode Seaborn `rc` blocks per-script.
  - Use `get_ordered_palette(n)` for color sets; default is `viridis` ordered.
  - Default look: serif (`DejaVu Serif`), `whitegrid` theme, `paper` context, white background, light grey grid, dark grey axes, publication-oriented (not slide-oriented).
  - No figure or panel titles by default; only add them when the comparison would otherwise be ambiguous.
  - Vertically stacked comparison figures: per-panel axis labels (not one shared), tight inter-panel spacing, compact shared colorbars close to the panels.
  - Integer-valued color quantities use a stepped colorbar, not a continuous ramp.
  - Aitoff guide overlays: black dotted longitude guides, black dashed Galactic-plane line, longitude labels with a thin white bezel.
  - Legend border is black; place legends outside the data area when possible.
  - Save publication figures as both PDF and PNG. Rasterize only dense scatter layers, not the whole figure.
- Update `doc/twirl_progress_log.md` after completing a milestone.

## PDO Rules (any `pdo*` host)

- MIT-adapted TGLC fork lives on MIT PDO machines (not in this repo).
- TGLC CLI stages: `catalogs → cutouts → epsfs → lightcurves`.
- **Read-only outside `/pdo/users/tehan/`.** Never modify, delete, or move files under `/pdo/qlp-data/`, `/sw/`, or other users' dirs. If another user's file is needed in our tree, symlink from `/pdo/users/tehan/` to the upstream path — do not copy or edit upstream.
- Prefer `pdogpu1` for CPU-only benchmark + Stage-2 work; `pdogpu6` for GPU ePSF runs.
- PDO fork requires a patch for pre-Sector-67 TICA headers (see `doc/mit_tglc_usage_guide.md`).
- **Cap BLAS/OpenMP threads before any multiprocessing pool** (`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS=NUMEXPR_NUM_THREADS=1`). Without this, each pool worker spawns BLAS threads equal to the full core count (128 on pdogpu1) → ~8× oversubscription, 10–100× slowdown. QLP detrend/HLSP wrappers bake these in at module load time.
- `tglc epsfs` GPU path: the installed `tglc 1.0.2` is GPU-enabled by default (`tglc/scripts/epsfs.py:78` gates on `if use_gpu and HAS_CUPY:`; `--no-gpu` is opt-out). But `cupy` is not installed in any shared PDO Python env, so by default `HAS_CUPY=False` and runs silently fall back to CPU. To get GPU ePSFs, use the TWIRL-local venv at `/pdo/users/tehan/twirl-gpu-venv` (built `2026-04-27` with `cupy-cuda11x`, `--system-site-packages`) and set `LD_LIBRARY_PATH=/sw/python-versions/python-3.11.9/lib` plus `PYTHONPATH=/pdo/users/tehan/tess-gaia-light-curve-twirl:/sw/qlp-environment/.venv/lib/python3.11/site-packages`. Benchmark on `pdogpu6`: S56 cam4/ccd1 finished in `30:42` (`-n 4`) vs `4:01:31` CPU baseline → `~7.9×` speedup; outputs match CPU within float roundoff. CUDA driver determines the wheel — `pdogpu6` is CUDA 11.8 ⇒ `cupy-cuda11x`.

## End-of-Day Wrap

When the user says **"wrap for the day"**:

1. Update `doc/twirl_plan.md` `## Current Status Snapshot` to current state.
2. Append dated notes to `doc/twirl_progress_log.md` for this session's work.
3. Commit all tracked changes with a concise message.
4. Push to remote.
5. Leave the next concrete step in the progress log as a `**Next:**` line.

## Stop and Ask Before

- Selecting a new benchmark sector
- Broadening science claims beyond the WD 1856-like regime
- Merging content from `TWIRL_proposal`
- Committing to a follow-up instrument
- Changing repo scope
- **Changing collaboration scope or paper-leadership commitments.** TWIRL is now a Schwamb-group collaboration (`2026-05-13`); ownership division and paper leadership are tracked in `doc/twirl_plan.md` §Collaboration & Ownership. Don't unilaterally promise contributions to or away from this division.
