# AGENTS.md

This file is the operating protocol for agents working in the standalone `TWIRL` repository.

## Read First

Before making substantive changes, read:

- `doc/twirl_plan.md`
- `doc/twirl_progress_log.md`
- `doc/a2v1_production_protocol.md` before planning or running Stage 1 production
- `doc/mit_tglc_usage_guide.md`
- `doc/ideas.md`
- `doc/science_background.md` before making scientific background or literature claims
- `doc/local_data.md`
- `doc/orcd_h200_usage.md` before planning or running ORCD/H200 jobs
- `doc/plotting_style.md` before making or revising publication-facing figures

Use `doc/twirl_plan.md` as the single forward-looking project plan, `doc/twirl_progress_log.md` for dated execution history and benchmark notes, `doc/a2v1_production_protocol.md` for the accepted Stage 1 production contract, `doc/ideas.md` for unresolved questions and brainstorming, `doc/local_data.md` for local-only data conventions, and `doc/plotting_style.md` for figure styling.
Use `doc/orcd_h200_usage.md` for ORCD/Engaging partition names, storage, Slurm snippets, and the PDO-vs-ORCD boundary.

Follow the Stage 1-5 pipeline structure defined in `doc/twirl_plan.md`.

If this file and those docs disagree, treat the docs as authoritative and update `AGENTS.md`.

## Repo Identity

- This repo is a standalone project and must remain separate from `TWIRL_proposal`.
- Do not pull proposal-only files, plots, or git history into this repo unless the user explicitly asks.
- Keep this repo focused on the executable survey pipeline, not proposal drafting artifacts.
- The active survey-paper manuscript is the separate sibling repo
  `/Users/tehan/PycharmProjects/twirl-survey-paper`, intended for Overleaf
  git sync. Do not vendor the paper repo into `TWIRL`.
- `TWIRL I: A Systematic TESS Search for Transiting Planetary Remnants
  around White Dwarfs` is the current framework/overview paper. It is not the
  final occurrence-rate paper and not an individual discovery paper.

## Working Style

- Plan first, then edit.
- Prefer minimal, local changes.
- Do not refactor unrelated files.
- Production code belongs in `src/twirl/` and `scripts/`, not notebooks.
- For structured catalog data, prefer `astropy.table.Table`; use FITS for catalogs, HDF5 for light curves, and JSON for metadata/provenance.

## Key Constants And Layout

- `FIRST_200S_SECTOR = 56`.
- High-confidence WD reference cut: `Pwd > 0.75`.
- Authoritative target identifier: Gaia DR3 `source_id`; TIC is operational metadata.
- `src/twirl/` is the importable production package:
  - `catalogs/` builds the WD master catalog and sector/orbit/camera/CCD coverage maps.
  - `io/hlsp.py` reads both `hlsp_qlp_*` and `hlsp_twirl_*` FITS schemas.
  - `lightcurves/` reads TGLC HDF5, detrends in flux space, and writes TWIRL HLSP FITS.
  - `plotting/style.py` is the plotting style authority.
  - `search/` contains BLS/search tooling and candidate consolidation.
  - `vetting/` contains heuristic vetting, LEO-Vetter adapters, and centroid checks.
- Add reusable `injections/`, `review/`, or `models/` package boundaries only
  as coherent workflows migrate out of overloaded Stage 5 drivers; do not
  create duplicate implementations solely to populate the target layout.
- `scripts/stage{1_lightcurves,2_search,3_injections,4_search,5_validation}/`
  is the target driver layout. Stage 4 is not yet built.
- Existing Stage 5 experiment drivers are migrated incrementally. Preserve their
  paths while PDO/ORCD workflows depend on them, but put new reusable logic in
  `src/twirl/` and new full-survey inference drivers in `stage4_search/`.
- The package is installable from `pyproject.toml`, and pytest resolves `src/`
  without an environment-specific `PYTHONPATH`. New local code should import
  the package normally. Preserve direct `sys.path` bootstraps only in active
  remote/compatibility wrappers until their callers migrate.

## Local Data Policy

- Large external inputs belong under `data_local/` or another user-configured local path, not in the git-managed repo root.
- Do not commit raw FITS catalogs or staged survey data.
- On ORCD, compact TWIRL exports and downstream results belong under `/orcd/data/mki_aryeh/001/twirl/`; use ORCD scratch only for active job-local temporary data.
- Never initiate Duo, password, or keyboard-interactive authentication from
  Codex, scripts, or automated probes. Use only non-interactive SSH checks; if
  ORCD access is unavailable, stop and ask the user to open the authenticated
  control socket from their own terminal.
- If code depends on a local catalog or local TGLC staging area, document the path convention and record provenance in outputs.
- Prefer concise, stable filenames for accepted master-catalog states, but do not duplicate large integrated FITS products solely to create cleaner names when the existing file already represents the accepted state.

## Scientific Scope

- TWIRL searches for any transiting or occulting objects around white dwarfs.
- The core survey uses `200 s` TESS FFIs only.
- That means production work should target `Sector >= 56`.
- Do not propose or hardcode sectors `1-55` as the primary TGLC survey input.
- First-year science is optimized for large, deep, short-duration events in the WD 1856-like regime.
- Do not overclaim Earth-size or habitable-zone occurrence constraints before completeness is demonstrated.
- For TWIRL I writing, describe the pipeline, survey framework, benchmark,
  validation plan, and completeness strategy. Do not write occurrence-rate
  constraints, final null-result claims, or discovery claims before the locked
  parent sample, end-to-end completeness, and candidate validation support them.

## Parent Sample Definition

- The current seed WD catalogue is the local Gentile Fusillo et al. (2021) main Gaia EDR3 catalogue, recommended at `data_local/catalogs/GaiaEDR3_WD_main.fits`.
- It is based on Gentile Fusillo et al. (2021), the Gaia EDR3 white dwarf catalogue.
- Treat Gaia DR3 as the authoritative TWIRL target identifier.
- Treat `Pwd > 0.75` as the default high-confidence reference sample for pilot work and comparisons.
- The exact TWIRL WD denominator is still to be determined later.
- Keep two concepts separate:
  - statistical survey samples derived from the locked high-confidence parent sample
  - broader exploratory searches that do not automatically feed occurrence-rate claims
- Record the exact parent-sample cuts in code and output metadata.

## Benchmark And Search Strategy

- `WD 1856+534` is the mandatory benchmark target.
- The current fixed first benchmark is `Sector 56`, orbits `119` and `120`, `cam4/ccd1`, chosen explicitly by the user.
- Do not silently change the benchmark sector/orbit/CCD unless the user asks.
- The first search version should be interpretable:
  - periodic short-duration box or trapezoid search
  - separate dip-search branch for non-periodic or weakly periodic events
- Machine learning is optional triage, not the first discovery engine.

## MIT TGLC Operating Assumptions

- Stage 1 uses the MIT-adapted TGLC fork on MIT PDO machines.
- ORCD/H200 is downstream compute for compact exports, injection-recovery, GPU search, feature extraction, and later ML triage; do not move primary Stage 1 TGLC/ePSF production there by default.
- On ORCD, default GPU jobs to `1` H200 and use at most `2` H200s unless the
  user explicitly approves a larger request for a specific run. Do not occupy
  the full H200 node for routine smokes, BLS/ranker jobs, or exploratory tests.
- Treat the MIT fork as an orbit/camera/CCD production pipeline, not a `quick_lc.py` replacement.
- Assume:
  - TICA FFIs are pre-staged on disk
  - `pyticdb` access to `tic_82` and `gaia3`
  - MIT `tglc-data` directory layout
- Use the MIT CLI stages: `catalogs`, `cutouts`, `epsfs`, `lightcurves`, or `all`.
- For WD survey production, do not impose a science target magnitude cap at TGLC catalog construction. The TWIRL wrappers pass an effectively unbounded catalog limit (`--max-magnitude 99`) so faint requested WD TICs are not dropped by the MIT fork's bright-star defaults. Smaller limits are only for controlled smoke or diagnostic runs.
- If `--max-magnitude 99` makes TIC catalog queries too slow or too large, solve that by adding requested-TIC inclusion to the catalog bridge, not by reintroducing a science magnitude cap.
- For the current CPU-only benchmark, use `--nprocs 16` as the default one-CCD `epsfs` worker count unless a new timing test motivates a change.
- Prefer `pdogpu1` for CPU-only benchmark and Stage 2 work; use `pdogpu6` for GPU ePSF production when the `cupy` environment is active.
- Preserve the MIT raw data layout. Add TWIRL metadata and indices on top rather than inventing a parallel raw layout.
- On PDO, keep TWIRL-run outputs, symlinks, staging, and scratch products under `/pdo/users/tehan/` rather than writing into shared PDO trees.
- The TGLC env for `tglc` stages is `/sw/qlp-environment/.venv/bin/python` with the TWIRL TGLC fork on `PYTHONPATH`; it has the newer Astropy layout needed by `tglc lightcurves`.
- The legacy QLP env for `lctools detrend` and `lctools hlsp` is `/pdo/app/qlp-environment/.venv/bin/python` with `qlp==0.13.2`.
- For Sector `< 67` including S56, HLSP uses `--flag-type spoc --flag-source fits`; SPOC flags live under `/pdo/qlp-data/spocflags/`.
- Cap BLAS/OpenMP threads before multiprocessing pools on PDO: `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS=NUMEXPR_NUM_THREADS=1`.
- The GPU ePSF path requires the TWIRL-local venv at `/pdo/users/tehan/twirl-gpu-venv` with `cupy-cuda11x`, `LD_LIBRARY_PATH=/sw/python-versions/python-3.11.9/lib`, and the TGLC fork resolving first on `PYTHONPATH`.

## Data Product Expectations

- MIT TGLC outputs HDF5 light curves, not the old FITS products.
- For TGLC extraction, treat the HDF5 `lightcurves` stage as the Stage 1 extraction checkpoint; public-style FITS products belong to the later QLP detrend + HLSP export path.
- For each completed A2v1 sector, the ADP/ADP015-only sector-level FITS export is required after both orbit HDF5 trees pass coverage and HDF5-openability validation. An HDF5-only state is an intermediate production checkpoint, not an accepted sector product.
- Downstream TWIRL code should expect decontaminated aperture photometry with `1x1`, `3x3`, and `5x5` apertures.
- Do not assume old columns such as `cal_psf_flux` or `cal_aper_flux`.
- TIC IDs are useful metadata and may still matter for the current MIT implementation, but they do not define the scientific sample.
- A major early engineering risk is TIC-oriented target emission:
  - audit how Gaia-selected targets propagate through the MIT fork
  - be prepared to recommend a Gaia-first extension if scientifically important WDs are missing
- QLP HLSP export operates at the sector level, not per orbit. All orbits in a sector must be detrended before legacy HLSP export runs.
- Before legacy QLP detrending, verify each TGLC `LightCurve/Cadence` array is a subset of the camera `camC_quat.txt` cadence column. If TGLC kept an orbit-boundary cadence missing from the quat file, trim the TGLC HDF5 with `scripts/stage1_lightcurves/fix_tglc_quat_cadence_mismatch.py`; the quat file is authoritative.

## Stage 1 Data Flow

```text
GaiaEDR3_WD_main.fits (external seed, data_local/)
  -> catalogs.master_catalog
  -> twirl_wd_master_catalog_v*.fits + JSON sidecar
  -> catalogs.tess_coverage
  -> TIC merge on PDO
  -> detector tables
  -> MIT TGLC HDF5 light curves
  -> TWIRL-FS flux-space detrend
  -> TWIRL-FS HLSP FITS (hlsp_twirlfs_*)
```

Current faint-end LC production planning should prioritize the adaptive TWIRL-FS branches used for human labeling and vetting: ADP (`DET_FLUX_ADP*`) and ADP015 (`DET_FLUX_ADP015*`). Keep canonical/default `DET_FLUX*` only where needed for compatibility or comparison until downstream readers and exports no longer require it.
The short name for the regenerated production product family is `A2v1`: no TIC magnitude cap (`--max-magnitude 99`), saturated-pixel ePSF masking, and ADP plus ADP015 saved for the `1x1`, `3x3`, and `5x5` apertures only. Use that name when referring to the S56 remake or future sectors following the same production rule.
To reuse existing prepared source pickles, prefer small `source_tic/source_X_Y.ecsv` TIC overlays plus a TGLC light-curve hook that assigns `source.tic` from the sidecar after unpickling. For requested TWIRL WD target emission, build those overlays from the TWIRL observation table detector coordinates (`--overlay-from-observations`) rather than regenerating max-99 TIC catalogs by default. Use max-99 TIC catalogs only when a broad all-TIC match table is required. Do not rewrite or regenerate large source pickles unless the overlay path fails validation.
For a queued A2v1 campaign, preserve sector order and require prepared-input preflight, edge-aware HDF5 coverage plus openability validation, FITS production, and full schema validation before advancing to the next sector. Stop the queue on any failed gate; do not silently skip or bypass a failed sector.

## Safety And Boundaries

- Do not modify raw input data.
- Do not change archive formats without updating the relevant schema docs.
- Do not introduce ML-first search changes unless the baseline search is already working.
- Do not change the survey denominator or `Pwd` threshold implicitly; make that an explicit documented decision.
- On PDO, never create, edit, overwrite, move, or delete anything outside `/pdo/users/tehan/`.
- Reading shared PDO locations such as `/pdo/qlp-data/` is allowed, but treat them as read-only inputs.
- If a PDO workflow needs files from a shared location, stage them under `/pdo/users/tehan/` with user-owned copies or symlinks rather than changing the shared source tree.

## Required QA And Validation Logic

- Run photometric QA before ML or large-scale search.
- The QA minimum includes:
  - RMS or MAD versus magnitude
  - sector-level failures
  - missing cadence statistics
  - completeness by target and orbit
  - WD 1856 recovery in the `200 s` products
  - aperture-to-aperture consistency
  - at least one independent extraction comparison on the benchmark set
- Injection-recovery must go through the real end-to-end stack:
  - search
  - classifier, if used
  - vetter
  - candidate merging

## Benchmark Acceptance

- Do not scale beyond the pilot stage until WD 1856 is a credible benchmark in at least one `Sector >= 56` product.
- Check the event timing against the published ephemeris.
- Check the event across the `1x1`, `3x3`, and `5x5` apertures.
- Investigate and document any major aperture disagreement.
- Compare against at least one independent extraction path on the benchmark set.

## Validation

- After code changes, run `make test-fast`.
- For detection changes, also run `make run-detection-sample`.
- After documentation or plan changes, run `make check-docs`.
- For catalog or index changes, validate the schema and at least one sample target.

## Outputs

- Put generated reports in `reports/`.
- Put reusable config in `configs/`.
- Record assumptions and provenance in output metadata.
- When progress is made against a milestone or deliverable, record it in the same turn if practical.
- Keep `doc/twirl_plan.md` as the forward-looking plan plus compact milestone status summaries only; do not let it accumulate command-level run logs.
- Record detailed dated execution history in `doc/twirl_progress_log.md` under the most relevant stage subsection.
- Do not add progress-log entries for purely aesthetic figure tweaks such as label nudges, tick colors, colorbar placement, or overlay opacity; record only scientifically meaningful figure changes such as changed inputs, model definitions, plotted quantities, smoothing assumptions, sample sizes, or recovery metrics.
- In `doc/twirl_plan.md`, keep only 1-3 milestone-level status bullets per active subsection and link to the relevant `doc/twirl_progress_log.md` section.
- Keep the prose part of each progress bullet readable; do not stuff raw long paths into the sentence.
- Keep file references inline when helpful, using short clickable Markdown link labels such as `[script]`, `[PNG]`, `[PDF]`, `[FITS]`, or `[CSV]` rather than raw long paths.
- Use plain Markdown that renders everywhere.
- Append a short status marker in parentheses to subsection headings when helpful: `(✓)` done, `(...)` in progress, `(?)` problematic, and nothing if not started.
- Treat these markers as manually tuned section status labels rather than auto-generated output.

## Current Priority Authority

Do not duplicate an ordered priority list in this protocol. The sole current
list is `doc/twirl_plan.md` under `## Immediate implementation priorities`.
Update that list when priorities change, and update this file only when the
change alters standing operating rules.

## Follow-Up Planning Assumptions

- Current default planning assumption is:
  - Magellan is the primary MIT-affiliated follow-up path
  - MMT is backup or collaborator-driven
- Julien is now part of TWIRL collaboration/follow-up planning; track details in `doc/twirl_plan.md` and `doc/twirl_progress_log.md`, and do not treat meeting-note instrument ideas as verified until access/proposal details are checked.
- Favor high-cadence follow-up planning for short predicted transit windows.
- If a task depends on a specific instrument or access policy, verify it before writing it into the repo as fact.
- Do not assume aperiodic events can be confirmed the same way as periodic candidates.
- Before follow-up is treated as ready, require independent re-extraction, contamination checks, and ephemeris adequacy for the planned facility.

## Code And Repo Conventions

- Prefer pipeline code, scripts, and structured metadata over notebook-only workflows.
- Treat `pyproject.toml` as the dependency and packaging authority;
  `requirements.txt` is only a compatibility entrypoint.
- Follow the staged repo layout in `doc/twirl_plan.md` as the codebase grows.
- Keep implementation decisions transparent and easy to audit.
- Prefer small smoke tests on pilot samples before scaling to full production.
- Avoid hidden scientific assumptions in code. Put important survey assumptions in docs or config.
- For publication-facing plots, use the shared style module in `src/twirl/plotting/style.py` rather than script-local Seaborn theme blocks.
- For publication-facing plots, do not add figure titles by default unless the user explicitly asks for them or the figure would otherwise be ambiguous.
- For publication-facing vertically stacked comparison figures, prefer per-panel axis labels, concise panel titles when needed for disambiguation, and compact nearby shared colorbars or legends rather than oversized figure-level decorations.
- For long-running scripts or batch jobs, add periodic progress reporting when practical so PDO runs are observable without attaching a debugger or guessing from zero-byte output files.
- When giving shell commands to the user, prefer simple multi-line fenced command blocks without prompt text or forced one-line compression unless the user explicitly asks for one-liners.
- When the user says to "wrap for the day", automatically do the end-of-day hygiene when the repo state is coherent: record milestone status in `doc/twirl_plan.md`, record dated execution notes in `doc/twirl_progress_log.md`, commit the current checkpoint, push it if possible, and leave the next concrete step clear.

## When To Stop And Ask

Ask the user before:

- selecting a benchmark sector on their behalf
- broadening the science claims beyond the current plan
- merging proposal-repo content into this repo
- committing to a follow-up instrument path that has not been verified
- changing the repo scope away from the standalone pipeline project
- changing collaboration scope or paper-leadership commitments. TWIRL is now a Schwamb-group collaboration; ownership and paper-leadership expectations are tracked in `doc/twirl_plan.md`.

## Maintenance Rule

When the project plan changes materially, update this file in the same turn if possible.
