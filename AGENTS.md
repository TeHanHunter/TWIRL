# AGENTS.md

This file is the operating protocol for agents working in the standalone `TWIRL` repository.

## Read First

Before making substantive changes, read:

- `doc/twirl_plan.md`
- `doc/twirl_progress_log.md`
- `doc/mit_tglc_usage_guide.md`
- `doc/ideas.md`
- `doc/local_data.md`
- `doc/plotting_style.md` before making or revising publication-facing figures

Use `doc/twirl_plan.md` as the single forward-looking project plan, `doc/twirl_progress_log.md` for dated execution history and benchmark notes, `doc/ideas.md` for unresolved questions and brainstorming, `doc/local_data.md` for local-only data conventions, and `doc/plotting_style.md` for figure styling.

Follow the Stage 1-5 pipeline structure defined in `doc/twirl_plan.md`.

If this file and those docs disagree, treat the docs as authoritative and update `AGENTS.md`.

## Repo Identity

- This repo is a standalone project and must remain separate from `TWIRL_proposal`.
- Do not pull proposal-only files, plots, or git history into this repo unless the user explicitly asks.
- Keep this repo focused on the executable survey pipeline, not proposal drafting artifacts.

## Working Style

- Plan first, then edit.
- Prefer minimal, local changes.
- Do not refactor unrelated files.
- Production code belongs in `src/twirl/` and `scripts/`, not notebooks.

## Local Data Policy

- Large external inputs belong under `data_local/` or another user-configured local path, not in the git-managed repo root.
- Do not commit raw FITS catalogs or staged survey data.
- If code depends on a local catalog or local TGLC staging area, document the path convention and record provenance in outputs.
- Prefer concise, stable filenames for accepted master-catalog states, but do not duplicate large integrated FITS products solely to create cleaner names when the existing file already represents the accepted state.

## Scientific Scope

- TWIRL searches for any transiting or occulting objects around white dwarfs.
- The core survey uses `200 s` TESS FFIs only.
- That means production work should target `Sector >= 56`.
- Do not propose or hardcode sectors `1-55` as the primary TGLC survey input.
- First-year science is optimized for large, deep, short-duration events in the WD 1856-like regime.
- Do not overclaim Earth-size or habitable-zone occurrence constraints before completeness is demonstrated.

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
- Treat the MIT fork as an orbit/camera/CCD production pipeline, not a `quick_lc.py` replacement.
- Assume:
  - TICA FFIs are pre-staged on disk
  - `pyticdb` access to `tic_82` and `gaia3`
  - MIT `tglc-data` directory layout
- Use the MIT CLI stages: `catalogs`, `cutouts`, `epsfs`, `lightcurves`, or `all`.
- For WD work, raise the target magnitude limit beyond the bright-star defaults, typically around `--max-magnitude 20`.
- For the current CPU-only benchmark, use `--nprocs 16` as the default one-CCD `epsfs` worker count unless a new timing test motivates a change.
- Prefer `pdogpu1` for the current CPU-only Stage 1 benchmark work; reserve `pdogpu6` for follow-up only if a working `cupy`/GPU TGLC environment is prepared there.
- Preserve the MIT raw data layout. Add TWIRL metadata and indices on top rather than inventing a parallel raw layout.
- On PDO, keep TWIRL-run outputs, symlinks, staging, and scratch products under `/pdo/users/tehan/` rather than writing into shared PDO trees.

## Data Product Expectations

- MIT TGLC outputs HDF5 light curves, not the old FITS products.
- For TWIRL Stage 1, treat `TGLC end-to-end` as complete at the HDF5 `lightcurves` stage; public-style FITS products belong to the later QLP detrend + HLSP export path.
- Downstream TWIRL code should expect decontaminated aperture photometry with `1x1`, `3x3`, and `5x5` apertures.
- Do not assume old columns such as `cal_psf_flux` or `cal_aper_flux`.
- TIC IDs are useful metadata and may still matter for the current MIT implementation, but they do not define the scientific sample.
- A major early engineering risk is TIC-oriented target emission:
  - audit how Gaia-selected targets propagate through the MIT fork
  - be prepared to recommend a Gaia-first extension if scientifically important WDs are missing

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

- After code changes, run `make test-fast` when that target exists. If it does not exist yet, say so explicitly.
- For detection changes, also run `make run-detection-sample` when that target exists.
- For catalog or index changes, validate the schema and at least one sample target.

## Outputs

- Put generated reports in `reports/`.
- Put reusable config in `configs/`.
- Record assumptions and provenance in output metadata.
- When progress is made against a milestone or deliverable, record it in the same turn if practical.
- Keep `doc/twirl_plan.md` as the forward-looking plan plus compact milestone status summaries only; do not let it accumulate command-level run logs.
- Record detailed dated execution history in `doc/twirl_progress_log.md` under the most relevant stage subsection.
- In `doc/twirl_plan.md`, keep only 1-3 milestone-level status bullets per active subsection and link to the relevant `doc/twirl_progress_log.md` section.
- Keep the prose part of each progress bullet readable; do not stuff raw long paths into the sentence.
- Keep file references inline when helpful, using short clickable Markdown link labels such as `[script]`, `[PNG]`, `[PDF]`, `[FITS]`, or `[CSV]` rather than raw long paths.
- Use plain Markdown that renders everywhere.
- Append a short status marker in parentheses to subsection headings when helpful: `(✓)` done, `(...)` in progress, `(?)` problematic, and nothing if not started.
- Treat these markers as manually tuned section status labels rather than auto-generated output.

## Current Priorities

When choosing what to implement next, prefer this order:

1. WD master catalog builder
2. orbit/camera/CCD mapper
3. PDO batch wrappers around the MIT TGLC CLI
4. consolidated HDF5-to-TWIRL index format
5. WD 1856 QA benchmark
6. Gaia-first target-support decision in the MIT fork
7. baseline periodic and dip-search pipelines
8. follow-up coordination support

## Follow-Up Planning Assumptions

- Current default planning assumption is:
  - Magellan is the primary MIT-affiliated follow-up path
  - MMT is backup or collaborator-driven
- Favor high-cadence follow-up planning for short predicted transit windows.
- If a task depends on a specific instrument or access policy, verify it before writing it into the repo as fact.
- Do not assume aperiodic events can be confirmed the same way as periodic candidates.
- Before follow-up is treated as ready, require independent re-extraction, contamination checks, and ephemeris adequacy for the planned facility.

## Code And Repo Conventions

- Prefer pipeline code, scripts, and structured metadata over notebook-only workflows.
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

## Maintenance Rule

When the project plan changes materially, update this file in the same turn if possible.
