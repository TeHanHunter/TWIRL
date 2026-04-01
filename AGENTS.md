# AGENTS.md

This file is the operating protocol for agents working in the standalone `TWIRL` repository.

## Read First

Before making substantive changes, read:

- `doc/twirl_plan.md`
- `doc/mit_tglc_usage_guide.md`
- `doc/ideas.md`
- `doc/local_data.md`
- `doc/plotting_style.md` before making or revising publication-facing figures

Use `doc/twirl_plan.md` as the single project plan, `doc/ideas.md` for unresolved questions and brainstorming, `doc/local_data.md` for local-only data conventions, and `doc/plotting_style.md` for figure styling.

Follow the Step 1-5 pipeline structure defined in `doc/twirl_plan.md`.

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
- Use canonical versioned master-catalog releases such as `twirl_wd_master_catalog_v1.fits` for accepted states; treat suffix-heavy stage filenames such as `*_ticmatched.fits` or `*_tesscoverage.fits` as temporary intermediates.

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
- The exact benchmark `Sector >= 56` containing WD 1856 is intentionally not fixed yet.
- Do not silently choose or hardcode the benchmark sector unless the user asks.
- The first search version should be interpretable:
  - periodic short-duration box or trapezoid search
  - separate dip-search branch for non-periodic or weakly periodic events
- Machine learning is optional triage, not the first discovery engine.

## MIT TGLC Operating Assumptions

- Step 1 uses the MIT-adapted TGLC fork on MIT PDO machines.
- Treat the MIT fork as an orbit/camera/CCD production pipeline, not a `quick_lc.py` replacement.
- Assume:
  - TICA FFIs are pre-staged on disk
  - `pyticdb` access to `tic_82` and `gaia3`
  - MIT `tglc-data` directory layout
- Use the MIT CLI stages: `catalogs`, `cutouts`, `epsfs`, `lightcurves`, or `all`.
- For WD work, raise the target magnitude limit beyond the bright-star defaults, typically around `--max-magnitude 20`.
- Preserve the MIT raw data layout. Add TWIRL metadata and indices on top rather than inventing a parallel raw layout.

## Data Product Expectations

- MIT TGLC outputs HDF5 light curves, not the old FITS products.
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
- When progress is made against a milestone or deliverable in `doc/twirl_plan.md`, record that progress in the same turn if practical.
- In `doc/twirl_plan.md`, record progress as concise dated bullets under the most relevant subsection rather than in a separate free-floating log when possible.
- Include the relevant code or output path(s) in each progress bullet so the artifact is easy to inspect.
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
- For long-running scripts or batch jobs, add periodic progress reporting when practical so PDO runs are observable without attaching a debugger or guessing from zero-byte output files.

## When To Stop And Ask

Ask the user before:

- selecting a benchmark sector on their behalf
- broadening the science claims beyond the current plan
- merging proposal-repo content into this repo
- committing to a follow-up instrument path that has not been verified
- changing the repo scope away from the standalone pipeline project

## Maintenance Rule

When the project plan changes materially, update this file in the same turn if possible.
