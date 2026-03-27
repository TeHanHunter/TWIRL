# plan.md

This file is the short operational summary of the current TWIRL repo state.

For the full pipeline specification, read `docs/twirl_pipeline_plan.md`.
For agent protocol, read `AGENTS.md`.
For brainstorming, referee risks, and unresolved questions, read `ideas.md`.

## Current Repo Role

- `plan.md`: current status, completed decisions, and next tasks
- `AGENTS.md`: execution protocol for future agent work
- `ideas.md`: open questions, referee-style concerns, and brainstorming
- `docs/twirl_pipeline_plan.md`: detailed scientific and pipeline plan
- `docs/mit_tglc_usage_guide.md`: MIT TGLC operating notes

## Completed Decisions

- TWIRL is a standalone repo and remains separate from `TWIRL_proposal`.
- The core survey uses `200 s` TESS FFIs only.
- That means the production survey should focus on `Sector >= 56`.
- The first benchmark target is `WD 1856+534`.
- The exact benchmark sector is intentionally not fixed yet.
- The first search should be interpretable:
  - periodic short-duration baseline search
  - separate dip-search branch for non-periodic or weakly periodic events
- ML is not the first discovery engine.
- The seed WD catalogue is the local Gentile Fusillo et al. (2021) Gaia EDR3 main catalogue.
- `Pwd > 0.75` is the default high-confidence reference sample for pilot work.
- The exact final TWIRL WD denominator is still to be determined later.

## Current Data Policy

- Large external files stay local under `data_local/` or another user-configured path.
- The Gaia EDR3 WD FITS file should live at `data_local/catalogs/GaiaEDR3_WD_main.fits` by default.
- Raw external data and staged MIT TGLC products should not be committed to git.

## Current Scientific Position

- The near-term survey should optimize for large, deep, short-duration events in the WD 1856-like regime.
- Earth-size or habitable-zone occurrence statements are a longer-term goal and should not be overclaimed early.
- Statistical survey products and exploratory searches should remain distinct.
- Any occurrence-rate claim requires a locked denominator and end-to-end completeness.

## Completed Documentation Work

- The detailed pipeline plan was updated to reflect:
  - `200 s`-only survey scope
  - local-data policy
  - deferred final WD denominator
  - benchmark acceptance criteria
  - follow-up readiness criteria
  - null-result requirements
- `AGENTS.md` was created and aligned with the current science and implementation rules.
- `data_local/README.md` was added to document local-only data handling.

## Current Priorities

1. WD master catalog builder
2. orbit/camera/CCD mapper
3. PDO batch wrapper for MIT TGLC
4. consolidated HDF5-to-TWIRL index
5. WD 1856 QA benchmark
6. Gaia-to-TIC completeness characterization
7. baseline periodic and dip-search implementation
8. follow-up coordination support

## Near-Term Deliverables

- a reproducible TWIRL master catalog built from the local Gaia EDR3 WD seed file
- a documented crossmatch and indexing scheme
- a pilot `Sector >= 56` benchmark run containing WD 1856
- QA outputs sufficient to decide whether full production should proceed
- a transparent baseline search spec that is ready for end-to-end injections

## Known Unresolved Items

- the exact TWIRL WD denominator
- the final statistical event class for occurrence-rate work
- the exact benchmark sector for WD 1856
- the final Gaia-to-TIC handling policy for missing or ambiguous matches
- the precise MIT follow-up instrument path and scheduling route
