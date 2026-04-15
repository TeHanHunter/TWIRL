# twirl_plan.md

This document is the executable software and survey plan for TWIRL.

## Current Status Snapshot

- TWIRL is a standalone repo and remains separate from `TWIRL_proposal`.
- The core survey is limited to `200 s` TESS FFIs, which means `Sector >= 56`.
- `WD 1856+534` is the benchmark target, and the first fixed benchmark is `Sector 56`, orbits `119` and `120`, `cam4/ccd1`.
- The seed WD catalog is a local external dependency, not a git-tracked repo asset.
- Light curves will be gathered for all WDs in the full parent catalog. The statistical occurrence-rate sample (the locked denominator) will be defined separately after light-curve collection and QA are complete.
- The search is interpretable-first: build a transparent periodic + dip baseline before considering ML triage.
- Stage 1 Sector 56 benchmark (orbits 119 + 120, all 16 CCDs): TGLC lightcurves, detrend, and HLSP FITS production have all run end-to-end. 19,086 detrended HDF5 lightcurves per orbit, **14,807 HLSP FITS** written to `/pdo/users/tehan/tglc-deep-catalogs/hlsp_s0056/` (78% of targets). The remaining 22% are entirely in `cam3` of orbit 120: every LC fails detrend with a 6137-vs-6136 shape mismatch between TGLC's BJD array and the QLP quaternion feature matrix — a genuine cadence-length disagreement that blocks `quatspline.fit_trend`. Next session: reconcile the cam3/o120 quaternion vs BJD cadence grid, then re-run detrend + HLSP for those CCDs. Tooling fix uncovered this session: BLAS/OpenMP thread caps (single-threaded per worker) turned detrend from ~16 h to ~5 min at `-n 32`; the fix is now baked into both wrappers.

## Project Goal

TWIRL will use 200 s TESS FFIs, TGLC-based light-curve extraction, transparent transit/dip searches, optional machine-learning triage, automated vetting, and injection-recovery tests to:

1. search for any transiting or occulting objects around white dwarfs,
2. measure the completeness of that search in the regimes where 200 s data are strongest,
3. derive occurrence-rate constraints or upper limits for those regimes, and
4. produce validated candidate lists and follow-up targets.

The actual execution sequence:

1. generate WD light curves with the MIT-adapted TGLC on MIT PDO machines (and via QLP for newer sectors — see Stage 1.8),
2. build a transparent baseline search for periodic and non-periodic transit-like events,
3. run injection-recovery tests to establish completeness,
4. search the full WD sample,
5. validate resulting candidates,
6. add ML triage only if it materially improves the baseline search ranking or false-positive rate.

## Design Drivers From The Proposal

- Primary input sample: the local external Gentile Fusillo et al. (2021) Gaia EDR3 white dwarf catalogue, recommended at `data_local/catalogs/GaiaEDR3_WD_main.fits`, with Gaia DR3 as the authoritative target identifier and TIC/TESS metadata added when available.
- Photometry source: 200 s TESS FFIs only, which means sectors/orbits from Sector 56 onward.
- Extraction engine: the MIT-adapted TGLC pipeline already used in the QLP environment.
- Search target: any transiting or occulting object around a WD, with first-year emphasis on large, deep, short-duration events in the WD 1856-like regime.
- Science outputs: validated candidates, a ranked follow-up list, and occurrence-rate limits first for large occulters; Earth-size/HZ occurrence work is a longer-term goal only after completeness is demonstrated.

## Recommended Repo Layout

The repo should be organized around the pipeline stages rather than papers or notebooks.

```text
doc/
  twirl_plan.md
  twirl_progress_log.md
  mit_tglc_usage_guide.md
  ideas.md
  local_data.md
  plotting_style.md
catalogs/
  wd_master_catalog/
  sector_orbit_maps/
configs/
  tglc/
  detection/
  injections/
scripts/
  stage1_lcs/
  stage2_detection/
  stage3_injections/
  stage4_search/
  stage5_validation/
src/
  twirl/
    catalogs/
    io/
    lightcurves/
    detection/
    injections/
    search/
    validation/
notebooks/
reports/
```

## External Data Policy

Large external inputs and staged survey products should remain local and outside normal git tracking.

Recommended local layout:

```text
data_local/
  catalogs/
    GaiaEDR3_WD_main.fits
    twirl_master_catalog/
  tglc-data/
```

Rules:

- do not commit large raw FITS catalogs to this repo
- do not treat local external data as part of the git-managed source tree
- record local input path, file size, and provenance in derived metadata products
- version-control the code, configs, schemas, and derived table definitions that act on those local files

## Stage 1: Generate WD Light Curves On MIT PDO Machines

This is the foundation for everything else. The MIT fork is no longer a single-target `quick_lc.py` workflow. It is an orbit/camera/CCD production pipeline with separate catalog, cutout, ePSF, and light-curve stages.

### 1.1 Build the WD target table (✓)

The current seed catalogue for this repo is:

- `data_local/catalogs/GaiaEDR3_WD_main.fits` or another user-configured local path
- based on Gentile Fusillo et al. (2021), "A catalogue of white dwarfs in Gaia EDR3"
- this is the main Gaia EDR3 catalogue, not the reduced-proper-motion extension
- this file is a local external dependency and should not be committed to git

Important catalogue facts to treat as part of the survey definition:

- the main catalogue contains about `1.28 million` sources that passed the authors' quality selection
- the paper identifies `Pwd > 0.75` as the high-confidence WD threshold, yielding about `359,073` candidates
- the paper quotes an upper-limit completeness of about `93%` for white dwarfs with `G <= 20`, `Teff > 7000 K`, and `|b| > 20 deg`

TWIRL should build a single master catalog from this file with, at minimum:

- Gaia DR3 designation/source ID
- TIC ID, if available
- RA, Dec
- Gaia `G`, `BP`, `RP`
- Gaia `Pwd`
- TESS magnitude estimate
- WD selection probability / class flag from the input WD catalog
- sector/orbit coverage in TESS
- camera, CCD, and cutout mapping once known

This table becomes the control plane for the whole survey.

The master catalog should also include provenance fields that make later survey definitions reproducible:

- local seed-catalog path used at build time
- seed-catalog file size and modification time, or a stronger file hash when practical
- TWIRL catalog build version
- crossmatch recipe version
- sample-selection recipe version

The parent-sample protocol should be explicit in code and metadata:

- `seed_catalog`: all rows in the local Gentile Fusillo et al. (2021) main catalogue
- `reference_highconf_sample`: rows with `Pwd > 0.75`, used as a default high-confidence comparison sample for QA and pilot work
- `lc_gather_sample`: all WDs in the full parent catalog — light curves are extracted for all of them; no sample cut is applied at this stage
- `twirl_parent_sample`: the final locked denominator for occurrence-rate work, to be defined after light-curve collection and QA are complete; this cut must be frozen before any injection-recovery or occurrence-rate analysis begins
- until `twirl_parent_sample` is frozen, any completeness or occurrence summaries must be labeled provisional
- any exploratory search outside the final locked denominator must be labeled exploratory and excluded from occurrence-rate denominators by default

#### Status

- First Gaia-first WD master catalog and conservative TIC merge path are in place.
- Current milestone facts: `1,280,266` seed rows, `359,073` rows with `Pwd > 0.75`, `1,236,467` unique TIC matches, `43,035` ambiguous TIC matches, and `764` no-bridge cases.
- Detailed dated notes live in the Stage 1.1 [log](twirl_progress_log.md#11-build-the-wd-target-table).

### 1.2 Decide the production sample (✓)

Use two nested samples:

- `search_sample_primary`: the highest-priority WD sample for the first end-to-end run
- `search_sample_full`: the broader sample used for the final occurrence-rate analysis

A practical starting point is:

- use `Pwd > 0.75` as the default high-confidence reference sample for pilot QA and crossmatch studies
- primary sample near `T <= 18`
- full sample extended toward `T <= 20` where recovery remains useful
- include WD 1856+534 as the mandatory smoke-test target
- choose one `Sector >= 56` containing WD 1856+534 as the first end-to-end benchmark; the exact sector can be selected later

This creates two linked but distinct products:

- a statistical survey sample with a locked denominator and documented cuts
- an exploratory event-search layer that can be broader, but does not feed occurrence-rate claims by default

The exact TWIRL WD cut should be frozen only after:

- the Gaia-first target definition and any required MIT TGLC implementation changes are characterized well enough to support production
- benchmark QA is complete
- the baseline search stack is defined well enough to support end-to-end injections

#### Status

- The Stage 1 benchmark and pilot production-sample policy are set: use the Gaia-first parent catalog, keep `Pwd > 0.75` as the default high-confidence reference sample, and preserve unresolved final denominator choices for later occurrence-rate stages.
- Current milestone facts: `1,231,702 / 1,280,266` WDs (`96.2%`) have at least one `Sector >= 56` footprint, and `1,157,242` (`90.4%`) have observed coverage before Sector `100`.
- Detailed dated notes live in the Stage 1.2 [log](twirl_progress_log.md#12-decide-the-production-sample).

### 1.3 Prepare the MIT PDO environment (✓)

The MIT-adapted TGLC expects:

- local TICA FFI files already staged on disk,
- local `pyticdb` access to `tic_82` and `gaia3`,
- a `tglc-data` directory tree organized by orbit/camera/CCD.

For TWIRL, the immediate operational tasks are:

1. freeze the orbit list to process,
2. identify which orbit/camera/CCD combinations contain WD targets,
3. pre-stage the required 200 s FFIs,
4. confirm the `pyticdb` databases resolve on PDO,
5. define batch-job wrappers for one orbit-camera-CCD unit per job.

#### Status

- PDO helper-script environment and `pyticdb` access to `tic_82` are validated for the Gaia-to-TIC export path.
- Detailed dated notes live in the Stage 1.3 [log](twirl_progress_log.md#13-prepare-the-mit-pdo-environment).

### 1.4 Generate catalogs and cutouts (...)

The MIT fork runs in four logical stages:

1. `catalogs`
2. `cutouts`
3. `epsfs`
4. `lightcurves`

For TWIRL, the key configuration change is to lift the TIC magnitude limit well beyond the default value. The MIT fork defaults to bright-star QLP use cases; WD work needs a much deeper target cut.

Initial production assumptions:

- run `--max-magnitude 20`
- keep Gaia catalogs for all relevant field stars
- run per orbit and per CCD, not per target

#### Status

- The fixed Sector `56`, orbit `119/120`, `cam4/ccd1` benchmark has been pushed through `catalogs`, `cutouts`, `epsfs`, and WD-only `lightcurves`; WD 1856 is present in both orbit trees as `267574918.h5`.
- The pre-Sector-67 TICA header fallback is patched in the user-owned TGLC fork, and the current default per-CCD `epsfs` worker count is `--nprocs 16` based on the first one-CCD benchmark.
- The first full-orbit orbit `119` benchmark on `pdogpu1` has completed through HDF5 `lightcurves`: the 15-CCD sweep outside the original `cam4/ccd1` benchmark finished in `33.66 h` wall time with outer CCD concurrency `3` and inner workers `16/16/16/16`, and the final orbit tree now contains `32` catalog files, `3136` cutouts, `3136` ePSFs, and `19086` WD-only `.h5` light curves.
- The clean stopwatch for current performance planning is therefore `33.66 h` for the 15 non-benchmark CCDs; the orbit is now complete across all 16 CCDs, but the original `cam4/ccd1` benchmark was produced earlier in separate runs with different tuning, so it should be treated as a separate seed benchmark rather than folded into one apples-to-apples full-sector wall-time claim.
- Detailed dated notes live in the Stage 1.4 [log](twirl_progress_log.md#14-generate-catalogs-and-cutouts).

### 1.5 Extract and consolidate WD light curves

The full Stage 1 pipeline per orbit/camera/CCD is:

```
tglc catalogs → tglc cutouts → tglc epsfs → tglc lightcurves
  → qlp lctools detrend → qlp lctools hlsp (→ FITS per TIC)
```

`tglc lightcurves` writes one HDF5 file per TIC target. The `detrend` step (run via `qlp lctools detrend`) annotates those HDF5 files with a `bestdmagkey` attribute that `hlsp.py` reads to select the detrended magnitude column. **Without detrending, HLSP generation fails.** The `hlsp` step then converts HDF5 → FITS, requiring the full QLP environment (`lightcurvedb` + `qlp` package, not just TGLC).

For Sector 56 specifically, HLSP generation uses SPOC quality flags (not TICA), because `get_ticaflags` in `hlsp.py` raises `ValueError` for `sector < 67`. These SPOC flags live at `/pdo/qlp-data/spocflags/`.

**Standard recovery rate check** — run after both HDF5 and FITS production:

For each orbit/camera/CCD:
1. Count TIC IDs requested (from TWIRL detector target table)
2. Count HDF5 files produced (`tglc lightcurves` output)
3. Count FITS files produced (`qlp lctools hlsp` output; missing TICs logged as warnings by `generate_qlp_hlsp_fits_file`)
4. Bin the shortfall by TESS magnitude to detect faint-end dropout
5. Cross-check against Gaia-only targets (no TIC bridge): these will always be 0% recovered until the MIT fork is extended to emit Gaia-selected targets directly

This 3-level audit (requested → h5 → FITS) is the standard health check for every production run. Log counts and dropout fractions in the progress log.

After production, build a TWIRL light-curve index with:

- TIC ID and Gaia DR3 ID
- orbit, sector, camera/CCD
- HDF5 path and FITS path
- production status flag (`h5_ok`, `fits_ok`)
- summary quality metrics (RMS, MAD, cadence count)

The index format should explicitly distinguish raw catalog columns, derived TWIRL metadata, crossmatch outputs, processing-status fields, and sample-membership flags.

### 1.6 Run photometric QA before any ML work

Before training or search, measure:

- per-target RMS / MAD versus magnitude
- sector-level failure rates
- missing cadence statistics
- light-curve completeness by target and by orbit
- recovery of WD 1856+534 b in the 200 s products
- agreement or disagreement across the 1x1, 3x3, and 5x5 apertures for benchmark targets
- basic comparison against at least one independent extraction path on the benchmark target set

This step should produce a small QA report and a list of sectors/orbits to reprocess if needed.

### 1.7 Benchmark acceptance criteria

Before scaling beyond the pilot stage, the WD 1856 benchmark should satisfy all of the following at the documentation level:

- the expected transit window is visible in at least one `Sector >= 56` benchmark product
- the event timing is consistent with the published WD 1856 ephemeris within the expected propagated uncertainty
- the signal is present in at least one aperture product and checked across the `1x1`, `3x3`, and `5x5` apertures
- any strong aperture disagreement is documented and investigated as a contamination or extraction issue
- an independent extraction path has been compared on the benchmark set
- the benchmark result is strong enough to serve as a smoke test for the first baseline search implementation

### 1.8 Gaps to close in the current MIT fork (?)

The MIT fork is close to what TWIRL needs, but not identical to the survey requirements.

Known gaps to track early:

- Gaia DR3 should remain the survey-defining target identifier even if the current MIT implementation is still partly TIC-oriented.
- Gaia-to-TIC matching should be audited for metadata and implementation reasons, but it should not define the scientific parent sample by itself.
- If important WD targets lack TIC IDs, the light-curve stage will need to be extended to support Gaia-selected targets directly.
- For current terminology, distinguish `TGLC end-to-end` from `QLP/MAST end-to-end`: Stage 1 extraction is complete at the HDF5 `lightcurves` output, while public-style FITS deliverables require the later QLP detrend + HLSP export path.
- After the full Sector 56 benchmark run, verify that all metadata needed to launch whole-sector TIC light-curve extraction with a magnitude cut is already present and discoverable from the generated products, not hidden in ad hoc scripts.
- After the full Sector 56 benchmark run, confirm that the `cutouts` and `epsf` products under `/pdo/users/tehan/tglc-deep-catalogs/orbit-*/ffi/cam*/ccd*/` are exactly in the MKI TGLC `Manifest` layout so they can be reused by standard `tglc lightcurves` commands without special-case path handling.
- The new output product is decontaminated aperture photometry in HDF5, not the old per-target FITS product with `cal_psf_flux` and `cal_aper_flux`.
- The old `prior`-based single-target workflow is not exposed in the new CLI, so any need for floating-field-star priors must be added explicitly.
- The catalogue completeness quoted by Gentile Fusillo et al. (2021) applies only in specific regimes; TWIRL must not generalize that completeness to the full search sample without its own end-to-end validation.

Questions to answer before mass production:

- What outer CCD concurrency should be used per stage once the first full-orbit benchmark finishes, especially for the DB-heavy `catalogs` stage versus the CPU-heavy `epsfs` stage?
- Why do the WD-only `lightcurves` runs currently write fewer `.h5` files than requested TIC IDs, and is that shortfall astrophysical, metadata-driven, or an extraction/control-plane issue?
- What is the official QLP/LCDB metadata path needed to run detrend and HLSP export without the current local benchmark shims?
- Once the full Sector `56` benchmark is done, should whole-sector production stay WD-only at `lightcurves`, or should TWIRL also preserve a full-TIC sector archive for direct QLP/MAST-style downstream products?
- **QLP ingestion for Sectors 94+**: TGLC has been incorporated into QLP (Petitpas et al. 2025, in prep.) and is producing light curves from Sector 94 onward. Should TWIRL ingest those QLP/MAST products directly for recent sectors rather than re-running MIT PDO extraction? Key questions: format compatibility with the HDF5 archive, whether the WD magnitude limit (`--max-magnitude 20`) was applied, and whether the Gaia-first target list is recoverable from QLP outputs. Resolving this could substantially reduce the PDO compute load for newer sectors.

### Stage 1 Deliverables

- WD master target catalog (✓)
  Current status: v0 catalog built and full Gaia DR3 to TIC export completed.
- orbit/camera/CCD job table
- automated PDO production scripts
- consolidated WD light-curve archive
- QA report and reprocessing list

## Stage 2: Build The Search Stack And Optional WD-Specific Classifier

The proposal points toward Astronet-Triage-like deep-learning search, but the first survey version should not be ML-first. TWIRL should begin with an interpretable search stack and add ML only where it clearly improves ranking or false-positive control.

### 2.1 Define the detection unit

Decide what the search stack scores:

- sector-level candidate windows,
- period-search peaks,
- folded light curves,
- local/global transit views derived from a separate search stage,
- and non-periodic dip clusters that may represent debris or irregular occultations.

For the first version, the cleanest path is:

1. produce a transparent periodic candidate list with a short-duration box/trapezoid search,
2. run a separate dip-search branch for non-periodic or weakly periodic events,
3. pass candidates into automated vetting with LEO-vetter (Kunimoto et al. 2025) — assess input format compatibility with the TWIRL HDF5 archive before committing to it,
4. add a WD-specific classifier later if it materially improves triage.

### 2.1b Minimum baseline search specification

Before adding ML triage, the baseline search should define and document:

- the detrending approach used for minute-scale events, with explicit injection checks that the detrending does not erase the target signals
- how `1x1`, `3x3`, and `5x5` apertures are searched and compared
- the rule for choosing a preferred aperture per candidate
- the candidate statistics recorded by the periodic branch
- the candidate statistics recorded by the non-periodic dip branch
- the multi-sector merge logic for periodic candidates
- the rule for keeping irregular or apparently aperiodic events separate from the periodic survey sample

### 2.2 Build the labeled dataset

Positive examples:

- injected WD transits in real TGLC light curves
- WD 1856+534 b as a benchmark
- synthetic variants covering duration, depth, cadence loss, and crowding
- large-occulting benchmark cases spanning giant planets, brown dwarfs, and compact stellar companions

Negative examples:

- quiet WDs
- eclipsing binaries and blends
- scattered-light / systematics artifacts
- cadence-gap and momentum-dump failures

### 2.3 Train the model

If a classifier is added, the training pipeline should include:

- train/validation/test splits by target, not random cadence
- class-balance control
- calibration of classifier scores
- explicit evaluation versus magnitude, period, depth, and number of sectors

### Stage 2 Deliverables

- reproducible periodic-search baseline
- reproducible dip-search baseline
- reproducible training set definition
- trained detector checkpoint(s), if justified
- evaluation report
- inference script that scores the full TWIRL archive

## Stage 3: Injection-Recovery Tests

This is the completeness backbone of the survey and should run through the same search stack used in Stage 4.

### 3.1 Injection design

Inject over a grid in:

- orbital period
- transit depth
- duration / impact parameter
- host magnitude
- number of observed sectors

The grid should be densest near:

- the WD 1856-like large-occulting regime,
- Roche-limit boundary cases,
- deep short-duration transits that are most astrophysically plausible for WD systems,
- and only secondarily the Earth-size/HZ regime, which should not drive the first-year completeness claims.

### 3.2 Recovery protocol

Each injection should pass through:

1. the search stage,
2. the classifier, if used,
3. the vetter,
4. candidate merging rules.

Do not estimate completeness from only one stage; use the true end-to-end recovery fraction.

### 3.3 Completeness products

Produce completeness surfaces as a function of:

- `Tmag`
- period
- depth
- number of sectors
- crowding / contamination metrics

Completeness products should be tagged with the exact sample definition and search-branch definition they apply to. Do not present a completeness surface as survey-wide if it only applies to a provisional cut or one search mode.

### Stage 3 Deliverables

- injection engine
- recovery summaries
- completeness grids for occurrence-rate work

## Stage 4: Search For WD Transiting Objects

After the light curves, detector, and completeness machinery are stable, run the production search.

### 4.1 Full survey inference

Run the periodic and non-periodic search branches across all WD light curves and record:

- candidate ephemerides
- detection scores
- per-sector evidence
- merged multi-sector candidates
- flags for apparently aperiodic or irregular events

### 4.2 Candidate catalog construction

For each candidate, store:

- target identifiers
- discovery sectors
- transit parameters
- detection statistics
- vetter outputs
- links to diagnostic plots and files

### Stage 4 Deliverables

- machine-generated candidate table
- diagnostic plot package
- ranked follow-up target list

## Stage 5: Validate The Targets

Validation should combine automated filtering and external follow-up.

### 5.1 Automated checks

- contamination and crowding review
- centroid and aperture behavior
- odd/even and depth consistency checks where applicable
- ephemeris consistency across sectors
- image-level inspection around predicted transits

### 5.2 Archival time-domain first pass (ZTF + ATLAS)

Before requesting any new telescope time, cross-match surviving candidates against the ZTF and ATLAS public archives. Both surveys provide multi-year optical light curves for most TWIRL targets with `g/r/i` or `c/o` coverage, which is sufficient to:

- rule out obvious stellar variables misclassified as transit candidates (pulsators, eclipsing binaries, cataclysmic variables),
- check for additional dips or flares at other epochs that argue for or against a transit-like interpretation,
- extend the baseline for period refinement on periodic candidates using a single combined WD+ZTF+ATLAS ephemeris fit.

This is a cheap, fully automated pre-filter; candidates that do not survive ZTF/ATLAS consistency do not advance to Stage 5.3 telescope follow-up.

### 5.3 External follow-up

For candidates that survive the archival first pass:

- high-cadence photometry at predicted transit windows
- archival imaging and catalog checks (beyond ZTF/ATLAS)
- spectroscopy or RV constraints where physically meaningful
- high-resolution imaging if blending is a concern

Current planning assumption for MIT-affiliated follow-up:

- treat Magellan as the primary institutional follow-up path to develop early
- treat MMT as a collaborator-driven or backup path rather than the main plan
- identify which MIT-access high-speed photometric instrument is actually schedulable for WD transit windows
- build a rapid-response workflow around short predicted windows and sub-minute cadence requirements

Follow-up readiness criteria should be satisfied before requesting telescope time for a periodic candidate:

- the event survives the baseline search and automated vetting
- the signal is independently re-extracted or cross-checked in a second reduction path
- image-level and aperture-level contamination checks are complete
- the ephemeris is precise enough for the planned instrument and observing window
- the expected cadence and depth are compatible with the planned facility
- the candidate status, assumptions, and main false-positive modes are documented

For apparently aperiodic or irregular events:

- do not assume repeat-photometry confirmation is possible
- define a separate follow-up path centered on classification, continued monitoring, or archival/context constraints

### 5.4 Final population analysis

Regardless of the number of validated planets, TWIRL should finish with:

- occurrence-rate posteriors or upper limits, first in the large-occulting regime where the 200 s survey is strongest
- candidate catalog publication
- validated discoveries, if present

A publishable null result still requires:

- a locked parent sample
- a documented search definition
- end-to-end completeness tied to that search definition
- a final vetted candidate table or null-candidate statement
- upper limits stated only for the declared sample and event class

### Stage 5 Deliverables

- validated-candidate list
- follow-up status table
- occurrence-rate paper inputs

## Suggested Milestones

### Year 1

- finalize WD sample
- stand up MIT PDO light-curve production
- select and run a first `Sector >= 56` benchmark containing WD 1856+534
- recover WD 1856+534 b with the new 200 s workflow
- stand up the baseline periodic and dip-search branches
- create first injection datasets for large occulters
- identify the concrete MIT-affiliated follow-up path and required collaborators/instruments

### Year 2

- run the full search stack and automated vetting
- publish candidate catalog
- execute first follow-up campaign

### Year 3

- complete validation
- compute occurrence rates with completeness corrections, starting with the large-occulting regime
- publish discoveries and/or upper limits

## Immediate Implementation Priorities For This Repo

1. Create the WD master catalog builder.
2. Create the orbit/camera/CCD mapper for the WD sample.
3. Wrap the MIT TGLC CLI in PDO batch scripts.
4. Define the consolidated HDF5-to-TWIRL index format.
5. Build a QA notebook/report around WD 1856+534 b and a small control sample.
6. Audit Gaia-first target support and decide what MIT fork changes are needed for targets without TIC IDs.
7. Build the first transparent periodic and dip-search baselines before committing to an ML-heavy workflow.
8. Lock down the MIT-affiliated follow-up path for short, high-cadence transit confirmation.
