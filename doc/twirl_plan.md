# TWIRL survey and software plan

This is the single forward-looking plan for the standalone TWIRL pipeline.
Detailed execution history belongs in [the progress log](twirl_progress_log.md),
unresolved questions belong in [ideas](ideas.md), and operational commands
belong in the relevant runbook. Report-level plans and status files are dated
evidence, not project authority.

Last reconciled: `2026-07-20`.

## Current status

- The survey uses `200 s` TESS FFIs from `Sector >= 56`, with Gaia DR3
  `source_id` as the scientific identifier and TIC as operational metadata.
- `A2v1` is the active Stage 1 product family: no TIC-magnitude science cap,
  saturated-pixel ePSF masking, and sector-level ADP/ADP015 FITS products for
  the `1x1`, `3x3`, and `5x5` apertures. S56 (`31,450` FITS), S57
  (`27,213` FITS), and S58-S63 pass their edge-aware HDF5/FITS product gates.
  The next source-only production batch is S64-S69; it will refit all ePSFs
  because no legacy ePSFs are available. See
  [Stage 1 history](twirl_progress_log.md#stage-1).
- S56 passes the current **Tier-0 integrity/benchmark QA**, including WD 1856
  recovery in both active ADP apertures. Tier 0 verifies product integrity and
  benchmark behavior; it is not the Tier-1 science QA needed for survey
  release. The `s56_harmonic_cnn_v1` teacher remains the active-learning
  baseline. Teacher v2 was completed as an exploratory comparison but missed
  its promotion gates, so it is not a production ranker or student-label
  generator. See [Stage 2 history](twirl_progress_log.md#stage-2).
- LC-level injection/recovery, two-aperture vetting, and pixel-injection smokes
  exist, but there is not yet a frozen survey-wide completeness product. The
  transparent non-periodic dip branch, multi-sector search, stable archive
  index, and representative pixel-level calibration remain open.
- The active framework manuscript is *TWIRL I: A Systematic TESS Search for
  Transiting Planetary Remnants around White Dwarfs* in the separate sibling
  repo `/Users/tehan/PycharmProjects/twirl-survey-paper`. This repo remains the
  executable survey pipeline and must not absorb the paper or proposal repo.

## Locked survey decisions

- Seed catalog: Gentile Fusillo et al. (2021) Gaia EDR3 main catalog, kept as
  a local external input.
- High-confidence reference sample: `Pwd > 0.75`.
- Light-curve gathering: all reachable Gaia-first WD targets; magnitude cuts
  may define QA slices or a later statistical denominator but must not gate
  TGLC production.
- Statistical denominator: not yet frozen. The final `twirl_parent_sample`
  must be fixed before publication-grade end-to-end completeness or occurrence
  inference.
- Mandatory benchmark: WD 1856+534 in S56, orbits `119` and `120`,
  `cam4/ccd1`. Do not substitute another benchmark without an explicit user
  decision.
- Discovery engine: transparent periodic and non-periodic searches first.
  Machine learning may rank or triage candidates only after its inputs,
  provenance, and held-out behavior are auditable.
- Compute boundary: primary TGLC/ePSF production stays on PDO. ORCD consumes
  compact exports for downstream search, injections, features, and models.
  Routine ORCD work uses CPU nodes or one H200; more than two H200s requires
  explicit approval for the specific run.

## Product and interface contracts

### Master catalog

- One row per Gaia DR3 source, preserving seed-catalog columns and explicit
  sample-membership flags.
- Observation table is one row per target/sector/orbit/camera/CCD hit.
- All releases record seed provenance, selection recipe, crossmatch version,
  and build version. Gaia-only/no-TIC targets remain visible rather than being
  silently removed.
- Before the TWIRL I survey run, freeze a machine-readable release manifest
  that names the last included sector, accepted sector products and checksums,
  parent-sample/catalog versions, QA state, and exact search/injection
  contracts. `Sector >= 56` is the eligibility floor, not an indefinitely
  expanding publication release.

### A2v1 light curves

The operational contract is defined in the
[A2v1 protocol](a2v1_production_protocol.md):

- requested TICs come from observation-table `source_tic` overlays;
- prepared source pickles are reused without rewriting raw pixel payloads;
- ePSFs are reused only for empty saturated-pixel masks and otherwise refit;
- HDF5 extraction is an intermediate checkpoint;
- accepted sectors require ADP/ADP015-only FITS plus edge-aware full validation;
- compact exports and derived indices retain Gaia/TIC IDs, sector/orbit/
  detector provenance, cadence/quality information, product/schema version,
  and source checksums.

### Candidate and label provenance

- Search tables record branch/config version, aperture, all retained peaks,
  cadence cleaning, ephemeris, and cross-aperture evidence.
- Keep injection truth, BLS recovery, automated-vetter output, human
  morphology, reviewed fold factor, pseudo-labels, and model scores as
  separate fields. None may silently overwrite another.
- `morphology_fold_factor` is an evidence-view choice, not automatically a
  physical orbital-period claim.

## Stage 1: catalog, extraction, products, and QA (...)

### Implemented

- Gaia-first master catalog, TIC bridge, TESS coverage products, detector job
  tables, reusable TGLC wrappers, A2v1 source overlays, masked ePSF path,
  ADP/ADP015 FITS writer, compact export, and edge-aware validator.
- S56 and S57 are complete through required FITS validation. S56 also passes
  Tier-0 integrity/benchmark QA; Tier-1 science QA remains open.

### Current gate

S58-S63 completed their gated HDF5/FITS production on `pdogpu5`. The S64-S69
source-only batch is active on `pdogpu5` in the generic queue's explicit
all-refit mode; partial legacy ePSF inputs remain a hard failure. Do not create
sector-specific production logic unless the sector is a documented exception.

### Exit criteria

- Every production sector has validated HDF5 and FITS coverage, with only
  explicit `edge_warn` omissions.
- A stable light-curve index joins Gaia/TIC identity, paths, product/schema
  versions, detector visits, cadence retention, scatter metrics, and QA state.
- A frozen release manifest records the sector cutoff and accepted product,
  catalog, QA, search, and injection versions used by TWIRL I.
- Gaia-only/no-TIC support and the S94+ QLP-ingestion choice are characterized.
- Tier-1 science QA covers RMS/MAD versus magnitude with regression limits,
  missing cadences, quality flags, aperture outliers, WD 1856 timing, a fixed
  injection-preservation test, and a genuinely independent extraction
  comparison. A Tier-0 pass alone does not promote a sector for science use.

## Stage 2: transparent search and candidate generation (...)

### Implemented

- Per-sector multi-aperture BLS, retained peak tables, cross-aperture
  consolidation, heuristic vetting, LEO comparison, centroid diagnostics,
  human-review tooling, and an ADP-only S56 active-learning teacher.

### Current gate

The focused S56 compact revisit is complete (`407/407` sheets, including `11`
new Planet-like labels from the `400` model-selected compact rows), but the
separate blinded S56 `1,000`-TIC enrichment batch remains only partially
labeled (`177/1,000` at the preserved checkpoint). A bounded Franklin handoff
now deliberately extends enrichment review to `3,000` fresh real TICs across
S57--S59 (`1,000` per sector). It uses the frozen teacher-v1 ranker only to
enrich review: every displayed ephemeris is ADP-small BLS rank one, and model
scores and selection buckets remain hidden. The six prior S57 labels are
preserved as premature experimental evidence and the entire earlier S57 queue
is excluded. S57 is therefore no longer a pristine external holdout.

The immediate parallel work is to harden the existing S56 periodic/enrichment
path, not to add search branches. Require the Tier-1 target pass mask, freeze a
candidate-level aperture rule that keeps the small ADP aperture as the search
channel and primary ADP aperture as contamination evidence unless injection
and real-data tests justify a change, complete the bounded S56 review, and
freeze the resulting candidate/label set. Defer the non-periodic dip detector,
multi-sector aggregation, and false-alarm/background calibration until this
path is robust. Those capabilities remain mandatory before the full survey
search or a science-ready candidate catalog.

### Model gate

- Use `s56_harmonic_cnn_v1` as the active-learning baseline only.
- Treat the small real training set as the limiting resource: improve the
  versioned data/label manifest, TIC-grouped splits, source-separated
  evaluation, probability calibration, and bootstrap uncertainty before
  changing model architecture.
- Reach at least `50` unique real Planet-like labels and pass locked grouped
  real-data performance/calibration gates before student pseudo-labeling.
- Teacher v2 is an exploratory completed comparison that missed its external-
  retention and morphology-promotion gates. Do not promote or iterate it on
  the critical path. Any future model iteration requires the rare-factor and
  predicted-factor-versus-oracle-factor design in [ideas](ideas.md), plus a
  predeclared advantage over the v1 baseline.

### Exit criteria

- Reproducible periodic and dip baselines with versioned configs.
- WD 1856 recovery and deterministic synthetic smokes pass.
- Multi-sector candidate objects and branch-aware false-alarm statistics are
  defined before the full survey search.
- Any model used beyond enrichment has target-grouped splits, source-separated
  metrics, calibration, provenance-safe inputs, and a documented baseline
  comparison.

## Stage 3: injection-recovery and completeness (...)

### Implemented

- Raw-flux pre-detrend BATMAN injections, balanced and all-host samples,
  aperture/detrending/BLS audits, compact ORCD products, peak-recall/ranking
  diagnostics, and a working pixel/source-pickle full-chain smoke.

### Current gate

Run frozen LC-level candidate-retention recovery against the accepted A2v1 ADP
pair and the same Stage 2 search/vetting contract intended for Stage 4. Add a
representative pixel-level calibration subset spanning magnitude, crowding,
aperture disagreement, and detector-edge behavior. Reserve "end-to-end
completeness" for the later measurement that traverses extraction calibration,
search, any production ranker, automated vetting, and candidate merging.

### Exit criteria

- Recovery passes through search, any production ranker, automated vetting,
  and candidate merging.
- Completeness surfaces are reported by magnitude, period, radius/depth,
  duration, sector count, and crowding, with support boundaries shown.
- The LC-level/pixel-level delta is quantified before occurrence claims.
- Every result names the parent sample and exact search/product versions.

## Stage 4: frozen full-survey inference

Stage 4 begins only after the Stage 1 archive/index, Stage 2 periodic+dip
contracts, and Stage 3 recovery gates are frozen.

### Required implementation

- Add `scripts/stage4_search/` thin drivers over reusable package APIs.
- Run periodic and dip branches across all accepted sectors, merge multi-sector
  evidence, apply only approved rankers, and emit one versioned candidate
  catalog with diagnostic references.
- Keep broader exploratory discoveries separate from the locked statistical
  sample and occurrence denominator.

## Stage 5: validation and follow-up

### Implemented foundation

- Heuristic/LEO checks, two-aperture sheets, human-label schema and
  adjudication, scalar centroid diagnostics, and a WD 1856 pixel-map
  diagnostic.

### Required implementation

- Promote pixel-level on-target/off-target evidence to a reproducible candidate
  check and add image/crowding, odd-even, secondary-event, multi-sector, and
  independent-reduction checks.
- Use ZTF/ATLAS archival time-domain checks before new telescope requests.
- For periodic candidates, require a follow-up-ready ephemeris, contamination
  review, independent re-extraction, and facility-compatible cadence/depth.
  Irregular events use a separate monitoring/classification path.
- Treat Magellan as the current primary MIT-affiliated planning path and MMT as
  backup/collaborator-driven. Other facilities mentioned in meetings remain
  unverified until access, cadence, depth, and lead time are checked.

## Sample freeze and occurrence analysis

Freeze `twirl_parent_sample` only after:

1. Stage 1 coverage/QA and Gaia-to-TIC losses are characterized;
2. the production periodic and dip definitions are fixed;
3. the full extraction-to-candidate injection contract is fixed; and
4. exclusions and exploratory-only targets are encoded in catalog metadata.

Occurrence posteriors or upper limits must be limited to the declared sample,
event class, support region, and validated end-to-end completeness surface. A
benchmark recovery, BLS-to-teacher retention fraction, or classifier score is
never a completeness measurement.

## Collaboration and publication boundaries

- TWIRL is a Schwamb-group collaboration. Te Han owns the pipeline/data-product
  stewardship; detailed contributor ownership and paper leadership must be
  confirmed in writing rather than inferred from code contributions.
- TWIRL I is a framework/overview paper, not a final occurrence-rate or
  discovery paper. Discovery, catalog, methods, and occurrence-paper leadership
  remain separate decisions where not already explicitly agreed.
- Julien participates in comparison and follow-up planning. Instrument ideas
  from meetings are hypotheses until independently verified.
- Keep only agreed roles and open decisions in this repository. Private
  negotiation tactics and career-planning notes do not belong in the live
  public plan.

## Immediate implementation priorities

1. Run the gated S64-S69 A2v1 HDF5/FITS queue in source-only, all-ePSF-refit
   mode; retain its stop-on-failure gates and do not bypass a partial-input
   preflight failure.
2. Complete the bounded S56 `active_search_pair` Tier-1 evidence: build the
   authoritative cadence/quality reference, produce the genuinely independent
   WD 1856 comparison, rerun the current Tier-0 report, and publish the target
   QA pass mask. This scope may qualify enrichment but never science release.
3. Complete and audit the bounded S56 and Franklin S57--S59 enrichment
   reviews, freeze the candidate/aperture contract, and merge only confidently
   adjudicated labels into a versioned training set. Do not expand beyond this
   bounded `3,000`-TIC handoff until its yield and label consistency are audited.
4. Retrain and evaluate teacher v1 on that frozen set with TIC-grouped,
   source-separated, calibrated metrics and uncertainty intervals. Do not put
   teacher v2, student pseudo-labels, or another model family on this path.
5. After the periodic/enrichment path is robust, add the dip branch,
   multi-sector merging, and branch-aware false-alarm calibration; then rerun
   frozen-chain candidate-retention and representative pixel-level recovery
   before survey-wide enrichment or science claims.
6. Freeze the compact-export/index schema, release cutoff/manifest, and parent-
   sample criteria; characterize the `764` no-TIC-bridge WDs and the S94+ QLP
   boundary before the survey release is locked.

When a priority completes, record details in the progress log and retain only
one to three milestone-level status bullets here.
