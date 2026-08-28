# TWIRL survey and software plan

This is the single forward-looking plan for the standalone TWIRL pipeline. It
contains current milestones, locked decisions, stage gates, and priorities.
Read [the documentation guide](README.md) to find task-specific protocols.
Detailed results belong in [the progress log](twirl_progress_log.md), and
unresolved choices belong in [ideas](ideas.md).

Last reconciled: `2026-08-28`.

## Milestone dashboard

| Workstream | Current milestone | Next gate |
| --- | --- | --- |
| Stage 1 light curves | A2v1 S56--S65 accepted; S66 partially returned; S67--S78 retained ORCD compute checkpoints; S79 computing and S80 prepared | Re-ingress S66's two cleaned source grids, return S66, then run the unchanged PDO acceptance gates and continue chronological returns |
| Stage 2 search | Transparent periodic BLS exists; Teacher v3 is frozen for enrichment | Complete the sealed S63 test while independently implementing the dip branch and multi-sector merging |
| Stage 3 completeness | LC-level and pixel-level injection pilots exist | Freeze one extraction-to-candidate recovery chain and a representative pixel calibration subset |
| Stage 4 inference | Not started | Wait for the Stage 1 index, Stage 2 search contracts, and Stage 3 recovery gates |
| Stage 5 validation | Human review, aperture checks, LEO/centroid pilots, and WD 1856 diagnostics exist | Turn candidate checks into reproducible, versioned validation products |
| TWIRL-FM0 | FM0.2.1 completed its four-point canary trajectory; rank repair is learned, but masked Huber remains `0.109%` worse than zero, cross-sector retrieval did not improve, and the formal gate is incomplete | Admit and inventory genuinely later sectors, then freeze a development-only temporal panel for zero-shot evaluation; do not advance the model ladder |

The exact current run state and historical metrics are recorded in the
[progress log](twirl_progress_log.md). Reports are evidence snapshots, not
project authority.

## Active workstreams

### A2v1 production

- S56--S65 are accepted A2v1 sectors. Do not replay them.
- S66 completed all 32 ORCD cells, but only the two cam4/ccd4 cells have
  returned to canonical PDO storage. Its two cam4/ccd3 science payloads remain
  intact on ORCD; their cleaned source grids require a bounded re-ingress before
  retention. The other 28 cells remain retained on ORCD.
  S67--S78 completed all 32 cells and are retained as deferred checkpoints.
  None is accepted until returned to PDO and passed through the unchanged
  HDF5, FITS, schema, and QA gates.
- The S66--S93 campaign is computing S79 while S80's 32 prepared cells enter
  the live controller. The campaign follows the single frozen
  scientific recipe and the campaign resource profile in the
  [A2v1 protocol](a2v1_production_protocol.md). Stop on any failed scientific
  or provenance gate. Maintain the active sector plus three complete future
  sectors in the resident ORCD source buffer whenever capacity permits.
- S94+ will use the same frozen product recipe through a separate PDO reuse
  lane: shared QLP source/ePSF trees remain read-only, reusable empty-mask
  ePSFs are copied into a user-owned attempt, and only masked cutouts are
  refit. A fail-closed, read-only-by-default pilot wrapper is locally tested;
  no remote run has started. Start with one S94 `orbit-195/cam1/ccd1` smoke;
  full scaling waits for its unchanged gates and a variable-orbit queue for
  four-orbit sectors.

### Teacher v3 prospective use

- Teacher v3 is the frozen supervised morphology-enrichment baseline. It may
  rank a score-hidden human-review queue; it is not an automatic classifier or
  discovery engine.
- S63 is the sealed prospective sector. Primary hosts absent from the training
  corpus and repeated hosts are reported separately. Pre-score BLS failures
  remain in cohort accounting and never receive fabricated candidates.
- Teacher v4-SSL is complete and archived as representation-learning evidence.
  Its frozen linear probe was poor and its fine-tuned model did not decisively
  improve on Teacher v3, so it is not used for S63 or pseudo-labeling.

### TWIRL-FM0.1 proof of concept

- FM0.1 starts a new model family rather than extending the Teacher/CNN naming
  sequence. Its purpose is to expose mistakes and compare representations, not
  to claim a finished foundation model.
- The byte-exact scientific contract is the
  [frozen design](foundation_model_design.md),
  [configuration](../configs/models/twirl_fm0_1_s56_s67_poc.yaml), and
  [freeze receipt](../reports/stage5_validation/twirl_fm0_1_design_freeze_v1/freeze.json).
- The S56--S67 proof of concept is BLS-free. Gaia DR3 `source_id` is the
  scientific identity, but the fixed `80/10/10` train/development/sealed-test
  assignment is made by the connected Gaia--TIC `leakage_component_id`;
  ambiguous multi-Gaia components are quarantined. The release stores six views:
  `{raw-relative, ADP, ADP015} x {1x1, 3x3}`. It compares a TCN and Conformer
  through the frozen FM0.1.1--FM0.1.5 ladder. Repeated-sector hosts are an
  explicit consistency and aggregation experiment; their visits are not
  averaged into one light curve or counted as independent stars. Local elapsed
  time and deltas in nominal `200 s` units enter the window encoder; cross-visit
  offsets and gaps are retained only for audit and higher-level aggregation.
  Absolute BJD and explicit sector identifiers are not encoder inputs.
- S66--S67 may enter only as separately reported diagnostic checkpoints after
  FM-specific input gates. The first multi-sector release currently contains
  S56--S64: `342,721` immutable BLS-free shards and `3,764,787,143` cadences,
  each bound to the quality-aware A2v1 adapter and the fixed Gaia-component
  split. Its independent CPU release validation, pinned-environment loader and
  restart gate, fresh Stage-1-safe admission, and eight-step one-H200 FP32
  synthetic `FM0.1.1` mechanics test all passed with strict checkpoint and
  post-validation checks. The separately authorized real `FM0.1.1` TCN and
  `FM0.1.2` Conformer runs then completed 20,000 BF16 steps with one H200 each.
  Their [matched development comparison](../reports/stage5_validation/twirl_fm0_1_s56_s64_development_comparison_v1/README.md)
  found an optimization tie but effective ranks of only `6.20` and `8.58`
  against the provisional gate of `26`. Neither architecture is promoted,
  FM0.1.3 remains blocked, and the sealed test stays closed. The
  [FM0.2 objective-canary contract](../configs/models/twirl_fm0_2_s56_s64_objective_canary.yaml)
  keeps the exact S56--S64 release, split, TCN, ADP `1x1+3x3` views, masks,
  context, and optimizer while adding a scale-matched same-window VICReg loss
  directly to `z_window`. It is frozen for the single `FM0.2.1` TCN canary
  through step 2,000 only. Every H200 invocation requires a fresh Stage-1-safe
  admission and uses one H200, four CPUs, and 32 GiB. The sealed test,
  FM0.2.2, a 20,000-step extension, production use, and a foundation-model
  claim remain unauthorized. The exact-code CPU gate, FP32 smoke, real step-64
  milestone, and independent restart validation passed. Step 500 also passed
  strict validation. The matched step-500/1000/2000 development trajectory
  then raised `z_window` effective rank from `6.45` to `19.43` to `39.40` and
  made both separation intervals positive, but the final masked Huber remained
  `0.109%` worse than the frozen zero predictor and cross-sector retrieval was
  not demonstrated. The [step-2000 report](../reports/stage5_validation/twirl_fm0_2_s56_s64_step2000_evaluation_v1/README.md)
  is development-only evidence: the formal go/no-go was not applied,
  development event retention was not run, the exact same-seed step-0
  evaluation was initially missing, and the sealed test stayed closed. Detailed
  evidence is in the step-2000 report and its
  [new-Mac handoff](../reports/stage5_validation/twirl_fm0_2_s56_s64_step2000_evaluation_v1/HANDOFF.md).
  The dedicated [step-0 evaluation](../reports/stage5_validation/twirl_fm0_2_s56_s64_step0_evaluation_v1/README.md)
  then completed on the identical development panel with zero sealed access.
  The exact same-seed trajectory confirms that `z_window` effective rank was
  learned (`3.03 -> 39.40`) and masked Huber fell by `80.05%`, but
  cross-sector retrieval changed descriptively from `0.03380` to `0.02083`
  and the exact initialization already exceeded the separate seed-0 random
  control on paired separation. This is a targeted mechanism success, not an
  overall canary pass or useful-embedding result. The current operational
  campaign spans 25 sectors through S80, but only S56--S64 form a frozen FM
  release and only S56--S65 are accepted Stage-1 sectors. The fail-closed
  [later-sector policy](../configs/models/twirl_fm0_2_later_sector_admission_v1.yaml)
  now binds the real release, corpus-selection, and alias authorities and keeps
  new and repeated hosts separate. A separate ORCD source gate now recognizes
  the checksum-bound S65 archive and the 32-cell retained S66+ layout without
  promoting deferred sectors to accepted A2v1; real source receipts now cover
  S65--S78. A separate uniform mission-quality source gate follows current QLP
  semantics: SPOC before S67 and TICA from S67 onward, always joined with the
  detector QLP qflag. Provenance-bound mission-quality receipts pass S66--S77.
  S65 is now a predeclared, model-outcome-independent exclusion from the first
  later-sector panel and next training release because its camera-4/CCD-4 SPOC
  authority lacks `1,462` quaternion cadences; the machine-readable
  [exclusion ledger](../configs/models/twirl_fm0_later_sector_exclusions_v1.yaml)
  records the prohibition on retroactive insertion after repair. The frozen
  admission-policy v1 still requires S65 and therefore cannot be reused for
  the revised panel. The separately hashed admission-policy v2 consumes the
  exclusion ledger and requires exact S66--S77 quality, source, HDF5, and
  six-view bundles before it can pass; it does not freeze a panel or authorize
  training. S78 remains blocked by TICA rows absent beyond the delivered
  FFI/HDF5 tail. The full read-only ORCD HDF5/cadence
  chain passes S66--S77: all `1,356,463` declared products opened and
  `7,545,237,439` light-curve cadences reconciled, with no unreadable product
  or unexplained authority gap. S66 uses SPOC with zero exclusions. S77's `32`
  declared TICA-only detector rows propagate to `57,994` explicitly masked
  light-curve cadence occurrences. No temporal panel is frozen.
  Once at least five genuinely later sectors pass the ORCD HDF5-openability,
  six-view, cadence-quality, and identity gates, use the label-blind inventory
  to freeze the still-unset detector/cadence/cohort-adequacy thresholds and
  then a development-only panel. Evaluate the existing step-2000 checkpoint
  zero-shot before considering a separately controlled data-scale experiment.
- If that evidence justifies another training campaign, name it
  `TWIRL-FM0.3`; do not reuse `FM0.2.2`, which already denotes the blocked
  parameter-matched Conformer under the frozen FM0.2 contract. The provisional
  first candidate, `TWIRL-FM0.3.1`, is a fixed-TCN data-scale baseline: retain
  the FM0.2.1 two-ADP-view encoder and same-window objective, use only fully
  gated sectors, and compare at the same window-draw milestones before any
  longer exposure-normalized extension. Before freezing a changed-context
  candidate, run a matched development-only loader diagnostic with `256`,
  `512`, `1,024`, and `2,048` contiguous cadence-slot crops from the same
  immutable full-visit shards, holding physical-source sampling and total
  model-visible cadence exposure fixed. If shorter context is then selected,
  treat it as the single structural change in the provisional
  `TWIRL-FM0.3.2` candidate, not a silent continuation of `.3.1`. Freeze a
  genuinely newer temporal holdout before training. Subsequent candidates may
  optimize both independent reconstruction masks and then test a separate
  stable-host cross-visit head,
  but each is a one-mechanism ablation; raw/ADP015 view additions and a
  Conformer revisit remain later experiments. Exact same-seed initialization,
  source-paired changes, at least two frozen training seeds after the canary,
  new/repeated-host reporting, and event-retention checks are required before
  promotion. This naming/design direction is not training authorization.
- S63 may appear as unlabeled FM pretraining data, so it is never prospective
  FM evidence. The Teacher-v3 S63 test remains prospective only while no FM
  embedding, score, queue, threshold, or review decision influences it.

## Locked survey decisions

- **Scope:** search `200 s` TESS FFIs from `Sector >= 56` for transiting or
  occulting objects around white dwarfs.
- **Scientific identity:** Gaia DR3 `source_id`. TIC is operational metadata;
  ambiguous Gaia--TIC components remain grouped for leakage control.
- **Seed catalog:** Gentile Fusillo et al. (2021) Gaia EDR3 main catalog,
  retained as a local external input.
- **Reference sample:** `Pwd > 0.75`. The final statistical denominator and
  release sector cutoff are not yet frozen.
- **Benchmark:** WD 1856+534 in S56, orbits 119 and 120, camera 4/CCD 1.
- **Stage 1 product:** A2v1: no science magnitude cap, saturated-pixel ePSF
  masking, and ADP/ADP015 `1x1`, `3x3`, and `5x5` products.
- **Discovery policy:** transparent periodic and non-periodic searches precede
  machine-learning triage. A model score is not a detection, label, physical
  class, or completeness measurement.
- **Compute boundary:** PDO owns authoritative Stage 1 inputs and accepted
  products. ORCD may run the authorized source-only A2v1 compute campaign and
  downstream compact-data workflows. Deferred ORCD outputs are checkpoints,
  not accepted sectors.
- **Repository boundary:** the survey paper remains in the sibling
  `twirl-survey-paper` repository; proposal material is not vendored here.

## Stable product and interface rules

### Catalog and sample

- The master catalog has one row per Gaia DR3 source and keeps the original
  seed columns plus explicit sample-membership flags.
- Observation products retain sector, orbit, camera, CCD, TIC aliases,
  product/schema version, and source provenance without redefining the
  scientific identity.
- Gaia-only and ambiguous bridge cases remain visible. They are characterized
  before the survey denominator is frozen rather than silently dropped.
- The final release manifest names the sector cutoff, accepted product hashes,
  catalog and parent-sample versions, QA state, search configuration, and
  injection contract.

### Light curves

- The [A2v1 protocol](a2v1_production_protocol.md) and
  [machine-readable lock](../configs/a2v1_production_lock_v1.json) are the
  product authorities.
- HDF5 extraction is an intermediate checkpoint. A production sector becomes
  accepted only after sector FITS and the full edge-aware validation gates.
- Changing target selection, hooks, masks, apertures, detrending, schema, or
  scientific gates requires a new named product version. Hardware and
  concurrency may change only under the protocol's parity rules.

### Search, labels, and models

- Search tables preserve branch/config version, aperture, cadence cleaning,
  retained peaks, ephemeris, and cross-aperture evidence.
- Injection truth, BLS recovery, automated vetting, human morphology,
  morphology fold factor, pseudo-labels, and model scores remain distinct
  fields with distinct provenance.
- The legacy Teacher v3 release is grouped by TIC; it does not certify that
  distinct TIC aliases of one Gaia source stayed together. FM0.1 and future
  releases group the Gaia physical source and any unresolved Gaia--TIC leakage
  component before visits, windows, or augmentations are created. Evaluation
  and uncertainty follow the declared scientific unit.
- Calibration, subgroup behavior, and a transparent baseline accompany any
  model result used beyond exploratory development.
- Teacher and FM experiments may proceed as isolated research lanes, but they
  cannot enter production triage, discovery, or completeness claims before the
  transparent periodic and dip contracts are frozen.

## Stage gates

### Stage 1 — catalog, extraction, products, and QA (...)

Implemented: Gaia-first catalogs, TIC bridge, detector tables, A2v1 production,
FITS export, compact export, edge-aware validation, WD 1856 benchmark paths,
and bounded S56 QA.

Exit requires:

- validated HDF5 and FITS coverage for every admitted sector, with only
  declared edge omissions;
- a stable index joining identity, visits, paths, versions, checksums, cadence
  retention, scatter metrics, and QA state;
- a frozen release manifest and explicit sector cutoff;
- characterization of Gaia-only/no-TIC targets and the S94+ data boundary;
- Tier-1 science QA for population scatter, flags and gaps, aperture behavior,
  signal preservation, benchmark timing, and an independent extraction.

### Stage 2 — transparent search and candidate generation (...)

Implemented: per-sector multi-aperture BLS, retained peaks, candidate
consolidation, review tools, heuristic/LEO/centroid diagnostics, and the frozen
Teacher v3 enrichment workflow. The separate close-DWD public-TGLC diagnostic
is complete and remains outside the survey denominator.

Exit requires:

- reproducible periodic and dip-search baselines with versioned configuration;
- WD 1856 and deterministic synthetic smokes;
- multi-sector candidate objects and branch-aware false-alarm statistics;
- any promoted ranker to pass grouped, calibrated, source-separated tests and
  a declared transparent baseline comparison.

### Stage 3 — injection-recovery and completeness (...)

Implemented: pre-detrend LC injections, aperture/detrending/BLS audits,
candidate-retention diagnostics, and a pixel/source-pickle smoke.

Exit requires:

- recovery through extraction calibration, search, any production ranker,
  automated vetting, and candidate merging;
- supported completeness surfaces over magnitude, period, event size/depth,
  duration, sector count, crowding, and other declared nuisance dimensions;
- a measured LC-level versus representative pixel-level difference;
- explicit parent sample, product, and search versions for every result.

### Stage 4 — frozen full-survey inference

Stage 4 begins only after the Stage 1 index/release, Stage 2 periodic and dip
contracts, and Stage 3 recovery chain are frozen. It will use thin drivers in
`scripts/stage4_search/` to run approved branches across accepted sectors,
merge multi-sector evidence, and publish one versioned candidate catalog.

### Stage 5 — validation and follow-up

Exit requires reproducible aperture/pixel localization, crowding and neighbor
checks, odd-even and secondary-event checks, multi-sector consistency,
independent re-extraction, and appropriate archival or new follow-up. Periodic
and irregular events follow different confirmation paths. Facility access,
cadence, sensitivity, and lead time must be verified before a route is called
ready.

## Sample freeze and claim boundary

Freeze `twirl_parent_sample` only after coverage and Gaia-to-TIC losses are
characterized, the periodic and dip searches are fixed, the full recovery
contract is fixed, and every exclusion is encoded. Occurrence posteriors or
upper limits must be restricted to that declared sample, event class, support
region, and validated completeness surface.

TWIRL I is a framework and survey-method paper, not yet an occurrence-rate,
null-result, or discovery paper. A benchmark recovery, classifier metric,
embedding, or candidate-retention fraction cannot support those claims alone.

## Collaboration boundary

TWIRL is a Schwamb-group collaboration. Te Han owns pipeline and data-product
stewardship. Contributor roles, paper leadership, and discovery-response
responsibilities must be confirmed explicitly rather than inferred from code
or telescope access. Keep agreed roles and open project decisions here; private
negotiation or career-planning notes do not belong in the public plan.

## Immediate implementation priorities

1. Continue the frozen S66--S93 A2v1 campaign without disrupting the active
   S79/S80 controller. Once S80 ingress reaches a stable terminal state,
   re-ingress only the two cleaned S66 cam4/ccd3 source grids, retain their
   intact payloads, return S66 one cell at a time, and run the unchanged PDO
   acceptance gates before promoting or skipping any later sector. Then return
   and gate S67 onward in chronological order. Preserve Stage 1 priority over
   exploratory GPU work.
   In parallel, prepare
   only the isolated low-priority S94 shared-QLP reuse smoke described in the
   protocol; do not release a full S94+ queue before its gates pass.
2. Run the two isolated model experiments: complete the sealed Teacher-v3 S63
   score-hidden review, and keep its isolation rule. For FM work, preserve the
   completed four-point FM0.2.1 evidence and perform only the label-free
   later-sector FM admission/repeated-host inventory needed to define a
   genuinely later development panel. Evaluate the existing checkpoint
   zero-shot only after that panel is frozen. Keep the formal gate unapplied
   and do not run development event retention without explicit authorization.
   Do not advance FM0.1.3 or FM0.2.2, extend beyond step 2000, rerun a second
   seed, start a data-scale training comparison, or open the sealed FM test.
3. Advance the transparent survey path independently: implement the dip branch
   and multi-sector merging, freeze the archive/index and release boundary
   (including no-TIC and S94+ decisions), then complete the recovery chain
   required before Stage 4.

When a milestone changes, update this page briefly and place the dated details
in the progress log. Do not add detector/job diaries or full metric tables here.
