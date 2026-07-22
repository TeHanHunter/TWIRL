# TWIRL progress log

This is the live dated execution log. Use [the survey plan](twirl_plan.md) for
forward-looking decisions and acceptance gates. The complete April-through-
July 13 history is preserved in the
[archived log](archive/twirl_progress_log_through_2026-07-13.md).

Each active subsection ends with one `**Next:**` pointer. Add command-level
details here; keep the plan limited to milestone status.

## Infrastructure

### ORCD downstream compute

- `2026-07-11`: The production harmonic teacher trained on one H200 using the
  separate Torch environment. PDO remains the Stage 1 home; ORCD remains the
  compact downstream search/injection/model environment.
- `2026-07-13`: The five selected S56 checkpoints were transferred to the
  isolated PDO teacher-search run with recorded SHA-256 values. No new ORCD
  authentication should be initiated by automation when the user-owned control
  socket is absent.

**Next:** Use ORCD only for bounded, versioned downstream jobs required by the
current Stage 2/3 gates; do not expand GPU allocation or move raw TGLC trees.

## Stage 1

### A2v1 sector production and QA

- `2026-07-10`: S56 completed with `31,450` ADP/ADP015-only FITS. Full schema
  validation found zero bad FITS, zero non-edge omissions, and only documented
  detector-edge exclusions.
- `2026-07-13`: S57 completed with `27,213` FITS, zero build failures, zero bad
  schemas, zero zero-byte HDF5s, and zero non-edge omissions. Its compact ADP
  export is complete.
- `2026-07-13`: S58 HDF5 production completed after orbit `124` resumed from
  intact partial products; orbit `123` had completed before the original tmux
  parent ended without a recorded CCD-stage error. Required sector FITS and
  full validation remain next.
- `2026-07-15`: Added the reusable gated [A2v1 sector queue](../scripts/stage1_lightcurves/run_a2v1_sector_queue_pdo.sh) and its [S58-S63 manifest](../configs/a2v1_production_s58_s63.txt). All ten S59-S63 prepared orbit trees passed source/ePSF preflight (`3,136` each), and the serial tmux queue `twirl-a2v1-s58-s63-queue` started at `11:42 EDT` on `pdogpu6`. It begins with S58 HDF5 validation, then builds FITS and requires a full HDF5-plus-FITS schema pass before advancing through S63; it stops on the first failed gate.
- `2026-07-16`: The S58-S63 queue stopped at the S58 FITS gate after producing `23,139` FITS with two failures: TIC `1551609509` and TIC `1718164244` had zero-byte HDF5 outputs in orbit `124`, `cam4/ccd1`; their orbit-`123` counterparts were valid. The HDF5 coverage gate had only checked presence and size, so it was strengthened to open every nonzero HDF5 before FITS production. At `01:31 EDT`, the two zero-byte files were preserved with forensic suffixes and the targeted `twirl-s58-a2v1-h5-repair-r1` tmux job began the orbit-`124`, `cam4/ccd1` light-curve stage. It will validate S58 again, then resume the unchanged serial manifest.
- `2026-07-16`: The targeted repair completed in `0.09 h` with `1,974` orbit-`124`, `cam4/ccd1` HDF5 products. Both repaired TIC files were nonzero and opened successfully. The strengthened per-orbit gates passed: orbit `123` had `23,156` present HDF5s with `119` documented edge exclusions, and orbit `124` had `23,153` present with `122` documented edge exclusions; both had zero non-edge omissions, zero zero-byte files, and zero unreadable files. At `10:28 EDT`, the unchanged S58-S63 manifest resumed as `twirl-a2v1-s58-s63-queue-r2`; it begins by rebuilding and fully validating the incomplete S58 FITS product.
- `2026-07-16`: The `r2` queue's S58 HDF5 validator completed, but its `tee` subprocess blocked writing the full JSON report to the detached tmux TTY. No extraction or FITS process was active. The queue driver now sends verbose stage output only to the persistent queue log. After syntax and full-test validation, the deadlocked shell/logger were stopped without changing products, and the same manifest restarted as `twirl-a2v1-s58-s63-queue-r3` at `11:30 EDT`.
- `2026-07-16`: S58 completed at `12:03 EDT`. Its rebuilt full-product validation reports `46,309` present HDF5 rows, zero non-edge HDF5 omissions, zero zero-byte or unreadable HDF5 files, zero non-edge missing FITS targets, and zero bad checked FITS schemas. The `r3` queue then advanced to S59; orbit `125` ePSF/light-curve production is active.
- `2026-07-17`: S59 completed at `19:05 EDT` after HDF5 extraction and the required FITS/full-product gates. The `r3` S60 log stopped updating at `00:03 EDT` while `pdogpu6` became unresponsive; the shared tree retained `27,048` orbit-`127` and `15,390` orbit-`128` HDF5 files. After a clean `pdogpu5` GPU/runtime preflight, the unchanged manifest restarted as `twirl-a2v1-s58-s63-queue-r4-pdogpu5` at `11:00 EDT` on GPUs `4,5,6,7`. It revalidates accepted sectors and resumes S60 without deleting partial products.
- `2026-07-18`: The `r4` queue completed S60 (`15:45 EDT`), S61 (`01:27 EDT`), S62 (`10:59 EDT`), and S63 (`21:35 EDT`) after their HDF5, FITS, and full A2v1 schema gates. The terminal log reports `A2v1 queue complete`; no worker remained active.
- `2026-07-20`: Confirmed S64-S65's complete prepared source trees (`3,136` pickles per orbit) with no legacy ePSFs. The generic queue now permits either a complete reusable ePSF tree or an absent tree that forces all saturated-mask ePSF refits, while rejecting partial ePSF preparation. At `09:38 EDT`, the clean-clone tmux queue `twirl-a2v1-s64-s69-queue-r1-pdogpu5` started [S64-S69](../configs/a2v1_production_s64_s69.txt) on GPUs `4,5,6,7`; S64 preflight accepted both orbits with `epsf_mode=refit-all` and entered HDF5 extraction.
- `2026-07-20`: Rendered representative S56-S63 A2v1 pixel-mask diagnostics with the [mask-map renderer](../scripts/stage1_lightcurves/plot_a2v1_pixel_mask_maps.py): first-orbit `cam1/ccd1`, first `20` staged FFIs per sector. The left panel is the TGLC threshold-mask proxy; the right panel reports the corresponding per-cutout masked fraction and outlines local A2v1 ePSF refits, which are the nonempty-static-mask cases under the prefill contract. All eight PNGs opened successfully; the contact sheet is in `reports/stage1_lightcurves/a2v1_pixel_mask_maps/`.
- `2026-07-16`: Reclassified the existing S56 A2v1 QA as Tier-0
  integrity/benchmark QA. It remains valid evidence for product coverage and
  WD 1856 recovery, but science promotion now requires a separate Tier-1 pass
  with scatter-versus-magnitude limits, cadence-loss and aperture-outlier
  tests, fixed injection preservation, and a genuinely independent extraction.
- `2026-07-16`: Implemented the fail-closed S56 Tier-1
  `active_search_pair` contract for bounded enrichment. It scans both active
  ADP channels, emits TIC-level pass/review/fail flags, binds the exact compact
  product and a frozen `2,000`-host/four-shard injection canary, and requires
  authoritative cadence plus independent WD 1856 evidence. The production
  configuration remains intentionally non-runnable until the real external
  artifact hashes are reviewed; this scope cannot set `science_ready=true`.
- `2026-07-22`: Hardened the bounded S56 Tier-1 path around one authoritative
  cadence policy. The custom A2v1 compact product carries the internal TGLC
  flag only; the WD 1856 audit found `5.4268%` internally flagged cadences and
  another `4.2887%` flagged only by the external SPOC/QLP authorities. Tier-1,
  BLS, candidate-metadata, and real/injected teacher-native builders and
  validators now require the effective internal-OR-external mask and retain
  the cadence-table/manifest checksums. Added builders for the 16-detector SPOC
  intermediate and 32 orbit/detector QLP mappings, exact `(sector, TIC)` target
  eligibility, pinned Tier-0/BLS prerequisites, and a genuinely independent
  official TESSCut WD 1856 extraction. The downloaded S56 TESSCut product has
  `11,775` rows and passed the standalone extraction smoke; production evidence
  still requires the accepted compact product, authority files, reviewed
  hashes, and a full CPU run. No Tier-1 pass or science-readiness claim has
  been made. Because the effective mask changes native channels and BLS
  periodograms, the S56 native-input contract is now v2; candidate artifacts
  bind that v2 input while retaining a separately versioned provenance
  envelope. The old native-v1 checkpoint is intentionally incompatible. Any
  S56 path-validation retrain uses a distinct native-v2 namespace; the planned
  seven-sector retrain first requires a new observation-keyed multi-sector
  native contract and cannot simply reuse the S56-only file.

**Next:** Run the gated S64-S69 source-only production queue while the real S56
cadence/TESSCut evidence is published; regenerate the quality-aware BLS table,
rerun and pin Tier 0, then execute the bounded Tier-1 CPU gate before
enrichment uses S56.

### Catalog, archive index, and sample control

- The master catalog contains `1,280,266` seed rows, `359,073` rows with
  `Pwd > 0.75`, `1,236,467` unique TIC matches, `43,035` ambiguous matches,
  and `764` no-bridge cases.
- A consolidated production light-curve index and final parent-sample freeze
  criteria are not yet complete.

**Next:** Define the versioned A2v1 compact-export/index schema and a frozen
TWIRL I release manifest with an explicit sector cutoff; then characterize the
no-TIC cases and S94+ QLP ingestion before freezing the statistical parent
sample.

## Stage 2

### Transparent search and active-learning triage

- `2026-07-11`: `s56_harmonic_cnn_v1` selected the Shape + periodogram/BLS
  profile. The locked test has balanced accuracy `0.750`, macro F1 `0.757`, and
  ECE `0.048`, but only two real Planet-like test examples; the model is an
  enrichment ranker, not a student-label generator.
- `2026-07-13`: S56 A2v1 Tier-0 QA passed and WD 1856 is rank one in both ADP
  apertures. Native input assembly/scoring continued in bounded batches. S57
  BLS/scoring was designated experimental only.
- `2026-07-16`: The frozen teacher-v1 injection diagnostic contains `20,000`
  injections: BLS recovered `3,458` (`17.29%`) within its top five peaks and
  teacher v1 retained `1,527` (`7.635%`) through the tested BLS-to-teacher
  chain. These are candidate-retention diagnostics, not end-to-end survey
  completeness.
- `2026-07-15`: Teacher v2 completed as an exploratory comparison and opened
  its locked tests once. On the host-disjoint S57 injection set, BLS recovered
  `2,777/14,706` signals and teacher v2 retained `2,144/2,777` (`77.21%`),
  versus `56.25%` for teacher v1 and `61.94%` for metadata-only. The relative
  gains passed, but the preregistered `>=80%` absolute gate did not. The
  same-row real holdout independently blocked promotion: macro F1 was `0.653`
  versus teacher v1's `0.725`, with Eclipse/contact recall `0.364`. Teacher v1
  remains the accepted active-learning baseline; teacher v2 is neither a
  production ranker nor a student-label generator.
- `2026-07-15`: The focused S56 compact revisit completed `407/407` sheets.
  Seven Planet-like rows were known positives and the `400` model-selected
  compact rows contributed `11` new Planet-like labels (`2.75%` yield). The
  separate S56 mixed `1,000`-TIC batch remains a partial `177/1,000`
  checkpoint and must be completed and audited before expanding review.
- `2026-07-16`: The exploratory teacher-v2 ensembles scored `136,060`
  candidates over `27,212` real S57 TICs and produced a verified blinded
  `1,000`-TIC queue with the frozen category quotas and `S56-ADP-HV2` evidence
  sheets. Queue generation is complete, but human review is not: the live file
  contains exactly six nonempty labels. Preserve all six as premature
  experimental evidence, pause further S57 holdout consumption, and do not
  treat S57 as a pristine external holdout.
- `2026-07-17`: The accepted frozen teacher-v1 ensemble produced a blinded
  Franklin handoff of `3,000` real sector-observation rows: `1,000` each from
  S57, S58, and S59. Every row uses the ADP-small BLS rank-one ephemeris; the
  hidden quotas are `1,200` compact-transit, `900` Eclipse/contact, `300`
  Smooth-variable, `300` disagreement, and `300` controls. The package has
  `3,000` current `S56-ADP-HV2` PNG sheets, `14` real reference examples,
  zero PDFs, and no exposed model scores or selection provenance. It is an
  active-learning review set, not a new teacher test set. The self-contained package passed
  its browser `--check-only` validation in place on PDO with all `3,000`
  sheets resolved and a group-writable label CSV for Franklin.
- `2026-07-21`: Franklin returned all `3,000` labels with complete internal
  row/key/provenance coverage (SHA-256
  `2bd4d86870c70091eb7291ced067c63bc908118fd730083bfb2d12d52c5a09bf`).
  The return has `15` Planet-like labels rather than the reported `14`, plus
  `121` Eclipse/contact and `106` Broad-isolated-dip labels. Nine of the `15`
  Planet-like rows carry confidence, EB, wide-event, or faintness caveats, so
  these are enrichment morphology labels rather than confirmed planets.
- `2026-07-21`: Cross-sector auditing found `133` repeated TICs against the
  active real S56 training table and `42` active-label differences, including
  `23` against explicitly final S56 adjudications. The comprehensive second-
  read set contains `347` unique rows after including planet/EB/wide labels,
  cross-taxonomy notes, unresolved/refold cases, nine nonblank `uncertain`
  rows, and S56 differences. Combined S56 plus S57--S59 has `4,702` labeled
  sector observations but `4,569` unique TICs; every retraining and evaluation
  split must therefore be TIC-grouped. The return's standalone-app key is not
  the repository `label_io.candidate_key`, so it must be joined only to the
  exact frozen `3,000`-row queue under the app-key contract.
- `2026-07-21`: The user inspected the `347` special-case vet sheets and
  accepted Franklin's S57--S59 return at the batch-level morphology layer. The
  exact public queue (SHA-256 `6431ceef...`), returned labels (`2bd4d868...`),
  and normalized `3,000`-row morphology table are frozen in the
  [label-return report](../reports/stage5_validation/franklin_s57_s59_label_return_20260721/README.md).
  Franklin remains the original
  labeler; TeHan's acceptance is recorded separately. All `3,000` period
  factors/statuses remain audit metadata and produce zero harmonic targets.
  The model will still ingest all seven folds around the unchanged rank-one
  BLS ephemeris.
- `2026-07-21`: The next S60--S62 handoff is active on PDO with `3,000` rows
  (`1,000` per sector), `3,000/3,000` PNG vet sheets, and queue SHA-256
  `01f1e627...`. ORCD access is healthy and the accepted five-fold teacher-v1
  checkpoints remain hash-verified, but no retraining job was launched. The
  estimated merged S56--S59 corpus has about `42` unique real Planet-like TICs,
  still below the `50`-source student/promotion gate.
- Transparent per-sector BLS exists; the non-periodic dip branch and
  multi-sector aggregation remain unimplemented production gates.

**Next:** Apply the Tier-1 target mask and complete the bounded S56 and active
S60--S62 reviews. In parallel, design, implement, and validate an observation-
keyed multi-sector native contract plus one immutable TIC-grouped split
registry; only then assemble the seven-sector inputs. Freeze that corpus and
retrain teacher v1 once with source-separated, calibrated evaluation and
bootstrap uncertainty. Add dip, multi-sector, and false-alarm branches only
after this path is robust.

### Human labels and harmonic review

- The strict S56 adjudication audit completed `343/343` decisions across `323`
  unique real sources, with repeat agreement `0.85` and Cohen's kappa `0.814`.
- The locked adjudicated training set contains only `12` Planet-like sources.
  The compact revisit adds `11` new Planet-like enrichment candidates, but they
  remain separate evidence until the bounded S56 audit merges and resolves
  them. Student pseudo-labeling requires at least `50` unique real Planet-like
  sources and the locked real-data gates.
- Rare harmonic supervision remains inadequate for promoting or iterating the
  exploratory teacher-v2 design: `P/3` has no supervised example and `3P` has
  only three in the current harmonic table.
- Franklin morphology returns do not supply harmonic truth. The app's default
  `P/resolved` state and any saved morphology click are insufficient evidence
  that a factor was reviewed. Returned factors remain raw audit fields and are
  masked unless an explicit factor-only review or injection truth verifies
  them; all seven folds remain model inputs regardless of that mask.

**Next:** Use teacher v1 to enrich real labels while preserving blinded score
provenance and cross-review overlap; require a predeclared rare-factor and
oracle-factor evaluation before any future model iteration.

## Stage 3

### Injections and completeness

- LC-level raw-flux pre-detrend BATMAN products, balanced and all-host grids,
  aperture/detrending/BLS audits, and ORCD compact workflows are established.
- Pixel/source-pickle injection smokes run through extraction, TWIRL-FS, and
  BLS, but no representative pixel-level calibration set or frozen
  extraction-to-candidate completeness surface exists.
- `2026-07-22`: Bound the existing A2v1 candidate-retention recovery launcher
  to the native-v2 selected-checkpoint manifest. Its scorer now rejects
  mixed-generation candidate/native/checkpoint inputs, stages outputs
  atomically, and records exact before/after hashes. This was infrastructure
  hardening only; no new recovery or H200 job was launched.

**Next:** Run candidate-retention recovery against the accepted A2v1 ADP pair
and frozen Stage 2 contract, then build a pixel-level calibration subset
spanning magnitude, crowding, detector edge, and aperture disagreement before
claiming end-to-end completeness.

## Stage 4

No frozen full-survey inference run exists.

**Next:** Begin only after the archive/index, periodic+dip contract, and full
extraction-to-candidate injection gates are frozen; add thin drivers under
`scripts/stage4_search/` rather than placing inference in Stage 5.

## Stage 5

### Candidate validation and follow-up

- Heuristic/LEO checks, two-aperture review, human adjudication, scalar centroid
  diagnostics, and a WD 1856 pixel-map diagnostic exist.
- A reusable pixel-level on-target test, archival ZTF/ATLAS pass, and frozen
  follow-up-readiness workflow remain open.

**Next:** Turn the pixel-map smoke into a reproducible Tier-2 candidate check
before a publication-ready catalog; verify any proposed facility before
recording it as an available follow-up path.

## Documentation and repository hygiene

- `2026-07-13`: Reconciled the live plan with validated S56/S57 A2v1 products,
  the active S58 run, and the teacher-v1 enrichment boundary. The prior detailed
  log was frozen under `doc/archive/` so obsolete `Next` pointers no longer
  compete with current work.
- `2026-07-13`: Refreshed onboarding, science framing, operational-guide
  lifecycle labels, and script ownership; made `pyproject.toml` authoritative;
  added repeatable documentation and detection smokes. Plain pytest now imports
  `twirl` without a manual `PYTHONPATH`, and the full suite passes.
- `2026-07-16`: Completed a repository-wide audit and cleanup pass, hardened
  A2v1 validation, HLSP array-length checks, Stage 2 configuration/provenance,
  BLS-coverage QA, and orchestration failure handling, and reconciled the live
  plan with the exploratory teacher-v2 and enrichment artifacts. The compact
  repository audit under [reports](../reports/README.md) records the evidence,
  figures, deferred risks, and adjusted critical path.
- `2026-07-17`: Consolidated ORCD to one fresh, clean Git checkout and retained
  every prior checkout payload in a dated, inventoried archive. Code deployment
  now pins a Git SHA and rejects dirty checkouts. A separate sparse clean PDO
  checkout is staged against the existing local data tree; the active queue and
  live review session remain on the preserved legacy checkout until they reach
  a safe handoff point.

**Next:** Keep one authoritative priority list in the plan, one current pointer
per live-log subsection, mark report-level status files as snapshots, and keep
large generated payloads outside git while retaining compact provenance.
