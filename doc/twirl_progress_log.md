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
- `2026-07-23`: The first full quality-aware BLS attempt failed closed on
  cadence `699957`, which is present in the frozen S56 compact product for all
  camera-3 detectors but absent from the authoritative camera-3 quaternion
  timeline. Cadence-reference builder v3 now publishes the exact four
  detector/cadence rows in a hash-bound `authority_exclusions` sub-contract.
  The v2 external overlay assigns only those declared rows reserved bit `62`,
  masks and counts them, and still rejects every undeclared coverage gap.
  This provenance propagates through real BLS, injection and independent
  evidence, native HDF5, shard merging, candidate metadata, and teacher-input
  validation. The combined suite passed `413` tests with `3` skips, plus the
  detection sample and documentation checks; an independent review also
  caught and fixed a shard-merge publish-order failure before deployment.
- `2026-07-23`: Rebuilt the S56 cadence-reference manifest on PDO while
  preserving the `188,396`-row table hash. The new manifest binds four exact
  exclusions (`cam3/ccd1` through `cam3/ccd4`, cadence `699957`) and the
  official TESSCut WD 1856 evidence was regenerated under the v2 overlay.
  Both active apertures retain the accepted recovery: the 1x1 reference has
  period `1.407960330 d` and depth S/N `11.87`; the WCS-defined 2x2 reference
  has period `1.407896973 d` and depth S/N `11.33`. The legacy v1 four-shard
  injection canary was independently rehashed after its ORCD cleanup
  relocation: all `2,000` injection IDs, immutable metadata, compact-source
  parity, and ordered shard hashes match that frozen contract. The restarted
  `31,450`-target BLS run and its replacement Tier-0 audit both completed on
  `pdogpu1`. The table contains all `62,900` expected target/aperture pairs,
  with five `too_few_cadences` results and no missing pairs; all six Tier-0
  gates pass. The reviewed Tier-0 report and BLS hashes are now locked before
  the final ORCD Tier-1 population audit. Exact v1 evaluator commit
  `4831713f` was pushed and deployed, but CPU job `18641848` failed closed
  before the population scan: the v1 retention predictor measured
  normalization transfer rather than normalized depth retention, and one
  frozen epoch had no effective-good in-transit samples. Correcting for the
  stored injection baseline and detrending scale gives median retention
  `0.9908` in both apertures with tenth percentiles near `0.966`; no population
  products or scientific Tier-1 conclusion were emitted by the failed run.
- `2026-07-23`: Tier-1 contract v2 and the quality-aware epoch policy were
  locked in producer commit `f1b8b53c`, pushed, and deployed as a clean
  detached ORCD worktree. CPU array job `18650874` regenerated shard indices
  `00`, `13`, `26`, and `39` from the unchanged frozen schedule. Independent
  audit job `18651439` passed all `2,000` injection IDs and `2,000` unique
  hosts with a minimum of two effective-good in-transit cadences and zero
  alignment, copied-flux, raw-equation, uncertainty-equation, stored-quality,
  or provenance failures. The ordered shard hashes are `fe999651…6525`,
  `7b166ddf…fa07`, `4d6f8f11…ce9b`, and `9e309d6b…e5a4`; thresholds were not
  relaxed. The first v2 evaluator preflight, CPU job `18651862`, then failed
  closed on the legacy metadata digest before population scanning. Exact
  old/new manifest comparison showed bit-for-bit identity for TIC, TESS
  magnitude, detector, period, duration, and shard; the effective-quality
  policy resampled all `2,000` epochs and changed the cadence-sampled
  `model_depth` for `1,591`. Replacing only that realization field transforms
  the legacy digest `54996889…1910` into the observed v2 digest
  `b56a0b84…f0b0`, so the v2 lock was updated without changing the sample or
  thresholds. A future contract should separate schedule/host and
  epoch-realization digests explicitly.
- `2026-07-23`: The exact full-population v2 evaluator at commit `5adae3b6`
  completed on ORCD as CPU job `18652943` in `4 min 58 s`. Seven of eight
  gates pass: cadence authority, source parity, Tier-0, population scatter,
  aperture behavior, fixed injections, and independent WD 1856 extraction.
  The overall status is `review`, with `passed=false`,
  `enrichment_ready=false`, and `science_ready=false`, because the
  detector-stratified cadence gate is review for cam1/ccd1. Its median usable
  fraction is `0.7434`, below the predeclared `0.80` pass threshold despite
  zero cadence loss and every cam1/ccd1 target passing its individual cadence
  check. The locked external mask isolates the cause to orbit 119
  (`65.5%` quality-zero versus `99.0%` in orbit 120), dominated by one
  contiguous `1,929`-cadence flagged interval rather than a cadence-join gap.
  The published target mask contains `30,115` pass, `1,303` review, and `32`
  fail rows. It is not authorized for enrichment while the overall gate
  remains review.
- `2026-07-23`: Reinterpreted Tier-1 as a conservative input-usability gate
  after the v2 result showed that detector-level flagged fraction could veto
  thousands of otherwise searchable light curves. Contract v3 made
  population, cadence-fraction, and aperture reviews nonblocking, retained
  their metrics for sensitivity calibration, and kept all 16 detector cells.
  ORCD job `18655927` completed successfully, but independent review found
  that its target count used finite normalized flux without finite time and
  could let flagged flux values influence normalization. The v3 mask was
  therefore not promoted.
- `2026-07-23`: Contract v4 at commit `2633e8b4` now shares the exact BLS input
  preparation between search and Tier-1. ORCD job `18657818` completed in
  `9 min 17 s` with `passed=true`, `enrichment_ready=true`, and
  `science_ready=false`; the sole warning is the retained cadence-fraction
  diagnostic. The paired teacher keeps `31,446/31,450` targets, while
  `31,449/31,450` have at least one searchable aperture. Three exceptions have
  an entirely nonfinite small aperture but a healthy primary aperture; only
  one has neither channel. None is a negative label. The complete
  `628,955`-row BLS binding and all published output hashes validate, and the
  local suite passes `420` tests with `3` skips plus detection and docs checks.

**Next:** Keep the gated S64-S69 source-only production queue active in
parallel. Apply the v4 paired-input mask when the S56--S62 label corpus is
frozen; retain the three primary-only cases for manual or single-aperture
review and do not spend the current critical path on additional Tier-1
threshold work.

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
- `2026-07-23`: A schema-aware scan of tracked artifacts and reachable git
  history found S63 only in Stage-1 production/status QA, with no S63 rows in
  candidate, vetting, label, training, scoring, or evaluation products. S63
  can therefore be reserved conditionally for the first prospective teacher
  evaluation after non-git PDO/ORCD/local exposure is also audited. Freeze the
  S56--S62 checkpoint, calibration, thresholds, TIC-disjoint cohort, and
  metrics before blind S63 labeling; report repeated hosts separately and do
  not tune from the unblinded result.
- `2026-07-23`: Preflighted one final local Planet-like/EB pass over the
  accepted labels available through S59. The current queue has `223` candidate
  observations across `207` TICs (`45` prior Planet-like and `178`
  Eclipse/contact), preserves `230` source decisions and all `65` verified S56
  harmonic targets, and resolves all `223/223` PNG sheets locally. Its buttons
  are prefilled but `0/223` rows are marked reviewed; wait for the accepted
  S60--S62 return, rebuild once, and then freeze TeHan's explicit row-level
  pass.
- `2026-07-23`: Added two bounded pre-training contracts. The split registry
  deterministically fixes whole TICs into a nearest-feasible `20%` test
  population plus five nonempty development folds, and training preparation
  now consumes those assignments without recomputing them. The native registry
  binds each `(sector, TIC)` observation to an explicit HDF5 file/group,
  contract version, and SHA-256, preventing repeated hosts from colliding in
  TIC-only storage. The full fast suite passes `446` tests with `3` skips.
  Actual seven-sector registry artifacts and trainer/checkpoint binding remain
  pending the frozen label corpus and regenerated per-sector native inputs.
- `2026-07-24`: Franklin returned the complete S60--S62 handoff. A strict
  exact-key join against the frozen `3,000`-row queue passed with `1,000` rows
  per sector, no missing or duplicate decisions, and returned-label SHA-256
  `852a3135...`. The return contains `27` Planet-like, `259` Eclipse/contact,
  `215` Broad-isolated-dip, `1,449` uncertain, `602` stellar-variability,
  `442` instrumental/systematic, and `6` skip labels. The exact queue, return,
  normalized labels, and handoff manifest are frozen in the
  [label-return report](../reports/stage5_validation/franklin_s60_s62_label_return_20260724/README.md).
  Acceptance is batch-level staging only; all returned factor choices remain
  audit metadata with zero new harmonic targets.
- `2026-07-24`: Rebuilt the final S56--S62 Planet-like/EB review set from all
  accepted inputs. It has `509` sector observations across `493` unique TICs:
  `72` prior Planet-like and `437` prior EB labels, corresponding to `69` and
  `426` class-specific unique TICs. The provisional Planet-like support is
  therefore `19` TICs above the `50`-source training gate. All `509/509` exact
  PNG sheets are present locally and the app preflight passes with `0`
  explicitly reviewed rows. Sixteen TICs repeat, two have cross-sector class
  disagreements, and one TIC has two distinct S56 ephemerides; these remain
  separate observations for the final pass.
- `2026-07-24`: Replaced the mixed handoff-era review images with one uniform,
  versioned `s56_s62_a2v1_current_adp_v1` set. PDO re-rendered all `509` rows
  from `DET_FLUX_ADP_SML + DET_FLUX_ADP`, each row's frozen ephemeris, and the
  current `S56-ADP-HV2` renderer with `20,000` trial periods and `10` retained
  peaks. All seven sector summaries report only `ok`, and the local
  identity/provenance gate verifies the input-export identity and full
  decodability of all `509/509` PNGs. The final launcher now uses TeHan's
  standard repository vetter with exact-name sheet binding, suggestions that
  remain unreviewed until explicitly saved, `P/3` and `3P` controls, and a
  global Planet-like-first then EB order.
- `2026-07-24`: TeHan completed the uniform `509`-row Planet-like/EB pass and
  declared the untouched suggestions accepted after inspecting the full
  queue. The immutable
  [decision snapshot](../reports/stage5_validation/s56_s62_morphology_corpus_teacher_v3_v1/)
  keeps the raw `22`-row click log separate from the materialized acceptance:
  `22` explicit overrides plus `487` accepted-unchanged decisions. The
  accepted table SHA-256 is
  `47ffb23909dd2aceb3055d33969a903a356d6a79617b287ff172a6f156694278`;
  final labels are `70` Planet-like, `419` Eclipse/contact, `10`
  Broad-isolated-dip, `7` stellar variability, `2` uncertain, and `1`
  instrumental/systematic. Exactly `64` S56 harmonic decisions remain
  explicitly verified.
- `2026-07-24`: Froze the seven-sector `teacher_v3` data release while keeping
  the architecture identifier `s56_harmonic_cnn_v1`. The
  [corpus and split](../reports/stage5_validation/teacher_v3_s56_s62_a2v1_current_adp/)
  contain `8,181` rows (`7,724` real plus `457` injections), `8,163` active
  training rows, and `7,903` TIC groups. The corpus SHA-256 is
  `7f05bd9e8b7c7b04c7cf6ec0b45df6300fb0cb08623acbfbec65fb85cd66ebd7`;
  the immutable `20%` test/five-development-fold registry SHA-256 is
  `551747a18175ecc3c4dbcc9d464e2f521f0550cb9e045c28c6a7580455922f7c`
  and has zero TIC leakage. All seven sectors and all accepted label sources
  contribute; Franklin's period factors remain audit-only. S63 remains
  reserved for a later prospective evaluation.
- `2026-07-25`: Deployed the exact `teacher_v3` launch checkpoint
  (`2f4a00a4`) to a clean detached ORCD checkout and restaged the frozen corpus,
  split, and split-bound table from that Git object. Their ORCD SHA-256 values
  match the local freeze (`7f05bd9e...`, `551747a1...`, and `04a565a6...`).
  ORCD job `18882429` (`twirl-tv3-meta`) is the native-independent,
  metadata-only bootstrap: it entered `RUNNING` on `node4900` in CPU-only mode
  and emitted the `teacher_v3_start` record for all `8,181` source rows. It
  uses the immutable TIC grouping (`6,534` active development rows and `1,629`
  active fixed-test rows), five development folds, fold-local feature
  normalization, and pooled development-OOF temperature fitting. The fixed
  test remains sealed in this bootstrap. No H200 was requested or displaced.
  The full shape/BLS model, clustered-bootstrap intervals, and a genuinely
  retrained uncertain-masked sensitivity remain pending.
- `2026-07-25`: ORCD job `18882429` completed the five-fold metadata-only
  control in `36m41s` with exit `0:0`. Development OOF balanced accuracy is
  `0.812`, macro F1 is `0.806`, real Planet-like recall is `0.625`, real EB
  recall is `0.638`, and ECE is `0.0070`. These are development-control
  diagnostics only: the fixed test remained sealed and the result is not the
  full Teacher v3 model.
- `2026-07-26`: Prepared the one-H200 Teacher v3 launch path without opening
  the fixed test. PDO audits pass for all S56--S62 cadence authorities and
  exact-host raw subsets; the `1.77 GB` handoff tar is SHA-256 matched on ORCD
  and extracted into an isolated staging tree. The remaining large-file hashes
  and compact-path map paused when the user-owned ORCD control socket expired.
  The [training code](../src/twirl/vetting/teacher_v3_training.py) now requires
  CUDA for the full run, rejects dirty output directories, verifies native
  inputs and all frozen checkpoint bytes before and after test evaluation,
  imports the sealed metadata control by hash, independently retrains the
  uncertain-masked sensitivity, and produces `2,000`-draw TIC-clustered
  intervals. The S56 `300 + 156 = 456` injection-pair inputs have an atomic
  exact-group promotion gate, native preparation defaults to one CPU shard,
  and plotting is a dependent CPU job. Independent review found no remaining
  launch blocker; the full suite passes `487` tests with `3` skips. No H200
  training job has yet been submitted.
- `2026-07-27`: Completed all seven observation-keyed native-v2 inputs and
  the merged registry from exact Git object `0823d49a`. Registry job
  `18918806` passed with `7,705` real observations, `7,571` unique TICs, and
  `134` multi-sector TICs; its SHA-256 is `81b258590c...`. S62 required one
  checksum-bound compatibility correction because both raw and compact inputs
  inherited an orbit-partition defect: `65,478` target-cadence rows (`0.606%`)
  were remapped from orbit 132 to the authoritative orbit 131 assignment.
  Times, fluxes, uncertainties, quality values, labels, and splits were
  unchanged. This is accepted for Teacher v3 only; a science-ready S62 product
  still requires an upstream rebuild.
- `2026-07-27`: ORCD job `18918808` completed the fixed primary and
  independently retrained uncertain-masked profiles on exactly one H200 in
  `02:17:59` with exit `0:0`; all ten newly trained folds are present. The
  fixed test opened once after pretest freeze
  `140ad51dcbbad66cc8fbb10a9fd138f383fdeb0698be47cc1f083a3e48ed2eb3`,
  and all metadata, primary, and sensitivity checkpoint/calibration hashes are
  identical before and after evaluation. On the common `528`-row,
  `519`-TIC real non-uncertain support, primary versus metadata-only balanced
  accuracy is `0.779` versus `0.720`, macro F1 is `0.809` versus `0.744`, EB
  recall is `0.812` versus `0.694`, and ECE is `0.025` versus `0.054`. Real
  Planet-like recall is `8/14 = 0.571` and remains descriptive. The
  uncertain-masked retraining reaches balanced accuracy `0.817` and macro F1
  `0.818` on that same support after retaining `3,333/8,163` active rows,
  making label policy the clearest near-term model-development issue. Student
  training remains blocked and automatic production promotion is false.
  Dependent CPU job `18918849` rendered and hash-verified the
  [PNG, PDF, metrics, and short report](../reports/stage5_validation/teacher_v3_s56_s62_a2v1_current_adp/performance/).
  Final local validation passes `488` tests with `3` skips plus the
  documentation checker.
- `2026-07-27`: Froze the
  [Teacher SSL v1 preregistration](../reports/stage5_validation/teacher_ssl_v1_s56_s62_a2v1_current_adp/)
  and deployed the corrected exact Git object `634e865f` to a new detached
  ORCD worktree.
  The primary fold-local VICReg pool contains `6,168` real development
  observations (`6,054` TICs), including `3,780` uncertain observations whose
  labels are ignored by the SSL objective. Fixed-test and S63 model tensors,
  injected SSL rows, and scalar metadata are forbidden; immutable containing
  files are read only for integrity validation. The model-facing fine-tuned
  name is **Teacher v4-SSL**, but Teacher v3 remains the frozen operational
  baseline and automatic promotion/student labeling remain blocked. Local
  validation passes `514` tests with `6` optional-dependency skips. ORCD
  preflight job `18999895` passed all `38` Torch/model assertions and exposed
  one test-harness-only failure because it invoked a bare `python`; the
  corrected remote asset suite then passed `16/16`, and the obsolete dependent
  training job was canceled. Replacement one-H200 job `19000602` entered
  `RUNNING` on `node4900`, confirmed an H200 through `nvidia-smi`, and emitted
  the checksum-bound `teacher_ssl_v1_start` record with CUDA required. The user
  approved up to four H200s for later parallel folds/seeds or a validated
  broader pool; the small single-process native pilot intentionally requests
  one. The initial runner emits embeddings, neighbors, matched development
  OOF, and the exact-support Teacher v3 comparison; the preregistered
  permutation, collapse, confound, review-budget, injection, and WD 1856 gate
  evaluator remains the next dependent analysis.
- `2026-07-29`: Completed and registry-validated all `112/112` fresh
  Teacher v4-SSL full-pool native-v2 shards for the exact `175,347` eligible
  observations, with the `19` locked model exclusions retained in the
  `175,366`-row audit. The initial five-fold run then failed before epoch 1.
  One-H200 diagnostic job `19198212` completed normally on `node4900` and
  isolated the failure to observation `s0060-tic0000000722078603`, cadence
  `738841`, whose `quality=1` sentinel remained model-active. Its extreme
  harmonic activation first became nonfinite in
  `model.harmonic_encoder.sequence_encoder.stem.1.norm`; the same failure in
  BF16 and FP32, with finite pre-forward model and optimizer state, confirms a
  model-input LayerNorm overflow rather than an H200 fault. A clean local
  worktree now implements quality-aware photometry/error masking while
  preserving phase and quality, plus a fail-closed numerical release gate over
  all `175,347` eligible rows. No corrected remote numerical gate, fresh smoke,
  or replacement five-fold result is claimed yet.
- `2026-07-29`: Corrected full-pool numerical diagnostic array `19204195`
  completed all `112` checksum-valid reports and failed closed: `106` Slurm
  tasks completed, `6` failed, and `175,338/175,347` observations passed.
  Exactly nine S60 observations exceeded the harmonic/local numerical
  envelope; their periodograms remained normal. Read-only reconstruction
  traced all nine to compact ADP03q splines fitted before the authoritative
  external-quality overlay, allowing `801` external-only bad cadences in the
  affected detector/orbit to contaminate later nominal-quality values.
  Recomputing both ADP apertures in memory from immutable raw-v1 flux and
  uncertainty with the final effective quality reduced every affected
  transformed maximum to at most `6.44`, without clipping, imputation, or a
  new exclusion. No numeric release, smoke, or fold was launched. A wholly
  fresh, uniform effective-quality-aware ADP native rebuild is the recommended
  next contract and remains pending explicit approval.
- Transparent per-sector BLS exists; the non-periodic dip branch and
  multi-sector aggregation remain unimplemented production gates.

**Next:** Complete the Teacher v4-SSL follow-up gates, freeze Teacher v3 as
enrichment only, publish the reviewed S56--S62 candidate table/TIC roll-up,
and define the sealed S63 queue, thresholds, random-control slice, and
prospective metrics before blind review. Do not schedule more corpus labeling
or start student pseudo-labeling. Add dip, multi-sector, and false-alarm
branches only after this prospective periodic path is robust.

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

**Next:** Use the frozen Teacher v3 release to enrich the sealed S63 review
while preserving blinded score provenance and a random/control slice. Require
the prospective result plus a predeclared rare-factor and oracle-factor
evaluation before any architecture iteration or student pseudo-labeling.

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
