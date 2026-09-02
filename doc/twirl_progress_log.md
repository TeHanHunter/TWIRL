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
- `2026-08-07`: Verified the new Blackhole/Globus transport route and installed
  the persistent Globus CLI. `tehan` has `tso`/`globus` membership and write
  access under `/globus/tso`; both `TESS TSO` and `MIT ORCD Engaging
  Collection` are CLI-accessible with refreshable consent. A checksum-verified
  4 GiB/eight-file transfer completed without faults but reached only
  `30.8 MB/s` (`246 Mb/s`) over 139 s, despite Globus-selected concurrency 4,
  parallelism 8, and pipelining 20. Eight simultaneous tasks did not improve
  aggregate throughput, while a Blackhole-local write reached `706 MB/s`.
  Temporary probe payloads were removed. The detailed [transport report]
  (../reports/infrastructure/blackhole_globus_transfer_probe_20260807.md)
  records task IDs, timings, access state, and the proposed source-pickle
  streaming boundary. A user-opened Blackhole-to-PDO control socket and a
  CuPy/TGLC H200 parity environment remain required before production.

**Next:** Ask MKI/ORCD to diagnose the approximately `250 Mb/s` managed-endpoint
path versus the expected near-`1 GB/s` rate. Then prove a single S65
orbit-138/cam4/ccd1 prepared-source A2v1 smoke on one H200 before enabling the
authorized two-H200, at-most-78-CPU streaming lane; PDO remains authoritative.

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
- `2026-08-06`: A live read-only PDO audit reconciled the later-sector queue.
  S64 is fully accepted: `159,781/160,656` expected HDF5 observation products
  and `79,842/80,254` expected sector FITS targets are present; all `875`
  HDF5 and `412` FITS omissions are declared `edge_warn`, with zero zero-byte,
  unreadable, non-edge-missing, or bad-schema products. S65 orbit 137 is
  complete and orbit 138 completed `15/16` detector cells. Its only missing
  cell is `cam4/ccd1`, where ePSF construction stopped on a CuPy out-of-memory
  error while requesting about `98.6 MiB`. No Stage-1 queue is active, and
  S66-S69 have not completed. Resume only that failed CCD with safer GPU
  concurrency before re-entering the unchanged sector gates.
- `2026-08-06`: Published reproducible S60-S63 target and HDF5 filename
  inventories for collaborator crossmatching. The broad observation-table
  target lists contain `27,165`, `41,403`, `40,158`, and `53,512` TICs; the
  corresponding nonzero-HDF5 filename unions contain `27,036`, `41,168`,
  `39,964`, and `53,249`. Every difference (`129`, `235`, `194`, and `263`)
  is an explicit observation-table `edge_warn`; these filename scans do not
  replace the accepted HDF5-openability/full-schema validation.
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
- `2026-08-26`: A strict read-only PDO/ORCD reconciliation confirmed that the
  operational campaign now reaches S80 without changing the acceptance
  boundary. S65 remains the only accepted sector after the S56--S64 FM release:
  its PDO validation has `228,917` readable nonzero HDF5 products, zero
  non-edge omissions, and `114,419` passing FITS. S66 completed 32 ORCD cells,
  but only its two cam4/ccd4 cells (`3,564` HDF5 total) returned to PDO. The two
  cam4/ccd3 attempt payloads remain complete on ORCD at `1,886` HDF5 and `196`
  ePSF files each; only their cleaned source grids must be re-ingressed before
  retention. The other 28 cells remain retained.
  S67--S78 each have 32 completion and retention receipts explicitly marking
  PDO return deferred and the sector unaccepted. S79 compute and S80 input
  preparation remain active. These checkpoints are not FM temporal evidence
  until the unchanged PDO return, FITS, full-schema, and QA gates pass.

**Next:** Let the live S79/S80 controller remain coherent. After S80 ingress
reaches a stable terminal state, re-ingress only S66's two cleaned cam4/ccd3
source grids, retain their intact payloads, and return/promote all 30
unpromoted S66 cells one at a time before running the unchanged PDO acceptance
gates. Then return and gate S67 onward in sector order. Keep this parallel to
the S63 prospective Teacher-v3 lane; later-sector gathering does not block that
sealed test.

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
  next contract. The user explicitly approved this fresh v3 contract later on
  `2026-07-29`. Its implementation and fail-closed canary/controller assets
  passed two independent final reviews with no remaining P0--P2 findings,
  `726` local tests with `29` optional-environment skips, documentation/static
  checks, and an exact ORCD Torch 2.11 real-VICReg diagnostic. The diagnostic
  proved the locked embedding-only optimizer partition exactly: `146` total
  encoder/projector parameters = `90` active + `56` excluded, with all `90`
  active gradients and AdamW states present, no excluded gradients, finite
  loss, and the deployed `decoupled_weight_decay=true` schema. At that
  checkpoint, no v3 remote artifact or job had yet been claimed.
- `2026-07-29`: Fresh-v3 diagnostic canaries `19218089` (S56 shard 3) and
  `19218090` (S60 shard 4) reached terminal `COMPLETED` without downstream
  submission. S56 produced `1,952` source rows = `1,951` eligible + the one
  locked exclusion, zero stderr, and zero recorded rank warnings. S60
  produced `1,655` eligible rows but emitted exactly `3,324` NumPy
  `RankWarning`s and `3,324` corresponding source lines. Read-only diagnosis
  found that v3 supplied absolute BJD to the ADP spline; the warning-free S56
  result does not make that coordinate contract valid. Both HDF5 files remain
  diagnostic-only and no v3 registry, numerical release, smoke, or fold was
  launched. The corrected fresh-v3r1 contract detrends on compact BTJD,
  publishes exact absolute BJD separately, rejects non-authorized warnings,
  and requires a zero RankWarning ledger through every downstream gate.
- `2026-08-06`: Prepared a
  [Julien-team meeting package](../reports/meeting_julien_team_20260806/README.md)
  from the frozen `509`-observation S56--S62 morphology review, including a
  browsable index of all reviewed rows, the full `80`-observation
  Planet-like/broad-wide atlas, Teacher v1--v3 evidence, and a time-stamped
  Teacher v4-SSL status. The presentation keeps morphology separate from
  confirmation and enrichment separate from automatic classification; it
  introduces no discovery, completeness, occurrence-rate, or final SSL
  transfer-performance claim.
- `2026-08-06`: Completed the corrected full-pool Teacher v4-SSL development
  evaluation without opening S63. On the exact matched development-OOF
  support, Teacher v3 versus fine-tuned SSL balanced accuracy/macro F1 is
  `0.801/0.789` versus `0.775/0.789` under uncertain-as-other and
  `0.822/0.802` versus `0.821/0.808` under uncertain-masked. The latter
  profile's Planet-like average precision is `0.763` for Teacher v3 versus
  `0.614` for SSL, while the frozen SSL linear probe is substantially worse
  in both label policies. The confidence intervals overlap on the aggregate
  metrics, so the result does not establish a reliable improvement. The
  [figures and exact tables](../reports/stage5_validation/teacher_ssl_fullpool_v4_development_performance/)
  are archived as representation-learning evidence; Teacher v4-SSL is not
  promoted, not used for S63, and not a pseudo-label generator. Teacher v3
  remains the frozen operational enrichment model.
- `2026-08-06`: The user authorized the sealed prospective S63 Teacher-v3 run
  and clarified the scientific label language: the growing positive sample is
  **human-confirmed Planet-like transit morphology**, not confirmed
  exoplanets. A read-only PDO/ORCD audit found accepted S63 Stage-1 inputs and
  the hash-verified five-fold Teacher-v3 release, but no S63 cadence, compact,
  BLS, native, score, or review artifact and no active job. The frozen
  `53,512`-TIC reservation has SHA-256 `0560961c...`. The preliminary nonzero-
  HDF5 inventory contains `53,249` represented TICs: `52,487` are disjoint
  from the Teacher-v3 corpus and define the prospective primary cohort, while
  `762` repeated hosts are a separately reported side cohort. Final counts
  must be frozen from the validated compact/model-ready export. The
  [prospective plan](../reports/stage5_validation/teacher_v3_s63_prospective_v1/README.md)
  pins the model/checkpoint hashes, search rule, hidden queue quotas, control
  sampling, metrics, and one-time unblinding policy. The primary run will
  preserve Teacher v3's original current-ADP/native input contract;
  an SSL corrected-input sensitivity would be a separate future experiment.
  Scores will select an annotation-withheld enrichment queue plus a
  preregistered control and will never be treated as labels.
- `2026-08-06`: Built the checksum-bound S63 preprocessing, inference, and
  review-queue lane without opening S63 model scores. The implementation keeps
  the accepted Stage-1 validation source immutable, makes an attested run-local
  receipt copy, requires the compact ADP pair to match its accepted FITS root,
  counts, zero-skip record inventory, and HDF5 target identities, and propagates
  one exact clean producer Git SHA through cadence, model-ready, cohort, BLS,
  candidate, raw-source, native-input, scoring, and queue artifacts. The
  [prospective assets](../reports/stage5_validation/teacher_v3_s63_prospective_v1/)
  now include a pre-score selection policy; its current pre-score errata-bound
  SHA-256 is `51f238c5...5e3e`;
  the original prospective-plan bytes remain unchanged at `c2fcc7c7...ec8f`.
  Its explicit erratum narrows the review boundary honestly: score, bucket,
  control-stratum, and explicit cohort annotations are withheld, but visible
  TIC identity makes cohort membership technically joinable. Cohort-wise
  analysis still waits for all `1,100` accepted morphology decisions. A first
  PDO compact export from the earlier planning commit completed with `53,249`
  exported targets, zero read/filter/duplicate/missing-flux skips, and a
  `7.18`-GB HDF5. It is retained only as a diagnostic data-path smoke; it
  cannot authorize scoring and must be superseded by the exact tested
  implementation commit.
- `2026-08-06`: Final adversarial validation caught two pre-launch contract
  gaps before any S63 score was opened. The first candidate path omitted ten
  Teacher-v3 checkpoint features: own-ephemeris odd/even depths, their delta
  and significance, and trend amplitude in each ADP aperture. The BLS path now
  measures those quantities with the training-time two-aperture implementation,
  the candidate table carries their exact prefixed columns, and launch plus
  scoring fail closed on missing or wholly nonfinite features. The second gap
  allowed a queue consumer to verify hashes without proving the complete
  frozen Teacher-v3 release and candidate evidence. Queue publication now
  requires the exact five checkpoints, release documents, score policy,
  initial/final clean-checkout audits, and full candidate-column equality;
  public and hidden queue products are both mode `0700/0600`, use an exact
  public-column allowlist, and require paired completion markers without
  exposing hidden hashes. Repository validation passed with `837` tests and
  `29` skips; the detection sample and documentation checks also passed.
- `2026-08-06`: The first exact-commit PDO preflight stopped before writing
  because the immutable accepted S63 validation (`619d70aa...a1`) predates the
  validator's embedded `expected_contract` field. The receipt still contains
  independent exact expectation evidence: positive counts for orbits `133`
  and `134`, all `32` orbit/camera/CCD cells, and totals equal to
  `n_expected_h5=107,104`. The compatibility gate now accepts the older schema
  only when that complete key set and both sums agree; missing, nonpositive, or
  inconsistent legacy counts fail. The original accepted receipt remains
  unchanged, and its real PDO JSON passes the corrected read-only preflight.
- `2026-08-06`: Exact commit `af83863e` completed the authoritative PDO
  preprocessing gates without opening S63 model scores. The attested accepted
  Stage-1 receipt has SHA-256 `f351529b...dbc`; the cadence authority contains
  `179,040` rows for orbits `133/134`; and the `7,183,065,848`-byte compact ADP
  pair exported exactly `53,249` targets with zero read, filter, duplicate, or
  missing-flux skips (SHA-256 `aa299a9a...db48`). Model readiness froze the
  same `53,249` TICs and documented `263` reserved edge omissions. The exact
  partition is `52,487` Teacher-v3-disjoint primary TICs plus `762` repeated
  hosts, with disjointness and exhaustive-union checks passed. All artifacts
  identify the positive decision as human-confirmed Planet-like transit
  morphology and explicitly disclaim confirmed-exoplanet status. The sealed
  inputs were transferred to user-owned ORCD storage, recomputed byte-for-byte,
  and atomically published beside a clean detached `af83863e` checkout. Slurm
  test-only admission passed; the locked 16-shard CPU BLS array is job
  `19793979`, with after-success merge job `19793980`.
- `2026-08-07`: All `16` S63 BLS shards and dependent merge completed with
  Slurm exit code `0:0`. The merged summary passed exact coverage of `53,249`
  unique TICs, strict cadence/orbit reconciliation with zero mismatches, both
  ADP apertures, `50,000` periods, and up to `10` peaks. It contains
  `1,063,549` rows (`1,063,390` valid peak rows and `159` explicit
  `too_few_cadences` rows); the merged Parquet SHA-256 is
  `6738a582...1075`. Six shard stderr files contain eight NumPy
  `Mean of empty slice` warnings from sparse odd/even metadata, while the
  attested finite-count ledger and merged `passed=true` gate remain intact.
  Preserve those warnings as diagnostic provenance and require the candidate
  builder's one-row-per-TIC and feature-presence gates before native export.
- `2026-08-07`: Candidate job `19868952` failed closed before publishing any
  candidate or score artifact because the compact-derived model-ready set was
  broader than the successful two-aperture rank-one population. Exact audit
  found `53,160/53,249` candidate-eligible TICs: `80` lack ADP-small rank one,
  `79` lack primary-ADP rank one, and their union is `89`. All `89` belong to
  the Teacher-v3-disjoint primary cohort; no repeated host is affected. The
  candidate contract is corrected pre-score to preserve all `89` in
  model-ready/cohort accounting as explicit BLS-ineligible observations while
  excluding them from native input, scoring, and queue selection. The
  selection policy records that transparent status-only rule with S63 Teacher
  scores still unopened; the immutable prospective-plan bytes remain
  unchanged. Launch provenance now binds the exact upstream preprocessing Git
  object separately from the clean descendant that produces candidate/native
  and score artifacts.
- `2026-08-07`: Corrected CPU candidate job `19869834` completed from clean
  exact commit `46161232` with Slurm exit `0:0`. The published rank-one table
  has exactly `53,160` unique, duplicate-free TICs (SHA-256
  `0a844238...d2c2`); its summary SHA-256 is `4b0f9069...042e`. Independent
  cohort audit confirms all `89` BLS-ineligible TICs are in the frozen primary
  cohort, zero are repeated hosts, and all `762` repeated hosts remain
  candidate-eligible. The candidate summary records `80` missing ADP-small
  rank-one results, `79` missing primary-ADP rank-one results, the exact
  `89`-TIC union, and `teacher_scores_opened=false`. The compact candidate
  artifacts are hash-verified in ORCD and local transfer staging; the next
  authoritative step is the S63 raw-source export on PDO.
- Transparent per-sector BLS exists; the non-periodic dip branch and
  multi-sector aggregation remain unimplemented production gates.

**Next:** Export the exact candidate population's raw sources on PDO, then
build and merge the native-input shards, publish the immutable launch manifest, run
Teacher v3 once, and produce the annotation-withheld enrichment-plus-control
review queue. Do not inspect S63 scores before the launch freeze, schedule
unrelated corpus labeling, or start student pseudo-labeling.

### Human labels and harmonic review

- The strict S56 adjudication audit completed `343/343` decisions across `323`
  unique real sources, with repeat agreement `0.85` and Cohen's kappa `0.814`.
- The final uniform S56--S62 review contains `70` Planet-like morphology
  observations across `65` unique TICs. These are human-confirmed morphology
  decisions, not confirmed planets or an unbiased occurrence sample. Meeting
  the old `50`-TIC count threshold does not by itself authorize student
  pseudo-labeling; grouped performance, calibration, and prospective-
  enrichment gates remain unmet.
- Rare harmonic supervision remains inadequate for promoting or iterating the
  exploratory teacher-v2 design: `P/3` has no supervised example and `3P` has
  only three in the current harmonic table.
- Franklin morphology returns do not supply harmonic truth. The app's default
  `P/resolved` state and any saved morphology click are insufficient evidence
  that a factor was reviewed. Returned factors remain raw audit fields and are
  masked unless an explicit factor-only review or injection truth verifies
  them; all seven folds remain model inputs regardless of that mask.

**Next:** Use the frozen Teacher v3 release to enrich the sealed S63 review
while preserving score-hidden provenance and a random/control slice. Require
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

### Foundation-model validation

- `2026-08-26`: The frozen seed-`560067` FM0.2.1 canary stopped at step
  `2000`; its strict CPU validation and development-only representation
  evaluation passed with zero sealed-test access. The [step-2000 report]
  (../reports/stage5_validation/twirl_fm0_2_s56_s64_step2000_evaluation_v1/README.md)
  records `z_window` effective rank `39.40` and positive paired and
  trained-minus-random separation intervals, while masked Huber remained
  `0.109%` worse than the frozen zero predictor and cross-sector retrieval was
  not demonstrated. The formal gate and development event-retention run remain
  unapplied.
- `2026-08-26`: Consolidated the new-Mac checkout onto the pushed FM handoff
  branch while retaining the dirty rsync state in a recoverable Git stash.
  Added a dedicated CPU-only step-0 evaluation path that pins the exact
  seed-`560067` initialization, step-2000 validation/evaluator receipts,
  unchanged evaluator revision and code hashes, development population, and
  observation-sector authority. It refuses output-directory reuse and maps the
  frozen evaluator's `trained_encoder` field explicitly to the exact
  initialization rather than its separate seed-0 random control. Focused tests
  and the full fast suite passed (`915` passed, `33` optional skips). No ORCD
  job was submitted because the new Mac has no user-opened control socket.
- `2026-08-26`: Reused the user-opened new-Mac ORCD control socket with all
  authentication fallbacks disabled, fetched only the pushed handoff branch,
  and deployed clean detached orchestration revision `393ee8f9`. The exact
  seed-`560067` initialization evaluator ran as CPU job `21340937` on
  `node4702` with four CPUs, `16 GiB`, no GPU, and empty stderr; it completed
  `0:0` in `3m46s`. The step-0 evaluation SHA-256 is
  `bc75b962...1fe6`, and its receipt SHA-256 is `5e489c8a...e29c`. The
  receipt binds the frozen training/evaluator revisions, initialization
  checkpoint, step-2000 receipts, identical `256`-component/`383`-visit
  development population, and observation-sector authority. It records zero
  sealed access, no event-retention run, no formal gate, and no production or
  foundation-model authorization.
- `2026-08-26`: The completed
  [four-point analysis](../reports/stage5_validation/twirl_fm0_2_s56_s64_step0_evaluation_v1/README.md)
  shows a genuine same-seed rank repair: `z_window` effective rank rises from
  `3.03` at initialization to `39.40` at step 2000, while masked Huber falls by
  `80.05%` but remains `0.109%` worse than zero. Cross-sector retrieval does
  not improve descriptively (`0.03380 -> 0.02083`), and the exact
  initialization already exceeds the separate seed-0 control on paired
  separation, so the trained-minus-random gate is seed-sensitive. A read-only
  sector audit reconciled the approximately 20-sector statement: the
  operational footprint is 19 sectors through S74, but only S56--S65 are
  accepted Stage-1 products and only S56--S64 form the frozen FM release.
  Later-sector admission, repeated-host inventory, and a frozen development-
  only temporal panel precede any data-scale training comparison.
- `2026-08-26`: Added a separate fail-closed
  [later-sector admission policy](../configs/models/twirl_fm0_2_later_sector_admission_v1.yaml),
  [inventory builder](../scripts/stage5_validation/build_twirl_fm0_later_sector_inventory.py),
  and reusable admission module. The contract pins the exact S56--S64 input
  manifest, corpus-selection, and Gaia--TIC alias authorities; rebuilds the
  frozen leakage components; rejects any later alias edge that touches or
  merges the baseline graph; and requires checksum-bound accepted-sector
  evidence for provenance, coverage/openability, cadence quality, omissions,
  FM numerical/channel integrity, and source-registry identity. Incremental
  inventories keep new and repeated hosts separate, exclude sealed identities,
  and cannot freeze the panel. Their five-sector and repeated/new-host
  requirements are count floors only; detector/cadence/cohort-adequacy
  thresholds remain deliberately unfrozen until a real label-blind inventory
  exists. The [preflight report]
  (../reports/stage5_validation/twirl_fm0_2_later_sector_admission_v1/README.md)
  records that only S65 is currently later and accepted; no inventory receipt,
  temporal freeze, checkpoint evaluation, training, event retention, formal
  gate, or sealed-test access was performed. The final adversarial admission
  suite passes `100` focused tests; the repository fast suite passes `1,017`
  tests with `33` optional-environment skips, and documentation/static checks
  pass. Independent final review found no remaining P1/P2 contract issue.
- `2026-08-27`: Reconciled the later-sector FM path with the ORCD-resident
  light curves rather than waiting on PDO product return. S65 already has a
  checksum-bound `228,917`-HDF5 sector archive on ORCD, and S66--S78 use the
  immutable 32-cell completion-plus-retention layout. Added a fail-closed
  [ORCD source verifier](../scripts/stage5_validation/build_twirl_fm0_orcd_source_receipt.py)
  that binds either source form while keeping `FM_ORCD_SOURCE_READY` distinct
  from HDF5-openability, six-view, temporal-panel, A2v1-acceptance, and
  occurrence-rate claims. S79 is actively computing on two H200s and S80 is
  dependency-queued; the FM preparation lane remains CPU-only and separate.
- `2026-08-27`: Published real immutable ORCD source-readiness receipts for
  every sector S65--S78. Each retained S66--S78 sector binds all `32` cells;
  the newly verified S70--S78 sectors declare respectively `36,963`, `38,685`,
  `51,665`, `49,220`, `81,173`, `65,791`, `59,012`, `57,992`, and `63,394`
  HDF5 products without claiming openability or six-view admission. A complete
  PDO cadence audit found exact QLP-qflag coverage for every S66--S78 detector.
  Legacy SPOC files are exact through S74, but S75 and S76 omit `29,446` and
  `20,770` detector-cadence rows, demonstrating that SPOC cannot be forced
  beyond its intended boundary. Inspection of the activated QLP 0.14.6 source
  confirms `auto` selects SPOC for S<67 and TICA for S>=67. Added a reusable
  fail-closed [mission-quality gate](../scripts/stage5_validation/build_twirl_fm0_mission_quality_receipt.py)
  that enforces that provider boundary plus exact quaternion and QLP-qflag
  coverage while keeping HDF5, shard, and panel claims false.
- `2026-08-27`: Materialized current-QLP TICA quality files in immutable
  user-owned releases and replaced the initial exact-only quality receipt with
  a provenance-bound v2 contract. S66--S76 pass with zero quaternion-authority
  exclusions. S77 passes after recording exactly `32` mission-only boundary
  rows (cadences `936740` and `936741` in all `16` detectors) as explicit
  exclusions; retained HDF5 is not mutated. S65 remains a negative control and
  publishes no receipt because cam4/CCD4 still lacks `1,462` SPOC rows. S78
  remains blocked because its TICA database/delivery contains `6,405` rows per
  detector and misses trailing quaternion cadences; a sampled orbit-164 HDF5
  ends immediately before the first missing cadence, but detector-wide omission
  proof is still required.
- `2026-08-27`: Built and transferred one checksum-bound compact S66--S77
  mission-quality release to ORCD. It contains `12` receipts and `694` source
  files, preserves the QLP orbit/run and mission layouts, and binds
  `226,416,794` logical bytes. The archive SHA-256 is
  `f3d8b532...beda514`; all `706` manifest-listed receipts/sources re-hashed
  successfully after extraction.
- `2026-08-27`: Added a read-only retained-HDF5/cadence-quality admission gate
  and launched a strict CPU-only S67--S77 `afterok` chain. S67 job `21406846`
  completed in `12m59s` with exit `0:0`: all `215,405` declared HDF5 opened,
  zero were unreadable, and `1,262,165,389` cadences reconciled against the
  transferred quaternion, QLP-qflag, and TICA authorities. The immutable
  receipt SHA-256 is `71afae2c...52e4a3`; it reports zero unexplained or
  authority-excluded S67 cadences and keeps checksum-provenance, six-view,
  identity, A2v1 acceptance, and temporal-panel claims closed. S68 job
  `21407076` then completed in `6m38s` with exit `0:0`: all `108,991` declared
  HDF5 opened, zero were unreadable, and `634,273,102` cadences reconciled with
  zero authority exclusions. Its receipt SHA-256 is `6521f1aa...b17abf7`.
  S69--S73 then completed consecutively with exit `0:0`. Across those five
  sectors, all `256,326` declared HDF5 opened, zero were unreadable, and
  `1,374,212,726` cadences reconciled with zero authority exclusions. Their
  receipt SHA-256 prefixes are respectively `b612ad07`, `e761a8e9`,
  `cf1088f3`, `94edd224`, and `b64cd095`. S74--S76 then passed with all
  `205,976` declared HDF5 opened, zero unreadable files, `1,179,099,318`
  reconciled cadences, and zero authority exclusions. Their receipt SHA-256
  prefixes are `7023d6b3`, `f1ae5b15`, and `128548fe`. Final S77 job
  `21407087` completed in `3m41s` with exit `0:0`: all `57,992` HDF5 opened,
  zero were unreadable, and `223,416,351` cadences reconciled. Its `32`
  predeclared mission-only detector rows recur across `57,994` light-curve
  cadence occurrences and are explicitly authority-excluded; no unexplained
  gap remains. The S77 receipt SHA-256 is `2d3f6cfa...b1f8dc` and binds
  exclusion-ledger SHA-256 `850722d9...f7811`.
- `2026-08-27`: The complete S67--S77 retained-HDF5/cadence chain passed. All
  `844,690` declared products opened, zero were unreadable, and
  `4,673,166,886` light-curve cadences reconciled against quaternion, QLP
  detector-qflag, and TICA mission-quality authorities. The only authority
  exclusions are S77's `57,994` expected light-curve occurrences of the `32`
  declared mission-only detector rows. These receipts close HDF5 openability
  and cadence-quality only; six-view, identity, product-acceptance, temporal-
  panel, training, and model-result gates remain closed.
- `2026-08-27`: Froze a model-outcome-independent
  [later-sector exclusion ledger](../configs/models/twirl_fm0_later_sector_exclusions_v1.yaml).
  S65 is excluded from the first later-sector zero-shot panel and first
  post-FM0.2 training release because the camera-4/CCD-4 SPOC authority lacks
  quaternion cadences `802061--803522` (`1,462` rows). This does not demote its
  accepted A2v1 product state. A future repair may enter only a newly versioned
  release after every FM gate passes. The frozen admission-policy v1 remains
  immutable and cannot express this skip, so a separately hashed v2 must bind
  the ledger before an inventory or panel is frozen. Reserved `TWIRL-FM0.3`
  for any justified post-zero-shot campaign: `.3.1` is provisionally a
  fixed-TCN, fixed-view, fixed-objective data-scale baseline, while dual-mask
  reconstruction and a stable-host cross-visit head remain separately
  controlled later ablations. No training or temporal-panel freeze was
  authorized.
- `2026-08-28`: Corrected CPU-only S66 HDF5/cadence audit job `21498338`
  completed on `node4702` in `1h32m11s` with exit `0:0`. All `511,773`
  declared HDF5 products opened, zero were unreadable, and `2,872,070,553`
  light-curve cadences reconciled against the SPOC plus QLP authorities with
  zero authority exclusions. The immutable receipt SHA-256 is
  `c1c4bb61...1c5e985`. The complete S66--S77 chain therefore covers
  `1,356,463` opened products and `7,545,237,439` reconciled cadences with no
  unreadable file or unexplained authority gap; only S77's predeclared
  `57,994` mission-only light-curve cadence occurrences are explicitly
  excluded. S65 remains outside the release.
- `2026-08-28`: Recorded Jeroen's short-event-retention concern as a matched
  context-window experiment rather than a photometry rebuild. FM0.2.1 uses a
  `2,048`-cadence (`4.74 d`) input crop, a theoretical `1,543`-cadence
  (`3.57 d`) TCN receptive field, and masked-mean window pooling. The next data
  products therefore retain complete visits with no model context bound. A
  later development-only loader diagnostic will compare contiguous `256`,
  `512`, `1,024`, and `2,048` cadence-slot crops from the same shards while
  holding physical-source sampling and total model-visible cadence exposure
  fixed. `TWIRL-FM0.3.1` remains the fixed-context data-scale baseline; a
  justified shorter-context-only candidate is provisionally
  `TWIRL-FM0.3.2`. Neither training nor a temporal-panel freeze was authorized.
- `2026-08-28`: Implemented the fail-closed S66--S77 full-visit preparation
  path. It builds provider-neutral cadence references using SPOC for S66 and
  TICA for S67+, creates label-blind source inventories without opening HDF5
  photometry, writes immutable six-view shards only for train/development
  identities, and keeps the complete visits independent of model context
  length. A checksum-bound [authority map]
  (../configs/orcd/twirl_fm0_s66_s77_full_visit_preparation_v1.json) with
  SHA-256 `10adea2e...fa7d7c` fixes every common and per-sector input including
  the completed S66 receipt. The dry-run-first ORCD controller launches two
  chronological CPU Phase-A chains, freezes their exact generated hashes,
  then permits one chronological six-view chain and one full CPU admission-v2
  pass. All long-output builders now publish read-only, clean only their exact
  attempt-owned final/partial paths after interruption, and receive a
  ten-minute pre-termination signal. The aggregate sealed-data contract states
  precisely that pre-existing label-blind QA read only identity attributes,
  Cadence, and QualityFlag for all identities; no sealed time, aperture-flux,
  flux-error, or derived shard enters preparation. The complete fast suite
  passed (`1,109` tests; `33` optional skips), focused independent launch audit
  found no P0/P1 blocker, and no H200 or training path exists in the
  controller.
- `2026-09-01`: Repaired the final S66--S77 admission provenance handoff by
  separating the immutable Phase-A producer revision from the current
  admission-validator revision and binding the recorded Phase-A SHA-256.
  CPU retry job `21749958` completed `0:0` in `17m15s`. The read-only
  admission receipt SHA-256 is `9461e2f3...740701`; it closes exactly
  `673,404` source rows to `606,387` nonsealed six-view shards while retaining
  `67,017` sealed identities as excluded and keeping S65 outside the pool.
  Added a separate identity-only temporal-panel freezer, a paired CPU
  step-0/step-2000 evaluator, and evaluation-only short-context model support.
  These paths do not authorize training, event retention, a formal gate, or
  sealed-test access.
- `2026-09-01`: CPU job `21751614` froze the S66--S77 identity-only
  development panel. Its receipt SHA-256 is
  `78c370e10c556472c5997c20cfe95207a0b334bafe7f024bf7ba4fc7ec4de624`.
  The panel contains `67,403` visits in `52,344` leakage components: `16,409`
  repeated-host visits in `7,030` components and `50,994` genuinely new-host
  visits in `45,314` components. It contains zero sealed rows and remains
  ineligible for training.
- `2026-09-01`: CPU job `21752135` completed the exact paired S66--S77
  temporal zero-shot evaluation; the immutable output SHA-256 is
  `d7cb2d9d81f53fbfd2f599841339b8e7734243e5c3a8c7acd3e379510d6a20b7`.
  Repeated-host `z_window` effective rank increased
  `2.8898 -> 32.5069`, with paired-separation delta `+0.006781` and
  source-clustered 95% interval `[0.000417, 0.013358]`, while retrieval changed
  `2.9106% -> 1.6632%`. New-host rank increased `2.1358 -> 18.1707`, but
  separation changed `-0.006670`, interval `[-0.013505, -0.000181]`, and
  retrieval remained `5% -> 5%`. Step-2000/zero reconstruction ratios were
  `1.00303` repeated and `1.00249` new. Rank repair transfers, but the later-
  sector evidence is mixed and does not promote the model.
- `2026-09-01`: CPU job `21752626` completed the fixed-exposure context
  diagnostic; the immutable output SHA-256 is
  `5376cdf2e9eea46cb6d1804f3ba9ec220278d352175533c574402570df78c708`.
  Results were nearly flat across `256/512/1,024/2,048`, but the test cannot answer
  event dilution because it averaged crop embeddings over the same fixed
  `2,048`-cadence exposure. It validates evaluation-only short-context
  mechanics, not a shorter-context choice.
- `2026-09-01`: Clarified FM0's target as a BLS-free general light-curve
  representation backbone for later supervised classification and triage,
  with localized event morphology retained. Transparent periodic BLS remains
  separate. Froze the next bounded diagnostic as one centered trapezoid event
  at direct `128/256/512/2,048` contexts, with no period, repeated synthetic
  transits, cross-crop averaging, or BLS input. It compares training-free
  representation and visible-event reconstruction response at exact step 0
  and step 2,000. A separately frozen classifier remains the eventual transfer
  gate.
- `2026-09-01`: CPU job `21759626` completed the direct centered single-event
  context diagnostic in `30m02s` with exit `0:0`; the immutable output SHA-256
  is `9ad6d817ef5c3e00b44c62e3a110940964256a780b84a92c7b446a7001ba3d1c`.
  It used one eligible development visit for each of `48` repeated and `48`
  new leakage components, applied all `12` one-event duration/depth conditions
  per source, and opened no labels, candidates, BLS/period features, sealed
  shards, optimizer, or training path. At step 2,000, paired `128/2,048`
  `h_window` robust-scaled response ratios are `16.324` (source-clustered 95%
  interval `[14.954, 17.702]`) for repeated hosts and `16.299`
  (`[14.728, 17.919]`) for new hosts. The corresponding pooling- and
  signal-corrected ratios are only `1.081` (`[1.040, 1.130]`) and `1.034`
  (`[0.993, 1.077]`), showing that the large raw gain is predominantly the
  expected masked-mean dilution effect. Training also suppresses short-event
  response: at `128` cadences, step-2,000 `h_window` robust-scaled response is
  about `0.50x` step 0 for repeated and `0.59x` for new hosts, while
  `z_window` is only about `0.044x` and `0.040x`. The next transfer test should
  therefore use the encoder-side `h_window`, not the learned projection
  `z_window`, with `128` cadences as the leading local-event context.
- `2026-09-01`: CPU mechanics job `21765394` completed `0:0` in `16 s` and
  verified both provisional FM0.3 candidates at exact `128/128` cadence-token
  geometry. TCN and Conformer both use stride one, identity patch/time paths,
  no temporal downsampling, and finite cadence-indexed outputs; the immutable
  result SHA-256 is
  `4ed8bd35ceb586f75194bd3e193e2f5d1c7ed738a67a26bf0779a2289895f464`.
  This is a training-free mechanics pass, not an architecture selection.
- `2026-09-01`: CPU transfer job `21765395` completed `0:0` in `11m06s`; the
  immutable result SHA-256 is
  `57762a5e94ad03d4d1244ae8b50af1d3baf8a3ba4e693da715c6967e517a1798`.
  Its primary step-2,000 `128`-cadence AUROC was `0.499783`, but the result is
  not valid model evidence: the raw control was exactly `0.5` in every one of
  the `12` duration/depth cells, including `30%` dips. The hard-max training
  operator passed gradient only through the current argmax, so paired
  clean/injected gradients cancelled whenever the event cadence did not
  already win. A separately versioned v2 repair retains exact cadence tokens
  and hard-max evaluation but uses event support only as training truth for
  four pair-balanced per-cadence BCE strata. It adds fail-closed raw-flux and
  quality-only mechanics gates; no FM training, BLS, GPU, or sealed access is
  authorized.
- `2026-09-01`: Hardened the provisional FM0.3 implementation around the
  no-cadence-merging rule. Both variants now fail closed unless the context is
  exactly `128` native `200 s` cadences at stride one; the Conformer exposes
  its true post-context token sequence while keeping its stem skip confined to
  reconstruction. Added the isolated
  `position_centered_cadence_vicreg_loss`, which preserves `[batch, cadence,
  feature]` geometry, centers only across windows at each fixed cadence, and
  pools only scalar loss sufficient statistics. The Torch-enabled objective
  tests passed `7/7` on ORCD; trainer replay, the cadence-collapse gate, and
  the training-eligible manifest remain incomplete, so this does not authorize
  a real FM0.3 run.
- `2026-09-01`: Event-transfer v2 CPU job `21766861` completed `0:0` in
  `1h21m56s`; its immutable result SHA-256 is
  `42e9625c8bf1261b24379426fb732bc7e4289b76bde832bad0847916186930f3`.
  The repaired pair-balanced cadence fit behaved sensibly: the loss fell from
  about `0.6932` to `0.5597`, both standardized ADP coefficients were negative,
  and the quality-only paired control was exactly `0.5`. The all-window
  reduction nevertheless invalidated the mechanics gate. Raw AUROC was only
  `0.50225`, paired ranking was `0.51667`, and the three `30%`-depth cells were
  `0.50125`, `0.50906`, and `0.52125` for durations `1`, `3`, and `9` cadences.
  The per-cadence fit learned the dip direction, but the maximum over all `128`
  cadences usually remained on an unrelated real outlier shared by the clean
  and injected pair. Every FM metric in this artifact is therefore
  non-interpretable; it is not evidence for or against the TCN, the shorter
  context, or training.
- `2026-09-01`: Implemented the bounded FM0.3 pre-run path without launching a
  model. The trainer now binds the explicit cadence-objective schema, optimizes
  first-mask reconstruction plus position-centered cadence VICReg on the true
  post-context `[batch, cadence, feature]` stream, and replays exact staged
  gradients without creating a temporal representation. The collapse runner
  checks native cadence differences from exact step-0/step-2,000 checkpoints,
  but is explicitly a non-authoritative preflight until a panel receipt and
  real-run contract exist. A compact identity-only composite freezer binds the
  existing S56--S64 and S66--S77 manifests without copying shards; all S77
  `poc_train` rows become the internal temporal holdout and every overlapping
  earlier component is quarantined. The pre-freeze audit found `789,588`
  candidate earlier training rows before that component quarantine and
  `23,162` S77 holdout rows containing `178,463,210` cadences; the final counts
  await the immutable remote receipt. A separate raw/quality-only v3 mechanics
  control now uses exact paired event-center cadences and cannot load an FM.
  The full local suite passed `1,213` tests with `43` optional skips. No shard
  was rewritten, no sealed payload was opened, and no FM training was enabled.
- `2026-09-01`: V3 mechanics job `21772041` completed `0:0` in `43 s`; its
  immutable result SHA-256 is
  `10b78bf083595ad8f4ecacb0428d99b939abc9dae5e3f361d02559659d1e2833`.
  The raw center-cadence paired ranking was `0.98542` with a deterministic
  `95%` interval of `0.97495--0.99583`, both ADP flux coefficients were
  negative, and the quality-only control was exactly `0.5`. The evaluator
  emitted all `128` native cadence logits and used no off-support values,
  cadence merging, averaging, pooling, or downsampling. This validates the
  probe mechanics only; it does not score an FM or select an architecture.
- `2026-09-01`: Composite freeze job `21772119` published its first pass with
  receipt SHA-256
  `cc5cc3bce4c24e74bef1fbf084f407855233de7893183e6bc3486e284a2f44d9`,
  `723,461` provisional training rows, `66,232` quarantined-overlap rows, and
  `23,162` S77 holdout rows. This corrects the earlier pre-freeze count upward
  by `105` rows to `789,693` before quarantine. The ORCD control socket expired
  while the intentional second full validation was still running, so neither
  these counts nor the receipt are accepted until terminal state and immutable
  closure are re-verified.
- `2026-09-01`: Implemented the actual fixed-band v4 transfer screen. It emits
  one scalar score for each of `128` native cadence tokens, then takes a maximum
  over the same literal indices `36--92` for every held-out sample. Synthetic
  support is restricted to the training target and is absent from the held-out
  scorer/reducer. The four exact arms are raw, quality-only, step-0
  `h_cadence_128`, and step-2,000 `h_cadence_128`; the old `2,048`-cadence and
  projected-`z` arms are excluded. Focused v1--v4 tests passed `33/33`, and the
  full local suite passed `1,220` tests with `43` optional skips. This synthetic
  candidate-centered screen explicitly cannot choose TCN versus Conformer.
- `2026-09-02`: Reverified composite freeze job `21772119` as `COMPLETED 0:0`
  after `15m59s`, with two identical validation passes and an empty error log.
  The immutable four-file identity release has receipt SHA-256
  `cc5cc3bce4c24e74bef1fbf084f407855233de7893183e6bc3486e284a2f44d9`,
  source-bindings SHA-256
  `8cbcac99409ab89fe2dee0c36687f50122e78376206495073698489ca5424f2b`,
  and role-index SHA-256
  `abe9c616523f2486bf1b7be69dfcfda6193d534b40b3f4afc49d8ebf3e40ce5a`.
  It closes on `723,461` training, `66,232` quarantined-overlap, and `23,162`
  S77 holdout identities; the two-ADP-view adapter exposes `676,359` training
  and `20,617` holdout observations. No shard was copied or rewritten.
- `2026-09-02`: Fixed-band job `21818153` completed `0:0` in `2m40s` with an
  empty error log. Its immutable result SHA-256 is
  `c4f770f753c36b35fe42ce4c2a364a1a29d92542c6bedfae998b809e2e8f687d`.
  Cadence geometry was correct, but the raw positive control failed: AUROC was
  `0.50392` and paired ranking `0.53125`; even the `30%`-depth duration cells
  reached only `0.50531`, `0.51594`, and `0.52781`. The step-0 and step-2,000
  AUROCs were `0.51522` and `0.50136`, respectively, but those FM values are
  non-interpretable because the broad maximum usually selected an unrelated
  natural outlier. Quality-only stayed exactly at chance. V4 therefore rejects
  only the fixed `36--92` maximum, not the short context or either future
  architecture.
- `2026-09-02`: Retired a proposed fixed-center v5 on the obsolete FM0.2
  checkpoint before implementation. It would have repaired the readout but
  could not test the new cadence-level objective. The next experiment is the
  matched one-seed FM0.3 canary itself. Its primary downstream task is
  predeclared candidate-centered classification at fixed index `64`, with all
  `128` native `200 s` tokens still emitted and no BLS, pooling, cadence merge,
  or truth-routed scoring. The existing panel supplies only S66--S71 probe
  fitting and S72--S74 validation; a new identity-only S77 test is being frozen
  from `20,617` ADP-usable holdout observations (`15,137` repeated and `5,480`
  new), with zero training-component overlap.
- `2026-09-02`: Completed the local matched-canary implementation without
  launching training. FM0.3.1 and FM0.3.2 now share exact native-`128`
  cadence geometry, stride one, the cadence-level VICReg objective, and a
  contextual reconstruction input; the legacy stride-four Conformer stem
  skip is disabled for this matched comparison. The two-stage evaluation
  freezer first selects identities without payload access, then opens exactly
  those `1,440` bound shards once and freezes deterministic, unpadded crops
  with at least `103/128` joint-valid cadences and fully valid indices
  `60--68`. Its receipt fixes the index-`64` injection/readout, full-batch
  linear probe, raw/quality controls, and paired bootstrap before any real
  checkpoint. Real training now requires and revalidates that immutable
  receipt and the same variant's exact step-`8` FP32 smoke; the composite
  sampler rejects every segment shorter than `128` rather than padding it, and
  step `2,000` must exactly resume the same run's step-`64` checkpoint. The
  controls-only evaluator can run before any checkpoint and freezes its own
  checksum-bound result. The full local suite passed `1,275` tests with `43`
  optional skips; documentation, focused lint, shell syntax, and diff checks
  passed.

**Next:** Publish the exact code revision to ORCD, freeze the identity-only and
payload-screened evaluation bundles, and run the matched eight-step synthetic
smokes. Report readiness before submitting either real one-H200 canary;
preserve one native `200 s` cadence per token throughout.

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
