# S56 A2v1 active-search-pair Tier-1 quality assessment

**Scope:** bounded enrichment QA for `DET_FLUX_ADP_SML` and
`DET_FLUX_ADP`

**Contract:** `a2v1_tier1_science_qa_v1`

**Assessment date:** 2026-07-23

**Promotion enabled:** no

**Science ready:** no, by contract

## Abstract

We assess whether the Sector 56 A2v1 active search pair is sufficiently
controlled for target-filtered candidate enrichment. The audit covers the full
compact population of 31,450 targets, a hash-bound 188,396-row external
cadence/quality reference with four declared authority exclusions, a
2,000-injection fixed canary, aperture behavior, and an independent
official-TESSCut recovery of WD 1856 b. The final bounded-gate outcome is
**[PENDING ORCD TIER-1: overall status, pass flag, and `enrichment_ready`]**.
Regardless of that result, this two-aperture assessment cannot promote the
six-channel A2v1 product and must retain `science_ready=false`.

## Methods

### Data and scope

The evaluated compact product contains 31,450 S56 targets. The search channel
is the small ADP aperture (`DET_FLUX_ADP_SML`); the primary ADP aperture
(`DET_FLUX_ADP`) supplies cross-aperture contamination and consistency
evidence. The audit deliberately excludes the ADP015 branches and the 5x5
aperture, so it is an enrichment gate rather than release QA.

The external-quality reference combines QLP camera-quaternion cadence
authority with SPOC quality flags. Its 188,396 rows cover the S56 detector
timeline after declaring exactly four missing-authority cases: cadence 699957
for camera 3, CCDs 1–4. Those four rows are masked with the reserved external
quality bit and remain visible in provenance.

### Gate design

Eight fail-closed gates test: cadence-reference provenance, compact-to-injection
source parity, the quality-aware Tier-0 prerequisite, scatter versus TESS
magnitude, cadence retention and finite flux, aperture outliers, injected-depth
preservation, and independent extraction. The fixed canary comprises four
reviewed shards and 2,000 unique injected hosts. Every evidence input is
checksum-bound before the population scan and rechecked before output
publication.

The independent benchmark uses official MAST TESSCut data rather than another
TGLC production tree. WD 1856 b is recovered in its 1x1 reference aperture at
period 1.407960330 d with BLS S/N 11.87 and an approximately 1.27 min epoch
residual. The WCS-defined 2x2 reference aperture gives period 1.407896973 d,
BLS S/N 11.33, and an approximately 0.23 min epoch residual. These
measurements establish the independent evidence input; the final gate status
still depends on the locked end-to-end audit.

## Results

The PDO quality-aware BLS and replacement Tier-0 products passed review and
are checksum-pinned before the ORCD population audit.

| Gate | Status | Principal result |
|---|---:|---|
| Cadence-reference prerequisite | **[PENDING ORCD TIER-1]** | 188,396 rows; 4 declared authority exclusions |
| Injection-source parity prerequisite | **[PENDING ORCD TIER-1]** | **[PENDING ORCD TIER-1: comparison count and mismatch count]** |
| Tier-0 prerequisite | **[PENDING ORCD TIER-1 HASH CHECK]** | Tier-0 passed all 6 nested gates; summary `1f7865b9…432ca`, BLS table `c3a7bd9f…10692` |
| Population scatter | **[PENDING ORCD TIER-1]** | **[PENDING ORCD TIER-1: supported magnitude bins, slopes, and outlier fractions]** |
| Cadence and finite data | **[PENDING ORCD TIER-1]** | **[PENDING ORCD TIER-1: loss, finite-flux, and usable-cadence quantiles]** |
| Aperture outliers | **[PENDING ORCD TIER-1]** | **[PENDING ORCD TIER-1: valid, ratio-outlier, correlation, and anticorrelation fractions]** |
| Fixed-injection preservation | **[PENDING ORCD TIER-1]** | 2,000 unique injection IDs; **[PENDING: retention quantiles and in-band fraction]** |
| Independent extraction | **[PENDING ORCD TIER-1]** | WD 1856 b recovered in both independent apertures; **[PENDING: common-cadence fraction and final gate status]** |

**Overall status:** **[PENDING ORCD TIER-1: `status` and `passed`]**

**Target eligibility:** **[PENDING ORCD TIER-1: pass/review/fail counts and
Gaia-identified count]**

**Enrichment ready:** **[PENDING ORCD TIER-1: `enrichment_ready`]**

**Science ready:** `false` (fixed by scope, not a pending measurement)

## Figures

> **Figure 1 slot — `tier1_qa_diagnostics.png`.** Full-population diagnostic
> panels for MAD and five-point RMS versus TESS magnitude, primary-to-small
> aperture MAD ratio versus aperture correlation with magnitude color coding,
> and injected-depth retention. Pass and review boundaries are the locked
> configuration thresholds. **[PENDING ORCD TIER-1: insert rendered figure and
> summarize the dominant visible trend.]**

> **Figure 2 slot — `tier1_detector_eligibility.png`.** Camera-by-CCD map of
> the fraction of S56 targets passing the final Tier-1 target filter.
> **[PENDING ORCD TIER-1: insert rendered figure and identify any
> detector-localized deficit.]**

## Limitations

This assessment covers two ADP apertures in one sector. It does not audit all
six ADP/ADP015 channels, establish a survey denominator, calibrate periodic or
dip-search false alarms, merge multi-sector detections, or measure full
search-to-vetting completeness. WD 1856 b is a necessary benchmark but one
favorable system cannot establish survey-wide sensitivity. Passing this gate
would therefore authorize only the use of target-filtered S56 inputs for
bounded enrichment; a review or failure would require remediation before that
use.

## Repository audit and reproducibility

The accompanying repository scan found 2,181 tracked files occupying
approximately 802 MB, dominated by `reports/` (630 MB) and `outputs/`
(167 MB); `.git` occupies approximately 14 GiB. Fifteen tracked files are at
least 10 MB, and 80 already-tracked files now match ignore rules
(approximately 47 MB). No secret-scan, broken-symlink, cache, syntax, or
untracked-file hygiene blocker was identified. The current validation suite
passes 405 tests with 3 skips, together with the detection sample and
documentation checks.

These figures argue for a later, explicit artifact-retention migration rather
than history rewriting during the QA run. In particular, authoritative label
artifacts must be separated from regenerable tensors and presentation
scratch products before large files are untracked. Reproducibility of this
assessment rests on the locked configuration, exact compact checksum, cadence
table and manifest checksums, injection shard and source-parity checksums,
independent-extraction checksums, the Tier-0 summary (`1f7865b9…432ca`) and
BLS table (`c3a7bd9f…10692`) checksums, and **[PENDING DEPLOYED PRODUCER
COMMIT: code revision]**.

## Verdict

**[PENDING ORCD TIER-1: state pass/review/fail and whether S56 may enter the
target-filtered enrichment workflow.]** The non-negotiable interpretation is
that `science_ready=false`: even a clean pass is evidence for bounded
enrichment, not a science-ready search pipeline or survey release.

The teacher-model bottleneck is presently data and evaluation infrastructure,
not a demonstrated lack of model capacity. Teacher v1 should remain the
seven-harmonic enrichment baseline. The next retrain should use the frozen
S56–S62 observation corpus only after quality-aware native inputs, immutable
TIC-grouped splits, source-separated evaluation, one out-of-fold calibration
transform, and complete checksum provenance exist. Teacher v2, student
pseudo-labeling, and architecture sweeps should remain off the critical path.

## Adjusted next steps

1. **Close this bounded gate.** Review and pin the PDO Tier-0/BLS outputs,
   execute the exact ORCD Tier-1 audit, inspect both figures, and publish the
   target pass/review/fail mask.
2. **Freeze the seven-sector enrichment corpus.** Finish S56 and S60–S62 human
   review, preserve observation-level `(sector, TIC, candidate_key)` labels,
   and publish the sibling candidate table/TIC roll-up without implying
   multi-sector confirmation.
3. **Build the training contract before retraining.** Regenerate per-sector
   quality-aware BLS/native inputs, route ephemeris-incompatible labels to
   re-review, freeze one TIC split registry, and report calibration plus
   TIC-bootstrap uncertainty and the uncertain-label sensitivity test.
4. **Use S63 only if it remains genuinely sealed.** Audit git and non-git
   exposure first; freeze teacher v1, thresholds, cohort, and metrics before
   blind S63 labeling; evaluate TIC-disjoint hosts as primary and repeated
   hosts separately; unblind once and do not tune on the result.
5. **Add broader search capabilities after the periodic path is robust.**
   Implement the dip branch, multi-sector merging, and branch-aware
   false-alarm calibration before survey-wide enrichment or science claims.
