# S56 A2v1 active-pair Tier-1 assessment

**Contract:** `a2v1_tier1_science_qa_v4`

**Scope:** `DET_FLUX_ADP_SML` + `DET_FLUX_ADP` bounded enrichment
**Assessment date:** 2026-07-23

## Abstract

We audited all 31,450 Sector 56 A2v1 targets using the exact preprocessing
shared with the locked BLS search, the authoritative internal-plus-external
quality mask, 2,000 fixed injections, and an independent TESSCut recovery of
WD 1856 b. Diagnostic cadence review remains visible, but flagged-cadence
fraction is not a veto. The gate returns `status=review`, `passed=true`, and
`enrichment_ready=true`; it remains `science_ready=false` because this bounded
test covers only the two active ADP channels.

The paired teacher retains 31,446/31,450 targets (99.987%). Of the four
paired-input exceptions, three retain a healthy primary aperture and remain
available for single-aperture or manual review. Only one target lacks a usable
series in both apertures. None is converted into a negative label.

## Methods

The v4 eligibility predicate calls the same BLS preparation used in production:

1. apply the exact internal-or-authoritative-external quality mask;
2. require finite time and raw flux;
3. retain at least 200 cadences per required aperture;
4. require finite, nonzero normalization;
5. apply the locked conditional upper-tail clip; and
6. require a nondegenerate post-cleaning period baseline.

Population scatter, flagged fraction, and aperture disagreement remain QA and
sensitivity metadata. They do not reject a target or detector by themselves.
The exact compact product, cadence authority, Tier-0/BLS evidence, injection
shards, and independent extraction are checksum-bound.

## Results

| Quantity | Result |
|---|---:|
| Full S56 population | 31,450 |
| Paired-teacher eligible | 31,446 (99.987%) |
| Searchable in at least one aperture | 31,449 (99.997%) |
| Searchable in neither aperture | 1 (0.003%) |
| Diagnostic target QA | 30,178 pass; 1,240 review; 32 fail |
| Fixed injections | 2,000 hosts on 16 detectors |
| Median depth retention | 0.9908 small; 0.9906 primary |
| 10th-percentile retention | 0.9666 small; 0.9661 primary |
| Quality-aware BLS binding | 628,955 rows; pass |

All prerequisite, population-scatter, aperture, injection, and independent
extraction gates pass. The cadence/finite-data gate remains a nonblocking
review because cam1/ccd1 has an unusually large flagged interval. Its targets
still retain thousands of search cadences, so excluding the detector would
remove valid information without improving the current teacher corpus.

### Paired-input exceptions

| TIC | Detector | Searchable aperture(s) | Action |
|---:|---|---|---|
| 1201200937 | cam4/ccd4 | ADP primary | Retain for single-aperture/manual review; omit from paired teacher |
| 1201317288 | cam4/ccd3 | ADP primary | Retain for single-aperture/manual review; omit from paired teacher |
| 1973584484 | cam3/ccd1 | ADP primary | Retain for single-aperture/manual review; omit from paired teacher |
| 2019898202 | cam3/ccd2 | none | Omit as unusable input; never label as a negative |

## Figures

![Full-population diagnostics](tier1_qa_diagnostics.png)

*Figure 1.* Scatter-versus-magnitude, aperture consistency, and fixed-injection
depth retention. The pass/review/fail annotation is diagnostic QA, not the
paired-teacher inclusion decision. [PDF](tier1_qa_diagnostics.pdf)

![Paired-input eligibility by detector](tier1_detector_eligibility.png)

*Figure 2.* Paired-input eligible/total counts by camera and CCD. All detector
cells retain essentially their full populations; the four exclusions are
missing-channel products, not high-flag-fraction vetoes.
[PDF](tier1_detector_eligibility.pdf)

## Verdict and next step

S56 is ready for the bounded enrichment and paired-teacher workflow. Tier-1
should now stop consuming project attention unless a new sector exposes a
structural failure. This result does not make the survey pipeline science
ready and does not replace later end-to-end completeness, dip-search,
multi-sector, or false-alarm calibration.

The immediate path is to ingest Franklin's final labels, re-review every
Planet-like and EB example across S56--S62, freeze the observation-level
corpus and TIC-grouped split, then retrain the existing seven-harmonic teacher
architecture once.

## Reproducibility

ORCD job `18657818` completed in 9 min 17 s from commit `2633e8b4`. All 14
published output hashes match the local copies. The Tier-1 gate and complete
BLS evidence validators both pass. Local validation completed with 420 tests
passed and 3 skipped, plus the detection smoke and documentation checks.
