# Teacher v3 prospective S63 enrichment plan

This directory freezes the analysis choices for the first prospective use of
the S56--S62 Teacher v3 release. The run is an enrichment experiment over
accepted S63 A2v1 light curves. It is not a blind estimate of sector-wide
classification performance, a confirmed-exoplanet search product, or an
automatic-labeling exercise.

## Scientific language

The positive human decision is **Planet-like transit morphology**. A completed
human review confirms that morphology decision only; it does not confirm that
the source is an exoplanet. Teacher probabilities are queue-selection metadata
and never become labels.

## Frozen design

- Use the existing five-fold Teacher v3 `shape_plus_periodogram_bls` ensemble,
  selected before its fixed test was opened, with its pooled development-OOF
  temperature calibration. Every checkpoint and controlling artifact is pinned
  in [the prospective plan](preregistered/prospective_plan_v1.json).
- Preserve Teacher v3's original `s56_adp_raw_pair_v2` model-input contract for
  the primary prospective test. The corrected Teacher v4-SSL tensor contract is
  a different input distribution and is not silently substituted.
- Revalidate the accepted S63 Stage-1 products, build the authoritative cadence
  overlay and compact `DET_FLUX_ADP_SML + DET_FLUX_ADP` pair, and freeze the
  exact model-ready TIC allowlist before scoring.
- Define the primary cohort as S63 model-ready TICs absent from the frozen
  Teacher v3 corpus. Put S63 observations of Teacher-v3-corpus TICs in a
  separate repeated-host cohort. Never combine the two in the primary metric.
- Use the locked transparent BLS configuration: `50,000` periods, `10` retained
  peaks, both active ADP apertures, strict orbit/cadence authority, and the
  rank-one ADP-small ephemeris for one candidate row per TIC. The primary ADP
  aperture is retained as contamination/context evidence.
- Build a deterministic `1,000`-TIC primary queue and a `100`-TIC repeated-host
  side queue. The hidden bucket quotas are fixed in the JSON plan. Reviewers see
  neither probabilities nor bucket names.
- Complete the frozen review once, archive the human decisions, and unblind the
  selection provenance once. Do not tune Teacher v3, thresholds, quotas, or the
  label policy from the result.

The current filename inventory is useful planning evidence but is not the
final cohort authority. It contains `53,249` S63 TICs with nonzero HDF5 files:
`52,487` are Teacher-v3-disjoint and `762` are repeated hosts. Final counts come
only from the validated compact/model-ready export and will be bound in a new
immutable launch manifest before the first S63 model-scoring job.

## Primary reporting

Report human Planet-like morphology yield for each hidden selection bucket,
the fixed top-K Planet-like lift relative to the stratified control, and the
predeclared enrichment gate. Report the repeated-host side cohort separately.
The enriched sample cannot support a sector-wide recall, balanced-accuracy, or
prevalence claim; the deterministic top-score rows do not have survey-sampling
inclusion probabilities. The random/control stratum sizes and realized
sampling fractions will be retained so the control itself remains auditable.

After the prospective report is frozen, accepted human morphology decisions
may enter the next versioned training corpus. Teacher scores and pseudo-labels
may not.

## Freeze sequence

1. Commit this plan without reading S63 light-curve tensors or model scores.
2. Finish and validate the cadence, compact, BLS, candidate, and native-input
   artifacts using synthetic and prior-sector smokes first.
3. Publish `preregistered/launch_manifest.json` with the exact Git object and
   every S63 input/output hash listed as pending in the JSON plan.
4. Score once, construct the hidden queue atomically, and render the two-
   aperture review sheets.
5. Human-review all frozen rows, then unblind once and publish the prospective
   result. Any follow-up model or input-contract sensitivity is a new version.
