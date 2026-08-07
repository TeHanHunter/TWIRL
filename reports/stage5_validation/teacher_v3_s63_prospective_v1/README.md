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

## Blinding boundary and plan erratum

The committed prospective plan's statement that reviewer-visible fields
exclude “cohort membership” means that the queue withholds an explicit cohort
annotation. It does not mean cohort identity is cryptographically hidden. TIC,
candidate identity, and vet-sheet identity remain visible for scientific
vetting, so a reviewer could technically join them to the frozen Teacher v3
corpus and infer repeated-host status. The prospective-plan bytes remain
unchanged; the pre-score selection policy records this interpretation as an
explicit erratum. Cohort-wise unblinding and analysis remain deferred until all
`1,100` frozen rows have accepted morphology decisions.

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
  side queue. The exact formulas, depletion order, tie-breaks, within-cohort
  percentile population, seeded control draw, and hidden quotas are frozen in
  [the selection policy](preregistered/selection_policy_v1.json). Deterministic
  enriched rows have no sampling inclusion probability; only the seeded control
  rows retain stratum inclusion fractions and weights. Reviewers see neither
  probabilities, bucket names, nor an explicit cohort annotation. Visible TIC
  and candidate identity nevertheless make cohort status technically joinable.
- Complete the frozen review once, archive the human decisions, and unblind the
  selection provenance once. Do not tune Teacher v3, thresholds, quotas, or the
  label policy from the result.

The current filename inventory is useful planning evidence but is not the
final cohort authority. It contains `53,249` S63 TICs with nonzero HDF5 files:
`52,487` are Teacher-v3-disjoint and `762` are repeated hosts. Final counts come
only from the validated compact/model-ready export and will be bound in a new
immutable launch manifest before the first S63 model-scoring job.

The transparent BLS candidate gate subsequently identified `89` model-ready
TICs without a successful rank-one result in both required ADP apertures. All
`89` are in the disjoint primary cohort; none are repeated hosts. They remain
in model-ready/cohort accounting as explicit pre-score BLS failures but cannot
enter candidate, native, scoring, or queue artifacts. The resulting candidate
population is `53,160` TICs. This availability rule was recorded before any
Teacher-v3 S63 score was opened; it does not depend on model output and does
not change the immutable prospective-plan bytes.

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
   every S63 input/output hash listed as pending in the JSON plan, including
   the byte-exact selection-policy hash.
4. Score once, then publish the reviewer-safe queue and hidden selection
   provenance into distinct fresh, non-nested directories. Both directories
   are unpublished science and remain private (`0700` directories, `0600`
   files) until an explicit release decision. Both are fully staged before
   commit, and consumers require the matching bundle ID in the final public and
   private `bundle_complete.json` markers; an interrupted split publication is
   invalid. The public marker binds public files only, while the private marker
   retains the full cross-directory hash audit.
   The queue also re-proves the frozen Teacher-v3 release identity and all five
   checkpoint hashes, and it compares every candidate-derived field to the
   score table before selection. Render review sheets only from the verified
   reviewer-safe bundle.
5. Human-review all frozen rows without cohort-wise analysis, then unblind once
   and publish the prospective result. Any follow-up model or input-contract
   sensitivity is a new version.
