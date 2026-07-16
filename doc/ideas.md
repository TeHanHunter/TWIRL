# TWIRL unresolved questions and ideas

This file contains only unresolved scientific or engineering questions.
Accepted decisions belong in [the plan](twirl_plan.md); execution history
belongs in [the progress log](twirl_progress_log.md). Earlier brainstorming is
preserved in the [archived snapshot](archive/ideas_through_2026-07-13.md).

Last reconciled: `2026-07-16`.

## Blocking decisions

### A2v1 archive and sample control

- What exact versioned schema should join Gaia/TIC identity, detector visits,
  HDF5/FITS/compact paths, cadence retention, scatter metrics, product version,
  checksums, and QA state?
- What measurable gates freeze `twirl_parent_sample`? At minimum, define how
  `Pwd`, magnitude, coverage, Gaia/TIC bridge state, crowding, edge exclusions,
  and failed products affect the statistical denominator.
- What is the frozen last sector for the TWIRL I release, and what manifest
  schema records every accepted sector/product checksum, QA tier, catalog,
  search, and injection-contract version? `Sector >= 56` defines eligibility,
  not the publication cutoff.
- Do the `764` no-TIC-bridge targets occupy a scientifically important
  magnitude or sky regime that requires Gaia-first TGLC emission?
- For S94+, can the QLP/MAST products meet A2v1 target-coverage, raw-flux,
  aperture, quality, and provenance requirements, or must TWIRL re-extract?

### Transparent search contract

- Which non-periodic statistic should be the first dip-search baseline, and
  how should it represent isolated, clustered, asymmetric, and variable-depth
  events without assuming a stable period?
- Should the periodic baseline remain dense BLS alone or add a short-duration
  matched-filter/trapezoid branch informed by WD radius and expected duration?
- How are sector-level peaks combined into a target-level multi-sector object,
  including harmonics, ephemeris drift, non-detections, and irregular events?
- How will false-alarm rates be estimated for both periodic and dip branches?

### Completeness boundary

- What pixel-level calibration design is sufficient to measure the delta from
  LC-level injections across magnitude, crowding, detector edge, aperture
  disagreement, saturated neighbors, and centroid displacement?
- How should extraction failures and targets without TIC bridges enter
  completeness rather than being treated as search non-recoveries?
- Which event classes will receive separate completeness surfaces: compact
  periodic occultations, broad stellar eclipses, and irregular debris dips?

## Post-hoc teacher-v2 harmonic ambiguity

Teacher v1 received unbinned views at `P/4`, `P/3`, `P/2`, `P`, `2P`, `3P`,
and `4P`, but view availability did not provide rare-factor supervision. The
current harmonic table is concentrated at `P`; `P/3` has no supervised example
and `3P` has only three. Teacher v2 was nevertheless completed as an
exploratory comparison and missed its external-retention and morphology-
promotion gates. Teacher v1 therefore remains the active-learning baseline;
the questions below apply only to a future, predeclared model iteration, not a
retroactive rationale for promoting teacher v2.

Before any future teacher iteration is approved:

- Measure how often genuine A2v1 candidates need `P/3` or `3P`, separately for
  Planet-like and eclipse/contact morphology.
- Store `morphology_fold_factor` separately from any claimed physical orbital
  period, and require reviewers to mark it confirmed, changed, or unresolved.
- Compare morphology performance using the predicted factor with an oracle
  reviewed/truth factor so harmonic selection and morphology encoding failures
  are distinguishable.
- Test whether training-only reference-anchor augmentation improves rare-factor
  recovery without creating injection or provenance shortcuts. Final evaluation
  must remain anchored to genuine ADP-small BLS ephemerides.
- Report top-1/top-2 factor accuracy by factor and source, period-ratio error,
  calibration, and morphology metrics under both predicted and oracle factors.
- Define the minimum advantage over teacher v1 and the external-sector
  retention gate before training, and keep S57 out of any "pristine holdout"
  claim because its premature review queue contained six human labels at audit
  close. Preserve all six and pause further S57 holdout consumption.

## Catalog and contextual science

- Which polluted-WD, infrared-excess, known PCEB/double-WD, and pulsator
  catalogs have adequate provenance and sky coverage for target flags?
- Should CVZ and multi-sector WDs be processed first once the multi-sector
  search exists, and what priority metric avoids biasing the statistical run?
- Is a parallel WD pulsation catalog worth supporting without distracting from
  the transit/dip completeness gates?

## Candidate validation

- What quantitative pixel-map likelihood or offset threshold defines
  on-target versus neighboring-source variability?
- How should per-pixel residual cubes be stored or reconstructed without
  moving the full PDO cutout tree to ORCD?
- Which ZTF/ATLAS checks can be automated robustly for faint WDs, and what is
  the correct outcome when archival coverage is absent?
- What confirmation strategy applies to irregular events that cannot be
  scheduled at a deterministic future transit window?

## Follow-up and collaboration

- Verify access, cadence, filters, sensitivity, and lead time before adopting
  SPECULOOS, MISCOT/Avi, LCO, EPRV, proto-Lightspeed, WINTER/SPRING, Magellan,
  or MMT as a concrete route.
- Determine the funding/proposal path for rapid high-cadence follow-up and the
  ephemeris precision required by each facility.
- Confirm paper leadership and discovery-response responsibilities in writing;
  do not infer them from pipeline ownership or telescope access.

## Referee-risk checklist

- Parent sample and exclusions are frozen and machine-readable.
- Search sensitivity is not inferred from WD 1856 alone.
- Completeness traverses extraction calibration, search, production ranker,
  vetting, and candidate merging.
- Gaia-to-TIC losses and heterogeneous sector coverage are explicit.
- Claims name an observable event class unless physical classification is
  independently supported.
- Null-result limits are restricted to the supported parameter space and
  declared sample.
- Debris-search sensitivity is calibrated to variable/asymmetric events rather
  than periodic boxes alone.

## Maintenance

When an item becomes a decision, move its concise outcome to
`doc/twirl_plan.md` and remove it here. Do not append dated execution notes to
this file.
