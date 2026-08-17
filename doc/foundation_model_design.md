# TWIRL-FM working design

This document is the frozen scientific design for the `TWIRL-FM0.1`
proof-of-concept campaign and the revisable long-term design for a new TWIRL
light-curve foundation-model candidate. It is subordinate to the
[TWIRL survey plan](twirl_plan.md), which remains the single authority for
forward priorities and promotion gates. The frozen FM0.1 section authorizes
input construction, local validation, and gated ORCD proof-of-concept jobs. It
does not authorize production use, change the A2v1 product contract, or replace
the transparent periodic and non-periodic searches.

Design version: `TWIRL-FM0.1`  
Model family: `TWIRL-FM0`  
Status: FM0.1 scientific design frozen; implementation authorized; no FM input
release, job, or checkpoint exists yet  
Proof-of-concept corpus: S56--S67, with product state retained explicitly  
Long-term target: at least 50 accepted `200 s` A2v1 sectors

The machine-readable companion contract is
[`configs/models/twirl_fm0_1_s56_s67_poc.yaml`](../configs/models/twirl_fm0_1_s56_s67_poc.yaml).
The document and configuration are bound by the freeze receipt in
[`reports/stage5_validation/twirl_fm0_1_design_freeze_v1/`](../reports/stage5_validation/twirl_fm0_1_design_freeze_v1/).

## Short working roadmap

1. Build one checksum-bound, BLS-free FM input release spanning S56--S67.
   S56--S65 are accepted A2v1 sectors. S66--S67 are allowed only as explicitly
   marked, hash-bound `ORCD_COMPLETE_DEFERRED` development inputs after the
   FM-specific cadence, quality, channel, and numerical gates pass. This does
   not promote them to accepted A2v1 sectors.
2. Split complete Gaia--TIC leakage components into 80% training, 10%
   iteration-development, and 10% sealed proof-of-concept test groups. All
   twelve sectors may occur in each source partition; all visits and aliases of
   one source remain together.
3. Run the five controlled `FM0.1.1`--`FM0.1.5` comparisons on ORCD. They test
   TCN versus Conformer structure, ADP versus ADP+ADP015 versus
   raw-relative+ADP+ADP015 inputs, and masked reconstruction with or without a
   same-window consistency term.
4. Choose among those variants using the development partition only. Open the
   sealed 10% source test once for the frozen finalist and publish the complete
   error ledger. Its job is to expose data mistakes, representation collapse,
   detector shortcuts, event loss, weak transfer, and excessive compute
   cost--not to establish a final model.
5. If the sealed result motivates `FM0.2`, freeze a new sealed evaluation before
   using the FM0.1 test failures to redesign the model. The 20-sector
   and 50-plus-sector studies begin only after this loop produces a stable,
   useful representation and a credible scaling hypothesis.

This is similar to gradient boosting as a development philosophy: each version
focuses on the residual failures of the preceding version. It is not a literal
boosted ensemble unless a later, separately named experiment explicitly combines
the predictions of several frozen models.

## Naming and scope

This project starts a new naming sequence so that it cannot be confused with
the earlier `s56_harmonic_cnn_v1`, Teacher v1--v3, or Teacher v4-SSL work.

- `TWIRL-FM0` is the exploratory model family.
- `TWIRL-FM0.1` is the frozen S56--S67 proof-of-concept campaign and data/split
  contract.
- `TWIRL-FM0.1.1`--`TWIRL-FM0.1.5` are controlled candidates inside that
  campaign. A suffix changes only the declared architecture, input view, or
  objective component in the frozen comparison ladder below.
- Seeds do not create new model versions. Exact runs use names such as
  `twirl_fm0_1_2_seed560068`.
- `TWIRL-FM0.2` is reserved for an error-directed redesign after FM0.1 is
  complete. It is not another point in the FM0.1 grid.
- `TWIRL-FM1` is reserved for the first evaluated architecture/objective family
  that passes the multi-task transfer, prospective-sector, injection-retention,
  confound, and provenance gates in this document. A later all-data retrain is
  explicitly a deployment candidate until its own genuinely new evaluation.

Until those gates pass, the scientifically accurate description is **TWIRL
foundation-model candidate** or **pretrained white-dwarf light-curve encoder**.
Large-scale self-supervised training alone does not establish that a model is a
foundation model; useful reuse across distinct tasks is the central test.

The initial domain is all eligible Gaia-selected TWIRL white-dwarf light
curves, not all TESS stars. Expanding to a general TESS encoder would require a
new population, data-release, and evaluation contract.

Two identity fields serve different purposes throughout this design.
`physical_source_id` is the canonical Gaia DR3 `source_id` for one scientific
source. `leakage_component_id` is an opaque identifier for the complete
connected component of the Gaia--TIC alias graph. The first field says which
star is being studied; the second is the conservative grouping key used for
splits, sampling weights, deduplication audits, and clustered uncertainty.

## Decisions made in this working design

1. The physical source is the scientific unit. The more conservative
   `leakage_component_id` is the independence, split, sampling-weight, and
   uncertainty unit whenever aliases are involved.
2. Gaia DR3 `source_id` is the canonical `physical_source_id`. TIC IDs remain
   operational aliases, and every physical source connected through a shared
   alias stays in one leakage component.
3. Sector observations and cadence windows remain separate observations below
   the host. They are not averaged into one light curve.
4. Pretraining is population-wide and candidate-blind. It does not run or
   consume BLS.
5. The encoder sees only chronological, quality-aware FM derivatives of
   immutable A2v1 raw aperture photometry. The frozen ladder begins with ADP
   `1x1+3x3`, then adds ADP015, then adds a normalized raw-relative view; `5x5`
   is not model-visible in FM0.1.
6. Local time and gaps retain real duration in nominal `200 s` cadence units.
   A separate host-relative visit-time coordinate preserves the true elapsed
   time between sector visits without exposing absolute BJD or a sector ID.
7. The model preserves window-, sector-, and host-level representations so
   later tasks can use the appropriate scientific unit.
8. Architecture size and training duration are selected by held-out transfer,
   robustness, calibration, and cost, never by the minimum training loss.
9. A final all-allowed-data checkpoint is trained only after the design is
   frozen and honestly evaluated. It then requires genuinely later or external
   data for its next prospective test.

## What counts as one data point

One star is one **scientific source unit**, but not one flat tensor. A clean
one-source leakage component is also the independence unit; unresolved
multi-source components are quarantined. A host may have one sector or dozens
of sectors, and each sector may contain events that do not occur in the others.
Averaging all sectors into
one series would erase transient information; treating every window as an
independent star would overstate the sample size and let heavily observed
continuous-viewing-zone stars dominate.

The registry therefore separates scientific identity, leakage grouping,
stable observation identity, and a versioned product instance:

```text
leakage_component_id      connected Gaia--TIC alias component; split unit
  physical_source_id      one canonical Gaia DR3 source_id
    observation_key       physical_source_id plus sector; stable scientific visit
      product_instance_id observation_key plus selected TIC product, version, hash
        segment_id        orbit/gap-bounded chronological segment
          window_id       segment, start index, length, and window-contract version
```

`observation_key` is independent of file layout and product version. Orbit,
camera, CCD, cadence, and detector coordinates remain provenance attached to
the observation or its segments; they do not create additional scientific
visits. If more than one TIC product represents the same
`(physical_source_id, sector)`, a frozen selection/deduplication rule must choose
one product or quarantine the visit. The aliases are never counted as two
stars. `product_instance_id` records the selected TIC alias, A2v1 and FM schema
versions, and source checksum so that a rebuilt product cannot masquerade as a
new observation.

The source partition is assigned to `leakage_component_id` before windows,
augmentations, normalization fits, or downstream labels are constructed. Every
physical source, alias, observation, window, augmentation, and later injection
in that component inherits the same **source partition**. A separate sector-era
field supplies the chronological role described below; these are two
orthogonal assignments rather than one overloaded split. Metrics and
uncertainty intervals resample whole leakage components.

The leakage-component builder must take the connected closure of Gaia--TIC
aliases. A clean component contains one Gaia DR3 source and one or more TIC
aliases. Components containing multiple plausible Gaia sources fail closed in
the first input release and remain in an audit table rather than being silently
assigned or treated as negatives. Gaia-only sources also remain visible even
if the current TIC-oriented extraction cannot yet emit them.

## Corpus and sector admission

The long-term pretraining population is every eligible, product-validated A2v1
light curve in a frozen release. It is not the set of successful BLS candidates
and it is not limited to labeled sources, the `Pwd > 0.75` statistical sample,
or a magnitude-selected subset. Sample-membership fields remain available for
auditing but are not model inputs.

The release is named abstractly as `S56--S_cut` until the survey cutoff is
frozen. A target of at least 50 accepted sectors is a scale goal, not a claim
about the current archive. S56--S93 contains only 38 sectors, and the S94+
input path remains an open survey decision. ORCD HDF5 compute checkpoints are
not accepted sectors until their required return and unchanged A2v1 product
gates pass. If S56--S105 were all accepted, that interval would contain exactly
50 sectors, but this design does not preselect S105 as the release cutoff.

An accepted long-term FM sector must have:

- checksum-bound A2v1 HDF5 and FITS provenance;
- edge-aware coverage and openability validation;
- the authoritative internal-plus-external cadence-quality evidence;
- explicit omissions rather than fabricated light curves;
- an FM-specific channel, mask, finite-value, and numerical-envelope report;
- a stable join to the physical-source registry.

A diagnostic sector may support representation research before it is promoted
to an accepted product or the statistical survey sample, but its product state
and QA state must be retained for sensitivity analysis. The foundation-model
corpus, the accepted A2v1 archive, and the occurrence-rate denominator are
different objects.

For the bounded FM0.1 proof of concept, S66 and S67 may enter as diagnostic-only
retained ORCD HDF5 checkpoints after their immutable cell receipts,
openability, cadence authority, external-quality join, and FM channel gates pass.
Their rows carry a non-model-visible product-state field and are reported
separately in all data summaries. They remain ineligible for production,
prospective, completeness, or accepted-sector claims until the unchanged PDO
return, sector-FITS, and full-schema gates pass.

## Frozen FM0.1 model-visible input contract

`TWIRL-FM0.1` uses a strict allowlist. The complete compact release stores the
following six aligned flux views:

```text
                    1x1 aperture             3x3 aperture
raw-relative        raw_relative_1x1          raw_relative_3x3
ADP (0.3 d)         adp_1x1                   adp_3x3
ADP015 (0.15 d)     adp015_1x1                adp015_3x3
```

The variants select subsets of that one release; they do not rebuild or
silently renormalize the same star. The `5x5` aperture remains in provenance
and photometric QA but is not model-visible in FM0.1. Absolute raw detector-unit
flux is also not model-visible.

For each aperture, let `m` be the median raw flux over its effective-good
cadences and let `s` be the positive robust subtractive scale returned by the
unchanged ADP scale rule on those same cadences. The stored flux views are

```text
raw_relative = (RawFlux - m) / s
adp          = DET_FLUX_ADP - 1
adp015       = DET_FLUX_ADP015 - 1
```

For an ordinary positive-baseline source, the first expression is the familiar
fractional raw variation. The robust `s` keeps the definition finite for faint
or background-subtracted sources whose median is weak or negative. The builder
fails the aperture if `m` or `s` is not finite or `s <= 0`; it never substitutes
a per-window standard deviation. The two model-visible measurement-error
proxies are `RawFluxError_1x1 / s_1x1` and
`RawFluxError_3x3 / s_3x3`. They describe measurement noise and do not include
uncertainty in either fitted cotrend.

The encoder additionally receives:

- local time since the first valid cadence and cadence-to-cadence time
  differences, both divided by `200 s` so that an ordinary cadence interval is
  approximately one;
- segment-boundary, time-validity, per-view flux-validity,
  error-availability, and view-presence masks;
- a view identity indicating aperture and whether the view is raw-relative,
  ADP, or ADP015.

The release also stores two host-level timing coordinates for future
multi-sector combination: the start of each visit relative to that host's first
represented visit and the gap from the preceding represented visit, again in
nominal `200 s` cadence units. They retain the real elapsed duration across
sector gaps while withholding absolute BJD and explicit sector identity. They
are audit/aggregation inputs, not window-encoder inputs in FM0.1.

The model receives chronological values, not phase-folded arrays. Flux remains
in fractional units so transit depth and noise scale are preserved. The
pipeline must not standardize every star, channel, or window to unit variance,
because that would make a shallow quiet-star dip artificially resemble a deep
noisy-star dip. Any population-level numerical transform is fitted on training
hosts only and frozen by checksum.

Invalid or missing cadences are replaced by a neutral numerical value only
alongside an explicit mask. A gap is never interpreted as measured flux and is
never imputed from another source or split. Detailed quality bits remain in the
audit table; the primary encoder sees the effective validity mask.

Time, segmentation, and windows follow one versioned rule. The builder orders
rows by the authoritative cadence number and time **before** quality filtering;
flagged rows remain in place and are masked rather than compressed out. A new
segment begins at an orbit change or an authoritative time gap greater than
`0.2 d`, matching the ADP/ADP015 gap boundary. A window is a contiguous slice
of one segment in that uncompressed cadence order. Its declared length, stride
or random-start policy, edge-padding rule, and random seed are part of the
window-contract version. A window never crosses a segment boundary; a short
edge window is padded only with neutral values and false masks. The frozen
FM0.1 context is `2,048` cadences, about `4.74 d` at nominal cadence. Training
uses deterministic-seeded random window starts; evaluation tiles each segment
with a `1,024`-cadence stride. The encoder receives the dimensionless local
timing coordinates defined above, not absolute BJD. Absolute time, cadence
number, orbit, sector, camera, and CCD remain in the audit/provenance table.

### Quality-aware derived input release

The current A2v1 FITS files contain the three-aperture detrended products, but
only a primary-aperture error column for each detrending branch. Those FITS error
columns are a broadcast residual-MAD scatter estimate, not propagated
cadence-level measurement uncertainty. In addition, the internal `QUALITY`
field is not the final external SPOC/QLP quality authority. Earlier full-pool
SSL work showed that external-only bad cadences could contaminate detrending
before being masked.

The FM input builder therefore reads immutable A2v1 HDF5 photometry, joins the
authoritative external quality table, and reproduces ADP/ADP015 for `1x1` and
`3x3` with final effective quality. For S56--S65 the source products are the
accepted A2v1 HDF5 archive; for diagnostic S66--S67 they are the hash-bound
retained ORCD HDF5 checkpoints admitted by the stricter exception above. It
then writes one new checksum-bound **FM input release** containing the six
derived views, two aperture-level measurement-error proxies, masks,
dimensionless timing, and opaque registry keys. This is a derived model input,
not an overwrite or silent revision of A2v1. Any difference from the accepted
Stage 1 scientific recipe has a new named contract and parity report.

For each cadence,

```text
external_quality = spoc_quality | (qlp_quality << 30)
effective_good = finite(time)
                 AND internal_quality == 0
                 AND external_quality == 0
                 AND NOT authority_excluded
aperture_good[a] = effective_good AND finite(RawFlux[a])
error_good[a] = aperture_good[a]
                AND finite(RawFluxError[a])
                AND RawFluxError[a] > 0
view_good[a,v] = aperture_good[a] AND finite(derived_flux[a,v])
```

The detrending fit uses `effective_good` together with the aperture's finite raw
flux and positive finite raw error. The model uses `view_good` separately for
every view; it must not discard a healthy aperture merely because another
aperture is missing. A variant that requires a wholly missing view is excluded
for that observation or evaluated as an explicitly named reduced-view
sensitivity. The error proxy uses the aperture scale `s` above. It does not
include uncertainty in the fitted spline, so empirical residual scatter or a
detrending-uncertainty estimate must be named and stored separately.

The builder consumes the cadence authority's declared exclusions. In
particular, it propagates the bounded S56 cadence-`699957` camera-3 exclusions
and masks those rows. Every undeclared missing detector cadence or authority
join gap fails closed; the S56 exception is not generalized to another sector.
The parity report binds source identities, cadence order, raw flux/error arrays,
and the unchanged ADP/ADP015 settings, while recording the expected numerical
difference caused by the declared effective-quality change. It does not require
the new derived flux to equal the older internal-quality-only FITS columns.

### Inputs excluded from the primary encoder

The following fields are stored separately for audit, subgroup evaluation, or
later preregistered ablations and are not passed to `forward()`:

- Gaia or TIC identity, RA/Dec, `Pwd`, Tmag, or catalog labels;
- absolute BJD and explicit sector, camera, CCD, or detector identifiers;
- absolute detector-unit raw flux, background, centroid, or crowding
  measurements; the frozen normalized raw-relative view is the only raw-flux
  derivative allowed;
- human labels, pseudo-labels, Teacher scores, or Teacher embeddings;
- injection truth or injected training rows;
- every search-derived field listed in the next section.

Instrument context is not inherently a label, but direct identifiers and
detector-unit signals offer strong shortcuts. FM0.1.4 therefore tests only the
normalized raw-relative views on the identical source split. Audit-only
sector/camera/CCD metadata is used to measure how much instrument context can
still be inferred from the light curves. Metadata late fusion or a dual
physical/instrument representation requires a later named design.

## Absolute BLS and supervision boundary

The pretraining builder must never run, join, filter, oversample, center, fold,
or normalize using:

- BLS period, epoch, duration, depth, power, SDE, peak, rank, or periodogram;
- any phase-folded or candidate-centered view;
- candidate tables, review queues, selection buckets, or model scores;
- morphology labels, period-factor decisions, pseudo-labels, or injection
  truth.

BLS computed from the same light curve is not future-information leakage in
the usual sense. It is nevertheless a strong task-conditioned transformation
that selects periodic candidates and gives the model the answer to a specific
search question. Excluding it is what distinguishes `TWIRL-FM0` from the
archived Teacher v4-SSL encoder, whose inputs were built around rank-one BLS
periods, harmonic folds, local event windows, and periodograms.

Transparent BLS and the future non-periodic dip search remain independent
discovery baselines. Their outputs can be evaluated beside a frozen encoder or
combined in a separately trained downstream ranker. They never enter core
pretraining, and the foundation model never substitutes for the real
search-to-vetting-to-merging completeness chain.

Synthetic events also have a strict development/evaluation boundary. A
checksum-bound **development injection canary** may use only pretraining-train
or pretraining-development source components and training-era or temporal-
validation observations to choose augmentations, loss weights, or checkpoints.
A separately generated **sealed injection evaluation set** uses source
components and/or temporal-test sectors absent from that tuning and is opened
only after the model and downstream evaluation contract are frozen. Injection
truth is available to the evaluator, never to the encoder input or core
self-supervised loss. WD 1856 remains a known engineering benchmark; because it
may guide development, it is never described as a prospective source test.

## Split and evaluation design

The development layout has two independent axes. `source_partition` is frozen
once for every `leakage_component_id`; `sector_era` is frozen once for every
accepted sector. A sample's role is the intersection of those fields. A later
visit to a pretraining-train source therefore remains in that source partition
while legitimately belonging to temporal validation or temporal test.

### FM0.1 proof-of-concept split

FM0.1 uses all available S56--S67 sector coverage but does not claim temporal
generalization. Complete leakage components are hashed once into 80%
`poc_train`, 10% `poc_development`, and 10% `poc_sealed_test`. The exact hash is
`sha256("twirl_fm0_1_source_split_v1:" + leakage_component_id)` interpreted as
one big-endian unsigned 256-bit integer modulo `10,000`: buckets `0--7,999`
train, `8,000--8,999`
development, and `9,000--9,999` sealed test. Encoder gradients, normalization,
augmentation fits, architecture choice, input-view choice, loss choice, and
checkpoint selection use no sealed-test row. The development partition
supplies the FM0.1.1--FM0.1.5 comparison and error ledger. The sealed test is
opened once for the frozen finalist. If its failures guide FM0.2, a new test
partition or genuinely later data must be frozen before FM0.2 evaluation.

All sectors may occur in every source partition, but every alias and sector
visit of one leakage component stays together. This means "use all sectors,"
not "train on every star." Repeated hosts provide a small realistic test of
cross-sector aggregation, while host-first sampling prevents them from
receiving more weight merely because they have more visits.

Because S63 enters the unlabeled FM corpus, S63 can never be described as a
prospective FM evaluation. This does not change the independently frozen
Teacher-v3 S63 experiment provided no FM embedding, score, queue, threshold, or
review decision influences it. Any later FM probe using S63 labels is explicitly
transductive transfer evidence.

### Long-term 50-plus-sector split

For a frozen release of `N >= 50` accepted sectors:

1. Hash every clean `leakage_component_id`, using one frozen salt, into 80%
   `pretrain_train`, 10% `pretrain_development`, and 10%
   `sealed_internal_host_test`. All aliases and physical sources in the
   component receive the same assignment. The internal-host-test components
   are excluded from encoder fitting, normalization, augmentation fitting,
   architecture selection, and temporal-validation tuning in every sector.
2. Order only accepted sectors chronologically. Freeze a leading
   `training_era`, a following `temporal_validation` block containing at least
   five accepted sectors, and a final `sealed_temporal_test` block containing
   at least five accepted sectors. Five is a minimum, not an automatic choice.
   Before freezing, an inventory must show adequate eligible-source counts,
   unseen and repeated hosts, detector/camera coverage, cadence retention, and
   channel availability in both later blocks. A block may be enlarged or the
   proposed cutoff moved only from this label-free inventory, never after model
   or task results are inspected.
3. Encoder gradients and corpus-level transform fits use only
   `pretrain_train × training_era`. Architecture, objective, context length,
   capacity, and stopping choices may use
   `pretrain_development × training_era` and the non-internal-test portions of
   `temporal_validation`, but temporal-validation rows never receive encoder
   gradient updates.
4. `sealed_internal_host_test × training_era` is the familiar-era unseen-host
   test. It remains sealed until all source-era development choices are frozen.
   The complete `sealed_temporal_test` block is also sealed and cannot set
   normalization, augmentation, thresholds, probe hyperparameters, or stopping
   choices.
5. In each later-sector result, sources absent from the complete training-era
   inventory are the primary new-host cohort. Sources with an earlier
   training-era visit form a separate repeated-host time-transfer cohort. The
   two are never mixed into one primary number, and temporal-validation results
   are never relabeled as prospective test results.

This produces three scientifically different evaluations:

| Cohort | What it tests |
| --- | --- |
| sealed internal-test component, training-era sector | generalization to another star under familiar survey conditions |
| repeated host, sealed later sector | stability across time and detector visits without claiming a new star |
| new host, sealed later sector | simultaneous source and temporal/instrument generalization; primary prospective result |

Downstream heads receive their own label splits. A head evaluated on hosts that
were seen unlabeled during pretraining is a legitimate transductive transfer
test, but it must be labeled as such. The strict inductive result uses hosts
absent from both encoder pretraining and head training.

The S56--S67 proof of concept must not compare a newly split FM encoder casually with the
frozen Teacher-v3 number. If the original Teacher-v3 fixed-test result is reused,
every leakage component in that fixed test remains forbidden to FM pretraining
and tuning. Otherwise, the `s56_harmonic_cnn_v1` architecture is retrained under
the new source/era registry and reported as a fresh architecture baseline, not
as the frozen Teacher-v3 checkpoint. Any probe on a Teacher-v3 test host seen
unlabeled by the FM encoder is explicitly transductive and cannot replace the
strict comparison.

After the architecture, objectives, input release, and all evaluation choices
are frozen and the sealed test is reported once, a release checkpoint may be
trained from scratch on every eligible unlabeled observation in the frozen
`S56--S_cut` release. That all-data checkpoint no longer has the old sectors as
an untouched test; genuinely later sectors or an external light-curve set must
be reserved for its next evaluation.

## Sampling without allowing repeated hosts to dominate

The primary sampler is hierarchical:

1. sample a physical source approximately uniformly;
2. sample one of its sector visits;
3. sample one uninterrupted chronological segment;
4. sample a fixed-length window within the segment;
5. construct safe self-supervised views only after the split and window are
   known.

This gives approximately equal expected weight to every host rather than every
available window. Batch and epoch reports must show effective exposure by
host, sector count, magnitude, detector, and valid-cadence fraction.

Uniform windows are the unbiased primary stream. A secondary, preregistered
stream may devote a minority of batches to locally high-deviation windows so
rare events are not almost never observed. That selector may use only a simple
time-local robust deviation statistic on allowed flux and masks; it cannot use
BLS, labels, injections, or candidate catalogs. Its fraction and sampling
weights must be frozen before evaluation, and uniform-window performance is
reported separately.

## Model hierarchy

The [model hierarchy figure](../reports/exploratory/twirl_fm0_working_plan_v0_1/twirl_fm0_hierarchy.png)
summarizes the intended data and representation levels. One physical host
contains sector observations; sector observations contain chronological
windows. A shared encoder produces local representations, a sector aggregator
combines windows from one visit, and a host aggregator combines visits without
assuming a fixed sector count.

Three embeddings remain available:

- `z_window` retains local morphology and supports event localization,
  reconstruction, dip classification, and short-timescale regression;
- `z_sector` summarizes a visit and supports sector-level classification,
  variability, recurrence, and ranking;
- `z_host` summarizes evidence across visits and supports stable source
  classification, retrieval, and multi-sector context.

The host aggregator does not force all sectors of a star to be identical.
Stable host information and time-varying local information remain separate.

## Architecture candidates

The preferred candidate is a hierarchical local-to-global encoder, but it must
earn that status against simpler baselines.

### Cadence stem

A small residual 1-D convolution or temporal-convolution stack operates before
subsampling. Kernels span roughly one to 15 cadences so that one-cadence
outliers, short white-dwarf occultations, ingress/egress, and somewhat broader
local structure remain resolvable. All convolutions are mask-aware.

### Window encoder

The leading design combines the local stem with time-aware Conformer blocks.
FM0.1 freezes the context at `2,048` cadences, about `4.74 d` at `200 s`
cadence. A local convolution acts before any patching or stride so that a
one-cadence event remains visible. The Conformer receives local time and gap
coordinates in nominal-cadence units; it never receives absolute time or a
sector ID.

A pure TCN with the same input and approximate parameter budget is the required
non-attention baseline. If the hybrid model cannot beat it outside uncertainty,
the simpler TCN remains preferred.

### Sector and host aggregation

The sector aggregator combines window embeddings with relative position,
coverage, and gap information. In FM0.1 it is compared with masked mean pooling
and the best single window. Repeated-host evaluation compares sector embeddings
and simple masked means, while retaining the true host-relative visit offsets
for interpretation. A learned host aggregator and any cross-sector consistency
loss are deferred to FM0.2 unless the repeated-host inventory and FM0.1
diagnostics justify them. A future host aggregator must use the stored elapsed
visit-time coordinates rather than treating sector visits as equally spaced.
Explicit sector and detector IDs are not encoder inputs; they are used only for
subgroup and shortcut probes.

## Frozen FM0.1 comparison ladder

Every adjacent candidate uses the same input-release rows, 80/10/10 source
split, `2,048`-cadence windows, host-first sampler, optimization budget, seeds,
and evaluation panel. Only the stated mechanism changes:

| Candidate | Architecture | Flux views exposed | Objective |
| --- | --- | --- | --- |
| `TWIRL-FM0.1.1` | TCN | ADP `1x1 + 3x3` | synchronized masked reconstruction |
| `TWIRL-FM0.1.2` | parameter-matched Conformer | ADP `1x1 + 3x3` | synchronized masked reconstruction |
| `TWIRL-FM0.1.3` | development winner of `.1`/`.2` | ADP + ADP015, both apertures | synchronized masked reconstruction |
| `TWIRL-FM0.1.4` | same as `.3` | raw-relative + ADP + ADP015, both apertures | synchronized masked reconstruction |
| `TWIRL-FM0.1.5` | same as `.4` | same six views | reconstruction + same-window VICReg |

If adding views reverses the apparent TCN/Conformer conclusion or the
difference is within source-clustered uncertainty, the other architecture is
rerun once with the winning view set before a finalist is named. That is a
declared interaction check, not a new model family.

Masked reconstruction hides a deterministic-seeded target fraction of `0.15`
using synchronized contiguous spans of `1--64` cadences. The same times are
hidden in every available raw-relative and detrended view so no view can copy
the answer from another. Invalid cadences are never unmasked and never enter
the target loss. The reconstruction loss is the mean of per-view Huber losses
over masked valid targets; views receive equal weight so adding channels does
not silently multiply their importance. The initial Huber transition is
`0.01` in fractional-flux units and is changed only through a new named
contract if a numerical smoke shows it is pathological.

FM0.1.5 adds two independently masked, otherwise safe views of the same window.
Its total objective is `L_reconstruction + 0.1 * L_VICReg`, with VICReg's
invariance, variance, and covariance terms weighted `25:25:1`. This is a frozen
proof-of-concept choice rather than an inherited claim that Teacher v4-SSL's
weights were optimal. The minimum pretraining loss is never a model-selection
metric.

- **Temporal-mask reconstruction** hides all flux views at the same selected
  cadences or contiguous spans. This prevents the model from solving the task
  by copying the simultaneous value from another aperture. Measurement-error
  proxies are model inputs and audit quantities in FM0.1; they do not weight
  this first reconstruction loss.
- **Same-window consistency** uses two safely perturbed views and the frozen
  variance/covariance anti-collapse term above.
- Whole-channel prediction, neighboring-window prediction, learned host
  aggregation, and stable-host cross-sector consistency are deliberately
  outside the frozen FM0.1 ladder. They are candidate responses to measured
  failures, not hidden extra objectives.

Injections are not part of the core self-supervised loss. Two disjoint,
checksum-bound injection suites are defined before objective selection. A
**development canary** may be inspected while choosing augmentations, masking,
and objectives. A source- and injection-parameter-disjoint **sealed canary** is
opened once with the FM0.1 sealed source test and cannot set thresholds or
stopping choices. That result is a proof-of-concept event-retention test, not a
prospective-sector result. Both suites traverse the same allowed input builder and
downstream head as real light curves, and every injected host inherits its real
host split. A later synthetic-event adaptation stage, if useful, receives a
new name, new training rows, and a separate evaluation contract.

## Safe and unsafe augmentations

Primary safe augmentations are uncertainty-scaled noise, additional mask-only
dropout, and alternative allowed aperture subsets. They may not alter valid
time coordinates or the physical fractional depth.

The primary contract forbids time reversal, arbitrary time warp, smoothing,
per-window variance normalization, depth rescaling, flux-sign inversion,
candidate recentering, and event insertion or removal. Cropping is permitted
only through the declared chronological window sampler, not as an
event-centered transformation.

Every augmentation is first tested on injected and real benchmark events to
show which morphology it preserves. Event preservation is evidence for the
augmentation contract, not training supervision.

## Capacity and data scaling

FM0.1 is not a capacity race. It uses one small proof-of-concept budget with a
`256`-dimensional window embedding. The TCN and Conformer implementations must
be within 10% of one another in trainable scalar parameters after the first
H200 compile/throughput smoke; the target range is approximately 8--12 million
parameters. Both see the same number of model-visible cadences, optimizer
steps, host sampler, and two frozen seeds (`560067`, `560068`). If either model
cannot fit the one-H200 smoke envelope without changing context or batch
semantics, that is a measured FM0.1 result rather than permission to give it a
different compute budget.

Nested 1%/10%/30%/100% data subsets and Tiny/Base/Large capacity scaling are
deferred until the architecture and input contract survive FM0.1. Those later
subsets are deterministic samples of training-only leakage components, never
samples of windows, and preserve fixed validation and test cohorts.

Training curves are used to diagnose optimization, collapse, and overfitting.
They do not determine model complexity by themselves. Checkpoints are selected
using a predeclared composite of `poc_development` transfer tasks,
representation health, event retention, and compute efficiency. The long-term
50-plus-sector design later adds a temporal-validation block.

## Downstream model family

The frozen encoder is evaluated through the same progression for every task:

1. a true linear probe;
2. a small nonlinear head with the encoder frozen;
3. partial unfreezing of the final encoder stages;
4. full fine-tuning;
5. the identical architecture trained from random initialization.

Frozen and fine-tuned label-efficiency curves use nested, source-grouped
`1%`, `5%`, `10%`, and `100%` subsets of each task's training labels. The
validation and test support is identical across fractions, and rare-class
minimums are reported; a nominal fraction that contains no supported positive
class is omitted rather than presented as a meaningful comparison.

Before a task may count toward foundation-model promotion, a versioned task
registry records its target authority, scientific unit, eligible population,
split registry, real/injected/transductive status, label and subgroup counts,
minimum supported-class count, head configuration, calibration set, primary
metric, and failure threshold. In particular, morphology truth remains
candidate-selected human morphology; localization, depth, and duration rows
identify whether their truth is injected or independently measured; period
truth records allowed harmonic equivalents; and anomaly evaluation names its
reference or blinded-review authority. Three synthetic tasks cannot by
themselves satisfy multi-task promotion: at least two promoted task families
must have real-data evaluation support.

Planned tasks are deliberately different enough to test reusable
representation rather than one classification problem:

- human morphology classification;
- periodic-candidate and non-periodic-dip ranking after the transparent search;
- per-cadence or per-window event localization;
- depth and duration regression with calibrated uncertainty intervals;
- visit- or host-level recurrence/period estimation;
- anomaly detection and same-host-excluded similarity retrieval;
- multi-sector source variability classification;
- sector/camera/CCD and magnitude probes used only to measure shortcuts.

Learned sector and host aggregators are compared with the best single
window/sector, masked mean pooling, and parameter-matched attention pooling
without host-level pretraining. These hierarchy baselines distinguish a better
representation from the simpler benefit of exposing a head to more visits.

Regression heads must predict uncertainty as well as a point value, for
example through quantiles or a conditional distribution. Calibration and
coverage are assessed on held-out hosts; an interval is not accepted merely
because the head outputs a scale parameter.

## Required baselines

Every task compares against the applicable row of a frozen task-to-baseline
matrix rather than claiming that one baseline applies to every problem:

| Task family | Trivial baseline | Transparent or scalar baseline | Learned baseline |
| --- | --- | --- | --- |
| candidate morphology/ranking | class prior and random review queue | versioned scalar logistic model with L1/L2 controls; BLS/dip heuristic where applicable | Teacher v3 on its exact common support; supervised CNN/TCN; random-initialized new backbone |
| localization and event regression | no-event predictor and training median | transparent BLS/dip event estimate where applicable | supervised CNN/TCN and random-initialized new backbone |
| period/recurrence | training-period prior | BLS and the declared non-periodic recurrence statistic | supervised sector/host model and random-initialized new backbone |
| variability classification | majority class | versioned scalar linear/logistic model | supervised CNN/TCN, simple pooling, and random-initialized new backbone |
| retrieval/anomaly | random retrieval/review | robust-deviation ranking and PCA/scalar distance | untrained embedding, objective ablations, and compatible public encoders |

Masked-reconstruction-only, VICReg-only, aperture-only, and host-objective
ablations accompany the relevant learned rows. Compatible public light-curve
or time-series foundation encoders are included only when their input and
evaluation assumptions permit a fair frozen benchmark.

Teacher v3 is mandatory only for the candidate-morphology task on identical
review IDs, native-input policy, label policy, and source grouping. It is not a
population-wide baseline and is not silently compared with chronology-only
tasks on different support.

The previous full-pool SSL result is especially important: end-to-end
fine-tuning nearly recovered Teacher v3 performance, but its frozen linear
probe was poor. `TWIRL-FM0` must demonstrate that useful information is
actually present in a reusable embedding rather than recoverable only through
complete supervised retraining.

## Metrics and uncertainty

### Representation health

- per-dimension variance, effective rank, leading-principal-component share,
  and duplicate/constant dimensions;
- similarity of safe augmented pairs versus unrelated hosts;
- same-host cross-sector retrieval with the query sector excluded;
- nearest-neighbor enrichment after excluding every observation of the query
  host;
- detector, sector, magnitude, missingness, and aperture predictability as
  confound probes.

### Downstream science and operations

- classification: balanced accuracy, macro F1, per-class recall, Planet-like
  average precision, precision/recall at fixed review budgets, Brier score,
  NLL, ECE, reliability curves, and the nested label-efficiency curves;
- localization: event precision-recall, event recall at fixed false events per
  source, time overlap, onset/center timing error, and false event windows per
  light curve;
- regression: MAE, robust fractional error, NLL or CRPS, bias by subgroup, and
  empirical 50%, 68%, 90%, and 95% interval coverage;
- period/recurrence: top-1 and top-K recovery, relative period error, explicit
  `P/2`, `P`, and `2P` outcomes, harmonic-aware recovery under the task's
  declared equivalence rule, and false alarms per source;
- retrieval: recall at K, mean average precision, class/source diversity, and
  same-sector-excluded results;
- anomaly review: known-reference and sealed-injection recovery plus blinded
  human-review yield/lift against a random and robust-deviation queue; an
  attractive embedding plot is not an anomaly-performance metric;
- event sensitivity: injected recovery by magnitude, depth/radius, duration,
  period, crowding, sector count, detector, and aperture disagreement, plus WD
  1856 behavior, with development and sealed injection canaries reported
  separately;
- hierarchy: every host-level metric is separated into one-sector and
  multi-sector hosts and compared with the single-visit and simple-pooling
  baselines;
- operations: examples per second, GPU-hours, peak memory, storage per host,
  inference latency, restart parity, and checkpoint reproducibility.

Confidence intervals and bootstrap comparisons cluster on
`leakage_component_id`. Adjacent windows and multiple sectors are never
counted as independent stars.

## Promotion and stopping gates

The final numeric gate table is stored in a checksum-bound configuration before
Phase 4 begins and cannot change after temporal-validation selection or sealed-
canary access. The following are the provisional minimums that Phase 3 must
either freeze or replace with a documented, scientifically justified value:

| Gate | Provisional minimum or tolerance |
| --- | --- |
| source/input isolation | zero split collisions and zero forbidden model-visible fields |
| representation rank | effective rank at least `max(16, ceil(0.10 * embedding_dim))`; zero constant dimensions |
| paired-view separation | source-clustered 95% lower confidence bound for paired-minus-unrelated similarity is greater than zero |
| seed stability | range across finalist seeds no larger than `0.03` in each primary normalized transfer metric |
| task breadth | at least three materially different task families, at least two with real-data evaluation, and at least one frozen or few-shot result |
| aggregate transfer | no primary balanced-accuracy or macro-F1 result more than `0.02` below its best applicable baseline; a promotion gain must have a source-clustered 95% interval excluding zero |
| supported-class retention | no supported-class recall decreases by more than `0.05`; the task registry fixes the minimum support before training |
| calibration retention | ECE no more than `0.02` worse than the calibrated applicable baseline, with Brier/NLL also reported |
| sealed injection retention | recovery difference at least `-0.05` overall and at least `-0.10` in every preregistered stratum with at least 20 injections |
| shortcut robustness | primary transfer metric decreases by no more than `0.05` under the declared same-sector-excluded or held-detector audit |
| WD 1856 | the benchmark receipt freezes the exact event, ephemeris, and aperture-recovery rule before Phase 4; provisionally, the periodic ranking head also places the unique host in its top 5% |

Thresholds that cannot be supported by the Phase-3 sample size remain blocking
open decisions; they are not relaxed after looking at Phase-4 results.

### Data-contract gate

- all source products, quality authorities, aliases, splits, transforms, and
  shards are checksum-bound;
- every leakage component occupies exactly one source partition;
- the model tensor allowlist contains no BLS, candidate, label, injection, ID,
  or score field;
- all six stored views, the two aperture-level error proxies, timing fields,
  and mask semantics pass the declared numerical gate; each FM0.1.x candidate
  receives only its allowlisted view subset. Missing samples within a present
  view retain explicit masks, while an observation missing a whole required
  view is excluded or routed to an explicitly named reduced-view sensitivity.

### Representation gate

- no collapse or constant embedding dimensions;
- safe paired views are distinguishable from unrelated hosts without merely
  encoding sector/detector identity;
- representation-health results are stable across seeds;
- injected short events remain detectable at the cadence/window
  representation level.

### Transfer gate

- useful transfer is demonstrated on at least three materially different
  tasks, including at least one frozen or few-shot task;
- the frozen representation beats both the random-initialized representation
  and simple baselines outside source-clustered uncertainty;
- fine-tuning provides a reproducible advantage over the same architecture
  trained from scratch;
- supported-class recall, calibration, and review-budget metrics do not hide a
  material regression behind one aggregate score.

### Prospective and science-integration gate

- the full contract is frozen before later-sector scores are opened;
- unseen-host and repeated-host later-sector results are reported separately;
- event-retention and WD 1856 checks pass the preregistered thresholds;
- a deployed head cannot bypass the transparent periodic/dip search or be
  interpreted as end-to-end completeness;
- only after all gates pass may the evaluated architecture/objective family be
  named `TWIRL-FM1`.

Stop and archive a run if a source crosses a split, a forbidden field enters
pretraining, the representation collapses, frozen probes remain near random,
fine-tuning cannot beat training the same architecture from scratch, detector
shortcuts dominate, event retention materially worsens, or the run misses a
frozen gate. Failure of a larger model to improve beyond uncertainty stops
scaling and selects the smaller Pareto-optimal candidate; it does not
invalidate or archive that smaller model.

## Compute and artifact contract

Development must respect the existing ORCD boundary:

- one H200 for smokes and most experiments;
- one H200 per FM0.1 smoke or training run; FM0.1 has no multi-GPU or GPU-array
  entitlement;
- no competition with the active Stage 1 production lane, which currently has
  priority and its own four-H200/78-CPU ceiling;
- FM0.1 input building and CPU validation may proceed in bounded CPU-only jobs.
  Every GPU submission requires a fresh admission receipt showing that it will
  not delay a runnable Stage 1 cell; otherwise the FM job remains unsubmitted
  until Stage 1 pauses or releases capacity;
- compact manifested input shards under `/orcd/data/mki_aryeh/001/twirl/`, not
  raw FFIs or a duplicate PDO source tree;
- BF16 only after a full-precision numerical smoke; gradient accumulation,
  activation checkpointing, and streaming shards as needed;
- every Slurm training allocation is at most `47:30:00`, below ORCD's 48-hour
  partition limit, and no scientific result depends on one allocation finishing
  the full run;
- atomic resumable checkpoints with model, optimizer, scheduler, scaler,
  sampler, RNG state, epoch/step, contract, and input hashes;
- write a rotating latest/previous recovery checkpoint at least every 30
  minutes and preserve milestone/epoch checkpoints; a deterministic smoke
  compares an uninterrupted run with a
  checkpoint/resume run at the same optimizer step and requires identical
  identities plus loss/state agreement under a frozen numerical tolerance;
- periodic progress, throughput, memory, and loss-component reporting.

Every run records the exact code revision, input-release manifest, split
registry, sampler, architecture, parameter count, objective weights,
augmentation contract, numerical precision, seeds, GPU type/count, total
tokens/windows/hosts seen, checkpoints, evaluation registry, and artifact
hashes.

Every `FM0.1.x` run also publishes an error ledger with five columns: observed
failure, evidence and affected subgroup, most likely experiment layer, proposed
single change, and the result expected if that diagnosis is correct. A version
may change one mechanism or one tightly coupled set. Input, split, objective,
architecture, and evaluation cannot all drift together because that would make
the improvement uninterpretable.

## Staged implementation plan

### Phase 0: design and registry

- preserve the frozen FM0.1 document/configuration receipt and maintain the
  hierarchy figure as a visual explanation;
- define the Gaia/TIC physical-source closure and ambiguity report;
- inventory accepted sectors, channel availability, external-quality evidence,
  repeated hosts, and expected storage;
- build and freeze the FM0.1 input manifest and source split under the frozen
  design while keeping long-term FM1 choices open for continued textbook
  review.

### Phase 1: ORCD input prototype

- implement a BLS-independent FM eligibility registry and strict tensor
  allowlist;
- build a small S56 derived input shard from quality-aware A2v1 HDF5;
- verify the six `{raw-relative, ADP, ADP015} x {1x1, 3x3}` views, two error
  proxies, masks, fractional-depth preservation, local and cross-sector timing,
  gaps,
  source grouping, deterministic windows, and restart parity;
- define the development and sealed injection canaries and run only the
  development canary during data-loader and tiny-model smokes;
- build compact allowlisted shards for S56--S65 near the accepted PDO data and
  transfer only those shards to ORCD; build S66--S67 shards from their retained
  ORCD checkpoints only after the same FM gate and external-quality authority
  are available.

### Phase 2: FM0.1 S56--S67 proof of concept

- use every eligible S56--S67 observation under the frozen 80/10/10
  leakage-component split, without a temporal-generalization claim;
- run `FM0.1.1`--`FM0.1.5` in order, selecting only on `poc_development` and
  opening `poc_sealed_test` once for the frozen finalist;
- compare parameter-matched TCN and Conformer encoders, then add ADP015,
  raw-relative, and same-window VICReg one declared mechanism at a time;
- create the first task registry, truth-authority audit, label-efficiency
  subsets, and hierarchy baselines;
- run the two frozen seeds on ORCD, measure throughput, report single- versus
  repeated-host behavior, and publish the FM0.1 error ledger.

### Phase 3: FM0.2 error-directed revision

- choose the smallest coherent change predicted to address the FM0.1 error
  ledger and freeze it before training;
- reuse the same training/development source split, input release, labels,
  baselines, and evaluation panel so the change is interpretable, but freeze a
  new sealed test if FM0.1's sealed result informed the redesign;
- add or retain sector/host aggregation only if the repeated-host subset can
  support the comparison;
- publish what improved, what regressed, and which diagnosis was unsupported.

### Phase 4: intermediate scaling

- repeat the selected frozen comparisons on approximately 20 accepted sectors;
- add cross-aperture learning and expanded confound probes only through named
  ablations;
- finalize and checksum the task-to-baseline matrix, metric composite, numeric
  promotion thresholds, and sealed-canary opening rule;
- estimate the exact time/storage budget for the 50-sector run.

### Phase 5: 50-plus-sector development

- freeze the accepted `S56--S_cut` manifest and the chronological split;
- complete the nested train-only 1%/10%/30%/100% data and Tiny/Base/Large
  capacity study under both fixed-compute and near-convergence comparisons;
- choose the architecture, objective, context, and checkpoint on development
  hosts plus the frozen temporal-validation block;
- run the sealed temporal prospective evaluation and sealed injection canary
  once under the frozen gate configuration.

### Phase 6: final release candidate

- if the gates pass, promote the evaluated architecture/objective family to
  `TWIRL-FM1` and train one encoder from scratch on all allowed unlabeled data
  in the frozen release;
- publish model/data cards, manifests, limitations, frozen embeddings for a
  bounded reference set, and the complete downstream benchmark;
- reserve genuinely later sectors or an external dataset for the next test;
- name the retrained checkpoint `TWIRL-FM1-all-data-candidate` until that exact
  checkpoint passes the genuinely later or external prospective test; do not
  transfer the earlier checkpoint's prospective claim to new weights.

### Phase 7: optional survey integration

- train versioned classification, regression, localization, retrieval, or
  ranking heads without altering the encoder record;
- evaluate any ranker inside the full search, vetting, candidate-merging, and
  injection-retention chain;
- keep all model scores separate from human labels and occurrence-rate
  evidence.

## Open design choices

These remain open for FM0.2 or later and are not silently changed inside the
frozen `TWIRL-FM0.1` campaign:

- the exact `S_cut` and the S94+ A2v1-equivalent input path;
- whether ambiguous multi-Gaia leakage components can be resolved safely
  rather than quarantined;
- whether the two aperture-level measurement-error proxies calibrate well
  enough for weighted losses and downstream uncertainty heads;
- window sizes other than the frozen 2,048-cadence FM0.1 context, the final
  patch width and local kernel bank, and a learned host aggregator;
- the minority high-deviation sampling fraction;
- whole-channel masking, neighboring-window prediction, cross-sector
  consistency, and aperture-residual objectives beyond the frozen FM0.1 loss;
- the Tiny/Base/Large widths after measured H200 throughput;
- the exact nested source subsets, fixed-compute cadence budget, and
  near-convergence stopping rule;
- whether absolute raw detector-unit flux, metadata late fusion, or a dual
  physical/instrument model earns inclusion in a later named design;
- the exact supported rows in the downstream task registry and the final
  numerical promotion thresholds, which must be frozen by Phase 3;
- the external or later-sector benchmark for the all-data checkpoint.

## Literature anchors

- Goswami et al. (2024),
  [*MOMENT: A Family of Open Time-series Foundation Models*](https://proceedings.mlr.press/v235/goswami24a.html),
  motivates multi-task transfer benchmarks rather than judging a model by its
  pretraining loss.
- Donoso-Oliva et al. (2023),
  [*ASTROMER: A transformer-based embedding for the representation of light curves*](https://doi.org/10.1051/0004-6361/202243928),
  is the closest published light-curve representation anchor.
- Nie et al. (2023),
  [*A Time Series is Worth 64 Words: Long-term Forecasting with Transformers*](https://openreview.net/forum?id=Jbdc0vTOcol),
  motivates masked patch pretraining, while TWIRL must retain much shorter
  event timescales before patching.
- Bardes, Ponce, and LeCun (2022),
  [*VICReg: Variance-Invariance-Covariance Regularization for Self-Supervised Learning*](https://openreview.net/forum?id=xm6YD62D1Ub),
  motivates an anti-collapse consistency baseline.
- Audenaert et al. (2025),
  [*Causal Foundation Models: Disentangling Physics from Instrument Properties*](https://arxiv.org/abs/2507.05333),
  is workshop evidence for using repeated targets under different observing
  conditions to audit or separate physical and instrumental information.
- Parker et al. (2025),
  [*AION-1: Omnimodal Foundation Model for Astronomical Sciences*](https://proceedings.neurips.cc/paper_files/paper/2025/hash/893df77404832e974b097b361ef49623-Abstract-Conference.html),
  provides a broader astronomy example in which a reusable encoder is tested
  through several downstream task families.
