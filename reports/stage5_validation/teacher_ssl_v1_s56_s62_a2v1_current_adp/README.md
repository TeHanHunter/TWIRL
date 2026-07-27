# Teacher SSL v1 / Teacher v4-SSL preregistration

**Status:** frozen development-pilot contract; no fixed-test or S63 model use.

## Question

Can self-supervised pretraining organize the existing S56--S62 periodic-search
candidates by astrophysically useful morphology, and does that representation
improve candidate enrichment after fine-tuning with the labels already in
hand?

This is a representation-learning experiment, not an assumption that every
class forms one clean cluster. Planet-like signals may occupy several sparse
neighborhoods, and nuisance variables can also produce apparent structure.

## Names and decision boundary

- `teacher_ssl_v1` is the self-supervised encoder/pretraining contract.
- **Teacher v4-SSL** is the model-facing name for the encoder after supervised
  fine-tuning.
- The pilot may earn only `companion_enrichment_candidate` status.
- It cannot replace Teacher v3, generate pseudo-labels, classify sources
  automatically, or support science claims without a later untouched labeled
  cohort.

Teacher v3 remains the frozen enrichment baseline.

## Data isolation

The immutable Teacher v3 table and TIC split registry remain authoritative.
The primary SSL pool is the frozen active candidate pool, not the full
unlabeled survey. It contains the `6,168` active, real S56--S62 development
observations (`6,054` TICs), including `3,780` real `uncertain` observations
whose labels are ignored by the SSL objective.

For each held development fold, the encoder sees only real observations from
the other four folds. Labels are not exposed to the SSL objective. The pilot
excludes:

- all fixed-test rows and TICs from constructed datasets, tensors, fitting,
  normalization, tuning, neighbors, and evaluation;
- the held development fold from its fold-specific encoder;
- all injections from the primary SSL loss;
- every Sector 63 product;
- scalar BLS/target metadata and raw chronology from SSL v1.

Injections remain available only for post-training signal-retention
diagnostics. The immutable containing files and split identities, including
HDF5 containers that also hold fixed-test groups, are read only for integrity
validation and checksum binding. Fixed-test model tensors are never
constructed.

## Representation and augmentation

SSL v1 reuses the `shape_plus_periodogram_bls` harmonic architecture without
its scalar-metadata branch during pretraining. Inputs are the seven unchanged
`P/4, P/3, P/2, P, 2P, 3P, 4P` folds, their primary/secondary local windows,
and the paired-aperture periodograms. The exported representation is the
pre-projection 128-dimensional fused embedding.

Two independent views use an event-preserving VICReg objective. Augmentations
are limited to small valid-sample noise and independent Bernoulli mask-only
dropout with probability `0.02` outside every repeated primary/secondary event
alias in the seven harmonic views. Local event windows retain their masks.
Metadata is zeroed. The implementation never unmasks rejected samples.

Cropping, phase/time warps, reversal, smoothing, flux inversion, depth
rescaling, aperture swapping, harmonic reassignment, or event
insertion/removal are forbidden.

VICReg is used instead of treating other batch members as known physical
negatives. Its variance and covariance terms provide an explicit
anti-collapse check without forcing two potentially similar astrophysical
signals apart.

## Development evaluation

Five TIC-grouped outer folds provide matched development estimates. The
fold-local SSL embedding is inductive with respect to its held fold. During
supervised fine-tuning, however, the held fold selects the best epoch, matching
the Teacher v3 protocol; those OOF metrics are therefore not an untouched
prospective estimate. The main fine-tuning policy masks `uncertain` labels
rather than mapping them to `Other`. The current Teacher v3 uncertain-masked
OOF predictions are the operational baseline on exactly matched review IDs.

The training runner exports fold-local pre-finetune embeddings, held-query
neighbors, fine-tuned OOF embeddings, and an exact-support Teacher v3
comparison. Before any pass/fail interpretation, a follow-up evaluator must
use the original embedding, never a UMAP projection, for:

- leave-TIC-out cosine-neighbor morphology enrichment;
- balanced-accuracy, macro-F1, per-class recall, and Planet-like average
  precision after supervised fine-tuning;
- unique-TIC recovery and known-positive yield at review budgets
  `50`, `100`, and `200`;
- sector, detector, magnitude, BLS-strength, missingness, and aperture-
  disagreement confound audits;
- effective-rank, leading-PC share, dimension-variance, and augmented-pair
  similarity checks.

The permutation null, same-sector exclusion, collapse diagnostics, confound
audit, review-budget yields, injection retention, WD 1856 rank, and whole-TIC
bootstrap are preregistered follow-up analyses; they are not emitted by the
initial training job. UMAP or PCA views may be made later for interpretation,
but they are not selection statistics.

## Predeclared gate

The representation must first pass all leakage and collapse checks. A
Teacher v4-SSL companion is considered worth a later replicated or
prospective evaluation only when:

1. morphology neighborhood structure exceeds a matched permutation null and
   remains present when same-sector reference neighbors are forbidden;
2. the effective embedding rank is at least 16 and augmented views are more
   similar than unrelated TICs;
3. common-support balanced accuracy and macro-F1 are no more than `0.02`
   below Teacher v3;
4. no supported class recall falls by more than `0.05`;
5. the primary `100`-TIC queue recovers at least as many known Planet-like and
   Planet-like-plus-EB examples as Teacher v3;
6. injection retention is within `0.05` overall, with no stratum of at least
   20 injections worse by more than `0.10`;
7. held-out WD 1856 remains in the top 5% of unique-TIC Planet-like scores.

If geometry passes but enrichment does not, the representation may be kept
for visualization or diversity sampling only. If leakage, collapse,
confounding, injection, or WD 1856 checks fail, it is not used for ranking.

## Scope limitation and next expansion

This native-only pilot learns morphology among already selected rank-one BLS
candidates. It is not a universal representation of TESS light curves.
Quality-aware compact ADP products for `209,077` fixed-test-excluded
S56--S62 observations (`147,430` TICs) are already staged on ORCD, but they
need a separate checksum-bound event/periodogram preprocessing contract before
they can safely broaden SSL training. That expansion follows only if the
native pilot is non-collapsed and scientifically informative.

The user approved up to four H200s for this SSL effort. The native integrity
pilot remains a one-H200 job because its small dataset and single-process input
path would not scale efficiently. Up to four H200s may be used later to run
predeclared folds or seeds in parallel, or for the broadened population after
its preprocessing contract is validated.
