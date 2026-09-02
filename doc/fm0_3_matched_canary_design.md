# TWIRL-FM0.3 matched canary design

Status: frozen matched canary, authorized but not started (`v1`,
`2026-09-02`).

This version is immutable: any scientific change requires a new versioned
document, config, and freeze receipt.

This document is the scientific contract for a bounded, parameter-matched
comparison of the cadence-preserving `TWIRL-FM0.3.1` TCN and
`TWIRL-FM0.3.2` Conformer. It is separate from
`foundation_model_design.md` so earlier model freezes remain byte-exact.

## Scientific question

The canary asks whether a Conformer or TCN better learns reusable local
light-curve representations when both receive the same short native-cadence
context and optimize the same self-supervised objective. It is not a BLS or
blind-search experiment. The intended downstream use is representation reuse
by later classifiers and triage models.

The architecture is the only intended difference between the two arms. Both
use ADP `1x1` and `3x3`, seed `560067`, and the same data roles, sampler,
masking, loss, optimizer, precision, batch size, and stopping schedule.

## Native-cadence contract

- Each example contains exactly `128` native TESS `200 s` cadences
  (`7.111111` hours).
- Each cadence remains one encoder token; patch stride is exactly one.
- Cadences may not be merged, averaged, pooled, resampled, or temporally
  downsampled before or inside either encoder.
- The optimized representation is cadence-indexed `h_cadence`; no pooled
  window representation enters the loss.
- In both arms, the reconstruction head consumes the contextual cadence
  sequence. The legacy Conformer stride-four stem skip is disabled at stride
  one, so the two arms do not optimize different reconstruction pathways.
- Mask spans are `1--4` native cadences and target `15%` of the window.

The small context is meant to retain minute-scale transit morphology with
local baseline on each side. It does not require multiple periodic events in
one window.

## Frozen objective and optimization

The objective is first-mask Huber reconstruction plus position-centered
cadence VICReg between two independently masked views. The second mask is used
only to form the paired cadence representation and is never a reconstruction
target. VICReg uses weights `25/25/1` and total weight `0.0002`.

Both arms use AdamW, learning rate `0.0003`, weight decay `0.01`, linear
warmup for `1,000` optimizer steps, a cosine schedule with invariant horizon
`20,000`, and an effective batch of `64` unique windows. Runtime stops do not
change that scientific contract.

## Frozen input

The input is the immutable identity-only S56--S64 plus S66--S77 composite
release. Gradients use only role `fm03_train`; S77 role `temporal_holdout` and
all sealed rows remain unavailable to training. S65 is absent. The composite
receipt and both identity tables are SHA-256-bound in the YAML config and
design-freeze receipt. Existing shards are read in place and are not rebuilt.

## Pre-checkpoint evaluation freeze

Before either arm sees real training data, one shared evaluation bundle must
be frozen in two stages. The first stage selects identities without opening
payloads. The second opens exactly those `1,440` bound shards once, verifies
their SHA-256 values, and freezes one deterministic, unpadded `128`-cadence
crop wholly inside one segment for every identity. Each crop must have at
least `103/128` cadences jointly valid in time and both ADP views, with all
indices `60--68` jointly valid. An ineligible preselected identity fails the
freeze; it is not silently replaced.

The identity stage must descend from temporal-panel receipt
`78c370e10c556472c5997c20cfe95207a0b334bafe7f024bf7ba4fc7ec4de624`.
That content-addressed receipt transitively binds the panel and sector-binding
tables. A different otherwise-valid temporal panel is not part of this
matched comparison.

The same bundle is used for both architectures. It fixes:

- probe training on S66--S71, diagnostic validation on S72--S74, and fresh
  S77 testing, each with `240` repeated-source and `240` new-source component
  pairs;
- one clean/injected pair per identity, with a symmetric cadence-center
  trapezoid centered at zero-based index `64`, durations `1`, `3`, and `9`
  cadences, and fractional depths `0.01`, `0.03`, `0.10`, and `0.30`;
- injection into ADP `1x1` and `3x3` flux only, using
  `(1 + flux) * (1 - profile) - 1`, leaving every other tensor and flux view
  bitwise unchanged;
- unmasked inference, all `128` cadence tokens retained, and a linear logit
  from `h_cadence[64]` only, with no token/window pooling and no duration,
  depth, support, period, BLS, or search-score feature;
- separate raw, quality-only, step-0, and step-2,000 feature arms; a full-batch
  Adam linear probe with train-only dimensionwise standardization, `400`
  epochs, learning rate `0.02`, L2 weight `0.001`, and seed `560306`, with no
  hyperparameter tuning; and
- sample AUROC with a `1,000`-replicate, seed-`560305`, component-pair
  clustered bootstrap using identical resamples for all arms and paired
  deltas.

The payload freeze itself injects no event, loads no checkpoint, fits no
probe, and computes no metric. Its immutable receipt and schedule hashes are
bound into every real-data run contract before the first checkpoint exists.

## Bounded execution

For each architecture, the exact-code sequence is:

1. an `8`-step synthetic FP32 smoke on one H200;
2. completion and validation of the shared payload-screened evaluation
   freeze;
3. a fresh real-data BF16 invocation stopping at step `64` on one H200; and
4. after the step-64 checkpoint validates, an exact resume stopping at step
   `2,000` on one H200.

Immutable milestones are steps `0`, `64`, and `2,000`. No invocation may pass
step `2,000` under this freeze. Before the first real-data submission for an
arm, the config, design, freeze receipt, clean Git revision, immutable
composite hashes, shared evaluation-plan receipt, and that arm's step-8 smoke
must pass. The two arms may run independently, but neither may change the
other's resources or contract.

## Claim boundary

This freeze authorizes only the two matched canaries through step `2,000`.
It does not authorize a longer run, a second seed, sealed-test access,
production use, an architecture promotion, or a foundation-model claim.
Architecture comparison and any downstream transfer decision require a
separately frozen, source-aware evaluation applied identically to both arms.
