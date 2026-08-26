# TWIRL-FM0.2.1 step-2000 development evaluation

## Result

The exact-revision, seed-`560067` canary reached step `2000`, passed strict
post-validation, and passed the frozen CPU-only evaluator's mechanical and
provenance checks. The development representation changed substantially from
steps `500` and `1000`: `z_window` effective rank reached `39.40`, its constant
dimension count remained zero, and both source-clustered separation lower
bounds are now positive.

This is a mixed result, not a canary pass. Four of the five predeclared
development representation criteria are met, while masked reconstruction
remains `0.109%` worse than the frozen zero predictor:

| Step-2000 criterion | Observed | Frozen criterion | Status |
| --- | ---: | ---: | --- |
| `z_window` effective rank | `39.3994` | `>= 26` | meets |
| `z_window` zero/constant dimensions | `0` | `<= 0` | meets |
| paired-minus-unrelated 95% lower bound | `0.05507` | `> 0` | meets |
| trained-minus-random 95% lower bound | `0.03311` | `> 0` | meets |
| development masked Huber | `0.004869893` | `<= 0.004864579` | **misses** |
| sealed-test access count | `0` | `== 0` | meets |

The evaluator receipt explicitly records
`scientific_go_no_go_applied: false`. Development event retention was not run,
the sealed test remained closed, and no production, science, or
foundation-model claim is authorized. The formal step-2000 canary gate is
therefore incomplete and must not be described as passed.

## Representation trajectory

All three milestones use the same `256` development leakage components,
`383` visits, observation-sector authority, evaluator revision, and source-
clustered bootstrap definitions.

| Diagnostic | Step 500 | Step 1000 | Step 2000 |
| --- | ---: | ---: | ---: |
| `h_window` effective rank | `6.28` | `18.46` | `33.80` |
| `z_window` effective rank | `6.45` | `19.43` | `39.40` |
| `h_window` constant dimensions | `1` | `1` | `1` |
| `z_window` constant dimensions | `0` | `0` | `0` |
| `z` paired-minus-unrelated | `0.01618` | `0.02774` | `0.05715` |
| `z` paired-minus-unrelated 95% interval | `[0.01504, 0.01725]` | `[0.02654, 0.02899]` | `[0.05507, 0.05924]` |
| `z` trained-minus-random | `-0.00397` | `0.00759` | `0.03700` |
| `z` trained-minus-random 95% interval | `[-0.00645, -0.00163]` | `[0.00485, 0.01035]` | `[0.03311, 0.04058]` |
| masked Huber / zero | `1.00484` | `1.00385` | `1.00109` |
| `z` cross-sector retrieval, clustered mean | `0.01597` | `0.02778` | `0.02083` |
| projection effective rank | `154.45` | `152.24` | `147.75` |
| projection condition number | `244.3` | `415.4` | `368.5` |

The main causal hypothesis is supported in a narrow sense: directly optimizing
`z_window` with the scale-matched VICReg term repaired the low-rank and
trained-versus-random separation failures by step `2000`. It did not produce a
clean win across the whole development panel:

- masked Huber improved monotonically but still misses the frozen zero-predictor
  criterion by `5.31e-6` in absolute loss;
- `h_window` retains one constant dimension;
- query-sector-excluded retrieval remains weak, with a step-2000 clustered
  mean of `0.02083` and a 95% interval `[0, 0.05556]`; its interval includes
  zero and its mean is below the matched-random (`0.03125`) and robust-scalar
  (`0.04167`) controls;
- paired separation improves, but the step-2000 `z` value (`0.05715`) remains
  below the robust-scalar (`0.09863`) and train-fit PCA (`0.60007`) baselines;
- the learned `256 x 256` projection remains full numerical rank, although its
  effective rank falls modestly across milestones. The representation result
  is therefore not explained by a singular projection matrix.

The [representation trajectory PNG](fm0_2_step2000_representation_trajectory.png)
and [PDF](fm0_2_step2000_representation_trajectory.pdf) show these comparisons.
The exact plotted values are in
[`step500_step1000_step2000_metrics.csv`](step500_step1000_step2000_metrics.csv),
and the frozen criteria are tabulated in
[`step2000_representation_criteria.csv`](step2000_representation_criteria.csv).

## Component-loss trajectory

The [training-components PNG](fm0_2_step2000_training_components.png) and
[PDF](fm0_2_step2000_training_components.pdf) plot every immutable history row
without smoothing. The optimized reconstruction is the first-mask term; the
second-mask and mean reconstructions are diagnostics. The final-100-step
medians are:

| Component | Through step 500 | Through step 1000 | Through step 2000 |
| --- | ---: | ---: | ---: |
| reconstruction first / optimized | `0.0050069` | `0.0048373` | `0.0047579` |
| reconstruction second / diagnostic | `0.0049801` | `0.0048281` | `0.0047516` |
| reconstruction mean / diagnostic | `0.0050072` | `0.0048377` | `0.0047550` |
| VICReg invariance | `7.65e-4` | `8.09e-4` | `1.12e-3` |
| VICReg variance | `0.6256` | `0.4886` | `0.3539` |
| VICReg covariance | `2.2144` | `3.2200` | `3.9518` |
| VICReg aggregate, unweighted | `17.9575` | `15.4545` | `12.9157` |
| VICReg weighted contribution | `0.0035915` | `0.0030909` | `0.0025831` |
| total loss | `0.0085859` | `0.0079343` | `0.0073228` |
| projection-gradient norm | `4.15e-4` | `6.74e-4` | `7.53e-4` |
| learning rate | `1.3545e-4` | `2.8545e-4` | `2.9815e-4` |

The projection gradient remains finite and nonzero throughout the preserved
history. The weighted VICReg contribution decreases while representation rank
rises, consistent with the anti-collapse objective acting on `z_window` rather
than the projection remaining untrained. The loss minimum is not a selection
metric.

The checksum-bound loss extract contains contiguous rows `1--2000`; rows
`1--1000` are value-identical to the earlier step-1000 extract. The
complete values are in [`training_trace.csv`](training_trace.csv).

## Provenance and access boundary

- training revision: `ddf442aafb8f62966e549e2287abad3474dd556a`
- evaluator revision: `83816d07975eebe3825d76dfe7096d22b70376f5`
- step-2000 checkpoint SHA-256:
  `976b5053c857c38b9fbf7a35c9d0605f0023318b2c9bd37a88995992e5aa7bd2`
- representation evaluation SHA-256:
  `3c802741744999fe553e71bc2dccfbef6c309a29a28c30dfaac693a798bb3718`
- post-validation receipt SHA-256:
  `80d678c4621a22258c5d0e8f0be6f7e08ffa5bf93841c315abf42e9bb006b110`
- loss-history extract SHA-256:
  `7d9eb81700f0185a725f814dc693ef8969eeb40953c4c42266b6d803d5909dfe`
- training job `21303060`: completed `0:0` in `1:30:23`, one H200,
  four CPUs, `32 GiB`, BF16;
- strict post-validation job `21327505`: completed `0:0`, CPU-only;
- representation evaluator job `21327655`: completed `0:0` in `3m19s`,
  four CPUs, `16 GiB`, no GPU;
- compact loss-history extraction job `21331169`: completed `0:0` in `1s`,
  four CPUs, `16 GiB`, explicitly excluding the H200 node. Attempts `21331054`
  and `21331101` failed before producing scientific output because of a
  node-local path and then a corrected checkpoint-key assumption.

All evaluation artifacts record zero sealed-test accesses. The matched random
encoder control is not the exact same-seed initialization checkpoint.

## Remaining evidence and recommendation

The immutable seed-`560067` initialization checkpoint remains preserved at
SHA-256
`92463070381486f0c6053c190d62da8d0c5c0d31be8072e2bdcd677329ac792c`,
but the frozen step-0 representation evaluation has not been executed. The
trajectory evidence bundle is therefore incomplete until a checksum-bound,
CPU-only step-0 evaluator path is reviewed and run, or the evidence requirement
is explicitly changed.

Do not extend this run, start `FM0.2.2`, use seed `560068`, open the sealed
test, or make a production/science/foundation-model claim from this result.
The useful conclusion is narrower: the added objective repairs the targeted
rank/separation pathology, but the current canary still misses its
reconstruction criterion and does not demonstrate cross-sector retrieval.

Continuation on another machine should start from the scoped
[new-Mac handoff](HANDOFF.md).
