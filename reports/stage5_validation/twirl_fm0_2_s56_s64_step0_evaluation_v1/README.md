# TWIRL-FM0.2.1 step-0 completion and four-point trajectory

## Outcome

The exact seed-`560067` initialization was evaluated on the unchanged
development panel after the step-2000 result. ORCD CPU job `21340937`
completed in `3m46s` on `node4702` with four CPUs, `16 GiB`, no GPU, and empty
stderr. The checksum-bound receipt passed and records zero sealed-test access,
no development event-retention run, no scientific go/no-go application, and
no production or foundation-model authorization.

The complete same-seed trajectory establishes that FM0.2.1 learned a much
higher-rank window representation. It does **not** establish a generally
useful embedding or cross-sector improvement:

| `z_window` diagnostic | Step 0 | Step 500 | Step 1000 | Step 2000 |
| --- | ---: | ---: | ---: | ---: |
| effective rank | `3.027` | `6.447` | `19.435` | `39.399` |
| paired-minus-unrelated mean | `0.02775` | `0.01618` | `0.02774` | `0.05715` |
| checkpoint-minus-seed-0-control mean | `0.00760` | `-0.00397` | `0.00759` | `0.03700` |
| masked Huber / zero predictor | `5.01856` | `1.00484` | `1.00385` | `1.00109` |
| cross-sector retrieval, clustered mean | `0.03380` | `0.01597` | `0.02778` | `0.02083` |

From exact initialization to step 2000, `z_window` effective rank rises by
`36.37`, or about `13.0x`, and paired-minus-unrelated separation rises by
`0.02940`. Masked Huber falls by `80.05%`, but its final value
`0.004869893` remains `5.31e-6`, or `0.109%`, worse than the frozen zero
predictor. These results support the narrow causal hypothesis that directly
optimizing `z_window` with same-window VICReg repairs the measured rank
pathology.

The [trajectory PNG](fm0_2_step0_step2000_representation_trajectory.png) and
[PDF](fm0_2_step0_step2000_representation_trajectory.pdf) show the complete
comparison. Exact plotted values are in
[the four-point metrics table](step0_step500_step1000_step2000_metrics.csv).

## What the initialization changes in the interpretation

The exact initialization already exceeds the evaluator's separate seed-`0`
random encoder on paired separation: the `z_window` difference is `0.00760`
with a source-clustered 95% interval `[0.00590, 0.00937]`. Therefore the
positive step-2000 checkpoint-minus-random interval is partly sensitive to
initialization seed. It cannot serve alone as evidence that training created
the separation. Future comparisons should use each checkpoint's exact
same-seed initialization as the primary causal baseline and add multiple seeds
for stability. A source-paired bootstrap of the step-0-to-checkpoint change
would be preferable to comparing two marginal intervals.

Cross-sector retrieval does not improve over exact initialization. The
`z_window` clustered point estimate changes from `0.03380` at step 0 to
`0.02083` at step 2000, a descriptive decrease of `38.4%`; the corresponding
intervals are `[0.01087, 0.06482]` and `[0, 0.05556]`. Because the evaluator
does not publish a checkpoint-paired interval for that change, this is not a
claim of a statistically resolved regression. At step 2000, the retrieval
point estimate remains below the seed-`0` random (`0.03125`) and robust-scalar
(`0.04167`) controls, while intervals overlap.

The final `z_window` paired separation (`0.05715`) also remains below the
robust-scalar (`0.09863`) and train-fit PCA (`0.60007`) baselines. Effective
rank values are dimension-dependent: `39.40/256` passes the frozen absolute
gate, while the PCA baseline is `10.81/32`. The learned projection stays full
numerical rank, but its effective rank changes from `155.33` at initialization
to `147.75` at step 2000. The result is evidence that encoder activations
changed, not that the linear projection became more diverse. `h_window` still
has one constant dimension at steps 500, 1000, and 2000.

## Would approximately 20 sectors help?

Potentially, but the current repository does not contain approximately 20
FM-ready sectors. The operational Stage-1 footprint through S74 is 19 sectors:
S56--S65 are accepted, S66--S72 are retained deferred checkpoints, S73 is
active, and S74 is resident/queued. Only S56--S64 form the frozen,
checksum-bound FM release used here. S65 still needs an FM-specific manifest
and admission gates; S66 and later products are not implicitly admissible.

More accepted sectors can improve sky/instrument diversity, enlarge the
repeated-host cohort, and support a genuine later-sector panel. They do not
change the present mechanism: FM0.2.1 pairs two masks of the **same window**
and has no cross-sector consistency loss or host aggregator. Adding independent
windows therefore does not directly train the retrieval property that failed
to improve. The step-2000 canary also saw `128,000` window draws against
`273,871` raw training-authority rows, about `0.47` draws per row before
host-first resampling effects. Expanding the corpus at the same step budget
would broaden coverage while reducing average exposure further.

The smallest defensible next sequence is:

1. Complete a label-free admission inventory for each later accepted sector,
   including FM channel/numerical gates, immutable manifests, unique Gaia
   components, repeated-host counts, sector multiplicity, and detector/camera
   transitions.
2. Once at least five genuinely later accepted sectors pass those gates,
   freeze them as a development-only temporal panel and evaluate the existing
   step-2000 checkpoint zero-shot. Report new hosts and repeated hosts
   separately and retain the random, robust-scalar, and PCA controls.
3. If a data-scaling training test is then justified, compare the current
   nine-sector release with an expanded-sector release under the same TCN,
   objective, exact initialization, split rule, window-draw budget, and held-
   out temporal panel. Do not change the objective in that comparison.
4. Test cross-sector consistency or a stable-host/visit head only as a later,
   separately named ablation if the repeated-host inventory supports it. Keep
   local event/morphology information separate so forced host invariance does
   not erase real time-variable signals.

This preserves the Phase-3 causal result and makes Phase-4 scaling an actual
data experiment. No further FM training, event-retention run, formal gate,
second seed, Conformer, or sealed-test access is authorized by this report.

## Provenance

- orchestration revision:
  `393ee8f9f1372ca2e017d4fea59f9593bc96e4fc`
- immutable training revision:
  `ddf442aafb8f62966e549e2287abad3474dd556a`
- immutable evaluator revision:
  `83816d07975eebe3825d76dfe7096d22b70376f5`
- initialization checkpoint SHA-256:
  `92463070381486f0c6053c190d62da8d0c5c0d31be8072e2bdcd677329ac792c`
- step-0 evaluation SHA-256:
  `bc75b96271487851162f542682501f30fba1bcaa40b425333ad1767c32da1fe6`
- step-0 receipt SHA-256:
  `5e489c8abf0ef7c70c815889701a8dfc1a7ad31eed7b5a4366158ec68d61e29c`
- step-2000 evaluation SHA-256:
  `3c802741744999fe553e71bc2dccfbef6c309a29a28c30dfaac693a798bb3718`
- development panel: `256` leakage components and `383` visits under
  observation-sector authority
  `bdd3e8039b312aeb21557662226a8a11b761dd1da031742be42ee9b0d1c6edc5`.

The compact [step-0 evaluation](receipts/step0_representation_health.json),
[receipt](receipts/step0_representation_health.receipt.json), and
[receipt sidecar](receipts/step0_representation_health.receipt.json.sha256)
are tracked here. The checkpoints remain on ORCD.
