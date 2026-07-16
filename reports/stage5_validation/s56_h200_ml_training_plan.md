# S56 H200 ML Training Plan

> **Historical snapshot (2026-07-01).** This report captures the compute and
> model-gating plan at that date. It does not override
> [`doc/twirl_plan.md`](../../doc/twirl_plan.md) or the
> [progress log](../../doc/twirl_progress_log.md).

Date: 2026-07-01

## Decision

Do not use H200s for the active BLS peak-table/ranker chain. The current
pipeline stage is CPU-bound and table-oriented:

1. build injected BLS peak truth table,
2. train/apply a lightweight injected-truth peak ranker,
3. render a ranker-selected real-candidate LEO queue on PDO,
4. collect human labels.

Use H200s only after those gates produce a labeled or pseudo-labeled training
set that actually needs tensor/light-curve modeling.

Routine H200 policy remains:

- default first GPU experiment: `1` H200,
- routine development cap: `2` H200s,
- no full `8`-H200 job without explicit approval for that specific run.

## Current Status

As of `2026-07-01`, the all-host real-candidate human-triage queue is ready on
PDO:

- source: the coverage-first all-host injection product and injected-truth peak
  ranker;
- real ephemerides selected: `57,204` rows across `19,068` S56 targets;
- human queue: `1,000` shuffled ranker-selected real candidates;
- LEO queue verification: `1,000` PDFs, `0` LEO metric errors, `0` LEO plot
  errors (`FA=994`, `FP=4`, `PC=2`);
- app: `pdogpu1:5007`, writing labels to
  `s56_allhost_ranker_selected_real_leo_queue_pdo/human_labels_vetted.csv`.

The ORCD balanced-grid `20k` CPU peak-table chain is still running as Slurm job
`16934436`; ranker/gate job `16934439` and real-candidate apply job `16934507`
are dependency-pending. No H200 job is currently needed for this chain.

## Current Compute Split

### CPU Comparison

We are using ORCD for this stage for operational reasons, not because its CPU
cores are known to be dramatically faster per core than PDO. It is better here
because the private partition provides dedicated Slurm-managed CPU nodes, large
memory, durable project storage, and a cleaner handoff path to the H200 node.
PDO remains the right home for Stage 1 production and LEO rendering because it
has the native TGLC/HLSP/LEO data environment.

### CPU Now

Use ORCD CPU nodes for:

- full `20k` injected BLS peak table,
- post-BLS recall/gate/failure-mode audits,
- injected-truth peak ranker,
- real S56 BLS peak-table scoring,
- ORCD-to-PDO selected-ephemeris handoff.

Use PDO for:

- raw Stage 1 / TGLC / HLSP operations,
- LEO rendering, because PDO has the S56 HLSP tree and working LEO-Vetter
  environment,
- browser-based human triage.

### H200 Later

Use H200s for:

- multi-aperture light-curve tensor models,
- representation learning / embeddings from compact S56 light curves,
- AstroNet-like folded-light-curve branches adapted to WD hosts,
- later GPU-native search or matched filters only after the CPU baseline is
  measured.

Do **not** use H200s for the current logistic peak ranker or the engineered
candidate-level self-training model. Those are intentionally dependency-light
CPU baselines.

## Gates Before First H200 Training

Do not submit a model-training H200 job until all of these are true:

1. ORCD CPU peak-table chain is complete and verified:
   - `peak_training_orcd_full/s56_20k_injection_bls_peaks_orcd_full_verification.json`
   - `peak_training_gate_orcd/host_coverage/summary.json`
   - `peak_training_gate_orcd/summary.json`
   - `peak_training_gate_orcd/failure_modes/summary.json`
   - `peak_ranker_orcd/summary.json`
2. Real S56 selected ephemerides are verified:
   - `s56_ranker_selected_real_candidates_orcd/selected_ephemerides_verification.json`
3. PDO LEO queue exists and passes verification:
   - `s56_allhost_ranker_selected_real_leo_queue_pdo/verification.json`
4. Human-label readiness audit passes at least one training gate:
   - visibility smoke is enough for a visibility-only diagnostic,
   - object-teacher training requires real-candidate false-positive classes,
   - do not train object taxonomy from injection-only labels.
5. Compact LC tensors are staged and documented:
   - `s56_twirlfs_v2_lc_export_pdo.h5`
   - `injected_lightcurves.h5`
   - manifest/label CSVs with injection truth and train/validation/test split.
   - data-shape smoke accepts candidate-centered windows from the compact export.

## First H200 Experiment

The first GPU experiment should be a small, bounded, one-H200 smoke:

**Goal:** learn whether a multi-aperture light-curve model adds information
beyond BLS peak features and LEO metrics.

**Inputs:**

- compact S56 LC export,
- `20k` injected light curves with hidden truth,
- verified BLS peak/ranker outputs,
- human-labeled real-candidate rows after LEO triage.

**Model shape:**

- three compact-export light-curve channels initially: `DET_FLUX_SML`,
  `DET_FLUX`, and `DET_FLUX_LAG`,
- candidate-centered/folded windows from ranker-selected period/t0/duration,
- optional scalar side features: Tmag, BLS SDE/depth/duration, aperture,
  cadence-cleaning diagnostics, LEO class/metrics when available,
- target labels kept separate:
  - injected signal visibility/recovery,
  - real human taxonomy,
  - pseudo-labels.

**Acceptance for continuing GPU work:**

- beats the engineered-feature baseline on held-out injected recovery *and*
  does not degrade bright/short-period regimes,
- reports metrics by Tmag, period, radius/depth, duration, and source kind,
- calibration and confusion matrix are documented,
- every pseudo-label carries model version, confidence, margin, and source.

## First Commands After Gates Pass

Current human-triage command:

```bash
ssh -L 5007:localhost:5007 pdogpu1
```

Open `http://127.0.0.1:5007/` and label the all-host real-candidate queue.

After labels exist on PDO, build the label summary, teacher table, priority
list, and readiness audit. A detached PDO monitor is already running as
`twirl-s56-allhost-label-audit`; it watches for label-count changes and runs
this audit automatically. A local monitor `twirl-pdo-label-sync` then pulls the
small label/audit products back to this checkout and refreshes
`s56_ml_handoff_readiness/`. The manual PDO audit command is:

```bash
cd /pdo/users/tehan/TWIRL
bash scripts/stage5_validation/run_s56_allhost_real_label_audit_pdo.sh
```

Use the readiness audit as the gate:

- `injection_visibility_smoke`: useful for checking whether humans see signals,
  but not enough for real object taxonomy.
- `binary_teacher_smoke`: enough for a first engineered-feature smoke.
- `object_teacher_training`: minimum gate before training a real taxonomy
  teacher/student model.

After ORCD apply finishes:

```bash
scripts/orcd/run_s56_orcd_pilot.sh --run sync-apply
```

Or wait for apply, sync automatically, and render a bounded first PDO LEO queue:

```bash
scripts/orcd/run_s56_orcd_pilot.sh --run monitor-apply-leo-smoke
```

On PDO, render a small first LEO queue from ORCD-selected ephemerides:

```bash
cd /pdo/users/tehan/TWIRL
TWIRL_SELECTED_EPHEMERIDES=reports/stage5_validation/s56_ranker_selected_real_candidates_orcd/selected_ephemerides.csv \
TWIRL_REVIEW_N_REAL=50 \
TWIRL_MAX_LEO_REPORTS=50 \
  bash scripts/stage5_validation/run_s56_ranker_selected_real_leo_pdo.sh
```

After human labels are added and readiness passes, run the existing
engineered-feature teacher/student baseline first:

```bash
PYTHONPATH=src python scripts/stage5_validation/run_self_training_triage.py \
  --candidate-table reports/stage5_validation/s56_allhost_ranker_selected_real_leo_queue_pdo/review_queue.csv \
  --labels reports/stage5_validation/s56_allhost_ranker_selected_real_leo_queue_pdo/human_training_table/teacher_labeled_rows.csv \
  --config configs/detection/self_training_s56_allhost_real.yaml \
  --out-dir reports/stage5_validation/self_training_s56_allhost_ranker_selected
```

The runner writes Parquet when the environment has a Parquet engine and falls
back to CSV otherwise, so the smoke can run in the PDO QLP environment as well
as on ORCD/local.

The tensor data-shape smoke now passes on PDO for the all-host real queue
(`20/20` candidates, shape `(20, 3, 65)`) and on ORCD/H200 (`128/128`
candidates, shape `(128, 3, 257)`, Slurm job `16939737`, one H200 for
`12 s`). The follow-up torch-required smoke also passes with the separate
`twirl-s56-torch` env: Slurm job `16942919` used one H200 for `47 s`, imported
`torch 2.11.0+cu128`, reported `cuda_available=true`, `device_count=1`, and
moved the tensor smoke batch onto `cuda:0`.

Do not treat this as model training. The GPU infrastructure is ready for a
small real-label training smoke, but the scientific gate remains human-label
readiness. The trainer [script](../../scripts/stage5_validation/train_candidate_tensor_classifier_smoke.py)
and Slurm wrapper [script](../../scripts/orcd/slurm_s56_tensor_train_smoke_h200.sbatch)
are already in place. Synthetic-label environment smoke `16944030` used one
H200 for `6 s`, trained on `128` tensor rows with train/validation/test split
`76/26/26`, and confirmed the training loop runs on `cuda:0`; this is not a
scientific model.

Use the combined local handoff audit to check all gates before submitting the
real-label H200 smoke:

```bash
python scripts/stage5_validation/audit_s56_ml_handoff_readiness.py
```

The current expected pre-label state is `ready_for_human_triage=true` but
`ready_to_submit_real_h200_training_smoke=false`.

After teacher labels exist and pass readiness, the real-label H200 path needs
one small staging step followed by the one-H200 training smoke:

```bash
scripts/orcd/run_s56_orcd_pilot.sh --run stage-labels
scripts/orcd/run_s56_orcd_pilot.sh --run h200-tensor-train-smoke
```

`stage-labels` copies only queue metadata, human labels, teacher rows, label
priorities, and readiness summaries into the ORCD checkout. It does not move
raw light curves or large injection HDF5 files, and it refuses real staging
unless the readiness audit has passed a teacher-training gate.

For the coverage-first all-host product, keep the same CPU-first discipline:
stage with `scripts/orcd/run_s56_orcd_pilot.sh --run stage-allhost`, run
`allhost-peak` only on CPU if needed, and run `allhost-ranker`/`allhost-apply`
only after the all-host peak verifier exists. Do not move this to H200 training
until the selected real queue and human-label readiness gates pass.

## What Not To Do

- Do not train on injection truth as if it were real object taxonomy.
- Do not treat LEO `PC/FP/FA` as ground truth.
- Do not train object taxonomy before real false-positive labels exist.
- Do not request all H200s for a first model.
- Do not move raw TGLC/TICA/source-pickle trees to ORCD for this first pass.
