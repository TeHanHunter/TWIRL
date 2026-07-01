# S56 Peak-Ranker Training Plan

Date: 2026-06-26

## Core Decision

Use injected truth to train the **search/ranking layer**, not the final object
taxonomy layer.

The current human-triage issue is structural: when BLS folds the light curve at
the wrong period, LEO and a human reviewer are vetting the wrong ephemeris. The
first ML problem should therefore be:

> Given all BLS peaks for an injected light curve, rank the peak that corresponds
> to the injected ephemeris above aliases, systematics, and unrelated high-SDE
> structure.

Human labels remain essential for the second problem:

> Given a candidate ephemeris from real data, decide whether the signal is a
> planet-like occultation, PCEB/eclipsing binary, variability, contaminant,
> instrumental/systematic, or uncertain.

Keep these labels separate:

- `is_injected_signal_peak`: injection-truth label for a BLS peak.
- `match_kind`: `exact`, `harmonic`, `mismatch`, `no_peak`, or `error`.
- `human_label`: reviewer taxonomy for the displayed candidate.
- `leo_class`: LEO-Vetter baseline output, not a ground-truth label.

The human-label taxonomy and teacher/audit/exclude semantics now live in
`src/twirl/vetting/label_schema.py`. Browser labels, joined human-training
tables, and readiness audits should consume that schema rather than duplicating
label lists.

## Immediate Artifact

New script:

```bash
scripts/stage5_validation/build_injection_peak_training_table.py
```

It reruns BLS on an injection HDF5 and writes one row per BLS peak per aperture,
with explicit exact/harmonic/mismatch truth labels. This is the missing table
for training and auditing a BLS peak ranker.

Training/evaluation script:

```bash
scripts/stage5_validation/train_injection_peak_ranker.py
```

This trains a lightweight logistic peak ranker using only BLS/light-curve
features available for real candidates. Truth columns are used only for the
binary target and recall@K audit. Outputs are a model, scored peak table,
selected ephemerides, coefficients, and summary JSON.

Application script:

```bash
scripts/stage5_validation/apply_injection_peak_ranker.py
```

This scores real BLS candidate tables and selects the top-N ephemerides per TIC
for LEO/human review. It tolerates older candidate tables that do not yet have
the cadence-cleaning breakdown, falling back to `n_cad_kept` and `dropout_frac`.

The BLS candidate schema now records cadence-cleaning diagnostics:

- `n_cad_quality`: cadences kept after `QUALITY==0` and finite-flux masking.
- `n_cad_edge_trimmed`: good cadences removed by optional orbit-edge trimming.
- `n_cad_sigma_clipped`: good cadences removed by upper-tail spike clipping.
- `quality_dropout_frac`: dropout from quality/finite-flux masking alone.
- `dropout_frac`: total dropout after all BLS cleaning.

These are available as ranker/audit features for real candidates and injected
rows. The current one-shot PDO peak-table job was launched before these columns
were added; use the restartable chunked launcher for a final cadence-diagnostic
peak table if the diagnostics matter for a training run.

## PDO Smoke / 20k Run

Run this on PDO because the current `20k` injection HDF5 already lives there.
Use the restartable chunked launcher so the output table is current-schema and
self-verified:

```bash
cd /pdo/users/tehan/TWIRL
export LD_LIBRARY_PATH=/pdo/app/anaconda/anaconda2-4.4.0/lib:/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}
export PYTHONPATH=src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

TWIRL_PEAK_CHUNK_JOBS=2 \
TWIRL_PEAK_WORKERS=8 \
TWIRL_PEAK_CHUNK_SIZE=250 \
  bash scripts/stage5_validation/run_s56_20k_peak_training_chunked.sh
```

The active PDO handoff uses tmux session `twirl-s56-20k-peak-training` and log:

```text
/pdo/users/tehan/TWIRL/logs/s56_20k_peak_training_chunked_pdo.log
```

First audit the self-verification JSON and summary JSON, especially
`recall_at_1`, `recall_at_2`, `recall_at_5`, and `recall_at_20`.

Interpretation:

- High `recall_at_20` but low `recall_at_1`: train a peak-ranker and send the
  best few ephemerides to LEO/humans.
- Low `recall_at_20`: BLS is not placing the injected signal in the candidate
  list, so tune/replace the search statistic before scaling human triage.
- If high-observability misses concentrate at short periods while the selected
  BLS peaks are long-period aliases/systematics, test a separate short-period
  branch or per-period-range peak quota. The partial current-schema smoke found
  this exact pattern: a `p_max=2 d` branch recovered `15/43` high-observability
  not-in-top-20 misses, including `11/43` at top-1.

The peak-table builder now supports branch provenance:

```bash
python scripts/stage5_validation/build_injection_peak_training_table.py \
  --search-branch short_pmax2 \
  --p-max-cap-d 2.0 \
  --n-periods 100000 \
  --out-table short_pmax2_peaks.csv
```

Rows include `search_branch` and BLS configuration columns. If standard and
short branch tables are merged for ranker training or recall audit, use the
branch-aware merge key:

```bash
python scripts/stage5_validation/merge_injection_peak_training_chunks.py \
  --chunk-root <branch-chunks> \
  --out-table merged_standard_short.csv \
  --id-columns search_branch,injection_id
```

Then compare the standard and short branch explicitly:

```bash
bash scripts/stage5_validation/run_s56_short_branch_compare.sh
```

This one-command runner performs the standard failure-mode audit, builds the
short branch on high-observability misses, runs the branch gate, and writes the
comparison outputs. If the short branch table has already been built, compare
without rerunning BLS:

```bash
TWIRL_SHORT_BRANCH_SKIP_BUILD=1 \
TWIRL_SHORT_BRANCH_PEAK_TABLE=short_pmax2_peaks.csv \
  bash scripts/stage5_validation/run_s56_short_branch_compare.sh
```

The lower-level comparison command is:

```bash
python scripts/stage5_validation/compare_injection_bls_branches.py \
  --baseline-peak-table standard_peaks.csv \
  --branch-peak-table short_pmax2_peaks.csv \
  --baseline-label standard \
  --branch-label short_pmax2 \
  --branch-search-branch short_pmax2 \
  --out-dir branch_compare_standard_vs_short_pmax2
```

The comparison is the acceptance gate for making the short-period branch part
of the candidate generator: it must increase union top-N recovery in the
high-observability miss population without introducing a large degradation in
previously recovered injections.
The comparison defaults to the overlap/intersection of evaluated injection IDs,
so a short branch run only on high-observability misses is not penalized for
missing the rest of the `20k` sample.

Do not deploy the branch-provenance builder onto the currently active PDO
`20k` chunked job; wait for that run to finish so the table schema remains
internally consistent.

## ORCD / H200 Role

ORCD should not run primary Stage 1 TGLC/ePSF production. Use it once compact
S56 products are staged under `/orcd/data/mki_aryeh/001/twirl/`.

Near-term ORCD jobs:

1. Copy compact S56 light-curve exports, injection HDF5s, BLS candidate tables,
   and peak-training tables to ORCD project storage.
2. Train the peak-ranker from injected peak rows.
3. Run larger injected grids and repeated ranker evaluations.
4. Later: run GPU-accelerated matched filters or neural feature extraction on
   compact light-curve tensors.

The current BLS implementation is CPU/Astropy based, so the H200s are not the
first bottleneck for this exact script. The H200s become useful once we train
models or build GPU-native search/feature extraction. Until then, use ORCD CPU
nodes for large table builds and reserve H200s for model experiments.

Concrete ORCD entrypoints now live in [runbook](s56_orcd_next_steps.md):

- [stage script](../../scripts/orcd/stage_s56_orcd_inputs.sh) syncs the current
  checkout plus compact PDO artifacts under the ORCD project tree.
- [CPU Slurm script](../../scripts/orcd/slurm_s56_peak_training_cpu.sbatch)
  runs smoke or full peak-table builds from the staged `20k` injections through
  the restartable [chunked launcher](../../scripts/stage5_validation/run_s56_20k_peak_training_chunked.sh).
- [ranker Slurm script](../../scripts/orcd/slurm_s56_peak_ranker_train_cpu.sbatch)
  trains/evaluates the peak ranker after the table exists.
- [application script](../../scripts/stage5_validation/apply_injection_peak_ranker.py)
  scores real S56 BLS peaks and writes top-N selected ephemerides per TIC.
- [H200 smoke](../../scripts/orcd/slurm_h200_smoke.sbatch) verifies CUDA/GPU
  visibility before model-training jobs are submitted.

## Training Stages

1. **Peak ranker from injections**
   - Positive: BLS peak exact/harmonic matches injected ephemeris.
   - Negative: other BLS peaks from the same injected light curves, including
     negative-only peak lists where BLS never put the injected signal in top-N.
   - Target metric: recall@K for recovering the injected signal.
   - Features: rank, SDE, depth SNR, duration, qtran, cadence-cleaning
     diagnostics, period aliases, aperture agreement, and Tmag. Injection
     depth/radius are audit strata only, not ranker features for real-data
     inference.

2. **Real-candidate object classifier**
   - Uses real S56 BLS candidates plus human labels.
   - Labels: planet-like, PCEB/EB, variable, contaminant, systematic, uncertain.
   - Must include real false positives; injections alone cannot teach this.

3. **Comparison baselines**
   - BLS top-1 only.
   - BLS plus peak-ranker.
   - LEO-only.
   - ML without LEO features.
   - ML with LEO features.

## Human Triage Guidance

The current `~100` labels are useful as a UI/search smoke test, but they should
not yet be treated as a final training set. Continue light triage only if it is
helpful for interface feedback. Do not scale to `10k` human labels until the
peak-level BLS recall audit is complete and the displayed ephemerides are more
often the right injected signal.

For the next real queue, show:

- the best ranker-selected ephemeris,
- whether an injected row is exact/harmonic/mismatch hidden from the reviewer,
- LEO report,
- optional top-2/top-3 alternate ephemerides when the ranker is uncertain.

## Acceptance Gates

- Peak-training table exists for the current `20k` injected set.
- Summary reports recall@1/2/3/5/10/20.
- Decision made from recall@20:
  - train ranker if the signal is usually present,
  - tune/replace BLS if not.
- Human-label queue is regenerated only after candidate ephemeris selection is
  improved.
