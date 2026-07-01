# S56 ORCD Next Steps

Date: 2026-07-01

## Decision

Use ORCD as downstream compute for compact S56 products, not as a replacement
for PDO Stage 1 extraction. The active S56 pre-detrend BATMAN recovery product
still depends on raw TGLC `RawFlux` arrays on PDO, but its compact outputs are
small enough to stage to ORCD:

- compact S56 light-curve export: `3.3 GB`
- `20k` injected light-curve HDF5: `11 GB`
- current BLS recovery table: `27 MB`
- real S56 BLS peak table for ranker application: about `1.1M` rows

That scale is appropriate for `/orcd/data/mki_aryeh/001/twirl/`.

## Current ORCD State

The ORCD control-socket protocol is working again. The project storage layout
now exists under `/orcd/data/mki_aryeh/001/twirl/`.

Live status on `2026-07-01`:

- balanced-grid peak-table job `16934436` is running CPU-only on `node4702`
  with `48` CPUs, `192 GB`, and no GPU/GRES;
- dependent CPU-only jobs `16934439` (ranker/gate) and `16934507`
  (real-candidate apply) are waiting on `afterok`;
- the peak-table log confirms manifest-derived chunk IDs and is advancing
  through the `80` restartable chunks;
- no H200s are currently allocated to TWIRL; the one-H200 tensor data smoke
  `16939737` completed cleanly.

The all-host PDO branch has already produced a verified real-candidate triage
queue, so human labeling can proceed while the ORCD balanced-grid CPU chain
finishes. Use `pdogpu1:5007` for triage and run
`scripts/stage5_validation/run_s56_allhost_real_label_audit_pdo.sh` after labels
exist. A detached PDO monitor `twirl-s56-allhost-label-audit` is also running
and will rerun that post-label audit whenever the label count changes. A local
monitor `twirl-pdo-label-sync` pulls those small label/readiness products back
to the local checkout and refreshes the combined ML handoff readiness report.

A detached local tmux monitor `twirl-orcd-apply-leo-smoke` is watching for the
balanced-grid ORCD apply outputs. When `s56_ranker_selected_real_candidates_orcd/`
is complete, it syncs the verified small result bundle back to both the local
checkout and PDO, then launches a bounded `50`-row LEO queue under
`s56_ranker_selected_real_leo_queue_orcd50_pdo/`. It exits rather than
initiating interactive auth if the ORCD control socket expires.

Staged checkout and compact data:

```text
/orcd/data/mki_aryeh/001/twirl/code/TWIRL/
```

The full `20k` injected light-curve product is staged at:

```text
/orcd/data/mki_aryeh/001/twirl/code/TWIRL/data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_20k_predetrend_batman_periodradius_bright_balanced/injected_lightcurves.h5
```

Verification on ORCD found the expected `11,283,081,701`-byte HDF5 plus
`20,000` manifest rows and `20,000` label rows. The reusable Python environment
is:

```text
/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python
```

Host-coverage audit on the staged ORCD manifests found `5,323 / 19,072`
unique S56 TIC hosts represented (`27.9%`) with `0` injected TICs missing from
the compact export. The injection grid has `2,500` observed period-radius cells
with exactly `8` injections per cell, so this is a deliberately balanced
recovery/ranker grid rather than literal full-S56 host coverage.

To build the coverage-first companion sample on PDO, use:

```bash
cd /pdo/users/tehan/TWIRL
N_INJECTIONS=19072 \
  bash scripts/stage5_validation/run_s56_predetrend_all_host_grid_pdo.sh
```

That launcher calls the pre-detrend BATMAN injector with
`--target-selection-mode shuffled_cycle`, which attempts every discovered raw
S56 TIC once before any repeats. Its summary and host-coverage audit record the
actual achieved coverage, so if the raw TGLC source set is smaller than the
HLSP export or some targets fail injection, the gap is explicit.

A 2-target PDO smoke passed under the QLP runtime after setting
`LD_LIBRARY_PATH=/sw/python-versions/python-3.11.9/lib:/pdo/app/anaconda/anaconda2-4.4.0/lib`.
The launcher now uses that default path.

For the full all-host product, prefer the restartable sharded launcher:

```bash
cd /pdo/users/tehan/TWIRL
CHUNK_JOBS=2 \
SHARD_SIZE=250 \
  bash scripts/stage5_validation/run_s56_predetrend_all_host_grid_sharded_pdo.sh
```

The sharded path groups both S56 orbit HDF5 files by TIC, runs chunk-local
pre-detrend BATMAN injections with global injection-index offsets, and merges
manifest/label metadata without trying to concatenate the per-chunk HDF5 files.
A real 2-target PDO smoke passed with `2` completed shards, `2` source targets,
`2` unique injected TICs, `all_source_targets_injected=true`, zero skip counts,
and `0` injected TICs missing from the `19,072`-target compact export.
`scripts/orcd/stage_s56_orcd_inputs.sh` now treats the sharded all-host output
directory as an optional compact artifact to stage to ORCD once it exists.

The full sharded all-host injection run completed on `pdogpu1`: `19,072`
injections across `77` TIC-grouped shards, covering `19,071 / 19,072` unique
exported TIC hosts. The one missing host is TIC `2019898202`, skipped for a
nonpositive baseline; TIC `2019893140` received the duplicate replacement
injection. The host audit records this explicitly.

Monitor the completed injection product with:

```bash
cd /pdo/users/tehan/TWIRL
LD_LIBRARY_PATH=/sw/python-versions/python-3.11.9/lib:/pdo/app/anaconda/anaconda2-4.4.0/lib \
PYTHONPATH=src \
  /sw/qlp-environment/.venv/bin/python \
  scripts/stage5_validation/summarize_s56_allhost_injection_progress.py \
  --out-dir reports/stage5_validation/s56_allhost_predetrend_batman_periodradius_grid_sharded_pdo/progress
```

Robust BLS peak labeling on the sharded product is active in tmux
`twirl-s56-allhost-peak-monitor`. The active run uses bounded concurrency:
`TWIRL_ALLHOST_PEAK_CHUNK_JOBS=4` and `TWIRL_PEAK_WORKERS=8` (`~32` CPU
workers on `pdogpu1`).

The equivalent manual command is:

```bash
cd /pdo/users/tehan/TWIRL
TWIRL_ALLHOST_PEAK_CHUNK_JOBS=4 \
TWIRL_PEAK_WORKERS=8 \
  bash scripts/stage5_validation/run_s56_allhost_peak_training_sharded_pdo.sh
```

A one-shard plumbing smoke passed with `2` injections, `80` peak rows, `2`
apertures, `20` positive peak ranks, and cadence-diagnostic columns present.
The first real all-host BLS wave also completed cleanly: `4/77` chunks,
`1,000` injections, `40,000` peak rows, and top-20 injected-signal recall
`504/1000`, with the expected bright-to-faint recovery gradient.

Audit the all-host handoff status from PDO with:

```bash
cd /pdo/users/tehan/TWIRL
LD_LIBRARY_PATH=/sw/python-versions/python-3.11.9/lib:/pdo/app/anaconda/anaconda2-4.4.0/lib \
PYTHONPATH=src \
  /sw/qlp-environment/.venv/bin/python \
  scripts/stage5_validation/audit_s56_peak_handoff_status.py \
  --layout pdo_allhost \
  --root reports/stage5_validation/s56_allhost_predetrend_batman_periodradius_grid_sharded_pdo \
  --out-dir reports/stage5_validation/s56_allhost_peak_handoff_audit_pdo
```

After the all-host BLS runner writes and verifies the merged peak table, run:

```bash
cd /pdo/users/tehan/TWIRL
bash scripts/stage5_validation/run_s56_allhost_peak_ranker_review_pdo.sh
```

That command reuses the existing injected-peak verifier, post-BLS recovery
gate, failure-mode audit, lightweight peak ranker, real-S56 peak verifier, and
ranker-selected review-queue builder with all-host paths and minimum counts.

A detached PDO monitor is also available and is currently running as
`twirl-s56-allhost-ranker-monitor`:

```bash
cd /pdo/users/tehan/TWIRL
WAIT_SECONDS=600 \
TWIRL_SKIP_LEO=1 \
  bash scripts/stage5_validation/monitor_s56_allhost_peak_to_ranker_pdo.sh
```

It waits for the merged all-host peak-table verifier before launching the
handoff. The default `TWIRL_SKIP_LEO=1` produces a verified ranker-selected real
queue first; run the LEO rendering step only after that queue passes.

Completed smokes on `2026-06-30`:

- CPU peak-table smoke `16869128`: no GPU request, ran on `node4701`, completed
  in `2:16`, wrote `4,000` peak rows for `100` injections, and passed the
  strict injected-peak verifier. Recall snapshot: top-1 `39/100`, top-20
  `45/100`.
- H200 visibility smoke `16869129`: requested exactly `1` H200, ran on
  `node4900`, completed in `1 s`, and confirmed `CUDA_VISIBLE_DEVICES=0` plus
  `nvidia-smi` visibility. The current baseline env does not include `torch` or
  `cupy`.
- H200 tensor data smoke `16939737`: requested exactly `1` H200, ran on
  `node4900`, completed in `12 s`, read the compact S56 LC export plus the
  staged all-host real queue metadata, and wrote `128/128` candidate-centered
  tensors with shape `(128, 3, 257)`. The tensor channels are `DET_FLUX_SML`,
  `DET_FLUX`, and `DET_FLUX_LAG`; the smoke confirmed again that the current
  `twirl-s56` env does not include `torch`.
- H200 torch tensor smoke `16942919`: requested exactly `1` H200, ran on
  `node4900`, completed in `47 s`, used
  `/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56-torch/bin/python`, imported
  `torch 2.11.0+cu128`, reported `cuda_available=true`, `device_count=1`, and
  moved the tensor smoke batch onto `cuda:0`.
- H200 synthetic training smoke `16944030`: requested exactly `1` H200, ran on
  `node4900`, completed in `6 s`, used synthetic labels only to validate the
  training loop, trained on `128` tensor rows with train/validation/test split
  `76/26/26`, and confirmed training on `cuda:0`.

After human labels exist and the readiness audit passes, stage the small
label/readiness products to ORCD and then submit the real-label training smoke
with:

```bash
scripts/orcd/run_s56_orcd_pilot.sh --run stage-labels
scripts/orcd/run_s56_orcd_pilot.sh --run h200-tensor-train-smoke
```

Current ORCD run state:

- Peak-table build `16934436`: CPU-only (`48` CPUs, `192 GB`, no GPU/GRES) on
  `pg_mki_aryeh`, running on `node4702` with manifest-derived chunk IDs.
- Ranker/gate `16934439`: dependency-pending on `afterok:16934436`, CPU-only.
- Real-candidate ranker application `16934507`: dependency-pending on
  `afterok:16934439`, CPU-only.

The active peak-table output root is:

```text
/orcd/data/mki_aryeh/001/twirl/code/TWIRL/reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training_orcd_full/
```

The full run uses `80` restartable chunks of `250` injections. The earlier
`100`-injection smoke outputs are separated under `peak_training_orcd_smoke100`
so they cannot poison restart/skip logic for the full run.

The staged bundles have since been repaired and verified:

- all-host: `19,072` manifest/label rows and the expected `19,071` unique TICs
- balanced: `3.3 GB` compact LC export, `11 GB` injected HDF5, `20,000`
  manifest/label rows, and HDF5 attr `n_injections=20000`

A subset stage that was started before the staging helper was fixed removed
excluded `data_local/` and `reports/` trees during checkout sync. The local
helper now preserves excluded trees; if these artifacts are missing again,
rerun:

```bash
scripts/orcd/run_s56_orcd_pilot.sh --run stage-balanced
```

Historical note: peak job `16869261` failed after chunk `038`, and restarted job
`16881607` stalled before chunk-manifest creation because the old runner
enumerated HDF5 group names. Those jobs and their dependent ranker/apply jobs
were cancelled. The runner now derives chunk IDs from `injection_manifest.csv`
and `h5_group`, falling back to HDF5 only if no manifest exists; the live job
`16934436` confirms this path with `key_source=manifest` and `n_chunks=80`.

The real S56 BLS table is staged as:

```text
/orcd/data/mki_aryeh/001/twirl/code/TWIRL/data_local/stage2/bls_first_pass_v2/sector_0056/candidates.csv
```

Monitor the chain from the local checkout with:

```bash
bash scripts/orcd/check_s56_orcd_pipeline_status.sh
```

## Immediate ORCD Pilot

Current access state: non-interactive access works after the user opens the
control socket. Agents must still never initiate Duo/password prompts.

Authentication rule: do not let Codex, scripts, or automated probes initiate
Duo, password, or keyboard-interactive login attempts. If the control socket is
missing, stop and ask the user to open it from their own terminal.

1. Open the ORCD SSH control socket from a user terminal, per
   [ORCD guide](../../doc/orcd_h200_usage.md).
2. Run the smoke pilot wrapper from the local checkout:

```bash
scripts/orcd/run_s56_orcd_pilot.sh smoke-pilot
```

This is a dry run by default. After the printed commands look right and the
ORCD control socket is open, execute:

```bash
scripts/orcd/run_s56_orcd_pilot.sh --run smoke-pilot
```

That probes ORCD, refreshes the staged checkout and compact PDO artifacts if
needed, submits a `100`-injection CPU peak-table smoke, and submits an H200
smoke job.

3. If the smoke jobs pass, submit the full/review sequence as separate,
   inspectable steps:

```bash
scripts/orcd/run_s56_orcd_pilot.sh --run peak-full
scripts/orcd/run_s56_orcd_pilot.sh --run ranker
scripts/orcd/run_s56_orcd_pilot.sh --run apply
scripts/orcd/run_s56_orcd_pilot.sh --run status
```

If the standard peak table shows high-observability not-in-top-N misses, run
the short-period branch comparison as a separate CPU job:

```bash
scripts/orcd/run_s56_orcd_pilot.sh --run short-branch
```

Lower-level submission:

```bash
ssh -o BatchMode=yes -o ControlMaster=auto \
  -o ControlPath=~/.ssh/cm/%r@%h:%p \
  tehan@orcd-login.mit.edu \
  'cd /orcd/data/mki_aryeh/001/twirl/code/TWIRL &&
   sbatch scripts/orcd/slurm_s56_short_branch_compare_cpu.sbatch'
```

The lower-level commands are preserved below for debugging.

### Lower-Level Commands

Stage the local checkout and compact PDO artifacts:

```bash
scripts/orcd/stage_s56_orcd_inputs.sh --run
```

After the all-host PDO product or all-host peak tables change, use the smaller
subset stage instead of restaging the balanced `20k` HDF5:

```bash
scripts/orcd/run_s56_orcd_pilot.sh --run stage-allhost
```

The staging helper preserves excluded `data_local/`, `reports/`, and `logs/`
trees during checkout sync. If an older helper was used, verify both the
balanced-grid and all-host artifacts before restarting Slurm jobs.

If the PDO all-host peak build remains the bottleneck and ORCD CPU capacity is
available, the all-host BLS peak build can also be run on ORCD after
`stage-allhost`:

```bash
scripts/orcd/run_s56_orcd_pilot.sh --run allhost-peak
```

Run the following only after the all-host peak table exists and passes
verification:

```bash
scripts/orcd/run_s56_orcd_pilot.sh --run allhost-ranker
scripts/orcd/run_s56_orcd_pilot.sh --run allhost-apply
scripts/orcd/run_s56_orcd_pilot.sh --run allhost-sync-apply
```

Use `allhost-monitor-apply` instead of `allhost-sync-apply` if the apply job is
still running and you want the helper to wait.

For the balanced-grid product only:

```bash
scripts/orcd/run_s56_orcd_pilot.sh --run stage-balanced
```

The script preserves repo-relative paths under:

```text
/orcd/data/mki_aryeh/001/twirl/code/TWIRL/
```

Submit a small CPU peak-training smoke:

```bash
ssh -o BatchMode=yes -o ControlMaster=auto \
  -o ControlPath=~/.ssh/cm/%r@%h:%p \
  tehan@orcd-login.mit.edu \
  'cd /orcd/data/mki_aryeh/001/twirl/code/TWIRL &&
   sbatch --export=ALL,TWIRL_N_INJECTIONS=100 scripts/orcd/slurm_s56_peak_training_cpu.sbatch'
```

If the smoke writes a summary with nonzero rows, submit the full CPU run:

```bash
ssh -o BatchMode=yes -o ControlMaster=auto \
  -o ControlPath=~/.ssh/cm/%r@%h:%p \
  tehan@orcd-login.mit.edu \
  'cd /orcd/data/mki_aryeh/001/twirl/code/TWIRL &&
   sbatch scripts/orcd/slurm_s56_peak_training_cpu.sbatch'
```

After the full peak table exists, train/evaluate the lightweight peak ranker:

```bash
ssh -o BatchMode=yes -o ControlMaster=auto \
  -o ControlPath=~/.ssh/cm/%r@%h:%p \
  tehan@orcd-login.mit.edu \
  'cd /orcd/data/mki_aryeh/001/twirl/code/TWIRL &&
   sbatch scripts/orcd/slurm_s56_peak_ranker_train_cpu.sbatch'
```

Apply the trained ranker to real S56 BLS peaks to produce candidate
ephemerides for LEO/human review:

```bash
ssh -o BatchMode=yes -o ControlMaster=auto \
  -o ControlPath=~/.ssh/cm/%r@%h:%p \
  tehan@orcd-login.mit.edu \
  'cd /orcd/data/mki_aryeh/001/twirl/code/TWIRL &&
   sbatch scripts/orcd/slurm_s56_peak_ranker_apply_cpu.sbatch'
```

Submit the H200 smoke independently:

```bash
ssh -o BatchMode=yes -o ControlMaster=auto \
  -o ControlPath=~/.ssh/cm/%r@%h:%p \
  tehan@orcd-login.mit.edu \
  'cd /orcd/data/mki_aryeh/001/twirl/code/TWIRL &&
   sbatch scripts/orcd/slurm_h200_smoke.sbatch'
```

## What Runs Where

PDO remains the correct place for raw-flux pre-detrend injection generation
until we intentionally export raw TGLC arrays or a smaller raw-aperture shard.

ORCD CPU nodes are the correct immediate place for:

- peak-level BLS truth table builds,
- ranker cross-validation over peak tables,
- larger table joins and queue construction without NFS pressure on PDO.

The CPU peak-table Slurm job uses the restartable chunked launcher
[script](../../scripts/stage5_validation/run_s56_20k_peak_training_chunked.sh),
so interrupted jobs can resume from completed chunk tables instead of losing the
whole `20k` run.

The peak ranker itself is trained with
[script](../../scripts/stage5_validation/train_injection_peak_ranker.py). It
uses only BLS/LC features available for real candidates and writes a model,
scored peaks, selected ephemerides, coefficients, and recall@K summaries.
Before training, the ORCD ranker wrapper now also runs
[script](../../scripts/stage5_validation/audit_s56_injection_host_coverage.py)
to quantify whether the `20k` injection hosts represent the full S56 export or
a deliberately balanced subset. This is a host-coverage audit only; it does not
change the injected peak labels.
The Slurm training wrapper also runs
[script](../../scripts/stage5_validation/summarize_injection_peak_gate.py),
which writes the explicit search-vs-ranker gate: top-1 recovery, top-N recovery,
ranking losses, and binned period-radius-Tmag empirical recovery cells.
For short-period/systematic peak-crowding diagnostics, use
[runner](../../scripts/stage5_validation/run_s56_short_branch_compare.sh) or
[Slurm wrapper](../../scripts/orcd/slurm_s56_short_branch_compare_cpu.sbatch).
Apply it to real BLS candidate tables with
[script](../../scripts/stage5_validation/apply_injection_peak_ranker.py); this
writes the top ranker-selected ephemerides per TIC for the next LEO/human queue.
The corresponding Slurm wrapper is
[script](../../scripts/orcd/slurm_s56_peak_ranker_apply_cpu.sbatch). The wrapper
uses the CSV real-peak table and disables the full scored parquet output, so it
does not depend on a parquet engine in the ORCD Python environment.

PDO remains the right place to render the ranker-selected LEO PDFs because it
has the full S56 HLSP tree and working LEO-Vetter environment. After ORCD
writes `s56_ranker_selected_real_candidates_orcd/selected_ephemerides.csv`,
sync the small verified ranker outputs back to PDO:

```bash
scripts/orcd/run_s56_orcd_pilot.sh --run sync-apply
```

Then render the browser-ready LEO queue on PDO from the ORCD-selected
ephemerides:

```bash
TWIRL_SELECTED_EPHEMERIDES=reports/stage5_validation/s56_ranker_selected_real_candidates_orcd/selected_ephemerides.csv \
  bash scripts/stage5_validation/run_s56_ranker_selected_real_leo_pdo.sh
```

For a small first render, set `TWIRL_REVIEW_N_REAL=50` and
`TWIRL_MAX_LEO_REPORTS=50`.

The guarded monitor combines the wait, sync, and bounded PDO LEO-smoke render:

```bash
scripts/orcd/run_s56_orcd_pilot.sh --run monitor-apply-leo-smoke
```

It treats missing ORCD apply outputs as waitable, but stops on socket or
transfer errors. The legacy `monitor-apply` command remains available for
sync-only operation.

While the PDO one-shot `20k` peak-table job is running, use
[monitor](../../scripts/stage5_validation/monitor_s56_peak_training_pdo.sh) to
avoid manual polling:

```bash
cd /pdo/users/tehan/TWIRL
bash scripts/stage5_validation/monitor_s56_peak_training_pdo.sh --wait
```

After the gate summary is acceptable, the same monitor can trigger the
ranker/review handoff. Use `--skip-leo` for a queue-only smoke before rendering
the full PDF set:

```bash
TWIRL_REVIEW_N_REAL=50 \
  bash scripts/stage5_validation/monitor_s56_peak_training_pdo.sh \
    --wait --run-review --skip-leo
```

H200s are useful after those tables exist:

- peak-ranker and teacher/student model training,
- raw-light-curve tensor experiments,
- GPU-native matched-filter or BLS replacement work.

Do not spend H200 time rerunning the current Astropy BLS script unless the CPU
queue becomes the bottleneck and a GPU-native search implementation exists.

## Gate Before Human-Label Scale-Up

Use the `20k` injected peak table to decide what problem we actually have:

- high recall@20 but low recall@1: train a peak ranker and show ranker-selected
  ephemerides to LEO/humans;
- low recall@20: improve the search statistic or detrending before asking for
  more human labels.
- incomplete host coverage: keep using this as a balanced recovery/ranker
  training set, but do not call it an all-S56 host-distribution model until the
  host-coverage audit says the exported S56 population is represented well
  enough for that claim.
- all-host companion sample ready: stage
  `pdo_allhost_predetrend_batman_periodradius_grid/` to ORCD, run the same
  robust peak-table build/ranker gate, and compare recall/gate summaries
  against the balanced `20k` grid.

For the final training table, prefer the restartable chunked peak build because
it writes the current robust-BLS cadence diagnostics (`n_cad_quality`,
`n_cad_edge_trimmed`, `n_cad_sigma_clipped`, `quality_dropout_frac`,
`dropout_frac`). The older one-shot PDO run is still useful for a fast
recall@K answer, but it was launched before those diagnostic columns were
added.

Only after that gate should the mixed real+injected queue become a training set.
Human labels, LEO classes, and injected-truth recovery labels should remain
separate columns.
