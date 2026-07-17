# ORCD H200 Usage For TWIRL

This note records the initial ORCD/Engaging reconnaissance for the MKI Aryeh
CPU/GPU nodes and how TWIRL should use them. It is operational guidance, not a
scientific pipeline decision.

Last updated: `2026-07-02`. ORCD access last verified: `2026-06-30`.

## Access And Scheduler Names

Use the ORCD login host:

```bash
ssh tehan@orcd-login.mit.edu
```

Login requires the account password plus Duo. Do not put passwords in scripts.
For repeated local probes, use an SSH control socket that you authenticate once
from your own terminal.

Important authentication rule: agents and scripts must never initiate Duo,
password, or keyboard-interactive login attempts. If a non-interactive probe
does not find an existing control socket, stop immediately and ask the user to
log in or re-open the socket from their own terminal. Do not retry in a way
that can send multiple Duo requests.

### User-Initiated ORCD Access Protocol

Use this protocol whenever Codex or another agent needs ORCD access. ORCD is
different from PDO because login requires a password plus Duo, so the user owns
the interactive authentication step and the agent only reuses the resulting SSH
control socket.

1. From a user terminal, open the shared ORCD control socket:

```bash
mkdir -p ~/.ssh/cm
chmod 700 ~/.ssh ~/.ssh/cm

ssh -MNf \
  -o ControlMaster=yes \
  -o ControlPath="$HOME/.ssh/cm/%r@%h:%p" \
  -o ControlPersist=8h \
  -o ServerAliveInterval=60 \
  -o ServerAliveCountMax=3 \
  tehan@orcd-login.mit.edu
```

Enter the password and approve Duo in that terminal only. After the command
returns, tell Codex that the ORCD socket is open.

If this prints:

```text
ControlSocket ... already exists, disabling multiplexing
```

do not start another master connection. It usually means a control socket is
already present at that path. Run the check command below; if it succeeds, the
socket is already ready for Codex. If it fails because the socket is stale,
close or remove that stale socket from the user terminal, then rerun the
`ssh -MNf ...` command.

The TWIRL ORCD helper scripts use this same path by default via
`ORCD_CONTROL_PATH=$HOME/.ssh/cm/%r@%h:%p`. If a different control path is used,
set `ORCD_CONTROL_PATH` consistently before running the helpers.

2. Optional user-side check, still from the same terminal:

```bash
ssh \
  -O check \
  -o BatchMode=yes \
  -o PasswordAuthentication=no \
  -o KbdInteractiveAuthentication=no \
  -o NumberOfPasswordPrompts=0 \
  -S "$HOME/.ssh/cm/%r@%h:%p" \
  tehan@orcd-login.mit.edu
```

3. Agent-side probe. Codex may run exactly this style of non-interactive check:

```bash
ssh \
  -o BatchMode=yes \
  -o PasswordAuthentication=no \
  -o KbdInteractiveAuthentication=no \
  -o NumberOfPasswordPrompts=0 \
  -o ControlMaster=auto \
  -o ControlPath="$HOME/.ssh/cm/%r@%h:%p" \
  tehan@orcd-login.mit.edu \
  'hostname; whoami; sinfo -h -p pg_mki_aryeh -o "%P %D %t %G" | head -5'
```

If that probe fails with `Permission denied`, `keyboard-interactive`,
`Connection closed`, or a missing-control-socket message, the agent must stop
and ask the user to re-open the socket. Do not retry ORCD login in a way that
can trigger another password or Duo request.

4. Agent-side TWIRL commands after the probe succeeds:

```bash
scripts/orcd/run_s56_orcd_pilot.sh --run probe
scripts/orcd/run_s56_orcd_pilot.sh --run smoke-pilot
scripts/orcd/run_s56_orcd_pilot.sh --run status
```

5. Close the ORCD socket when done, from the user terminal:

```bash
ssh \
  -O exit \
  -o BatchMode=yes \
  -o PasswordAuthentication=no \
  -o KbdInteractiveAuthentication=no \
  -o NumberOfPasswordPrompts=0 \
  -S "$HOME/.ssh/cm/%r@%h:%p" \
  tehan@orcd-login.mit.edu
```

The administrative wording from Paul was:

```text
orcd_ug_pg_mki_aryeh_all
```

That is not the runnable Slurm partition name observed on ORCD. The runnable
partition is:

```bash
-p pg_mki_aryeh
```

Verification from `2026-06-10` and re-check on `2026-06-30`:

- `sinfo` lists `pg_mki_aryeh`.
- `scontrol show partition orcd_ug_pg_mki_aryeh_all` reports the partition is
  not found.
- `sbatch --test-only -p orcd_ug_pg_mki_aryeh_all ...` fails with an invalid
  partition.
- `sbatch --test-only -p pg_mki_aryeh ...` succeeds.
- `2026-06-30`: `sbatch --test-only -p pg_mki_aryeh --gres=gpu:h200:1 ...`
  and `--gres=gpu:h200:8` both succeed and target `node4900`.
- `2026-06-30`: `scontrol show partition pg_mki_aryeh` reports
  `AllowGroups=orcd_rg_par_pg_mki_aryeh,slurm_admin`, `AllowAccounts=ALL`,
  `QoS=pg_mki_aryeh`, `MaxTime=2-00:00:00`, nodes `node[4701-4702,4900]`,
  and TRES `cpu=312`, `gres/gpu:h200=8`.
- The observed Unix access group is `orcd_rg_par_pg_mki_aryeh`.

If ORCD staff later says a different Slurm name should be used, update this
file before changing any scripts.

Other partitions accepted `sbatch --test-only` for this account on
`2026-06-30`, including `mit_normal`, `mit_preemptable`, `mit_quicktest`,
`mit_normal_gpu`, and `mit_data_transfer`. Some `ou_mki*` partitions also
accepted test-only submissions, but TWIRL should not use other group partitions
without explicit policy confirmation. The default TWIRL target remains
`pg_mki_aryeh`; use MIT general/preemptable partitions only for clearly
appropriate fallback or short diagnostic work.

## Hardware

The private Slurm partition currently exposes three nodes:

| Node | Role | Slurm resources |
| --- | --- | --- |
| `node4701` | CPU | `96` CPUs, about `1.5 TB` RAM |
| `node4702` | CPU | `96` CPUs, about `1.5 TB` RAM |
| `node4900` | GPU | `120` CPUs, about `2 TB` RAM, `gpu:h200:8` |

Request H200 GPUs with Slurm GRES syntax:

```bash
--gres=gpu:h200:N
```

where `N` is the number of H200s requested. ORCD public documentation also
shows the shorter `-G h200:N` form, but the GRES form is explicit and matches
what `scontrol show node node4900` reports.

TWIRL usage rule: default to `1` H200 for smoke tests and early GPU prototypes,
and use at most `2` H200s for routine development runs. Do not request all `8`
H200s unless the user explicitly approves that specific job after the smaller
run has passed and the scaling need is clear. The standard ORCD submit helper
refuses Slurm scripts requesting more than `2` H200s unless
`TWIRL_ALLOW_MANY_H200=1` is set for a specifically approved submit.

The `2026-06-30` control-socket probe succeeds non-interactively from Codex
after the user opens the socket. `pg_mki_aryeh` currently exposes two CPU nodes
and one H200 node; test-only submissions for both CPU and H200 jobs succeed.
TWIRL can use ORCD operationally once compact S56 light-curve exports are
staged there; do not move raw PDO TGLC, TICA, ePSF, or source-pickle trees as
the first step.

The private partition has:

- max wall time: `2-00:00:00` (`48 h`)
- nodes: `node[4701-4702,4900]`
- total TRES: `cpu=312`, `gres/gpu:h200=8`

Paul's email also noted a MGHPCC shutdown from `2026-06-15` to `2026-06-18`;
avoid starting fragile long jobs across that window.

## Storage

Project storage:

```text
/orcd/data/mki_aryeh/001
```

Observed state on `2026-06-10`:

- size: `200 TB`
- writable by this account
- group: `orcd_rg_hstor003_pg_mki_aryeh`
- empty at first probe

Recommended TWIRL layout:

```text
/orcd/data/mki_aryeh/001/twirl/
  exports/
  results/
  logs/
  envs/
```

Personal ORCD locations:

| Path | Quota | Use |
| --- | --- | --- |
| `/home/tehan` | `200 GB` | code, small configs, environment metadata |
| `/home/tehan/orcd/pool` | `1 TB` | personal staged data, not backed up |
| `/home/tehan/orcd/scratch` | `1 TB` | active job scratch, not backed up |

Use `/orcd/data/mki_aryeh/001/twirl/` for shared TWIRL products. Use scratch
for temporary job-local expansions and delete or compact outputs afterward.

## TWIRL Policy

PDO remains the Stage 1 TGLC/ePSF production home because it has the TICA/QLP
file layout, pyticdb/catalog access, and existing TGLC production tree. Do not
try to re-home the primary raw extraction pipeline onto ORCD unless a separate
data-staging design is approved.

ORCD should consume compact downstream products:

- TWIRL-FS HLSP FITS exports or, preferably, sector/camera/CCD HDF5 or Parquet
  shards built from those products
- candidate tables
- injection-recovery input grids and output summaries
- vetting feature tables
- compact manifests with input paths, checksums, code versions, and sample cuts

Operational boundary as of `2026-07-02`: use PDO only for Stage 1/TGLC/HLSP
production and for building compact exports from the PDO-resident light-curve
tree. Run downstream testing, vetting, injection-recovery, BLS branch audits,
and model-training experiments on ORCD after those compact products are staged.
This includes two classes of CPU-bound Astropy-BLS work: production-level
raw-flux detrending-strength audits, and post-ADP vetting diagnostics such as
ADP+/two-aperture sheets. Run these on ORCD CPU nodes first, not on the H200
node and not as long PDO tmux jobs.

Current branch name: the July 2026 raw-flux audit promotes
`twirl-fs-v2-adp015q` (`bkspace_d=0.15 d`, `gap_split_d=0.2 d`, quantile
knots) as the next S56 search/vetting candidate. Build its compare FITS on PDO
with `scripts/stage1_lightcurves/run_s56_twirlfs_adp015q_compare_pdo.sh`, then
stage compact `DET_FLUX_ADP015*` exports to ORCD. Do not emulate this by adding
a second ADP+ layer to older `DET_FLUX_ADP*` columns.

As of `2026-07-04`, the S56 ADP015 compare FITS tree and compact two-aperture
export have both been built on PDO. The next ORCD action is to stage the compact
ADP015 export and run full two-aperture BLS/vet-sheet production from that
export; keep PDO for product generation and short smokes only.

The two-aperture renderer supports compact-HDF5 input via `--lc-export-h5`, so
ORCD does not need the full ADP015 FITS tree. Until the ADP015 code path is
merged/pushed, run `sync-code` before staging/submitting:

```bash
scripts/orcd/run_s56_orcd_pilot.sh --run sync-code stage-balanced twoap-smoke
```

After the smoke Slurm job completes, sync the smoke products for inspection:

```bash
scripts/orcd/run_s56_orcd_pilot.sh --run sync-twoap-smoke
```

If the smoke passes, submit the full CPU render:

```bash
scripts/orcd/run_s56_orcd_pilot.sh --run twoap-full
```

After the full render completes, sync the full sheet directory plus metrics back
to local and PDO:

```bash
scripts/orcd/run_s56_orcd_pilot.sh --run sync-twoap
```

Do not transfer raw TGLC cutout, source-pickle, ePSF, or FFI trees by default.
Those products are large, PDO-layout-specific, and not the right first target
for ORCD.

Best TWIRL uses for the 8xH200 node:

1. LC-level Stage 3 injection-recovery grids.
2. Stage 4 GPU search acceleration, including matched filters, dip searches,
   and a cuvarbase-style BLS path.
3. Multi-sector feature extraction and candidate merging.
4. Pixel-level injection or pixel-map vetting only on calibrated subsets.
5. ML/anomaly triage only after transparent search and vetting baselines are
   working and labeled false positives exist.

## Software Stack

Observed modules/tools:

- CUDA modules: `cuda/12.9.1`, `cuda/13.0.1`, `cuda/13.1.0`
- cuDNN module: `cudnn/9.8.0.87-cuda12`
- Miniforge: `miniforge/25.11.0-0`, Python `3.12.12`
- GCC `12.2.0`
- HDF5 and OpenMPI modules
- `apptainer` / `singularity`

No prebuilt PyTorch, CuPy, or RAPIDS module was found during the first probe.
Current TWIRL S56 environment:

```text
/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56
```

Build/update it from the staged checkout with:

```bash
cd /orcd/data/mki_aryeh/001/twirl/code/TWIRL
bash scripts/orcd/build_orcd_env.sh
```

The environment spec is [env](../configs/orcd/twirl-s56-env.yml). The builder
writes `twirl-s56-build-info.txt` inside the environment prefix. It currently
provides Python `3.11`, the core scientific stack, `pyarrow`, `batman-package`,
`tess-point`, and editable TWIRL. The CPU Slurm wrappers default to:

```bash
TWIRL_ORCD_PYTHON=/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python
```

The optional PyTorch/H200 environment is separate from the CPU/ranker env:

```text
/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56-torch
```

Build/update it from the staged checkout with:

```bash
cd /orcd/data/mki_aryeh/001/twirl/code/TWIRL
bash scripts/orcd/build_orcd_torch_env.sh
```

The environment spec is [env](../configs/orcd/twirl-s56-torch-env.yml). The
builder installs the same scientific stack plus the official PyTorch CUDA
`cu128` wheel set and writes `twirl-s56-torch-build-info.txt` inside the
environment prefix. The first successful build installed `torch 2.11.0+cu128`;
login-node import smoke reports `torch.version.cuda=12.8`, while the H200 Slurm
smoke below verifies actual CUDA visibility.

Keep environment manifests and install logs with the ORCD results they produce.

The first H200 visibility smoke on `2026-06-30` requested exactly `1` H200 and
confirmed `nvidia-smi` visibility on `node4900`. The baseline `twirl-s56`
environment is CPU/scientific-Python oriented and does not yet include
`torch` or `cupy`.

The torch-required tensor smoke on `2026-07-01` requested exactly `1` H200 and
completed as Slurm job `16942919` on `node4900`. It read `128/128`
ranker-selected real-candidate windows from the compact S56 LC export, wrote
shape `(128, 3, 257)` tensors, and confirmed `torch.cuda.is_available()`,
`device_count=1`, and `device='cuda:0'`.

The bounded tensor training smoke wrapper is:

```bash
scripts/orcd/run_s56_orcd_pilot.sh --run h200-tensor-train-smoke
```

It consumes real teacher rows from the human-label audit and refuses to train if
they are absent or too sparse. The explicit environment-only synthetic smoke is:

```bash
scripts/orcd/run_s56_orcd_pilot.sh --run h200-tensor-train-synthetic-smoke
```

Synthetic smoke job `16944030` completed on `node4900` in `6 s` with one H200
and confirmed the trainer runs on `cuda:0`; it is not a scientific model.

## Starter Slurm Snippets

CPU smoke test:

```bash
#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-cpu-smoke
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 00:05:00

set -euo pipefail
hostname
python3 --version
```

Single-H200 smoke test:

```bash
#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-h200-smoke
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 00:10:00

set -euo pipefail
module load cuda/13.1.0
hostname
nvidia-smi
```

Two-H200 development job:

```bash
#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-h200-2gpu
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=256G
#SBATCH --gres=gpu:h200:2
#SBATCH -t 2-00:00:00

set -euo pipefail
module load cuda/13.1.0
hostname
nvidia-smi

# Replace with a TWIRL runner after the ORCD environment is pinned.
# Example shape:
# python scripts/stage3_injections/run_lc_injection_grid.py --config ...
```

Full-node `--gres=gpu:h200:8` jobs are intentionally not a default TWIRL
template. Use them only after explicit user approval for a specific production
run.

Use `sbatch --test-only` before submitting a new script:

```bash
sbatch --test-only -p pg_mki_aryeh -t 00:05:00 -c 1 --mem=1G --wrap="hostname"
```

## First TWIRL Pilot

1. Build compact S56 TWIRL-FS v2 export shards on PDO or locally from synced
   products.
2. Stage those shards and injection products under the ORCD checkout's
   `data_local/` tree so existing TWIRL scripts can run unchanged.
3. Run CPU and 1xH200 smoke tests.
4. Run S56 GPU search equivalence against the existing CPU BLS results.
5. Confirm WD 1856 recovery before expanding to multi-sector or injection
   grids.

Current staged S56 inputs on ORCD:

```text
/orcd/data/mki_aryeh/001/twirl/code/TWIRL/data_local/stage3_injections/s56_twirlfs_v2_lc_export/
/orcd/data/mki_aryeh/001/twirl/code/TWIRL/data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_20k_predetrend_batman_periodradius_bright_balanced/
```

The `20k` injection HDF5 was verified on `2026-06-30` with size
`11,283,081,701` bytes and `20,000` manifest/label rows.

Acceptance checks:

- source and export TIC/Gaia IDs match
- cadence counts match for representative targets

Current helper scripts:

- `scripts/orcd/stage_s56_orcd_inputs.sh` stages compact PDO S56 products into
  `/orcd/data/mki_aryeh/001/twirl/code/TWIRL/`. The ORCD code tree is now a
  sparse Git checkout tracking `origin/main`; staging updates code through Git
  by default and preserves local `data_local/`, `reports/`, and `logs/`
  artifacts.
- `scripts/orcd/bootstrap_orcd_git_checkout.sh` initializes or refreshes the
  ORCD checkout as a sparse Git worktree over HTTPS, because ORCD currently has
  no GitHub SSH deploy key. It uses non-interactive SSH to ORCD and
  non-prompting Git fetches.
- `scripts/orcd/build_orcd_env.sh` builds/updates the reusable ORCD Python
  environment under `/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56`.
- `scripts/orcd/check_s56_orcd_pipeline_status.sh` is a read-only local helper
  for monitoring the current S56 ORCD peak/ranker dependency chain through the
  user-opened control socket; it reports tracked-job GPU/GRES use, chunk
  completion fraction, active chunks, expected downstream artifacts, and the
  ORCD-layout handoff audit.
- `scripts/orcd/run_s56_orcd_pilot.sh` submits named ORCD stages through the
  user-opened control socket and enforces the default `<=2` H200 cap before
  submitting Slurm scripts.
- `scripts/orcd/slurm_s56_peak_training_cpu.sbatch` runs the injected
  peak-level BLS truth table on ORCD CPU nodes; use `TWIRL_N_INJECTIONS=100`
  for the first smoke. Smoke outputs go under `peak_training_orcd_smoke100`,
  while the full run goes under `peak_training_orcd_full`.
- `scripts/orcd/slurm_s56_peak_ranker_train_cpu.sbatch` trains/evaluates the
  lightweight BLS peak ranker after a peak table exists.
- `scripts/orcd/slurm_s56_peak_ranker_apply_cpu.sbatch` applies the trained
  peak ranker to the staged real S56 BLS peak table, verifies the real peak
  table, verifies the selected ephemerides, and writes top-N ephemerides per
  TIC for downstream LEO/human review on PDO.
- `scripts/orcd/slurm_s56_predetrend_detrending_bls_audit_cpu.sbatch` runs the
  production-level raw-flux detrending-strength audit on ORCD CPU nodes. This
  is the appropriate branch-selection test because it starts from stored
  original/injected raw aperture flux, applies alternative `FluxDetrendConfig`
  strengths, then reruns BLS and compares against injected truth. Its default
  method set now uses `twirl_fs_v2_adp015q` as the named candidate branch.
- `scripts/orcd/slurm_s56_adpplus_bls_audit_cpu.sbatch` runs the ADP+
  post-processing and two-aperture BLS diagnostic on ORCD CPU nodes. Submit via
  `scripts/orcd/run_s56_orcd_pilot.sh --run adpplus-smoke` first, then
  `adpplus-full` after the smoke summary passes. This job should use the
  compact ADP export HDF5 rather than walking PDO HLSP FITS. Treat its output
  as a vetting/display diagnostic, not the production detrending choice.
- `scripts/orcd/sync_s56_orcd_ranker_outputs_to_pdo.sh` copies only the small
  verified ORCD ranker-selected outputs back to the PDO checkout for LEO
  rendering; it refuses to sync until the expected selected-ephemerides and
  summary artifacts exist.
- `scripts/orcd/stage_s56_human_label_products_to_orcd.sh` stages only the
  small all-host human-label, teacher-row, and readiness products back to the
  ORCD checkout after the PDO post-label audit passes. It refuses a real run
  unless the readiness summary recommends a teacher smoke or object-teacher
  training.
- `scripts/orcd/slurm_h200_smoke.sbatch` checks CUDA/H200 visibility before any
  model-training job is submitted.
- checksums and provenance manifests are written
- WD 1856 period, epoch, duration, and ranking behavior are consistent with
  the existing S56 baseline

## References

- ORCD docs: <https://orcd-docs.mit.edu>
- ORCD resource requests:
  <https://orcd-docs.mit.edu/running-jobs/requesting-resources/>
- ORCD available resources:
  <https://orcd-docs.mit.edu/running-jobs/available-resources/>
- ORCD filesystems:
  <https://orcd-docs.mit.edu/filesystems-file-transfer/filesystems/>
- ORCD transfers:
  <https://orcd-docs.mit.edu/filesystems-file-transfer/transferring-files/>
