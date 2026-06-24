# ORCD H200 Usage For TWIRL

This note records the initial ORCD/Engaging reconnaissance for the MKI Aryeh
CPU/GPU nodes and how TWIRL should use them. It is operational guidance, not a
scientific pipeline decision.

Last verified: `2026-06-10`.

## Access And Scheduler Names

Use the ORCD login host:

```bash
ssh tehan@orcd-login.mit.edu
```

Login requires the account password plus Duo. Do not put passwords in scripts.
For repeated local probes, use an SSH control socket that you authenticate once
from your own terminal.

The administrative wording from Paul was:

```text
orcd_ug_pg_mki_aryeh_all
```

That is not the runnable Slurm partition name observed on ORCD. The runnable
partition is:

```bash
-p pg_mki_aryeh
```

Verification from `2026-06-10`:

- `sinfo` lists `pg_mki_aryeh`.
- `scontrol show partition orcd_ug_pg_mki_aryeh_all` reports the partition is
  not found.
- `sbatch --test-only -p orcd_ug_pg_mki_aryeh_all ...` fails with an invalid
  partition.
- `sbatch --test-only -p pg_mki_aryeh ...` succeeds.
- The observed Unix access group is `orcd_rg_par_pg_mki_aryeh`.

If ORCD staff later says a different Slurm name should be used, update this
file before changing any scripts.

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

At probe time, `node4900` was also visible through `mit_preemptable`, and all
8 GPUs were allocated by other users. Treat `pg_mki_aryeh` as the correct
partition, but do not assume the nodes are idle or reserved exclusively for
TWIRL.

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
Plan to build either:

- a project Conda environment under `/orcd/data/mki_aryeh/001/twirl/envs/`, or
- an Apptainer image with pinned CUDA/PyTorch/CuPy dependencies.

Keep environment manifests and install logs with the ORCD results they produce.

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

Full-node H200 job:

```bash
#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-h200-8gpu
#SBATCH -N 1
#SBATCH -c 120
#SBATCH --mem=0
#SBATCH --gres=gpu:h200:8
#SBATCH -t 2-00:00:00

set -euo pipefail
module load cuda/13.1.0
hostname
nvidia-smi

# Replace with a TWIRL runner after the ORCD environment is pinned.
# Example shape:
# python scripts/stage3_injections/run_lc_injection_grid.py --config ...
```

Use `sbatch --test-only` before submitting a new script:

```bash
sbatch --test-only -p pg_mki_aryeh -t 00:05:00 -c 1 --mem=1G --wrap="hostname"
```

## First TWIRL Pilot

1. Build compact S56 TWIRL-FS v2 export shards on PDO or locally from synced
   products.
2. Transfer those shards to `/orcd/data/mki_aryeh/001/twirl/exports/`.
3. Run CPU and 1xH200 smoke tests.
4. Run S56 GPU search equivalence against the existing CPU BLS results.
5. Confirm WD 1856 recovery before expanding to multi-sector or injection
   grids.

Acceptance checks:

- source and export TIC/Gaia IDs match
- cadence counts match for representative targets
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
