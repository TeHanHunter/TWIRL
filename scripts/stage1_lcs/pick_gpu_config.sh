#!/bin/bash
# pick_gpu_config.sh -- emit run_tglc_orbit_pipeline.py args based on the
# currently-free GPU memory on this host. Lets a launcher scale automatically:
# if other users grab GPUs, we cut concurrency; if the box is quiet, we use
# all 8 GPUs.
#
# Usage:
#   pick_gpu_config.sh [min_free_gb] [max_parallel_cap] [epsfs_nprocs]
#
# Defaults: 8 GB free, max-parallel cap 4, epsfs-nprocs 4.
#
# Output (to stdout, one line, space-separated; designed to be inserted via
# $(...) in the launcher):
#
#   --max-parallel-ccd-jobs N --gpu-list X,Y,Z --epsfs-nprocs M
#
# Strategy:
#   - Pick GPUs whose free memory >= min_free_gb GB.
#   - max_parallel = min(max_parallel_cap, num_free_gpus).
#   - --gpu-list = first max_parallel GPU indices, comma-separated.
#     Each CCD subprocess pins to one GPU via CUDA_VISIBLE_DEVICES (driver
#     handles this), and runs --epsfs-nprocs workers all on that one GPU.
#   - --epsfs-nprocs = epsfs_nprocs (worker count per pinned GPU).
#
# Diagnostic info goes to stderr; the stdout line is purely the args.

set -euo pipefail

MIN_FREE_GB="${1:-8}"
MAX_PARALLEL_CAP="${2:-4}"
EPSFS_NPROCS="${3:-4}"

if ! command -v nvidia-smi >/dev/null; then
  echo "[pick-gpu-config] no nvidia-smi available; emitting --no-gpu config" >&2
  echo "--max-parallel-ccd-jobs 1 --epsfs-nprocs 16 --no-gpu"
  exit 0
fi

# Per-GPU "index, free_mib" lines.
free_list=$(nvidia-smi --query-gpu=index,memory.free --format=csv,noheader,nounits)
min_free_mib=$(awk -v g="$MIN_FREE_GB" 'BEGIN{printf "%d", g*1024}')

free_idxs=()
while IFS=',' read -r idx free_mib; do
  idx=$(echo "$idx" | tr -d ' ')
  free_mib=$(echo "$free_mib" | tr -d ' ')
  if (( free_mib >= min_free_mib )); then
    free_idxs+=("$idx")
  fi
done <<< "$free_list"

n_free=${#free_idxs[@]}
if (( n_free == 0 )); then
  echo "[pick-gpu-config] no GPUs with >= ${MIN_FREE_GB} GB free; falling back to GPU 0 with --epsfs-nprocs 2" >&2
  echo "--max-parallel-ccd-jobs 1 --gpu-list 0 --epsfs-nprocs 2"
  exit 0
fi

max_parallel=$(( n_free < MAX_PARALLEL_CAP ? n_free : MAX_PARALLEL_CAP ))
gpu_list=$(IFS=,; echo "${free_idxs[*]:0:$max_parallel}")

echo "[pick-gpu-config] free_gpus=$n_free (idx=${free_idxs[*]}) max_parallel=$max_parallel pinned=$gpu_list epsfs_nprocs=$EPSFS_NPROCS" >&2
echo "--max-parallel-ccd-jobs $max_parallel --gpu-list $gpu_list --epsfs-nprocs $EPSFS_NPROCS"
