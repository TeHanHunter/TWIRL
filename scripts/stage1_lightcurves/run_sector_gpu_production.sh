#!/bin/bash
# Generic per-orbit GPU production launcher. Drop-in successor to
# run_s56_gpu_production.sh, parameterized for any sector/orbit.
#
# Usage:
#   run_sector_gpu_production.sh <sector> <orbit> <orbit_tag> [extra args]
# Example:
#   run_sector_gpu_production.sh 57 121 o1            # first orbit of S57
#   run_sector_gpu_production.sh 57 122 o2            # second orbit of S57
#   run_sector_gpu_production.sh 57 122 o2 --reuse-catalogs-from-orbit 121
#
# Output root: /pdo/users/tehan/tglc-gpu-production/orbit-<orbit>/...
# Outer until-retry handles host reboots and transient driver crashes; the
# driver auto-skips CCDs whose epsf+LC outputs already exist.
#
# GPU concurrency is auto-tuned via pick_gpu_config.sh on every attempt.

set -uo pipefail

SECTOR="${1:?usage: $0 <sector> <orbit> <orbit_tag> [driver args]}"
ORBIT="${2:?usage: $0 <sector> <orbit> <orbit_tag> [driver args]}"
ORBIT_TAG="${3:?usage: $0 <sector> <orbit> <orbit_tag> [driver args]}"
shift 3 || true

REPO=/pdo/users/tehan/TWIRL
ROOT=/pdo/users/tehan/tglc-gpu-production
LOG=$ROOT/twirl_logs
RUN_LABEL="s$(printf '%02d' "$SECTOR")-gpu"

mkdir -p "$ROOT" "$LOG"

# QLP venv site-packages so non-epsf stages resolve their deps; libpython
# from /sw/python-versions for the GPU venv.
export PYTHONPATH=/sw/qlp-environment/.venv/lib/python3.11/site-packages
export LD_LIBRARY_PATH=/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}

cd "$REPO"

MAX_RETRIES=10
attempt=0
while (( attempt < MAX_RETRIES )); do
  attempt=$(( attempt + 1 ))
  echo "[launcher] sector=$SECTOR orbit=$ORBIT attempt=$attempt/$MAX_RETRIES at $(date -Iseconds)"

  GPU_ARGS=$(bash "$REPO/scripts/stage1_lightcurves/pick_gpu_config.sh" 8 4 4)
  echo "[launcher] gpu_args: $GPU_ARGS"

  /pdo/users/tehan/twirl-gpu-venv/bin/python \
    scripts/stage1_lightcurves/run_tglc_orbit_pipeline.py \
    --orbit "$ORBIT" --sector "$SECTOR" --orbit-tag "$ORBIT_TAG" \
    --tglc-data-dir "$ROOT" \
    --log-dir "$LOG" \
    --python-bin /pdo/users/tehan/twirl-gpu-venv/bin/python \
    --ld-library-prefix /sw/python-versions/python-3.11.9/lib \
    --catalogs-nprocs 16 \
    --cutouts-nprocs 16 \
    --lightcurves-nprocs 16 \
    --run-label "$RUN_LABEL" \
    $GPU_ARGS \
    "$@"
  rc=$?

  if (( rc == 0 )); then
    echo "[launcher] sector=$SECTOR orbit=$ORBIT completed cleanly on attempt $attempt"
    exit 0
  fi

  echo "[launcher] driver exited with rc=$rc; auto-resume will skip done CCDs. Sleeping 60s before retry."
  sleep 60
done

echo "[launcher] giving up after $MAX_RETRIES attempts" >&2
exit 1
