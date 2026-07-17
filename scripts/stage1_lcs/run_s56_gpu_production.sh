#!/bin/bash
# Production GPU pipeline driver for one orbit (sector 56). Designed to be
# survivable: an outer retry loop auto-resumes after host reboots or driver
# errors, and the driver itself auto-skips CCDs whose epsf+LC outputs already
# exist (so reruns are cheap).
#
# Output root: /pdo/users/tehan/tglc-gpu-production/  (parallel to tglc-deep-catalogs/).
# The existing /pdo/users/tehan/tglc-deep-catalogs/hlsp_s0056/ tree is untouched.
#
# Pipeline per orbit (16 CCDs each):
#   tglc catalogs -> cutouts -> epsfs (GPU) -> lightcurves
# After both orbits finish (separate launches):
#   qlp lctools detrend (sector-level; see scripts/detrend_wrapper.py)
#   scripts/stage1_lcs/fix_tglc_quat_cadence_mismatch.py (cam3/orbit-120 trim)
#   qlp lctools hlsp (sector-level; see scripts/hlsp_wrapper.py, pass --basedir)
#
# Per-stage CLI overrides for the GPU path:
#   --python-bin /pdo/users/tehan/twirl-gpu-venv/bin/python   (HAS_CUPY=True)
#   --ld-library-prefix /sw/python-versions/python-3.11.9/lib (libpython3.11)
# Pre-set PYTHONPATH so non-epsf stages still find their deps inside the QLP venv.
#
# GPU concurrency is auto-tuned via pick_gpu_config.sh, which queries
# nvidia-smi for free GPU memory and emits --max-parallel-ccd-jobs / --gpu-list
# / --epsfs-nprocs accordingly. If other users are running heavy GPU jobs,
# we shrink to whatever's free; if the box is quiet, we use up to 4 CCDs in
# parallel (one GPU pinned per CCD with -n 4 ePSF workers).
#
# Usage: run_s56_gpu_production.sh <orbit:119|120> [extra args passed to driver]

set -uo pipefail

ORBIT="${1:?usage: $0 <orbit:119|120> [extra driver args]}"
shift || true
case "$ORBIT" in
  119) ORBIT_TAG=o1 ;;
  120) ORBIT_TAG=o2 ;;
  *)   echo "orbit must be 119 or 120"; exit 2 ;;
esac

REPO=/pdo/users/tehan/TWIRL
ROOT=/pdo/users/tehan/tglc-gpu-production
LOG=$ROOT/twirl_logs

mkdir -p "$ROOT" "$LOG"

# Pre-seed PYTHONPATH so the GPU venv (--system-site-packages off /sw/python-versions/)
# can resolve tglc's runtime deps (pyticdb, astroquery, etc.) from the QLP venv.
export PYTHONPATH=/sw/qlp-environment/.venv/lib/python3.11/site-packages
# The GPU venv's python needs libpython3.11.so.1.0 from /sw/python-versions/.
export LD_LIBRARY_PATH=/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}

cd "$REPO"

MAX_RETRIES=10
attempt=0
while (( attempt < MAX_RETRIES )); do
  attempt=$(( attempt + 1 ))
  echo "[launcher] orbit=$ORBIT attempt=$attempt/$MAX_RETRIES at $(date -Iseconds)"

  # Re-pick GPU concurrency on every attempt: lets the next retry adapt to
  # whatever's free now (handy after a host reboot or when other users join).
  GPU_ARGS=$(bash "$REPO/scripts/stage1_lcs/pick_gpu_config.sh" 8 4 4)
  echo "[launcher] gpu_args: $GPU_ARGS"

  /pdo/users/tehan/twirl-gpu-venv/bin/python \
    scripts/stage1_lcs/run_tglc_orbit_pipeline.py \
    --orbit "$ORBIT" --sector 56 --orbit-tag "$ORBIT_TAG" \
    --tglc-data-dir "$ROOT" \
    --log-dir "$LOG" \
    --python-bin /pdo/users/tehan/twirl-gpu-venv/bin/python \
    --ld-library-prefix /sw/python-versions/python-3.11.9/lib \
    --catalogs-nprocs 16 \
    --cutouts-nprocs 16 \
    --lightcurves-nprocs 16 \
    --max-magnitude 20 \
    --run-label "s56-gpu-rerun" \
    $GPU_ARGS \
    "$@"
  rc=$?

  if (( rc == 0 )); then
    echo "[launcher] orbit=$ORBIT completed cleanly on attempt $attempt"
    exit 0
  fi

  echo "[launcher] driver exited with rc=$rc; auto-resume will skip done CCDs. Sleeping 60s before retry."
  sleep 60
done

echo "[launcher] giving up after $MAX_RETRIES attempts" >&2
exit 1
