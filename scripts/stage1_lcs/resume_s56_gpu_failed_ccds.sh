#!/bin/bash
# Resume the 3 orbit-119 CCDs that OOMed at the ePSF stage on 2026-04-28
# (cam3/ccd2, cam3/ccd4, cam4/ccd1). Source pickles already exist; we only
# need ePSF + lightcurves. Uses -n 4 (matches the proven 2026-04-27 benchmark)
# and runs serially to avoid the multi-GPU contention that triggered the
# original cupy.cuda.memory.OutOfMemoryError under the 4-parallel × 8-procs
# config.

set -euo pipefail

ROOT=/pdo/users/tehan/tglc-gpu-production
TICDIR=$ROOT/twirl_logs/s56-gpu-rerun/tic_lists
LOG=$ROOT/twirl_logs/s56-gpu-rerun

export PYTHONPATH=/sw/qlp-environment/.venv/lib/python3.11/site-packages
export LD_LIBRARY_PATH=/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}
PY=/pdo/users/tehan/twirl-gpu-venv/bin/python

ORBIT=119
for ccd in "3,2" "3,4" "4,1"; do
  cam="${ccd%,*}"; col="${ccd#*,}"
  label="cam${cam}_ccd${col}"
  ticlist="$TICDIR/orbit-${ORBIT}_${label}_wd_tics.txt"
  echo "=== $(date -Iseconds) $label epsfs ==="
  $PY -m tglc epsfs --orbit "$ORBIT" --ccd "$ccd" --nprocs 4 \
    --tglc-data-dir "$ROOT" \
    2>&1 | tee -a "$LOG/orbit-${ORBIT}_${label}_epsfs_resume.log"
  echo "=== $(date -Iseconds) $label lightcurves ==="
  mapfile -t tics < "$ticlist"
  $PY -m tglc lightcurves --orbit "$ORBIT" --ccd "$ccd" --nprocs 16 \
    --tglc-data-dir "$ROOT" --replace \
    --tic "${tics[@]}" \
    2>&1 | tee -a "$LOG/orbit-${ORBIT}_${label}_lightcurves_resume.log"
  echo "=== $(date -Iseconds) $label DONE ==="
done

echo "=== $(date -Iseconds) ALL 3 RESUMED CCDs DONE ==="
