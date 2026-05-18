#!/bin/bash
# CPU-only prep step: run catalogs + cutouts for both orbits of a sector.
# Designed to run on pdogpu1 (no GPU drivers) so pdogpu6 stays free for ePSF.
# After both orbits' cutouts succeed, writes
#   /pdo/users/tehan/tglc-gpu-production/markers/s<NN>_cutouts_done.flag
# which the finalize side polls for.
#
# Usage:
#   prep_sector_cpu.sh <sector> <orbit_1> <orbit_2>
# Example:
#   prep_sector_cpu.sh 58 123 124

set -uo pipefail

SECTOR="${1:?usage: $0 <sector> <orbit_1> <orbit_2>}"
ORBIT_1="${2:?usage: $0 <sector> <orbit_1> <orbit_2>}"
ORBIT_2="${3:?usage: $0 <sector> <orbit_1> <orbit_2>}"

REPO=/pdo/users/tehan/TWIRL
ROOT=/pdo/users/tehan/tglc-gpu-production
LOG=$ROOT/twirl_logs
MARKERS=$ROOT/markers
mkdir -p "$ROOT" "$LOG" "$MARKERS"

# Same env machinery as the GPU launcher. The TGLC venv has all deps via
# --system-site-packages from /sw/qlp-environment/.venv (where catalogs deps
# like pyticdb live). On a node without GPU drivers, ePSFs would fail at
# import — but we only request --stages catalogs,cutouts here.
export PYTHONPATH=/sw/qlp-environment/.venv/lib/python3.11/site-packages
# /pdo/app/anaconda.../lib provides libffi.so.6 needed on pdogpu1 (CLAUDE.md).
# Harmless on pdogpu6.
export LD_LIBRARY_PATH=/sw/python-versions/python-3.11.9/lib:/pdo/app/anaconda/anaconda2-4.4.0/lib:${LD_LIBRARY_PATH:-}

cd "$REPO"

run_orbit() {
  local ORBIT=$1
  local ORBIT_TAG=$2
  local EXTRA=$3
  echo "[prep] sector=$SECTOR orbit=$ORBIT tag=$ORBIT_TAG at $(date -Iseconds)"
  /pdo/users/tehan/twirl-gpu-venv/bin/python \
    scripts/stage1_lightcurves/run_tglc_orbit_pipeline.py \
    --orbit "$ORBIT" --sector "$SECTOR" --orbit-tag "$ORBIT_TAG" \
    --tglc-data-dir "$ROOT" \
    --log-dir "$LOG" \
    --python-bin /pdo/users/tehan/twirl-gpu-venv/bin/python \
    --ld-library-prefix /sw/python-versions/python-3.11.9/lib \
    --max-parallel-ccd-jobs 3 \
    --catalogs-nprocs 16 \
    --cutouts-nprocs 16 \
    --max-magnitude 20 \
    --run-label "s$(printf '%02d' "$SECTOR")-prep" \
    --stages catalogs,cutouts \
    $EXTRA
}

MAX_RETRIES=5
attempt=0
while (( attempt < MAX_RETRIES )); do
  attempt=$(( attempt + 1 ))
  echo "[prep-launcher] sector=$SECTOR attempt=$attempt/$MAX_RETRIES at $(date -Iseconds)"

  if ! run_orbit "$ORBIT_1" o1 ""; then
    echo "[prep-launcher] orbit-$ORBIT_1 failed; retry in 60s"
    sleep 60; continue
  fi
  if ! run_orbit "$ORBIT_2" o2 "--reuse-catalogs-from-orbit $ORBIT_1"; then
    echo "[prep-launcher] orbit-$ORBIT_2 failed; retry in 60s"
    sleep 60; continue
  fi

  # Verify each CCD has 196 source pkls before declaring sector ready.
  bad=0
  for orbit in "$ORBIT_1" "$ORBIT_2"; do
    for cam in 1 2 3 4; do
      for ccd in 1 2 3 4; do
        n=$(ls "$ROOT/orbit-$orbit/ffi/cam$cam/ccd$ccd/source/" 2>/dev/null | wc -l)
        if (( n != 196 )); then
          echo "[prep-launcher] WARN orbit-$orbit cam$cam/ccd$ccd source=$n (expect 196)"
          bad=1
        fi
      done
    done
  done
  if (( bad == 1 )); then
    echo "[prep-launcher] some CCDs missing source pkls; retry"
    continue
  fi

  marker=$MARKERS/s$(printf '%02d' "$SECTOR")_cutouts_done.flag
  date -Iseconds > "$marker"
  echo "[prep-launcher] sector=$SECTOR cutouts complete; wrote $marker"
  exit 0
done

echo "[prep-launcher] giving up on sector=$SECTOR after $MAX_RETRIES attempts" >&2
exit 1
