#!/bin/bash
# Reproduce S56 A2v1 TGLC HDF5 light curves on PDO from existing source pickles.
#
# This runner intentionally skips catalogs/cutouts. It builds source_tic overlays
# from the TWIRL observation table, symlinks the existing source_*.pkl payloads
# into a user-owned A2v1 root, then reruns ePSFs and lightcurves.

set -uo pipefail

REPO=${TWIRL_REPO:-/pdo/users/tehan/TWIRL}
SOURCE_ROOT=${TWIRL_SOURCE_ROOT:-/pdo/users/tehan/tglc-gpu-production}
A2V1_ROOT=${TWIRL_A2V1_ROOT:-/pdo/users/tehan/tglc-gpu-production-A2v1}
TGLC_FORK=${TWIRL_TGLC_FORK:-/pdo/users/tehan/tess-gaia-light-curve-twirl}
QLP_PY=${TWIRL_QLP_PY:-/sw/qlp-environment/.venv/bin/python}
TGLC_PY=${TWIRL_TGLC_PY:-/pdo/users/tehan/twirl-gpu-venv/bin/python}
SCRIPT_DIR=${TWIRL_SCRIPT_DIR:-$REPO/scripts/stage1_lightcurves}
LOG_DIR=${TWIRL_A2V1_LOG_DIR:-$A2V1_ROOT/twirl_logs}
RUN_LABEL=${TWIRL_A2V1_RUN_LABEL:-s56-a2v1-gpu}
SECTOR=56
ORBIT_1=119
ORBIT_2=120

export HDF5_USE_FILE_LOCKING=${HDF5_USE_FILE_LOCKING:-FALSE}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export VECLIB_MAXIMUM_THREADS=${VECLIB_MAXIMUM_THREADS:-1}
export NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS:-1}
export LD_LIBRARY_PATH=/sw/python-versions/python-3.11.9/lib:/pdo/app/anaconda/anaconda2-4.4.0/lib:${LD_LIBRARY_PATH:-}
export PYTHONPATH=/pdo/users/tehan/TWIRL/src:/sw/qlp-environment/.venv/lib/python3.11/site-packages:${PYTHONPATH:-}

mkdir -p "$A2V1_ROOT" "$LOG_DIR"

stamp() { date '+%Y-%m-%d %H:%M:%S %Z'; }
log()   { echo "[$(stamp)] $*" | tee -a "$LOG_DIR/s56_a2v1_reproduction.log"; }
abort() { log "ABORT: $1"; exit 1; }

require_hook() {
  local file=$1
  local pattern=$2
  local label=$3
  if ! grep -q "$pattern" "$file"; then
    abort "missing $label in $file"
  fi
}

clean_smoke_outputs() {
  log "cleaning previous A2v1 epsf/LC outputs under $A2V1_ROOT"
  A2V1_ROOT="$A2V1_ROOT" "$QLP_PY" - <<'PY'
import os
from pathlib import Path
root = Path(os.environ["A2V1_ROOT"])
patterns = [
    "orbit-*/ffi/cam*/ccd*/epsf/epsf_*.npy",
    "orbit-*/ffi/cam*/ccd*/LC/*.h5",
    "twirl_logs/s56-a2v1-gpu/orbit-*_summary.json",
    "twirl_logs/s56-a2v1-gpu/run_summary.json",
]
removed = 0
for pattern in patterns:
    for path in root.glob(pattern):
        if path.exists() or path.is_symlink():
            path.unlink()
            removed += 1
print(f"removed={removed}")
PY
}

build_overlays() {
  local orbit=$1
  log "building observation-table source_tic overlays for orbit-$orbit"
  PYTHONPATH="$REPO/src:$TGLC_FORK:$PYTHONPATH" \
  "$QLP_PY" "$SCRIPT_DIR/build_source_tic_overlays.py" \
    --source-tglc-data-dir "$SOURCE_ROOT" \
    --output-tglc-data-dir "$A2V1_ROOT" \
    --orbit "$orbit" \
    --overlay-from-observations \
    --link-sources \
    --apply \
    --overwrite \
    --summary-json "$LOG_DIR/s56_o${orbit}_source_tic_overlay_observations.json" \
    2>&1 | tee -a "$LOG_DIR/s56_a2v1_reproduction.log"
}

run_orbit() {
  local orbit=$1
  local orbit_tag=$2
  log "starting orbit-$orbit epsfs+lightcurves"
  local gpu_args
  if [ -x "$SCRIPT_DIR/pick_gpu_config.sh" ]; then
    gpu_args=$(bash "$SCRIPT_DIR/pick_gpu_config.sh" 8 4 4)
  else
    gpu_args=""
  fi
  if [ -z "$gpu_args" ]; then
    gpu_args="--max-parallel-ccd-jobs 4 --gpu-list 0,1,2,3 --epsfs-nprocs 4"
  fi
  log "orbit-$orbit gpu_args: $gpu_args"
  HDF5_USE_FILE_LOCKING=FALSE \
  PYTHONPATH=/sw/qlp-environment/.venv/lib/python3.11/site-packages \
  "$TGLC_PY" "$SCRIPT_DIR/run_tglc_orbit_pipeline.py" \
    --orbit "$orbit" \
    --sector "$SECTOR" \
    --orbit-tag "$orbit_tag" \
    --tglc-data-dir "$A2V1_ROOT" \
    --log-dir "$LOG_DIR" \
    --python-bin "$TGLC_PY" \
    --fork-path "$TGLC_FORK" \
    --ld-library-prefix /sw/python-versions/python-3.11.9/lib \
    --catalogs-nprocs 16 \
    --cutouts-nprocs 16 \
    --lightcurves-nprocs 16 \
    --run-label "$RUN_LABEL" \
    --stages epsfs,lightcurves \
    $gpu_args \
    2>&1 | tee -a "$LOG_DIR/s56_a2v1_orbit_${orbit}.log"
}

log "S56 A2v1 reproduction start"
log "source_root=$SOURCE_ROOT"
log "a2v1_root=$A2V1_ROOT"

require_hook "$TGLC_FORK/tglc/scripts/light_curves.py" "_source_tic_overlay_path" "source_tic lightcurve hook"
require_hook "$TGLC_FORK/tglc/scripts/epsfs.py" "source.mask.mask" "saturated-pixel ePSF mask"

if [ "${TWIRL_A2V1_CLEAN_OUTPUTS:-0}" = "1" ]; then
  clean_smoke_outputs
fi

build_overlays "$ORBIT_1"
build_overlays "$ORBIT_2"

run_orbit "$ORBIT_1" o1 || abort "orbit-$ORBIT_1 failed"
run_orbit "$ORBIT_2" o2 || abort "orbit-$ORBIT_2 failed"

log "S56 A2v1 HDF5 reproduction complete"
