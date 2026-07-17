#!/bin/bash
# Reproduce A2v1 TGLC HDF5 light curves on PDO from existing source pickles.
#
# Usage:
#   run_a2v1_reproduction_pdo.sh <sector> <orbit:tag> [<orbit:tag> ...]
#
# Example:
#   run_a2v1_reproduction_pdo.sh 56 119:o1 120:o2
#   run_a2v1_reproduction_pdo.sh 94 195:o1 196:o2
#
# The output root stays sector-agnostic and follows the MIT TGLC orbit layout:
#   /pdo/users/tehan/tglc-gpu-production-A2v1/orbit-<orbit>/...
#
# This runner intentionally skips catalogs/cutouts. It builds source_tic overlays
# from the TWIRL observation table, symlinks existing source_*.pkl payloads into
# the A2v1 root, links old ePSFs only when the source cutout has an empty mask,
# then reruns masked ePSFs and lightcurves.

set -euo pipefail

REPO=${TWIRL_REPO:-/pdo/users/tehan/TWIRL}
export TWIRL_REPO="$REPO"
SOURCE_ROOT=${TWIRL_SOURCE_ROOT:-/pdo/users/tehan/tglc-gpu-production}
A2V1_ROOT=${TWIRL_A2V1_ROOT:-/pdo/users/tehan/tglc-gpu-production-A2v1}
TGLC_FORK=${TWIRL_TGLC_FORK:-/pdo/users/tehan/tess-gaia-light-curve-twirl}
QLP_PY=${TWIRL_QLP_PY:-/sw/qlp-environment/.venv/bin/python}
TGLC_PY=${TWIRL_TGLC_PY:-/pdo/users/tehan/twirl-gpu-venv/bin/python}
SCRIPT_DIR=${TWIRL_SCRIPT_DIR:-$REPO/scripts/stage1_lightcurves}
LOG_DIR=${TWIRL_A2V1_LOG_DIR:-$A2V1_ROOT/twirl_logs}

SECTOR=${1:-${TWIRL_A2V1_SECTOR:-56}}
if [ "$#" -gt 0 ]; then
  shift
fi

if [ "$#" -gt 0 ]; then
  ORBIT_SPECS=("$@")
elif [ -n "${TWIRL_A2V1_ORBIT_SPECS:-}" ]; then
  read -r -a ORBIT_SPECS <<< "$TWIRL_A2V1_ORBIT_SPECS"
elif [ "$SECTOR" = "56" ]; then
  ORBIT_SPECS=("119:o1" "120:o2")
else
  echo "ABORT: orbit specs are required for sector $SECTOR; use '<orbit>:<o1|o2>'" >&2
  exit 1
fi

SECTOR_TAG="s$(printf '%02d' "$SECTOR")"
RUN_LABEL=${TWIRL_A2V1_RUN_LABEL:-${SECTOR_TAG}-a2v1-gpu}
REPRO_LOG="$LOG_DIR/${SECTOR_TAG}_a2v1_reproduction.log"

export HDF5_USE_FILE_LOCKING=${HDF5_USE_FILE_LOCKING:-FALSE}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export VECLIB_MAXIMUM_THREADS=${VECLIB_MAXIMUM_THREADS:-1}
export NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS:-1}
export LD_LIBRARY_PATH=/sw/python-versions/python-3.11.9/lib:/pdo/app/anaconda/anaconda2-4.4.0/lib:${LD_LIBRARY_PATH:-}
export PYTHONPATH=$REPO/src:/sw/qlp-environment/.venv/lib/python3.11/site-packages:${PYTHONPATH:-}

"$REPO/scripts/assert_clean_checkout.sh" "$REPO"

mkdir -p "$A2V1_ROOT" "$LOG_DIR"

stamp() { date '+%Y-%m-%d %H:%M:%S %Z'; }
log()   { echo "[$(stamp)] $*" | tee -a "$REPRO_LOG"; }
abort() { log "ABORT: $1"; exit 1; }

require_hook() {
  local file=$1
  local pattern=$2
  local label=$3
  if ! grep -q "$pattern" "$file"; then
    abort "missing $label in $file"
  fi
}

parse_orbit_spec() {
  local spec=$1
  if [[ "$spec" != *:* ]]; then
    abort "orbit spec must have format '<orbit>:<orbit_tag>'; got '$spec'"
  fi
  PARSED_ORBIT="${spec%%:*}"
  PARSED_ORBIT_TAG="${spec#*:}"
  if ! [[ "$PARSED_ORBIT" =~ ^[0-9]+$ ]]; then
    abort "orbit must be numeric in spec '$spec'"
  fi
  if [ -z "$PARSED_ORBIT_TAG" ]; then
    abort "orbit tag is empty in spec '$spec'"
  fi
}

clean_selected_outputs() {
  log "cleaning selected A2v1 epsf/LC outputs under $A2V1_ROOT"
  local orbits=()
  local spec
  for spec in "${ORBIT_SPECS[@]}"; do
    parse_orbit_spec "$spec"
    orbits+=("$PARSED_ORBIT")
  done

  ORBIT_LIST="${orbits[*]}" A2V1_ROOT="$A2V1_ROOT" "$QLP_PY" - <<'PY'
import os
from pathlib import Path

root = Path(os.environ["A2V1_ROOT"])
orbits = os.environ["ORBIT_LIST"].split()
removed = 0
for orbit in orbits:
    patterns = [
        f"orbit-{orbit}/ffi/cam*/ccd*/epsf/epsf_*.npy",
        f"orbit-{orbit}/ffi/cam*/ccd*/LC/*.h5",
    ]
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
    --summary-json "$LOG_DIR/${SECTOR_TAG}_o${orbit}_source_tic_overlay_observations.json" \
    2>&1 | tee -a "$REPRO_LOG"
}

prefill_empty_mask_epsfs() {
  local orbit=$1
  if [ "${TWIRL_A2V1_PREFILL_EMPTY_MASK_EPSFS:-1}" = "0" ]; then
    log "skipping empty-mask ePSF prefill for orbit-$orbit"
    return 0
  fi
  log "prefilling empty-mask ePSFs for orbit-$orbit from $SOURCE_ROOT"
  PYTHONPATH="$REPO/src:$TGLC_FORK:$PYTHONPATH" \
  "$QLP_PY" "$SCRIPT_DIR/prefill_epsfs_for_empty_masks.py" \
    --source-tglc-data-dir "$SOURCE_ROOT" \
    --output-tglc-data-dir "$A2V1_ROOT" \
    --orbit "$orbit" \
    --apply \
    --nprocs "${TWIRL_A2V1_PREFILL_NPROCS:-8}" \
    --summary-json "$LOG_DIR/${SECTOR_TAG}_o${orbit}_empty_mask_epsf_prefill.json" \
    2>&1 | tee -a "$REPRO_LOG"
}

pick_gpu_args() {
  local gpu_args=""
  if [ -n "${TWIRL_A2V1_GPU_LIST:-}" ]; then
    gpu_args="--max-parallel-ccd-jobs ${TWIRL_A2V1_GPU_MAX_PARALLEL:-4} --gpu-list ${TWIRL_A2V1_GPU_LIST} --epsfs-nprocs ${TWIRL_A2V1_EPSFS_NPROCS:-4}"
  elif [ -x "$SCRIPT_DIR/pick_gpu_config.sh" ]; then
    gpu_args=$(bash "$SCRIPT_DIR/pick_gpu_config.sh" \
      "${TWIRL_A2V1_GPU_MIN_FREE_GB:-8}" \
      "${TWIRL_A2V1_GPU_MAX_PARALLEL:-4}" \
      "${TWIRL_A2V1_EPSFS_NPROCS:-4}")
  fi
  if [ -z "$gpu_args" ]; then
    gpu_args="--max-parallel-ccd-jobs ${TWIRL_A2V1_GPU_MAX_PARALLEL:-4} --gpu-list ${TWIRL_A2V1_GPU_LIST:-0,1,2,3} --epsfs-nprocs ${TWIRL_A2V1_EPSFS_NPROCS:-4}"
  fi
  echo "$gpu_args"
}

run_orbit() {
  local orbit=$1
  local orbit_tag=$2
  log "starting orbit-$orbit epsfs+lightcurves"
  local gpu_args
  gpu_args=$(pick_gpu_args)
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
    --catalogs-nprocs "${TWIRL_A2V1_CATALOGS_NPROCS:-16}" \
    --cutouts-nprocs "${TWIRL_A2V1_CUTOUTS_NPROCS:-16}" \
    --lightcurves-nprocs "${TWIRL_A2V1_LIGHTCURVES_NPROCS:-16}" \
    --run-label "$RUN_LABEL" \
    --stages epsfs,lightcurves \
    $gpu_args \
    2>&1 | tee -a "$LOG_DIR/${SECTOR_TAG}_a2v1_orbit_${orbit}.log"
}

log "$SECTOR_TAG A2v1 reproduction start"
log "source_root=$SOURCE_ROOT"
log "a2v1_root=$A2V1_ROOT"
log "orbit_specs=${ORBIT_SPECS[*]}"

require_hook "$TGLC_FORK/tglc/scripts/light_curves.py" "_source_tic_overlay_path" "source_tic lightcurve hook"
require_hook "$TGLC_FORK/tglc/scripts/epsfs.py" "source.mask.mask" "saturated-pixel ePSF mask"

if [ "${TWIRL_A2V1_CLEAN_OUTPUTS:-0}" = "1" ]; then
  clean_selected_outputs || abort "failed to clean selected A2v1 outputs"
fi

for spec in "${ORBIT_SPECS[@]}"; do
  parse_orbit_spec "$spec"
  build_overlays "$PARSED_ORBIT" || abort "orbit-$PARSED_ORBIT source_tic overlay failed"
done

for spec in "${ORBIT_SPECS[@]}"; do
  parse_orbit_spec "$spec"
  prefill_empty_mask_epsfs "$PARSED_ORBIT" || abort "orbit-$PARSED_ORBIT ePSF prefill failed"
done

for spec in "${ORBIT_SPECS[@]}"; do
  parse_orbit_spec "$spec"
  run_orbit "$PARSED_ORBIT" "$PARSED_ORBIT_TAG" || abort "orbit-$PARSED_ORBIT failed"
done

log "$SECTOR_TAG A2v1 HDF5 reproduction complete"
