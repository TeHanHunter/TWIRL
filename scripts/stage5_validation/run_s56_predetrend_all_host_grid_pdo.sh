#!/usr/bin/env bash
# Build a coverage-first S56 pre-detrend BATMAN injection set on PDO.
#
# This is the "all hosts first" companion to the balanced 20k recovery grid:
# it shuffles the discovered raw S56 TIC hosts and attempts one injection per
# host before repeating any host. The summary records achieved unique-host
# coverage, so failures are visible instead of being hidden by random resampling.
set -euo pipefail

log() {
  printf '[s56-allhost-injections] %s\n' "$*" >&2
}

PYTHON_BIN="${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}"
export LD_LIBRARY_PATH="/sw/python-versions/python-3.11.9/lib:/pdo/app/anaconda/anaconda2-4.4.0/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${PYTHONPATH:-src}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

OUT_DIR="${OUT_DIR:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_allhost_predetrend_batman_periodradius_grid}"
N_INJECTIONS="${N_INJECTIONS:-19072}"
REBUILD_INJECTIONS="${REBUILD_INJECTIONS:-0}"
RANDOM_STATE="${RANDOM_STATE:-5656}"
INJECTION_INDEX_OFFSET="${INJECTION_INDEX_OFFSET:-0}"
INJECTION_ID_PREFIX="${INJECTION_ID_PREFIX:-predet}"
MAX_ATTEMPTS_PER_INJECTION="${MAX_ATTEMPTS_PER_INJECTION:-5}"
MIN_IN_TRANSIT="${MIN_IN_TRANSIT:-2}"
BASELINE_SOURCE="${BASELINE_SOURCE:-raw_median}"

GRID_PERIOD_RANGE_D="${GRID_PERIOD_RANGE_D:-0.12,13.0}"
GRID_RADIUS_RANGE_REARTH="${GRID_RADIUS_RANGE_REARTH:-0.18,18.0}"
GRID_PERIOD_BINS="${GRID_PERIOD_BINS:-50}"
GRID_RADIUS_BINS="${GRID_RADIUS_BINS:-50}"

APERTURES=(${APERTURES:-DET_FLUX_ADP_SML DET_FLUX_SML})
ORBIT_ROOTS=(${ORBIT_ROOTS:-/pdo/users/tehan/tglc-gpu-production/orbit-119/ffi /pdo/users/tehan/tglc-gpu-production/orbit-120/ffi})
TARGET_H5_LIST="${TARGET_H5_LIST:-}"
EXPORT_MANIFEST="${EXPORT_MANIFEST:-data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_lc_export_pdo.manifest.json}"
HOST_AUDIT_DIR="${HOST_AUDIT_DIR:-reports/stage5_validation/s56_allhost_predetrend_batman_periodradius_grid_pdo/host_coverage}"
RUN_HOST_AUDIT="${RUN_HOST_AUDIT:-1}"

log "python=${PYTHON_BIN}"
log "out_dir=${OUT_DIR}"
log "n_injections=${N_INJECTIONS} target_selection_mode=shuffled_cycle"
log "grid_period_range=${GRID_PERIOD_RANGE_D} bins=${GRID_PERIOD_BINS}"
log "grid_radius_range=${GRID_RADIUS_RANGE_REARTH} bins=${GRID_RADIUS_BINS}"
log "apertures=${APERTURES[*]}"
log "orbit_roots=${ORBIT_ROOTS[*]}"
if [[ -n "${TARGET_H5_LIST}" ]]; then
  log "target_h5_list=${TARGET_H5_LIST}"
fi

if [[ "${REBUILD_INJECTIONS}" == "1" || ! -s "${OUT_DIR}/injected_lightcurves.h5" ]]; then
  log "building coverage-first pre-detrend injection set"
  args=(
    scripts/stage3_injections/make_s56_predetrend_injection_set.py
    --out-dir "${OUT_DIR}"
    --orbit-roots "${ORBIT_ROOTS[@]}"
    --apertures "${APERTURES[@]}"
    --n-injections "${N_INJECTIONS}"
    --sampling-mode period_radius_grid
    --grid-period-range-d "${GRID_PERIOD_RANGE_D}"
    --grid-radius-range-rearth "${GRID_RADIUS_RANGE_REARTH}"
    --grid-period-bins "${GRID_PERIOD_BINS}"
    --grid-radius-bins "${GRID_RADIUS_BINS}"
    --target-selection-mode shuffled_cycle
    --random-state "${RANDOM_STATE}"
    --injection-index-offset "${INJECTION_INDEX_OFFSET}"
    --injection-id-prefix "${INJECTION_ID_PREFIX}"
    --min-in-transit "${MIN_IN_TRANSIT}"
    --max-attempts-per-injection "${MAX_ATTEMPTS_PER_INJECTION}"
    --baseline-source "${BASELINE_SOURCE}"
    --progress-every 250
    --overwrite
  )
  if [[ -n "${TARGET_H5_LIST}" ]]; then
    args+=(--target-h5-list "${TARGET_H5_LIST}")
  fi
  "${PYTHON_BIN}" "${args[@]}"
else
  log "reusing existing ${OUT_DIR}/injected_lightcurves.h5"
fi

if [[ "${RUN_HOST_AUDIT}" == "1" ]]; then
  log "running host-coverage audit"
  "${PYTHON_BIN}" scripts/stage5_validation/audit_s56_injection_host_coverage.py \
    --export-manifest "${EXPORT_MANIFEST}" \
    --injection-manifest "${OUT_DIR}/injection_manifest.csv" \
    --injection-summary "${OUT_DIR}/summary.json" \
    --out-dir "${HOST_AUDIT_DIR}" \
    --expected-total-targets 19072
fi

log "complete"
log "injections=${OUT_DIR}/injected_lightcurves.h5"
log "manifest=${OUT_DIR}/injection_manifest.csv"
log "host_audit=${HOST_AUDIT_DIR}/summary.md"
