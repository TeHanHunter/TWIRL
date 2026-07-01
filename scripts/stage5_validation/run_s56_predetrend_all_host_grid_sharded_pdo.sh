#!/usr/bin/env bash
# Build all-host S56 pre-detrend injections as restartable TIC-grouped shards.
set -euo pipefail

log() {
  printf '[s56-allhost-sharded] %s\n' "$*" >&2
}

PYTHON_BIN="${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}"
export LD_LIBRARY_PATH="/sw/python-versions/python-3.11.9/lib:/pdo/app/anaconda/anaconda2-4.4.0/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${PYTHONPATH:-src}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

OUT_ROOT="${OUT_ROOT:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_allhost_predetrend_batman_periodradius_grid_sharded}"
SHARD_DIR="${SHARD_DIR:-${OUT_ROOT}/target_shards}"
CHUNK_ROOT="${CHUNK_ROOT:-${OUT_ROOT}/chunks}"
REPORT_DIR="${REPORT_DIR:-reports/stage5_validation/s56_allhost_predetrend_batman_periodradius_grid_sharded_pdo}"
SHARD_SIZE="${SHARD_SIZE:-250}"
CHUNK_JOBS="${CHUNK_JOBS:-2}"
REBUILD_SHARDS="${REBUILD_SHARDS:-0}"
REBUILD_CHUNKS="${REBUILD_CHUNKS:-0}"
LIMIT_TARGETS="${LIMIT_TARGETS:-}"
RANDOM_STATE="${RANDOM_STATE:-5656}"
GRID_PERIOD_RANGE_D="${GRID_PERIOD_RANGE_D:-0.12,13.0}"
GRID_RADIUS_RANGE_REARTH="${GRID_RADIUS_RANGE_REARTH:-0.18,18.0}"
GRID_PERIOD_BINS="${GRID_PERIOD_BINS:-50}"
GRID_RADIUS_BINS="${GRID_RADIUS_BINS:-50}"
APERTURES="${APERTURES:-DET_FLUX_ADP_SML DET_FLUX_SML}"
ORBIT_ROOTS="${ORBIT_ROOTS:-/pdo/users/tehan/tglc-gpu-production/orbit-119/ffi /pdo/users/tehan/tglc-gpu-production/orbit-120/ffi}"
EXPORT_MANIFEST="${EXPORT_MANIFEST:-data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_lc_export_pdo.manifest.json}"

mkdir -p "${OUT_ROOT}" "${CHUNK_ROOT}" "${REPORT_DIR}"

log "out_root=${OUT_ROOT}"
log "shard_size=${SHARD_SIZE} chunk_jobs=${CHUNK_JOBS}"

if [[ "${REBUILD_SHARDS}" == "1" || ! -s "${SHARD_DIR}/manifest.json" ]]; then
  log "building TIC-grouped target shards"
  shard_args=(
    scripts/stage5_validation/build_s56_predetrend_target_shards.py
    --out-dir "${SHARD_DIR}"
    --shard-size "${SHARD_SIZE}"
    --overwrite
  )
  # shellcheck disable=SC2206
  orbit_roots_array=(${ORBIT_ROOTS})
  shard_args+=(--orbit-roots "${orbit_roots_array[@]}")
  if [[ -n "${LIMIT_TARGETS}" ]]; then
    shard_args+=(--limit-targets "${LIMIT_TARGETS}")
  fi
  "${PYTHON_BIN}" "${shard_args[@]}"
else
  log "reusing target shards: ${SHARD_DIR}/manifest.json"
fi

pending_list="${OUT_ROOT}/pending_shards.txt"
: > "${pending_list}"
for meta in "${SHARD_DIR}"/chunk_*.meta.json; do
  [[ -e "${meta}" ]] || continue
  chunk_name="$(basename "${meta}" .meta.json)"
  chunk_out="${CHUNK_ROOT}/${chunk_name}"
  if [[ "${REBUILD_CHUNKS}" != "1" && -s "${chunk_out}/summary.json" && -s "${chunk_out}/injected_lightcurves.h5" ]]; then
    log "skip ${chunk_name}: completed"
    continue
  fi
  printf '%s\n' "${meta}" >> "${pending_list}"
done

n_pending="$(wc -l < "${pending_list}" | tr -d ' ')"
log "pending shards=${n_pending}"
if [[ "${n_pending}" != "0" ]]; then
  export PYTHON_BIN LD_LIBRARY_PATH PYTHONPATH
  export CHUNK_ROOT RANDOM_STATE GRID_PERIOD_RANGE_D GRID_RADIUS_RANGE_REARTH GRID_PERIOD_BINS GRID_RADIUS_BINS APERTURES ORBIT_ROOTS EXPORT_MANIFEST
  xargs -n 1 -P "${CHUNK_JOBS}" bash -c '
    set -euo pipefail
    meta="$1"
    chunk_name="$(basename "${meta}" .meta.json)"
    list_path="${meta%.meta.json}.txt"
    chunk_out="${CHUNK_ROOT}/${chunk_name}"
    n_tics="$("${PYTHON_BIN}" -c "import json,sys; print(json.load(open(sys.argv[1]))[\"n_tics\"])" "${meta}")"
    start_index="$("${PYTHON_BIN}" -c "import json,sys; print(json.load(open(sys.argv[1]))[\"start_index\"])" "${meta}")"
    seed=$((RANDOM_STATE + start_index))
    echo "[$(date -Is)] [s56-allhost-sharded] run ${chunk_name}: n_tics=${n_tics} start_index=${start_index}"
    OUT_DIR="${chunk_out}" \
    TARGET_H5_LIST="${list_path}" \
    N_INJECTIONS="${n_tics}" \
    RANDOM_STATE="${seed}" \
    INJECTION_INDEX_OFFSET="${start_index}" \
    INJECTION_ID_PREFIX="predet" \
    RUN_HOST_AUDIT=0 \
    REBUILD_INJECTIONS=1 \
    GRID_PERIOD_RANGE_D="${GRID_PERIOD_RANGE_D}" \
    GRID_RADIUS_RANGE_REARTH="${GRID_RADIUS_RANGE_REARTH}" \
    GRID_PERIOD_BINS="${GRID_PERIOD_BINS}" \
    GRID_RADIUS_BINS="${GRID_RADIUS_BINS}" \
    APERTURES="${APERTURES}" \
    ORBIT_ROOTS="${ORBIT_ROOTS}" \
      bash scripts/stage5_validation/run_s56_predetrend_all_host_grid_pdo.sh
  ' _ < "${pending_list}"
fi

log "merging shard metadata"
"${PYTHON_BIN}" scripts/stage5_validation/merge_s56_predetrend_injection_shards.py \
  --chunk-root "${CHUNK_ROOT}" \
  --out-dir "${OUT_ROOT}"

log "running combined host-coverage audit"
"${PYTHON_BIN}" scripts/stage5_validation/audit_s56_injection_host_coverage.py \
  --export-manifest "${EXPORT_MANIFEST}" \
  --injection-manifest "${OUT_ROOT}/injection_manifest.csv" \
  --injection-summary "${OUT_ROOT}/summary.json" \
  --out-dir "${REPORT_DIR}/host_coverage" \
  --expected-total-targets 19072

log "complete"
log "combined_manifest=${OUT_ROOT}/injection_manifest.csv"
log "combined_summary=${OUT_ROOT}/summary.json"
log "host_audit=${REPORT_DIR}/host_coverage/summary.md"
