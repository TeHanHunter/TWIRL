#!/usr/bin/env bash
# Build a dense S56 pre-detrend BATMAN injection sample and run the canonical
# small-aperture BLS recovery pass. This is for visualizing the BLS sensitivity
# boundary in period/depth space; it intentionally skips LEO and human-queue
# rendering.
set -euo pipefail

cd /pdo/users/tehan/TWIRL

export LD_LIBRARY_PATH="/pdo/app/anaconda/anaconda2-4.4.0/lib:/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH=src
export MPLCONFIGDIR=/pdo/users/tehan/.cache/matplotlib
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

PYTHON="${PYTHON:-/sw/qlp-environment/.venv/bin/python}"
ORBIT_ROOTS=(
  "${ORBIT_119_ROOT:-/pdo/users/tehan/tglc-gpu-production/orbit-119/ffi}"
  "${ORBIT_120_ROOT:-/pdo/users/tehan/tglc-gpu-production/orbit-120/ffi}"
)
N_INJECTIONS="${N_INJECTIONS:-10000}"
RANDOM_STATE="${RANDOM_STATE:-5623}"
SAMPLING_MODE="${SAMPLING_MODE:-period_depth_grid}"
GRID_PERIOD_RANGE_D="${GRID_PERIOD_RANGE_D:-0.12,13.0}"
GRID_RADIUS_RANGE_REARTH="${GRID_RADIUS_RANGE_REARTH:-0.18,18.0}"
GRID_PERIOD_BINS="${GRID_PERIOD_BINS:-50}"
GRID_RADIUS_BINS="${GRID_RADIUS_BINS:-50}"
GRID_DEPTH_BINS="${GRID_DEPTH_BINS:-50}"
GRID_DEPTH_RANGE="${GRID_DEPTH_RANGE:-0.003,0.995}"
GRID_DEPTH_SPACING="${GRID_DEPTH_SPACING:-linear}"
TARGET_TMAG_BIN_EDGES="${TARGET_TMAG_BIN_EDGES:-}"
TARGET_TMAG_BIN_WEIGHTS="${TARGET_TMAG_BIN_WEIGHTS:-}"
TARGET_TMAG_TABLE="${TARGET_TMAG_TABLE:-}"
TARGET_TMAG_COLUMN="${TARGET_TMAG_COLUMN:-tessmag}"
TARGET_H5_LIST="${TARGET_H5_LIST:-}"
BASELINE_SOURCE="${BASELINE_SOURCE:-raw_median}"
INJECTION_DIR="${INJECTION_DIR:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_10k_predetrend_batman_depthgrid_dense_bls_map}"
OUT_DIR="${OUT_DIR:-reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo}"
CHUNK_DIR="${OUT_DIR}/chunk_ids"
RUN_CHUNK_DIR="${OUT_DIR}/chunks"
CHUNK_SIZE="${TWIRL_CHUNK_SIZE:-100}"
WORKERS="${TWIRL_CHUNK_WORKERS:-16}"
CHUNK_JOBS="${TWIRL_CHUNK_JOBS:-1}"
SWEEP_NAME="${SWEEP_NAME:-small_pair_200k}"
SWEEP_SPEC="${SWEEP_NAME}:DET_FLUX_ADP_SML+DET_FLUX_SML:200000"
REBUILD_INJECTIONS="${REBUILD_INJECTIONS:-0}"

log() { echo "[$(date -Is)] [dense-bls-map] $*"; }

mkdir -p "${CHUNK_DIR}" "${RUN_CHUNK_DIR}" logs

log "python=${PYTHON}"
log "n_injections=${N_INJECTIONS} sampling=${SAMPLING_MODE}"
log "period_grid=${GRID_PERIOD_BINS} range=${GRID_PERIOD_RANGE_D}"
log "radius_grid=${GRID_RADIUS_BINS} range=${GRID_RADIUS_RANGE_REARTH}"
log "depth_grid=${GRID_DEPTH_BINS} range=${GRID_DEPTH_RANGE} spacing=${GRID_DEPTH_SPACING}"
if [[ -n "${TARGET_TMAG_BIN_EDGES}" ]]; then
  log "target_tmag_bin_edges=${TARGET_TMAG_BIN_EDGES} weights=${TARGET_TMAG_BIN_WEIGHTS:-equal}"
fi
if [[ -n "${TARGET_TMAG_TABLE}" ]]; then
  log "target_tmag_table=${TARGET_TMAG_TABLE} column=${TARGET_TMAG_COLUMN}"
fi
if [[ -n "${TARGET_H5_LIST}" ]]; then
  log "target_h5_list=${TARGET_H5_LIST}"
fi
log "injection_dir=${INJECTION_DIR}"
log "out_dir=${OUT_DIR}"
log "chunk_size=${CHUNK_SIZE} chunk_jobs=${CHUNK_JOBS} chunk_workers=${WORKERS}"

if [[ "${REBUILD_INJECTIONS}" == "1" || ! -s "${INJECTION_DIR}/injected_lightcurves.h5" ]]; then
  log "building dense raw-flux pre-detrend BATMAN injections"
  INJECTION_ARGS=(
    --orbit-roots "${ORBIT_ROOTS[@]}"
    --out-dir "${INJECTION_DIR}"
    --n-injections "${N_INJECTIONS}"
    --apertures DET_FLUX_ADP_SML DET_FLUX_SML
    --sampling-mode "${SAMPLING_MODE}"
    --grid-period-range-d "${GRID_PERIOD_RANGE_D}"
    --grid-radius-range-rearth "${GRID_RADIUS_RANGE_REARTH}"
    --grid-depth-range "${GRID_DEPTH_RANGE}"
    --grid-period-bins "${GRID_PERIOD_BINS}"
    --grid-radius-bins "${GRID_RADIUS_BINS}"
    --grid-depth-bins "${GRID_DEPTH_BINS}"
    --grid-depth-spacing "${GRID_DEPTH_SPACING}"
    --baseline-source "${BASELINE_SOURCE}"
    --random-state "${RANDOM_STATE}"
    --min-in-transit 2
    --max-attempts-per-injection 100
    --progress-every 100
    --overwrite
  )
  if [[ -n "${TARGET_TMAG_BIN_EDGES}" ]]; then
    INJECTION_ARGS+=(--target-tmag-bin-edges "${TARGET_TMAG_BIN_EDGES}")
  fi
  if [[ -n "${TARGET_TMAG_BIN_WEIGHTS}" ]]; then
    INJECTION_ARGS+=(--target-tmag-bin-weights "${TARGET_TMAG_BIN_WEIGHTS}")
  fi
  if [[ -n "${TARGET_TMAG_TABLE}" ]]; then
    INJECTION_ARGS+=(--target-tmag-table "${TARGET_TMAG_TABLE}" --target-tmag-column "${TARGET_TMAG_COLUMN}")
  fi
  if [[ -n "${TARGET_H5_LIST}" ]]; then
    INJECTION_ARGS+=(--target-h5-list "${TARGET_H5_LIST}")
  fi
  "${PYTHON}" scripts/stage3_injections/make_s56_predetrend_injection_set.py \
    "${INJECTION_ARGS[@]}"
else
  log "reusing injections: ${INJECTION_DIR}/injected_lightcurves.h5"
fi

log "writing restartable BLS chunk id files"
"${PYTHON}" - "${INJECTION_DIR}/injected_lightcurves.h5" "${CHUNK_DIR}" "${CHUNK_SIZE}" <<'PY'
from pathlib import Path
import sys

import h5py

injection_h5 = Path(sys.argv[1])
chunk_dir = Path(sys.argv[2])
chunk_size = int(sys.argv[3])

with h5py.File(injection_h5, "r") as h5:
    keys = sorted(h5["injections"].keys())

for old in chunk_dir.glob("chunk_*.txt"):
    old.unlink()

for start in range(0, len(keys), chunk_size):
    chunk = keys[start:start + chunk_size]
    path = chunk_dir / f"chunk_{start // chunk_size:03d}.txt"
    path.write_text("\n".join(chunk) + "\n", encoding="utf-8")

(chunk_dir / "manifest.txt").write_text(
    f"n_injections={len(keys)}\nchunk_size={chunk_size}\nn_chunks={(len(keys) + chunk_size - 1) // chunk_size}\n",
    encoding="utf-8",
)
print((chunk_dir / "manifest.txt").read_text(), end="")
PY

if [[ "${CHUNK_JOBS}" -le 1 ]]; then
  for id_file in "${CHUNK_DIR}"/chunk_*.txt; do
    chunk_name="$(basename "${id_file}" .txt)"
    chunk_out="${RUN_CHUNK_DIR}/${chunk_name}"
    summary_path="${chunk_out}/${SWEEP_NAME}/recovery_mode_summary/summary.json"
    if [[ -f "${summary_path}" ]]; then
      log "skip ${chunk_name}: ${summary_path} exists"
      continue
    fi
    log "run ${chunk_name} with ${WORKERS} workers"
    "${PYTHON}" scripts/stage5_validation/run_injection_recovery_mode_sweep.py \
      --injection-h5 "${INJECTION_DIR}/injected_lightcurves.h5" \
      --out-dir "${chunk_out}" \
      --n-injections 0 \
      --workers "${WORKERS}" \
      --injection-id-file "${id_file}" \
      --sweep "${SWEEP_SPEC}"
  done
else
  pending_chunks="${OUT_DIR}/pending_chunk_files.txt"
  : > "${pending_chunks}"
  for id_file in "${CHUNK_DIR}"/chunk_*.txt; do
    chunk_name="$(basename "${id_file}" .txt)"
    chunk_out="${RUN_CHUNK_DIR}/${chunk_name}"
    summary_path="${chunk_out}/${SWEEP_NAME}/recovery_mode_summary/summary.json"
    if [[ -f "${summary_path}" ]]; then
      log "skip ${chunk_name}: ${summary_path} exists"
      continue
    fi
    printf '%s\n' "${id_file}" >> "${pending_chunks}"
  done

  n_pending="$(wc -l < "${pending_chunks}" | tr -d ' ')"
  log "run ${n_pending} pending chunks with ${CHUNK_JOBS} concurrent jobs and ${WORKERS} workers each"
  if [[ "${n_pending}" -gt 0 ]]; then
    export PYTHON INJECTION_DIR RUN_CHUNK_DIR SWEEP_NAME SWEEP_SPEC WORKERS
    xargs -n 1 -P "${CHUNK_JOBS}" bash -c '
      set -euo pipefail
      id_file="$1"
      chunk_name="$(basename "${id_file}" .txt)"
      chunk_out="${RUN_CHUNK_DIR}/${chunk_name}"
      summary_path="${chunk_out}/${SWEEP_NAME}/recovery_mode_summary/summary.json"
      if [[ -f "${summary_path}" ]]; then
        echo "[$(date -Is)] [dense-bls-map] skip ${chunk_name}: ${summary_path} exists"
        exit 0
      fi
      echo "[$(date -Is)] [dense-bls-map] run ${chunk_name} with ${WORKERS} workers"
      "${PYTHON}" scripts/stage5_validation/run_injection_recovery_mode_sweep.py \
        --injection-h5 "${INJECTION_DIR}/injected_lightcurves.h5" \
        --out-dir "${chunk_out}" \
        --n-injections 0 \
        --workers "${WORKERS}" \
        --injection-id-file "${id_file}" \
        --sweep "${SWEEP_SPEC}"
    ' _ < "${pending_chunks}"
  fi
fi

log "merging BLS chunks"
"${PYTHON}" scripts/stage5_validation/merge_recovery_sweep_chunks.py \
  --chunk-root "${RUN_CHUNK_DIR}" \
  --out-dir "${OUT_DIR}" \
  --sweep-name "${SWEEP_NAME}" \
  --apertures DET_FLUX_ADP_SML DET_FLUX_SML \
  --n-periods 200000 \
  --summary-aperture DET_FLUX_ADP_SML \
  --expect-n "${N_INJECTIONS}"

log "complete"
log "recovery_csv=${OUT_DIR}/${SWEEP_NAME}/injection_bls_recoveries.csv"
