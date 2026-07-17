#!/usr/bin/env bash
# Run a non-overlapping range of dense-map BLS chunks on PDO.
#
# The main dense-map runner creates the injection HDF5 and chunk id files. This
# helper lets us fan out independent chunk ranges in additional tmux sessions;
# the final merge is still handled by run_s56_predetrend_dense_bls_map_pdo.sh
# or by merge_recovery_sweep_chunks.py.
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
INJECTION_DIR="${INJECTION_DIR:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_10k_predetrend_batman_depthgrid_dense_bls_map}"
OUT_DIR="${OUT_DIR:-reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo}"
CHUNK_DIR="${OUT_DIR}/chunk_ids"
RUN_CHUNK_DIR="${OUT_DIR}/chunks"
WORKERS="${TWIRL_CHUNK_WORKERS:-16}"
SWEEP_NAME="${SWEEP_NAME:-small_pair_200k}"
SWEEP_SPEC="${SWEEP_NAME}:DET_FLUX_ADP_SML+DET_FLUX_SML:200000"
START_CHUNK="${START_CHUNK:?set START_CHUNK, e.g. 20}"
END_CHUNK="${END_CHUNK:?set END_CHUNK, e.g. 39}"

log() { echo "[$(date -Is)] [dense-bls-range] $*"; }

if [[ ! -s "${INJECTION_DIR}/injected_lightcurves.h5" ]]; then
  echo "[dense-bls-range] missing injection HDF5: ${INJECTION_DIR}/injected_lightcurves.h5" >&2
  exit 2
fi

log "range=${START_CHUNK}-${END_CHUNK} workers=${WORKERS}"
for idx in $(seq "${START_CHUNK}" "${END_CHUNK}"); do
  chunk_name="$(printf 'chunk_%03d' "${idx}")"
  id_file="${CHUNK_DIR}/${chunk_name}.txt"
  chunk_out="${RUN_CHUNK_DIR}/${chunk_name}"
  summary_path="${chunk_out}/${SWEEP_NAME}/recovery_mode_summary/summary.json"
  if [[ ! -s "${id_file}" ]]; then
    log "skip ${chunk_name}: missing ${id_file}"
    continue
  fi
  if [[ -f "${summary_path}" ]]; then
    log "skip ${chunk_name}: ${summary_path} exists"
    continue
  fi
  log "run ${chunk_name}"
  "${PYTHON}" scripts/stage5_validation/run_injection_recovery_mode_sweep.py \
    --injection-h5 "${INJECTION_DIR}/injected_lightcurves.h5" \
    --out-dir "${chunk_out}" \
    --n-injections 0 \
    --workers "${WORKERS}" \
    --injection-id-file "${id_file}" \
    --sweep "${SWEEP_SPEC}"
done
log "range complete"
