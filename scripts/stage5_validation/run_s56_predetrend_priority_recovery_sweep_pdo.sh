#!/usr/bin/env bash
# Run a targeted recovery-mode sweep on the pre-detrend priority debug set.
#
# This is the cheap first pass before rerunning the full 1k or scaling to 10k.
# It uses the priority_injection_ids.txt generated locally from the current
# BLS-failure and aperture-survival diagnostics.
set -euo pipefail

cd "${TWIRL_ROOT:-/pdo/users/tehan/TWIRL}"

PYTHON_BIN="${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}"
export LD_LIBRARY_PATH="/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="src"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/pdo/users/tehan/.cache/matplotlib}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

mkdir -p "${MPLCONFIGDIR}"

INJECTION_H5="${INJECTION_H5:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_1k_predetrend_batman_depthgrid_adp_compare/injected_lightcurves.h5}"
SURVIVAL_CSV="${SURVIVAL_CSV:-reports/stage5_validation/s56_1k_predetrend_batman_adp_only_review_queue_pdo/predet_signal_survival_snr_diagnostics.csv}"
INJECTION_ID_FILE="${INJECTION_ID_FILE:-reports/stage5_validation/s56_1k_predetrend_debug_priority_queues/priority_injection_ids.txt}"
PRIORITY_QUEUE_CSV="${PRIORITY_QUEUE_CSV:-reports/stage5_validation/s56_1k_predetrend_debug_priority_queues/priority_debug_queue.csv}"
OUT_DIR="${OUT_DIR:-reports/stage5_validation/s56_1k_predetrend_priority_recovery_sweep_pdo}"
N_INJECTIONS="${N_INJECTIONS:-0}"
WORKERS="${WORKERS:-8}"
REUSE="${REUSE:-0}"
EXPECTED_ID_COUNT="${EXPECTED_ID_COUNT:-94}"

sweep_args=(
  --sweep "adp_sml_priority_5k:DET_FLUX_ADP_SML:5000"
  --sweep "det_sml_priority_5k:DET_FLUX_SML:5000"
  --sweep "adp_priority_5k:DET_FLUX_ADP:5000"
  --sweep "adp_lag_priority_5k:DET_FLUX_ADP_LAG:5000"
  --sweep "small_pair_priority_5k:DET_FLUX_ADP_SML+DET_FLUX_SML:5000"
  --sweep "small_pair_priority_200k:DET_FLUX_ADP_SML+DET_FLUX_SML:200000"
)

cmd=(
  "${PYTHON_BIN}" scripts/stage5_validation/run_injection_recovery_mode_sweep.py
  --injection-h5 "${INJECTION_H5}"
  --survival-csv "${SURVIVAL_CSV}"
  --injection-id-file "${INJECTION_ID_FILE}"
  --out-dir "${OUT_DIR}"
  --n-injections "${N_INJECTIONS}"
  --workers "${WORKERS}"
  "${sweep_args[@]}"
)
if [[ "${REUSE}" == "1" ]]; then
  cmd+=(--reuse)
fi

echo "[pdo-priority-sweep] start $(date -Is)"
echo "[pdo-priority-sweep] python=${PYTHON_BIN}"
echo "[pdo-priority-sweep] injection_h5=${INJECTION_H5}"
echo "[pdo-priority-sweep] survival_csv=${SURVIVAL_CSV}"
echo "[pdo-priority-sweep] injection_id_file=${INJECTION_ID_FILE}"
echo "[pdo-priority-sweep] priority_queue_csv=${PRIORITY_QUEUE_CSV}"
echo "[pdo-priority-sweep] out_dir=${OUT_DIR}"
"${PYTHON_BIN}" scripts/stage5_validation/check_s56_predetrend_priority_sweep_inputs.py \
  --injection-h5 "${INJECTION_H5}" \
  --survival-csv "${SURVIVAL_CSV}" \
  --injection-id-file "${INJECTION_ID_FILE}" \
  --priority-queue-csv "${PRIORITY_QUEUE_CSV}" \
  --expected-id-count "${EXPECTED_ID_COUNT}" \
  --required-aperture DET_FLUX_ADP_SML \
  --required-aperture DET_FLUX_SML \
  --required-aperture DET_FLUX_ADP \
  --required-aperture DET_FLUX_ADP_LAG \
  --out-json "${OUT_DIR}/preflight_summary.json"
"${cmd[@]}"
"${PYTHON_BIN}" scripts/stage5_validation/summarize_priority_recovery_sweep.py \
  --sweep-dir "${OUT_DIR}" \
  --priority-queue-csv "${PRIORITY_QUEUE_CSV}" \
  --out-dir "${OUT_DIR}/priority_sweep_summary" \
  --baseline-sweep "adp_priority_5k"
echo "[pdo-priority-sweep] complete $(date -Is)"
echo "[pdo-priority-sweep] overview: ${OUT_DIR}/sweep_overview.csv"
echo "[pdo-priority-sweep] priority summary: ${OUT_DIR}/priority_sweep_summary/summary.json"
