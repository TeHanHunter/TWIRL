#!/usr/bin/env bash
# Sync the PDO dense S56 injection/BLS recovery output and render the local
# parameter-space sensitivity plots.
set -euo pipefail

REMOTE_HOST="${REMOTE_HOST:-pdogpu6}"
REMOTE_ROOT="${REMOTE_ROOT:-/pdo/users/tehan/TWIRL}"
REMOTE_DIR="${REMOTE_DIR:-${REMOTE_ROOT}/reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo}"
LOCAL_DIR="${LOCAL_DIR:-reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo}"
SWEEP_NAME="${SWEEP_NAME:-small_pair_200k}"
PERIOD_BINS="${PERIOD_BINS:-36}"
DEPTH_BINS="${DEPTH_BINS:-40}"
PYTHON_BIN="${PYTHON_BIN:-.venv/bin/python}"
PLOT_PARAMETER_SPACE="${PLOT_PARAMETER_SPACE:-1}"
PLOT_DURATION_AWARE="${PLOT_DURATION_AWARE:-1}"
LEO_QUEUE_CSV="${LEO_QUEUE_CSV:-reports/stage5_validation/s56_10k_predetrend_small_pair_200k_review_queue_pdo/review_queue.csv}"
LEO_METRICS_CSV="${LEO_METRICS_CSV:-reports/stage5_validation/s56_10k_predetrend_small_pair_200k_review_queue_pdo/leo_metrics.csv}"

if [[ ! -x "${PYTHON_BIN}" ]]; then
  echo "[sync-dense-map] missing Python executable: ${PYTHON_BIN}" >&2
  exit 2
fi

mkdir -p "${LOCAL_DIR}"
echo "[sync-dense-map] syncing ${REMOTE_HOST}:${REMOTE_DIR}/ -> ${LOCAL_DIR}/"
rsync -avhu --progress \
  --exclude '/chunks/' \
  --exclude '/chunk_ids/' \
  --exclude '/pending_chunk_files.txt' \
  "${REMOTE_HOST}:${REMOTE_DIR}/" "${LOCAL_DIR}/"

input_csv="${LOCAL_DIR}/${SWEEP_NAME}/injection_bls_recoveries.csv"
if [[ ! -s "${input_csv}" ]]; then
  echo "[sync-dense-map] missing recovery CSV after sync: ${input_csv}" >&2
  exit 3
fi

if [[ "${PLOT_PARAMETER_SPACE}" == "1" ]]; then
  echo "[sync-dense-map] plotting parameter-space maps from ${input_csv}"
  MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/twirl_mpl}" PYTHONPATH=src "${PYTHON_BIN}" \
    scripts/stage5_validation/plot_s56_injection_recovery_parameter_space.py \
    --input-csv "${input_csv}" \
    --out-dir "${LOCAL_DIR}/parameter_space" \
    --period-bins "${PERIOD_BINS}" \
    --depth-bins "${DEPTH_BINS}"
fi

if [[ "${PLOT_DURATION_AWARE}" == "1" ]]; then
  echo "[sync-dense-map] plotting duration-aware period-radius maps from ${input_csv}"
  MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/twirl_mpl}" PYTHONPATH=src "${PYTHON_BIN}" \
    scripts/stage5_validation/plot_s56_duration_aware_recovery.py \
    --bls-csv "${input_csv}" \
    --leo-queue-csv "${LEO_QUEUE_CSV}" \
    --leo-metrics-csv "${LEO_METRICS_CSV}" \
    --out-dir "${LOCAL_DIR}/duration_aware"
fi
