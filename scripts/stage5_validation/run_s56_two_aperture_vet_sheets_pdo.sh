#!/usr/bin/env bash
# Render TWIRL-native two-aperture vet sheets for a review queue on PDO.
set -euo pipefail

cd "${TWIRL_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"

PYTHON_BIN="${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}"
export PYTHONPATH="${PYTHONPATH:-src}"
export LD_LIBRARY_PATH="/pdo/app/anaconda/anaconda2-4.4.0/lib:/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

QUEUE_DIR="${QUEUE_DIR:-reports/stage5_validation/s56_mixed_teacher_queue_pdo}"
QUEUE_CSV="${QUEUE_CSV:-${QUEUE_DIR}/review_queue_1k.csv}"
OUT_DIR="${OUT_DIR:-${QUEUE_DIR}/twirl_vet_sheets}"
METRICS_CSV="${METRICS_CSV:-${QUEUE_DIR}/twirl_vet_metrics.csv}"
SUMMARY_JSON="${SUMMARY_JSON:-${QUEUE_DIR}/twirl_two_aperture_vet_summary.json}"
HLSP_ROOT="${HLSP_ROOT:-/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl_fs_v2_adp015q_compare}"
LC_EXPORT_H5="${LC_EXPORT_H5:-data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp015_lc_export_pdo.h5}"
BRANCH_NAME="${BRANCH_NAME:-twirl_fs_v2_adp015q}"
APERTURE_TAG="${APERTURE_TAG:-ADP015}"
LIMIT="${LIMIT:-0}"
WORKERS="${WORKERS:-8}"
N_PERIODS="${N_PERIODS:-20000}"
N_PEAKS="${N_PEAKS:-10}"

args=(
  --queue-csv "${QUEUE_CSV}"
  --hlsp-root "${HLSP_ROOT}"
  --lc-export-h5 "${LC_EXPORT_H5}"
  --out-dir "${OUT_DIR}"
  --metrics-csv "${METRICS_CSV}"
  --summary-json "${SUMMARY_JSON}"
  --branch-name "${BRANCH_NAME}"
  --aperture-tag "${APERTURE_TAG}"
  --limit "${LIMIT}"
  --workers "${WORKERS}"
  --n-periods "${N_PERIODS}"
  --n-peaks "${N_PEAKS}"
)
if [[ "${OVERWRITE:-0}" == "1" ]]; then
  args+=(--overwrite)
fi

echo "[two-aperture-pdo] queue=${QUEUE_CSV}"
echo "[two-aperture-pdo] out_dir=${OUT_DIR}"
echo "[two-aperture-pdo] summary_json=${SUMMARY_JSON}"
echo "[two-aperture-pdo] lc_export_h5=${LC_EXPORT_H5}"
echo "[two-aperture-pdo] branch=${BRANCH_NAME} aperture_tag=${APERTURE_TAG} n_periods=${N_PERIODS} workers=${WORKERS}"
"${PYTHON_BIN}" scripts/stage5_validation/render_two_aperture_vet_sheets.py "${args[@]}"
