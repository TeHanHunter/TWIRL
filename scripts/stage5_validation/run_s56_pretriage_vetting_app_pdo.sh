#!/usr/bin/env bash
set -euo pipefail

cd /pdo/users/tehan/TWIRL

export LD_LIBRARY_PATH="/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="src"

PYTHON_BIN="${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}"
HOST="${HOST:-127.0.0.1}"
PORT="${PORT:-5000}"
HLSP_ROOT="${HLSP_ROOT:-/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl_fs_v2_compare}"
REVIEW_DIR="${REVIEW_DIR:-reports/stage5_validation/s56_pretriage_review_queue_pdo}"
CANDIDATES="${CANDIDATES:-${REVIEW_DIR}/review_queue.csv}"
LABELS_OUT="${LABELS_OUT:-${REVIEW_DIR}/human_labels_vetted.csv}"
LEO_REPORT_ROOT="${LEO_REPORT_ROOT:-${REVIEW_DIR}/vet_reports}"

echo "[vet-app] candidates=${CANDIDATES}"
echo "[vet-app] labels_out=${LABELS_OUT}"
echo "[vet-app] leo_reports=${LEO_REPORT_ROOT}"

"${PYTHON_BIN}" scripts/stage5_validation/run_lightcurve_vetting_app.py \
  --host "${HOST}" \
  --port "${PORT}" \
  --candidates "${CANDIDATES}" \
  --labels-out "${LABELS_OUT}" \
  --hlsp-root "${HLSP_ROOT}" \
  --leo-report-root "${LEO_REPORT_ROOT}" \
  --shuffle-order \
  --shuffle-seed 56
