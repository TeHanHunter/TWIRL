#!/usr/bin/env bash
# Serve the S56 EB/PCEB active-learning queue in the standard TWIRL vetter.
set -euo pipefail

cd "${TWIRL_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"

PYTHON_BIN="${PYTHON_BIN:-.venv/bin/python}"
export PYTHONPATH="${PYTHONPATH:-src}"

N_REVIEW="${N_REVIEW:-100}"
QUEUE_DIR="${QUEUE_DIR:-reports/stage5_validation/s56_eb_miner_adp_only/review_queue_eb_priority_pilot100}"
CANDIDATES="${CANDIDATES:-${QUEUE_DIR}/review_queue_eb_priority_${N_REVIEW}.csv}"
LABELS_OUT="${LABELS_OUT:-${QUEUE_DIR}/human_labels_vetted.csv}"
TWIRL_VET_ROOT="${TWIRL_VET_ROOT:-${QUEUE_DIR}/twirl_vet_sheets}"
HLSP_ROOT="${HLSP_ROOT:-data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare}"
HOST="${HOST:-127.0.0.1}"
PORT="${PORT:-5004}"
APERTURE="${APERTURE:-DET_FLUX_ADP_SML}"

args=(
  --host "${HOST}"
  --port "${PORT}"
  --candidates "${CANDIDATES}"
  --labels-out "${LABELS_OUT}"
  --hlsp-root "${HLSP_ROOT}"
  --twirl-vet-root "${TWIRL_VET_ROOT}"
  --aperture "${APERTURE}"
  --labeler tehan
)

echo "[eb-priority-vet-local] queue=${CANDIDATES}"
echo "[eb-priority-vet-local] labels=${LABELS_OUT}"
echo "[eb-priority-vet-local] twirl_vet_sheets=${TWIRL_VET_ROOT}"
echo "[eb-priority-vet-local] url=http://${HOST}:${PORT}/"
"${PYTHON_BIN}" scripts/stage5_validation/run_lightcurve_vetting_app.py "${args[@]}"
