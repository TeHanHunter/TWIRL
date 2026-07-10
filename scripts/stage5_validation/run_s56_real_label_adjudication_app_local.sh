#!/usr/bin/env bash
# Serve the fixed S56 real-label adjudication queue in the standard vetter.
set -euo pipefail

cd "${TWIRL_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"

PYTHON_BIN="${PYTHON_BIN:-.venv/bin/python}"
export PYTHONPATH="${PYTHONPATH:-src}"

QUEUE_DIR="${QUEUE_DIR:-reports/stage5_validation/s56_label_adjudication_real343}"
CANDIDATES="${CANDIDATES:-${QUEUE_DIR}/review_queue_real343.csv}"
LABELS_OUT="${LABELS_OUT:-${QUEUE_DIR}/human_labels_adjudicated.csv}"
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

echo "[adjudication-vet-local] queue=${CANDIDATES}"
echo "[adjudication-vet-local] labels=${LABELS_OUT}"
echo "[adjudication-vet-local] sheets=${TWIRL_VET_ROOT}"
echo "[adjudication-vet-local] url=http://${HOST}:${PORT}/"
"${PYTHON_BIN}" scripts/stage5_validation/run_lightcurve_vetting_app.py "${args[@]}"
