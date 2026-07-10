#!/usr/bin/env bash
# Serve only discordant blind repeats for a third-pass final decision.
set -euo pipefail

cd "${TWIRL_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"

PYTHON_BIN="${PYTHON_BIN:-.venv/bin/python}"
export PYTHONPATH="${PYTHONPATH:-src}"
ROOT="${ROOT:-reports/stage5_validation/s56_label_adjudication_real343}"
QUEUE_DIR="${QUEUE_DIR:-${ROOT}/repeat_resolution}"
HOST="${HOST:-127.0.0.1}"
PORT="${PORT:-5004}"

args=(
  --host "${HOST}"
  --port "${PORT}"
  --candidates "${QUEUE_DIR}/review_queue_repeat_resolution.csv"
  --labels-out "${QUEUE_DIR}/human_labels_repeat_resolution.csv"
  --hlsp-root data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare
  --twirl-vet-root "${ROOT}/twirl_vet_sheets"
  --aperture DET_FLUX_ADP_SML
  --labeler tehan
)

echo "[repeat-resolution] queue=${QUEUE_DIR}/review_queue_repeat_resolution.csv"
echo "[repeat-resolution] labels=${QUEUE_DIR}/human_labels_repeat_resolution.csv"
echo "[repeat-resolution] url=http://${HOST}:${PORT}/"
"${PYTHON_BIN}" scripts/stage5_validation/run_lightcurve_vetting_app.py "${args[@]}"
