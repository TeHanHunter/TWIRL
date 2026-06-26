#!/usr/bin/env bash
# Serve the active S56 1k LEO-only injection queue for human triage.
#
# This is the current quick visual-check queue for the raw-flux pre-detrend,
# small-aperture, dense-BLS path. It is injection-only by design; labels are
# written beside the queue and can be joined to hidden truth metadata later.
set -euo pipefail

cd "${TWIRL_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"

PYTHON_BIN="${PYTHON_BIN:-.venv/bin/python}"
export PYTHONPATH="${PYTHONPATH:-src}"

QUEUE_DIR="${QUEUE_DIR:-reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_pdo}"
CANDIDATES="${CANDIDATES:-${QUEUE_DIR}/review_queue.csv}"
LABELS_OUT="${LABELS_OUT:-${QUEUE_DIR}/human_labels_vetted.csv}"
LEO_REPORT_ROOT="${LEO_REPORT_ROOT:-${QUEUE_DIR}/vet_reports}"
HLSP_ROOT="${HLSP_ROOT:-data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare}"
HOST="${HOST:-127.0.0.1}"
PORT="${PORT:-5006}"
APERTURE="${APERTURE:-DET_FLUX_ADP_SML}"
SHUFFLE_SEED="${SHUFFLE_SEED:-5606}"

args=(
  --host "${HOST}"
  --port "${PORT}"
  --candidates "${CANDIDATES}"
  --labels-out "${LABELS_OUT}"
  --hlsp-root "${HLSP_ROOT}"
  --leo-report-root "${LEO_REPORT_ROOT}"
  --aperture "${APERTURE}"
  --shuffle-order
  --shuffle-seed "${SHUFFLE_SEED}"
)

if [[ -n "${SMOKE_PLOT:-}" ]]; then
  args+=(--smoke-plot "${SMOKE_PLOT}" --smoke-index "${SMOKE_INDEX:-0}")
fi

echo "[small-pair-vet-app] queue=${CANDIDATES}"
echo "[small-pair-vet-app] labels=${LABELS_OUT}"
echo "[small-pair-vet-app] leo_reports=${LEO_REPORT_ROOT}"
echo "[small-pair-vet-app] aperture=${APERTURE}"
echo "[small-pair-vet-app] url=http://${HOST}:${PORT}/"
"${PYTHON_BIN}" scripts/stage5_validation/run_lightcurve_vetting_app.py "${args[@]}"
