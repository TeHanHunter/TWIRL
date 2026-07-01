#!/usr/bin/env bash
# Serve the all-host ranker-selected real S56 queue for human triage.
#
# This is the first real-candidate visual review queue after the coverage-first
# all-host injection run, robust BLS peak labeling, and injected-truth peak
# ranker. It should be launched only after WD-tuned LEO sheets have rendered
# and the queue verifier has passed.
set -euo pipefail

REPO_ROOT="${TWIRL_REPO_ROOT:-/pdo/users/tehan/TWIRL}"
cd "${REPO_ROOT}"

PYTHON_BIN="${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}"
export LD_LIBRARY_PATH="/pdo/app/anaconda/anaconda2-4.4.0/lib:/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${PYTHONPATH:-src}"

QUEUE_DIR="${QUEUE_DIR:-reports/stage5_validation/s56_allhost_ranker_selected_real_leo_queue_pdo}"
CANDIDATES="${CANDIDATES:-${QUEUE_DIR}/review_queue.csv}"
LABELS_OUT="${LABELS_OUT:-${QUEUE_DIR}/human_labels_vetted.csv}"
LEO_REPORT_ROOT="${LEO_REPORT_ROOT:-${QUEUE_DIR}/vet_reports}"
VERIFY_JSON="${VERIFY_JSON:-${QUEUE_DIR}/verification.json}"
HLSP_ROOT="${HLSP_ROOT:-/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl_fs_v2_compare}"
HOST="${HOST:-127.0.0.1}"
PORT="${PORT:-5007}"
TUNNEL_HOST="${TUNNEL_HOST:-pdogpu1}"
APERTURE="${APERTURE:-DET_FLUX_ADP_SML}"
SHUFFLE_SEED="${SHUFFLE_SEED:-5607}"

if [[ ! -s "${VERIFY_JSON}" ]]; then
  echo "[allhost-real-vet-app] missing verifier: ${VERIFY_JSON}" >&2
  echo "[allhost-real-vet-app] run run_s56_allhost_ranker_selected_real_leo_pdo.sh first" >&2
  exit 2
fi
if ! "${PYTHON_BIN}" - "${VERIFY_JSON}" <<'PY'
import json
import sys
from pathlib import Path

payload = json.loads(Path(sys.argv[1]).read_text())
raise SystemExit(0 if payload.get("passed") else 1)
PY
then
  echo "[allhost-real-vet-app] verifier did not pass: ${VERIFY_JSON}" >&2
  exit 2
fi

args=(
  --host "${HOST}"
  --port "${PORT}"
  --candidates "${CANDIDATES}"
  --labels-out "${LABELS_OUT}"
  --hlsp-root "${HLSP_ROOT}"
  --leo-report-root "${LEO_REPORT_ROOT}"
  --aperture "${APERTURE}"
  --tunnel-host "${TUNNEL_HOST}"
  --shuffle-order
  --shuffle-seed "${SHUFFLE_SEED}"
)

if [[ -n "${SMOKE_PLOT:-}" ]]; then
  args+=(--smoke-plot "${SMOKE_PLOT}" --smoke-index "${SMOKE_INDEX:-0}")
fi

echo "[allhost-real-vet-app] queue=${CANDIDATES}"
echo "[allhost-real-vet-app] labels=${LABELS_OUT}"
echo "[allhost-real-vet-app] leo_reports=${LEO_REPORT_ROOT}"
echo "[allhost-real-vet-app] aperture=${APERTURE}"
echo "[allhost-real-vet-app] url=http://${HOST}:${PORT}/"
"${PYTHON_BIN}" scripts/stage5_validation/run_lightcurve_vetting_app.py "${args[@]}"
