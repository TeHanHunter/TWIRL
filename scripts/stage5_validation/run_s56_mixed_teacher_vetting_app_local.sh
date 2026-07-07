#!/usr/bin/env bash
# Serve the synced S56 mixed-teacher 1k first-pass queue locally.
set -euo pipefail

cd "${TWIRL_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"

PYTHON_BIN="${PYTHON_BIN:-.venv/bin/python}"
export PYTHONPATH="${PYTHONPATH:-src}"

QUEUE_DIR="${QUEUE_DIR:-reports/stage5_validation/s56_mixed_teacher_queue_pdo}"
CANDIDATES="${CANDIDATES:-${QUEUE_DIR}/review_queue_1k.csv}"
LABELS_OUT="${LABELS_OUT:-${QUEUE_DIR}/human_labels_vetted.csv}"
LEO_REPORT_ROOT="${LEO_REPORT_ROOT:-${QUEUE_DIR}/vet_reports}"
DEFAULT_TWIRL_VET_ROOT="${QUEUE_DIR}/twirl_vet_sheets"
if [[ -d "${QUEUE_DIR}/twirl_vet_sheets_adp015q_orcd" ]]; then
  DEFAULT_TWIRL_VET_ROOT="${QUEUE_DIR}/twirl_vet_sheets_adp015q_orcd"
fi
TWIRL_VET_ROOT="${TWIRL_VET_ROOT:-${DEFAULT_TWIRL_VET_ROOT}}"
HLSP_ROOT="${HLSP_ROOT:-data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_adp015q_compare}"
VERIFY_JSON="${VERIFY_JSON:-${QUEUE_DIR}/verification.json}"
HOST="${HOST:-127.0.0.1}"
PORT="${PORT:-5006}"
APERTURE="${APERTURE:-DET_FLUX_ADP015_SML}"

if [[ ! -s "${VERIFY_JSON}" ]]; then
  echo "[mixed-teacher-vet-local] missing verifier: ${VERIFY_JSON}" >&2
  echo "[mixed-teacher-vet-local] run sync_s56_mixed_teacher_queue_pdo.sh first" >&2
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
  echo "[mixed-teacher-vet-local] verifier did not pass: ${VERIFY_JSON}" >&2
  exit 2
fi

args=(
  --host "${HOST}"
  --port "${PORT}"
  --candidates "${CANDIDATES}"
  --labels-out "${LABELS_OUT}"
  --hlsp-root "${HLSP_ROOT}"
  --leo-report-root "${LEO_REPORT_ROOT}"
  --twirl-vet-root "${TWIRL_VET_ROOT}"
  --aperture "${APERTURE}"
)

if [[ "${SHUFFLE_ORDER:-0}" == "1" ]]; then
  args+=(--shuffle-order --shuffle-seed "${SHUFFLE_SEED:-5608}")
fi
if [[ -n "${SMOKE_PLOT:-}" ]]; then
  args+=(--smoke-plot "${SMOKE_PLOT}" --smoke-index "${SMOKE_INDEX:-0}")
fi

echo "[mixed-teacher-vet-local] queue=${CANDIDATES}"
echo "[mixed-teacher-vet-local] labels=${LABELS_OUT}"
echo "[mixed-teacher-vet-local] leo_reports=${LEO_REPORT_ROOT}"
echo "[mixed-teacher-vet-local] twirl_vet_sheets=${TWIRL_VET_ROOT}"
echo "[mixed-teacher-vet-local] aperture=${APERTURE}"
echo "[mixed-teacher-vet-local] url=http://${HOST}:${PORT}/"
"${PYTHON_BIN}" scripts/stage5_validation/run_lightcurve_vetting_app.py "${args[@]}"
