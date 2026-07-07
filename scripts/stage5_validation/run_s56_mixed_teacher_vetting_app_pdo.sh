#!/usr/bin/env bash
# Serve the S56 mixed-teacher 1k first-pass queue on PDO.
set -euo pipefail

REPO_ROOT="${TWIRL_REPO_ROOT:-/pdo/users/tehan/TWIRL}"
cd "${REPO_ROOT}"

PYTHON_BIN="${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}"
export LD_LIBRARY_PATH="/pdo/app/anaconda/anaconda2-4.4.0/lib:/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
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
VERIFY_JSON="${VERIFY_JSON:-${QUEUE_DIR}/verification.json}"
HLSP_ROOT="${HLSP_ROOT:-/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl_fs_v2_adp015q_compare}"
HOST="${HOST:-127.0.0.1}"
PORT="${PORT:-5008}"
TUNNEL_HOST="${TUNNEL_HOST:-pdogpu1}"
APERTURE="${APERTURE:-DET_FLUX_ADP015_SML}"

if [[ ! -s "${VERIFY_JSON}" ]]; then
  echo "[mixed-teacher-vet-pdo] missing verifier: ${VERIFY_JSON}" >&2
  echo "[mixed-teacher-vet-pdo] run run_s56_mixed_teacher_queue_pdo.sh first" >&2
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
  echo "[mixed-teacher-vet-pdo] verifier did not pass: ${VERIFY_JSON}" >&2
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
  --tunnel-host "${TUNNEL_HOST}"
)

if [[ "${SHUFFLE_ORDER:-0}" == "1" ]]; then
  args+=(--shuffle-order --shuffle-seed "${SHUFFLE_SEED:-5608}")
fi
if [[ -n "${SMOKE_PLOT:-}" ]]; then
  args+=(--smoke-plot "${SMOKE_PLOT}" --smoke-index "${SMOKE_INDEX:-0}")
fi

echo "[mixed-teacher-vet-pdo] queue=${CANDIDATES}"
echo "[mixed-teacher-vet-pdo] labels=${LABELS_OUT}"
echo "[mixed-teacher-vet-pdo] leo_reports=${LEO_REPORT_ROOT}"
echo "[mixed-teacher-vet-pdo] twirl_vet_sheets=${TWIRL_VET_ROOT}"
echo "[mixed-teacher-vet-pdo] aperture=${APERTURE}"
echo "[mixed-teacher-vet-pdo] url=http://${HOST}:${PORT}/"
"${PYTHON_BIN}" scripts/stage5_validation/run_lightcurve_vetting_app.py "${args[@]}"
