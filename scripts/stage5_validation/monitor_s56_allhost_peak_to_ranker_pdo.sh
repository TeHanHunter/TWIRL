#!/usr/bin/env bash
# Wait for the all-host injected BLS peak table, then run the ranker handoff.
#
# Default behavior is deliberately conservative: build and verify the
# ranker-selected real queue with LEO rendering skipped. Full LEO report
# rendering can then be launched from the verified selected ephemerides.
set -euo pipefail

REPO_ROOT="${TWIRL_REPO_ROOT:-/pdo/users/tehan/TWIRL}"
cd "${REPO_ROOT}"

PYTHON_BIN="${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}"
export LD_LIBRARY_PATH="/sw/python-versions/python-3.11.9/lib:/pdo/app/anaconda/anaconda2-4.4.0/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${PYTHONPATH:-src}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

ROOT="${TWIRL_ALLHOST_ROOT:-reports/stage5_validation/s56_allhost_predetrend_batman_periodradius_grid_sharded_pdo}"
PEAK_TABLE="${TWIRL_PEAK_TABLE:-${ROOT}/peak_training/s56_allhost_injection_bls_peaks.csv}"
VERIFY_JSON="${TWIRL_VERIFY_JSON:-${PEAK_TABLE%.*}_verification.json}"
AUDIT_OUT_DIR="${TWIRL_AUDIT_OUT_DIR:-reports/stage5_validation/s56_allhost_peak_handoff_audit_pdo}"
HANDOFF_LOG_DIR="${TWIRL_HANDOFF_LOG_DIR:-${ROOT}/ranker_handoff}"
REVIEW_VERIFY_JSON="${TWIRL_REVIEW_VERIFY_JSON:-reports/stage5_validation/s56_allhost_ranker_selected_real_review_queue_pdo/verification.json}"
WAIT_SECONDS="${WAIT_SECONDS:-600}"
SKIP_LEO="${TWIRL_SKIP_LEO:-1}"

log() {
  printf '[%s] [allhost-ranker-monitor] %s\n' "$(date -Is)" "$*" >&2
}

json_passed() {
  local path="$1"
  "${PYTHON_BIN}" - "$path" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
if not path.exists() or path.stat().st_size == 0:
    raise SystemExit(1)
payload = json.loads(path.read_text())
if "passed" in payload:
    raise SystemExit(0 if payload["passed"] else 1)
failures = payload.get("failures")
if isinstance(failures, list):
    raise SystemExit(0 if not failures else 1)
raise SystemExit(0)
PY
}

mkdir -p "${HANDOFF_LOG_DIR}" "${AUDIT_OUT_DIR}"

if json_passed "${REVIEW_VERIFY_JSON}" 2>/dev/null; then
  log "review queue already verified; nothing to do"
  log "review_verify_json=${REVIEW_VERIFY_JSON}"
  exit 0
fi

while true; do
  "${PYTHON_BIN}" scripts/stage5_validation/audit_s56_peak_handoff_status.py \
    --layout pdo_allhost \
    --root "${ROOT}" \
    --out-dir "${AUDIT_OUT_DIR}" >/dev/null || true

  if [[ -s "${PEAK_TABLE}" ]] && json_passed "${VERIFY_JSON}" 2>/dev/null; then
    log "verified all-host injected peak table is ready"
    break
  fi
  log "waiting for verified all-host peak table"
  log "peak_table=${PEAK_TABLE}"
  log "verify_json=${VERIFY_JSON}"
  sleep "${WAIT_SECONDS}"
done

log "starting all-host ranker/review handoff; skip_leo=${SKIP_LEO}"
TWIRL_SKIP_LEO="${SKIP_LEO}" \
  bash scripts/stage5_validation/run_s56_allhost_peak_ranker_review_pdo.sh \
  >"${HANDOFF_LOG_DIR}/run.log" 2>&1

log "handoff complete"
log "run_log=${HANDOFF_LOG_DIR}/run.log"
log "review_verify_json=${REVIEW_VERIFY_JSON}"
