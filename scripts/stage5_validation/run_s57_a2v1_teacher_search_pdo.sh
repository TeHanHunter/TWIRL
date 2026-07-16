#!/usr/bin/env bash
# Resume S57 A2v1 validation, ADP-only BLS, and experimental teacher scoring.
set -euo pipefail

REPO="${TWIRL_PDO_TEACHER_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
A2V1_ROOT="${TWIRL_A2V1_ROOT:-/pdo/users/tehan/tglc-gpu-production-A2v1}"
PRODUCTION_REPO="${TWIRL_PDO_PRODUCTION_REPO:-/pdo/users/tehan/TWIRL}"
OBSERVATIONS="${TWIRL_OBSERVATIONS:-${PRODUCTION_REPO}/data_local/catalogs/twirl_master_catalog/twirl_wd_tess_observations_v0.fits}"
HLSP_ROOT="${A2V1_ROOT}/hlsp_s0057_A2v1"
HLSP_LOG="${TWIRL_S57_HLSP_LOG:-${A2V1_ROOT}/twirl_logs/s57_a2v1_hlsp_r1.log}"
SCHEMA_SUMMARY="${TWIRL_S57_SCHEMA_SUMMARY:-${A2V1_ROOT}/twirl_logs/s57_A2v1_validation_full_schema.json}"
OUT_ROOT="${TWIRL_S57_TEACHER_ROOT:-/pdo/users/tehan/twirl_stage5/s57_A2v1_teacher_search_v1}"
CHAIN_LOG="${OUT_ROOT}/logs/s57_teacher_chain.log"
BASE_PYTHON="${TWIRL_PDO_BASE_PYTHON:-/sw/qlp-environment/.venv/bin/python}"
POLL_SECONDS="${TWIRL_S57_POLL_SECONDS:-60}"
RUNNER="${REPO}/scripts/stage5_validation/run_s56_a2v1_teacher_search_pdo.sh"

export LD_LIBRARY_PATH="/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

mkdir -p "${OUT_ROOT}/logs"
exec > >(tee -a "${CHAIN_LOG}") 2>&1

log() {
  printf '[s57-teacher-chain] %s %s\n' "$(date -Iseconds)" "$*"
}

json_passed() {
  local path="$1"
  [[ -s "${path}" ]] || return 1
  "${BASE_PYTHON}" - "${path}" <<'PY'
import json, sys
payload = json.load(open(sys.argv[1]))
raise SystemExit(0 if payload.get("ok", payload.get("passed", False)) else 1)
PY
}

wait_for_hlsp() {
  while true; do
    if grep -Eq '\[build-twirl-hlsp\] done: [0-9,]+ ok / 0 fail' "${HLSP_LOG}" 2>/dev/null; then
      log "S57 A2v1 FITS build reports zero failures"
      return
    fi
    if grep -q '\[build-twirl-hlsp\] failure threshold exceeded' "${HLSP_LOG}" 2>/dev/null; then
      log "S57 A2v1 FITS build failed its production threshold"
      exit 20
    fi
    if ! pgrep -u "$(id -u)" -f 'run_a2v1_hlsp_pdo.sh 57 121 122' >/dev/null; then
      log "S57 FITS process ended without a zero-failure completion marker"
      exit 21
    fi
    log "waiting for the existing S57 A2v1 FITS build"
    sleep "${POLL_SECONDS}"
  done
}

validate_sector() {
  if json_passed "${SCHEMA_SUMMARY}"; then
    log "reuse passed S57 schema/coverage validation"
    return
  fi
  log "running full S57 A2v1 schema, finite-value, and coverage validation"
  "${BASE_PYTHON}" "${REPO}/scripts/stage1_lightcurves/validate_a2v1_product.py" \
    --a2v1-root "${A2V1_ROOT}" \
    --hlsp-root "${HLSP_ROOT}" \
    --observations "${OBSERVATIONS}" \
    --sector 57 \
    --orbits 121 122 \
    --allow-edge-warn-missing \
    --fits-workers "${TWIRL_S57_VALIDATE_WORKERS:-16}" \
    --progress-every 500 \
    --summary-json "${SCHEMA_SUMMARY}"
  json_passed "${SCHEMA_SUMMARY}" || {
    log "S57 validation did not pass"
    exit 22
  }
}

run_teacher() {
  log "starting S57 compact export, ADP-only BLS, QA, and native input preparation"
  TWIRL_TEACHER_SECTOR=57 \
  TWIRL_PDO_TEACHER_REPO="${REPO}" \
  TWIRL_PDO_NATIVE_PARALLEL="${TWIRL_PDO_NATIVE_PARALLEL:-4}" \
    bash "${RUNNER}" prep
  log "starting S57 teacher-v1 experimental scoring"
  TWIRL_TEACHER_SECTOR=57 \
  TWIRL_PDO_TEACHER_REPO="${REPO}" \
    bash "${RUNNER}" score
  log "S57 experimental score table and one-row-per-TIC Planet ranking are complete"
}

log "chain started; no human-review queue will be emitted"
wait_for_hlsp
validate_sector
run_teacher
