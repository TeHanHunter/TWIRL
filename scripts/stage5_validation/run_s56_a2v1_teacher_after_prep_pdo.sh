#!/usr/bin/env bash
# Wait for active S56 preparation, then score and run the A2v1 transfer audit.
set -euo pipefail

REPO="${TWIRL_PDO_TEACHER_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
OUT_ROOT="${TWIRL_TEACHER_SEARCH_ROOT:-/pdo/users/tehan/twirl_stage5/s56_A2v1_teacher_search_v1}"
NATIVE_H5="${OUT_ROOT}/inputs/s56_A2v1_adp_raw_pair_v1.h5"
TORCH_PYTHON="${TWIRL_PDO_TORCH_PYTHON:-/pdo/users/tehan/envs/twirl-teacher-pdo-v2/bin/python}"
RUNNER="${REPO}/scripts/stage5_validation/run_s56_a2v1_teacher_search_pdo.sh"
LOG="${OUT_ROOT}/logs/s56_teacher_after_prep.log"
POLL_SECONDS="${TWIRL_S56_POLL_SECONDS:-300}"

export LD_LIBRARY_PATH="/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${REPO}/src"
mkdir -p "${OUT_ROOT}/logs"
exec > >(tee -a "${LOG}") 2>&1

log() {
  printf '[s56-teacher-chain] %s %s\n' "$(date -Iseconds)" "$*"
}

wait_for_native() {
  while [[ ! -s "${NATIVE_H5}" ]]; do
    if ! pgrep -u "$(id -u)" -f 'run_s56_a2v1_teacher_search_pdo.sh prep' >/dev/null; then
      log "S56 preparation ended without the merged native HDF5"
      exit 30
    fi
    local completed
    completed=$(find "${OUT_ROOT}/native_shards" -maxdepth 1 -type f -name 'native_*.h5' 2>/dev/null | wc -l | tr -d ' ')
    log "waiting for S56 native assembly; completed shards=${completed}/16"
    sleep "${POLL_SECONDS}"
  done
  log "S56 merged native HDF5 is ready"
}

wait_for_torch() {
  while true; do
    if [[ -x "${TORCH_PYTHON}" ]] && "${TORCH_PYTHON}" - <<'PY'
import torch
raise SystemExit(0 if torch.cuda.is_available() else 1)
PY
    then
      log "PDO Torch environment sees CUDA"
      return
    fi
    if ! pgrep -u "$(id -u)" -f 'bootstrap_pdo_teacher_torch_env.sh' >/dev/null; then
      log "Torch bootstrap ended but the CUDA import check did not pass"
      exit 31
    fi
    log "waiting for the PDO Torch/CUDA environment"
    sleep "${POLL_SECONDS}"
  done
}

log "continuation started"
wait_for_native
wait_for_torch
log "scoring all S56 A2v1 real candidates with the five-fold teacher ensemble"
TWIRL_TEACHER_SECTOR=56 TWIRL_PDO_TEACHER_REPO="${REPO}" bash "${RUNNER}" score
log "running the labeled-overlap transfer audit"
TWIRL_TEACHER_SECTOR=56 TWIRL_PDO_TEACHER_REPO="${REPO}" bash "${RUNNER}" transfer
log "S56 scores, Planet ranking, and transfer audit are complete"
