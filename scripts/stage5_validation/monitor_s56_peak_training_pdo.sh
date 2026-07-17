#!/usr/bin/env bash
# Monitor the S56 20k peak table and run the post-BLS gate when it is ready.
#
# Intended PDO usage:
#   bash scripts/stage5_validation/monitor_s56_peak_training_pdo.sh --wait
#   bash scripts/stage5_validation/monitor_s56_peak_training_pdo.sh --wait --run-review --skip-leo
set -euo pipefail

REPO_ROOT="${TWIRL_REPO_ROOT:-/pdo/users/tehan/TWIRL}"
cd "${REPO_ROOT}"

PYTHON="${PYTHON:-/sw/qlp-environment/.venv/bin/python}"
export LD_LIBRARY_PATH="/pdo/app/anaconda/anaconda2-4.4.0/lib:/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

ROOT="${TWIRL_20K_ROOT:-reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo}"
PEAK_ROOT="${TWIRL_PEAK_ROOT:-${ROOT}/peak_training}"
CHUNK_DIR="${TWIRL_PEAK_CHUNK_DIR:-${PEAK_ROOT}/chunk_ids}"
RUN_CHUNK_DIR="${TWIRL_PEAK_RUN_CHUNK_DIR:-${PEAK_ROOT}/chunks}"
CHUNKED_TABLE="${TWIRL_CHUNKED_PEAK_TABLE:-${PEAK_ROOT}/s56_20k_injection_bls_peaks_chunked.csv}"
ONESHOT_TABLE="${TWIRL_ONESHOT_PEAK_TABLE:-${PEAK_ROOT}/s56_20k_injection_bls_peaks.csv}"
GATE_OUT_DIR="${TWIRL_GATE_OUT_DIR:-${ROOT}/peak_training_gate_pdo}"
RANKER_OUT_DIR="${TWIRL_RANKER_OUT_DIR:-${ROOT}/peak_ranker_pdo}"
LOG_PATH="${TWIRL_PEAK_LOG:-logs/s56_20k_peak_training_pdo.log}"
CHUNK_PROGRESS_OUT_DIR="${TWIRL_CHUNK_PROGRESS_OUT_DIR:-${PEAK_ROOT}/chunk_progress_live}"
WRITE_CHUNK_PROGRESS="${TWIRL_MONITOR_WRITE_CHUNK_PROGRESS:-1}"
SLEEP_SECONDS="${TWIRL_MONITOR_SLEEP_SECONDS:-300}"
MAX_WAIT_MINUTES="${TWIRL_MONITOR_MAX_WAIT_MINUTES:-0}"
WAIT=0
RUN_REVIEW=0
SKIP_LEO=0

usage() {
  cat <<'EOF'
Usage: monitor_s56_peak_training_pdo.sh [--wait] [--run-review] [--skip-leo]

Options:
  --wait        Poll until a peak table exists. Default is one status check.
  --run-review Run run_s56_peak_ranker_review_pdo.sh after the gate report.
  --skip-leo   When used with --run-review, build the review queue without PDFs.

Environment:
  TWIRL_MONITOR_SLEEP_SECONDS    Poll interval for --wait. Default: 300.
  TWIRL_MONITOR_MAX_WAIT_MINUTES 0 means no timeout. Default: 0.
  TWIRL_MONITOR_WRITE_CHUNK_PROGRESS
                                  Write read-only partial recall reports while
                                  waiting. Default: 1.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --wait)
      WAIT=1
      shift
      ;;
    --run-review)
      RUN_REVIEW=1
      shift
      ;;
    --skip-leo)
      SKIP_LEO=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "[peak-monitor] unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

log() { echo "[$(date -Is)] [peak-monitor] $*"; }

find_peak_table() {
  if [[ -f "${CHUNKED_TABLE}" ]]; then
    printf '%s\n' "${CHUNKED_TABLE}"
  elif [[ -f "${ONESHOT_TABLE}" ]]; then
    printf '%s\n' "${ONESHOT_TABLE}"
  fi
}

show_chunk_status() {
  local manifest n_chunks done_count active_chunks active_logs
  manifest="${CHUNK_DIR}/manifest.txt"
  if [[ ! -f "${manifest}" && ! -d "${RUN_CHUNK_DIR}" ]]; then
    return 0
  fi
  n_chunks="?"
  if [[ -f "${manifest}" ]]; then
    n_chunks="$(awk -F= '$1 == "n_chunks" {print $2}' "${manifest}")"
  fi
  done_count=0
  if [[ -d "${RUN_CHUNK_DIR}" ]]; then
    done_count="$(find "${RUN_CHUNK_DIR}" -mindepth 2 -maxdepth 2 -name injection_bls_peaks_summary.json 2>/dev/null | wc -l | tr -d ' ')"
  fi
  active_chunks=""
  if [[ -d "${RUN_CHUNK_DIR}" ]]; then
    active_chunks="$(
      find "${RUN_CHUNK_DIR}" -mindepth 1 -maxdepth 1 -type d -name 'chunk_*' 2>/dev/null |
        while read -r chunk_path; do
          [[ -f "${chunk_path}/injection_bls_peaks_summary.json" ]] || basename "${chunk_path}"
        done |
        sort |
        head -n 8 |
        paste -sd, -
    )"
  fi
  log "chunk progress: ${done_count}/${n_chunks} complete; active=${active_chunks:-none}"
  active_logs="$(
    find "${RUN_CHUNK_DIR}" -mindepth 2 -maxdepth 2 -name run.log 2>/dev/null |
      sort |
      tail -n 2
  )"
  if [[ -n "${active_logs}" ]]; then
    while read -r run_log; do
      [[ -n "${run_log}" ]] || continue
      log "tail ${run_log}"
      tail -n 6 "${run_log}" || true
    done <<< "${active_logs}"
  fi
}

show_worker_status() {
  local user_name worker_lines worker_count
  user_name="$(id -un)"
  worker_lines="$(
    ps -u "${user_name}" -o pid=,ppid=,stat=,pcpu=,etime=,command= 2>/dev/null |
      awk '/build_injection_peak_training_table.py/ && !/awk/ {print}' |
      sort -k4,4nr
  )"
  if [[ -z "${worker_lines}" ]]; then
    log "active peak-training workers: none"
    return 0
  fi
  worker_count="$(printf '%s\n' "${worker_lines}" | wc -l | tr -d ' ')"
  log "active peak-training workers: ${worker_count}"
  printf '%s\n' "${worker_lines}" | head -n 6 | sed 's/^/[peak-monitor] worker /'
}

show_chunk_progress_report() {
  if [[ "${WRITE_CHUNK_PROGRESS}" != "1" ]]; then
    return 0
  fi
  if [[ ! -d "${RUN_CHUNK_DIR}" ]]; then
    return 0
  fi
  if [[ ! -f scripts/stage5_validation/summarize_injection_peak_chunk_progress.py ]]; then
    return 0
  fi
  log "updating chunk progress report: ${CHUNK_PROGRESS_OUT_DIR}"
  "${PYTHON}" scripts/stage5_validation/summarize_injection_peak_chunk_progress.py \
    --chunk-root "${RUN_CHUNK_DIR}" \
    --out-dir "${CHUNK_PROGRESS_OUT_DIR}" |
    sed 's/^/[peak-monitor] progress /'
}

show_status() {
  log "repo=${REPO_ROOT}"
  log "chunked_table=${CHUNKED_TABLE}"
  log "oneshot_table=${ONESHOT_TABLE}"
  show_chunk_status
  show_chunk_progress_report
  show_worker_status
  if command -v tmux >/dev/null 2>&1 && tmux has-session -t twirl-s56-20k-peak-training 2>/dev/null; then
    log "tmux twirl-s56-20k-peak-training: running"
  else
    log "tmux twirl-s56-20k-peak-training: not running"
  fi
  if [[ -f "${LOG_PATH}" ]]; then
    log "tail ${LOG_PATH}"
    tail -n 12 "${LOG_PATH}" || true
  else
    log "missing log: ${LOG_PATH}"
  fi
}

start_epoch="$(date +%s)"
peak_table="$(find_peak_table || true)"
while [[ -z "${peak_table}" ]]; do
  show_status
  if [[ "${WAIT}" != "1" ]]; then
    log "peak table not ready"
    exit 1
  fi
  if [[ "${MAX_WAIT_MINUTES}" -gt 0 ]]; then
    now_epoch="$(date +%s)"
    elapsed_min="$(( (now_epoch - start_epoch) / 60 ))"
    if [[ "${elapsed_min}" -ge "${MAX_WAIT_MINUTES}" ]]; then
      log "timed out after ${elapsed_min} minutes"
      exit 3
    fi
  fi
  log "sleeping ${SLEEP_SECONDS}s"
  sleep "${SLEEP_SECONDS}"
  peak_table="$(find_peak_table || true)"
done

log "peak table ready: ${peak_table}"
mkdir -p "${GATE_OUT_DIR}"
"${PYTHON}" scripts/stage5_validation/summarize_injection_peak_gate.py \
  --peak-table "${peak_table}" \
  --out-dir "${GATE_OUT_DIR}" \
  --top-k 20

if [[ "${RUN_REVIEW}" == "1" ]]; then
  log "running ranker/review handoff"
  export TWIRL_PEAK_TABLE="${peak_table}"
  export TWIRL_GATE_OUT_DIR="${GATE_OUT_DIR}"
  export TWIRL_RANKER_OUT_DIR="${RANKER_OUT_DIR}"
  if [[ "${SKIP_LEO}" == "1" ]]; then
    export TWIRL_SKIP_LEO=1
  fi
  if ! bash scripts/stage5_validation/run_s56_peak_ranker_review_pdo.sh; then
    log "ranker/review handoff failed"
    log "search-gate outputs remain in ${GATE_OUT_DIR}"
    log "if the injected peak-table schema verification failed, rerun with the restartable chunked builder"
    exit 4
  fi
fi

log "complete"
