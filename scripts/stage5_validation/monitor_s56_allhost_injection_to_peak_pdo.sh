#!/usr/bin/env bash
# Wait for all-host injection shards, then run the sharded robust-BLS peak table.
set -euo pipefail

REPO_ROOT="${TWIRL_REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
cd "${REPO_ROOT}"

PYTHON_BIN="${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}"
export LD_LIBRARY_PATH="/sw/python-versions/python-3.11.9/lib:/pdo/app/anaconda/anaconda2-4.4.0/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${PYTHONPATH:-src}"
WAIT_SECONDS="${WAIT_SECONDS:-120}"
INJECTION_TMUX_SESSION="${INJECTION_TMUX_SESSION:-twirl-s56-allhost-injections}"
PROGRESS_DIR="${PROGRESS_DIR:-reports/stage5_validation/s56_allhost_predetrend_batman_periodradius_grid_sharded_pdo/progress}"
PEAK_LOG="${PEAK_LOG:-reports/stage5_validation/s56_allhost_predetrend_batman_periodradius_grid_sharded_pdo/peak_training/run.log}"
PEAK_CHUNK_JOBS="${TWIRL_ALLHOST_PEAK_CHUNK_JOBS:-1}"
PEAK_WORKERS="${TWIRL_PEAK_WORKERS:-8}"

log() {
  printf '[%s] [allhost-monitor] %s\n' "$(date -Is)" "$*" >&2
}

is_complete() {
  "${PYTHON_BIN}" - "${PROGRESS_DIR}/progress_summary.json" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
if not path.exists():
    raise SystemExit(1)
summary = json.loads(path.read_text())
complete = (
    summary.get("n_shards", 0) > 0
    and summary.get("n_completed_shards") == summary.get("n_shards")
    and summary.get("combined_summary_exists")
    and summary.get("combined_manifest_exists")
)
raise SystemExit(0 if complete else 1)
PY
}

mkdir -p "${PROGRESS_DIR}" "$(dirname "${PEAK_LOG}")"

while true; do
  log "refreshing all-host injection progress"
  "${PYTHON_BIN}" scripts/stage5_validation/summarize_s56_allhost_injection_progress.py \
    --out-dir "${PROGRESS_DIR}" >/dev/null
  if is_complete; then
    log "all-host injection product is complete"
    break
  fi
  if ! tmux has-session -t "${INJECTION_TMUX_SESSION}" 2>/dev/null; then
    log "injection tmux session ended before combined summary/manifest were complete"
    cat "${PROGRESS_DIR}/progress_summary.md" >&2 || true
    exit 1
  fi
  cat "${PROGRESS_DIR}/progress_summary.md" >&2 || true
  sleep "${WAIT_SECONDS}"
done

log "starting sharded all-host BLS peak table"
TWIRL_ALLHOST_PEAK_CHUNK_JOBS="${PEAK_CHUNK_JOBS}" \
TWIRL_PEAK_WORKERS="${PEAK_WORKERS}" \
  bash scripts/stage5_validation/run_s56_allhost_peak_training_sharded_pdo.sh >"${PEAK_LOG}" 2>&1

log "complete"
log "peak_log=${PEAK_LOG}"
