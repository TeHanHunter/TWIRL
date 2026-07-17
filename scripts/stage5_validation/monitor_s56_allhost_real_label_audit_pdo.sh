#!/usr/bin/env bash
# Watch the all-host real-candidate human-label CSV and rerun readiness audits.
#
# This is intentionally cheap and idempotent. It never edits labels or the
# review queue; it only calls run_s56_allhost_real_label_audit_pdo.sh after the
# label CSV exists and the labeled-row count changes.
set -euo pipefail

REPO_ROOT="${TWIRL_REPO_ROOT:-/pdo/users/tehan/TWIRL}"
cd "${REPO_ROOT}"

QUEUE_DIR="${QUEUE_DIR:-reports/stage5_validation/s56_allhost_ranker_selected_real_leo_queue_pdo}"
LABELS_CSV="${LABELS_CSV:-${QUEUE_DIR}/human_labels_vetted.csv}"
READINESS_JSON="${READINESS_JSON:-${QUEUE_DIR}/human_training_readiness/summary.json}"
WAIT_SECONDS="${WAIT_SECONDS:-300}"
MIN_LABEL_ROWS="${MIN_LABEL_ROWS:-1}"
RUN_ONCE="${RUN_ONCE:-0}"

log() {
  printf '[%s] [allhost-label-audit-monitor] %s\n' "$(date -Is)" "$*" >&2
}

label_count() {
  if [[ ! -s "${LABELS_CSV}" ]]; then
    printf '0\n'
    return
  fi
  awk 'END {print NR > 0 ? NR - 1 : 0}' "${LABELS_CSV}"
}

audited_count() {
  if [[ ! -s "${READINESS_JSON}" ]]; then
    printf -- '-1\n'
    return
  fi
  python - "${READINESS_JSON}" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
try:
    payload = json.loads(path.read_text())
except Exception:
    print(-1)
else:
    print(int(payload.get("n_labeled", -1)))
PY
}

log "labels=${LABELS_CSV}"
log "readiness=${READINESS_JSON}"
log "wait_seconds=${WAIT_SECONDS}; min_label_rows=${MIN_LABEL_ROWS}; run_once=${RUN_ONCE}"

while true; do
  current_n="$(label_count)"
  previous_n="$(audited_count)"

  if (( current_n < MIN_LABEL_ROWS )); then
    log "waiting for labels: current=${current_n}, required=${MIN_LABEL_ROWS}"
  elif (( current_n == previous_n )); then
    log "audit is current for ${current_n} labels"
    if [[ "${RUN_ONCE}" == "1" ]]; then
      exit 0
    fi
  else
    log "running post-label audit: previous=${previous_n}, current=${current_n}"
    bash scripts/stage5_validation/run_s56_allhost_real_label_audit_pdo.sh
    log "post-label audit complete"
    if [[ "${RUN_ONCE}" == "1" ]]; then
      exit 0
    fi
  fi

  sleep "${WAIT_SECONDS}"
done
