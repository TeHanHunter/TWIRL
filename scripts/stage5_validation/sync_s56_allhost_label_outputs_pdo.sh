#!/usr/bin/env bash
# Pull small all-host human-label audit products from PDO to the local checkout.
#
# This is the local companion to monitor_s56_allhost_real_label_audit_pdo.sh.
# It does not edit remote labels or queues. It only syncs the label CSV plus
# derived audit products and refreshes the combined ML handoff readiness report.
set -euo pipefail

LOCAL_REPO="${LOCAL_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
PDO_HOST="${PDO_HOST:-pdogpu1}"
PDO_REPO="${PDO_REPO:-/pdo/users/tehan/TWIRL}"
QUEUE_DIR="${QUEUE_DIR:-reports/stage5_validation/s56_allhost_ranker_selected_real_leo_queue_pdo}"
WAIT_SECONDS="${WAIT_SECONDS:-300}"
WATCH=0
RUN_AUDIT="${RUN_AUDIT:-1}"

usage() {
  cat <<'EOF'
Usage: sync_s56_allhost_label_outputs_pdo.sh [--watch]

Copies small all-host human-label products from PDO to the local checkout:

  human_labels_vetted.csv
  human_label_summary/
  human_training_table/
  human_label_priority_next/
  human_training_readiness/

If RUN_AUDIT=1, refreshes reports/stage5_validation/s56_ml_handoff_readiness/
after each sync attempt.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --watch)
      WATCH=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "[sync-pdo-labels] unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

PDO_SSH=(
  ssh
  -o BatchMode=yes
  -o ConnectTimeout=20
  -o ControlMaster=no
)

log() {
  printf '[%s] [sync-pdo-labels] %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" "$*" >&2
}

remote_exists() {
  local rel="$1"
  "${PDO_SSH[@]}" "${PDO_HOST}" "test -e '${PDO_REPO}/${rel}'" >/dev/null 2>&1
}

sync_relpath() {
  local rel="$1"
  local local_parent
  local_parent="$(dirname "${LOCAL_REPO}/${rel}")"
  mkdir -p "${local_parent}"
  rsync -avhu --progress "${PDO_HOST}:${PDO_REPO}/${rel}" "${local_parent}/"
}

run_once() {
  local synced=0
  local rel
  local paths=(
    "${QUEUE_DIR}/human_labels_vetted.csv"
    "${QUEUE_DIR}/human_label_summary"
    "${QUEUE_DIR}/human_training_table"
    "${QUEUE_DIR}/human_label_priority_next"
    "${QUEUE_DIR}/human_training_readiness"
  )

  for rel in "${paths[@]}"; do
    if remote_exists "${rel}"; then
      log "sync ${rel}"
      sync_relpath "${rel}"
      synced=1
    else
      log "pending ${rel}"
    fi
  done

  if [[ "${RUN_AUDIT}" == "1" ]]; then
    log "refresh local ML handoff readiness"
    (cd "${LOCAL_REPO}" && python scripts/stage5_validation/audit_s56_ml_handoff_readiness.py >/dev/null)
  fi

  if (( synced > 0 )); then
    return 0
  fi
  return 1
}

log "pdo=${PDO_HOST}:${PDO_REPO}/${QUEUE_DIR}"
log "local=${LOCAL_REPO}/${QUEUE_DIR}"
log "watch=${WATCH}; wait_seconds=${WAIT_SECONDS}; run_audit=${RUN_AUDIT}"

while true; do
  if run_once; then
    log "sync pass copied at least one artifact"
  else
    log "sync pass found no available label artifacts yet"
  fi
  if [[ "${WATCH}" != "1" ]]; then
    exit 0
  fi
  sleep "${WAIT_SECONDS}"
done
