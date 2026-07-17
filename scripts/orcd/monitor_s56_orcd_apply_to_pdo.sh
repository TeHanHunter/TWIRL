#!/usr/bin/env bash
# Wait for ORCD ranker-apply outputs, then sync them back to PDO.
#
# This does not initiate ORCD authentication. It reuses the user-opened ORCD
# control socket via sync_s56_orcd_ranker_outputs_to_pdo.sh and exits if that
# socket is unavailable. By default it only syncs the verified selected
# ephemerides and summaries; set TWIRL_RUN_PDO_LEO=1 to launch PDO LEO rendering
# after the sync succeeds.
set -euo pipefail

LOCAL_REPO="${LOCAL_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
PDO_HOST="${PDO_HOST:-pdogpu1}"
PDO_REPO="${PDO_REPO:-/pdo/users/tehan/TWIRL}"
DEST_DIR="${TWIRL_PDO_SELECTED_DIR:-reports/stage5_validation/s56_ranker_selected_real_candidates_orcd}"
WAIT_SECONDS="${WAIT_SECONDS:-600}"
MAX_ATTEMPTS="${MAX_ATTEMPTS:-0}"
RUN_PDO_LEO="${TWIRL_RUN_PDO_LEO:-0}"
LEO_OUT_DIR="${TWIRL_PDO_LEO_OUT_DIR:-reports/stage5_validation/s56_ranker_selected_real_leo_queue_orcd_pdo}"

log() {
  printf '[%s] [orcd-apply-monitor] %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" "$*" >&2
}

attempt=0
while true; do
  attempt=$((attempt + 1))
  log "checking ORCD selected-output readiness; attempt=${attempt}"
  set +e
  bash "${LOCAL_REPO}/scripts/orcd/sync_s56_orcd_ranker_outputs_to_pdo.sh" --run
  rc=$?
  set -e
  if [[ "${rc}" -eq 0 ]]; then
    log "ORCD selected outputs synced to PDO"
    break
  fi
  if [[ "${rc}" -ne 5 ]]; then
    log "sync helper failed with rc=${rc}; stopping"
    exit "${rc}"
  fi
  if [[ "${MAX_ATTEMPTS}" -gt 0 && "${attempt}" -ge "${MAX_ATTEMPTS}" ]]; then
    log "reached MAX_ATTEMPTS=${MAX_ATTEMPTS}; selected outputs are not ready"
    exit 5
  fi
  log "selected outputs are not complete yet; sleeping ${WAIT_SECONDS}s"
  sleep "${WAIT_SECONDS}"
done

if [[ "${RUN_PDO_LEO}" != "1" ]]; then
  log "PDO LEO rendering not requested"
  log "next command:"
  cat >&2 <<EOF
ssh ${PDO_HOST} 'cd ${PDO_REPO} && TWIRL_SELECTED_EPHEMERIDES=${DEST_DIR}/selected_ephemerides.csv TWIRL_RANKER_LEO_OUT_DIR=${LEO_OUT_DIR} bash scripts/stage5_validation/run_s56_ranker_selected_real_leo_pdo.sh'
EOF
  exit 0
fi

log "launching PDO LEO rendering"
ssh -o BatchMode=yes -o ConnectTimeout=20 "${PDO_HOST}" \
  "cd '${PDO_REPO}' && TWIRL_SELECTED_EPHEMERIDES='${DEST_DIR}/selected_ephemerides.csv' TWIRL_RANKER_LEO_OUT_DIR='${LEO_OUT_DIR}' bash scripts/stage5_validation/run_s56_ranker_selected_real_leo_pdo.sh"
log "complete"
