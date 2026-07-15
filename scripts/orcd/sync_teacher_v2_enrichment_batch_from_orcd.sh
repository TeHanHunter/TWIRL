#!/bin/bash
set -euo pipefail

SECTOR="${1:-56}"
BATCH_INDEX="${2:-0}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DEST_ROOT="${3:-${REPO_ROOT}/reports/stage5_validation/s56_s64_existing_teacher_enrichment/sector_$(printf '%04d' "${SECTOR}")/batch_$(printf '%02d' "${BATCH_INDEX}")}"
CONTROL_PATH="${ORCD_CONTROL_PATH:-${HOME}/.ssh/cm/%r@%h:%p}"
REMOTE_HOST="${ORCD_HOST:-tehan@orcd-login.mit.edu}"
REMOTE_REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_teacher_v2}"
REMOTE_ROOT="${REMOTE_REPO}/reports/stage5_validation/s56_s64_existing_teacher_enrichment/sector_$(printf '%04d' "${SECTOR}")/batch_$(printf '%02d' "${BATCH_INDEX}")"

expanded_control_path="${CONTROL_PATH//\%r/tehan}"
expanded_control_path="${expanded_control_path//\%h/orcd-login.mit.edu}"
expanded_control_path="${expanded_control_path//\%p/22}"
if [[ ! -S "${expanded_control_path}" ]]; then
  echo "Missing authenticated ORCD socket: ${expanded_control_path}" >&2
  exit 2
fi

mkdir -p "${DEST_ROOT}"
rsync -av --partial \
  --exclude hidden_selection_provenance.parquet \
  -e "ssh -o BatchMode=yes -o PasswordAuthentication=no -o KbdInteractiveAuthentication=no -o NumberOfPasswordPrompts=0 -o ControlMaster=no -o ControlPath=${CONTROL_PATH}" \
  "${REMOTE_HOST}:${REMOTE_ROOT}/" "${DEST_ROOT}/"
echo "[enrichment-sync] ${DEST_ROOT}"
