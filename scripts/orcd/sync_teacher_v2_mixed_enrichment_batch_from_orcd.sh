#!/bin/bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
REMOTE_HOST="${TWIRL_ORCD_HOST:-tehan@orcd-login.mit.edu}"
REMOTE_REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_teacher_v2}"
REMOTE_ROOT="${REMOTE_REPO}/reports/stage5_validation/s56_s57_existing_teacher_enrichment/batch_00"
DEST_ROOT="${1:-${REPO_ROOT}/reports/stage5_validation/s56_s57_existing_teacher_enrichment/batch_00}"
CONTROL_PATH="${TWIRL_ORCD_CONTROL_PATH:-${HOME}/.ssh/cm/%r@%h:%p}"

mkdir -p "${DEST_ROOT}/twirl_vet_sheets"
SSH_COMMAND="ssh -o BatchMode=yes -o PasswordAuthentication=no -o KbdInteractiveAuthentication=no -o NumberOfPasswordPrompts=0 -o ControlMaster=no -o ControlPath=${CONTROL_PATH}"
rsync -a --delete \
  --exclude hidden_selection_provenance.parquet \
  --exclude '*.pdf' \
  -e "${SSH_COMMAND}" \
  "${REMOTE_HOST}:${REMOTE_ROOT}/" "${DEST_ROOT}/"
echo "[sync] ${DEST_ROOT}"
