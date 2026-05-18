#!/bin/bash
# Pull per-sector QC PDFs from pdogpu6 to the local reports/stage1_lightcurves/qc_pdf/.
# Designed to be run on the laptop (where this repo lives) on demand or via
# cron. Only PDFs are pulled — the companion .npz files stay on PDO unless
# explicitly requested with --with-npz.
#
# Usage:
#   bash scripts/sync_qc_pdfs.sh                # PDFs only
#   bash scripts/sync_qc_pdfs.sh --with-npz     # PDFs + raw .npz
#
# Add to crontab for fully automatic sync, e.g. every 30 min:
#   */30 * * * * cd $HOME/PycharmProjects/TWIRL && bash scripts/sync_qc_pdfs.sh >> /tmp/twirl_qc_sync.log 2>&1

set -uo pipefail

REMOTE_HOST=${TWIRL_PDO_HOST:-pdogpu6}
REMOTE_QC_DIR=/pdo/users/tehan/tglc-gpu-production/qc_pdf
LOCAL_QC_DIR=$(cd "$(dirname "$0")"/.. && pwd)/reports/stage1_lightcurves/qc_pdf

mkdir -p "$LOCAL_QC_DIR"

INCLUDES=("--include=s*_qc.pdf")
if [[ "${1:-}" == "--with-npz" ]]; then
  INCLUDES+=("--include=s*_qc.npz")
fi

echo "[sync] $REMOTE_HOST:$REMOTE_QC_DIR -> $LOCAL_QC_DIR"
rsync -av --no-perms --no-owner --no-group \
  "${INCLUDES[@]}" --exclude='*' \
  "$REMOTE_HOST:$REMOTE_QC_DIR/" "$LOCAL_QC_DIR/"
echo "[sync] done"
