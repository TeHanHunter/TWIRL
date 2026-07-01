#!/usr/bin/env bash
# Sync the PDO mixed-teacher queue, LEO sheets, labels, and summaries locally.
set -euo pipefail

cd "${TWIRL_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"

REMOTE_HOST="${REMOTE_HOST:-pdogpu1}"
REMOTE_REPO="${REMOTE_REPO:-/pdo/users/tehan/TWIRL}"
REMOTE_DIR="${REMOTE_DIR:-${REMOTE_REPO}/reports/stage5_validation/s56_mixed_teacher_queue_pdo}"
LOCAL_DIR="${LOCAL_DIR:-reports/stage5_validation/s56_mixed_teacher_queue_pdo}"

mkdir -p "${LOCAL_DIR}"

echo "[sync-mixed-teacher] remote=${REMOTE_HOST}:${REMOTE_DIR}/"
echo "[sync-mixed-teacher] local=${LOCAL_DIR}/"
rsync -avhu --progress \
  --include='*/' \
  --include='*.csv' \
  --include='*.json' \
  --include='*.md' \
  --include='*.pdf' \
  --exclude='*' \
  "${REMOTE_HOST}:${REMOTE_DIR}/" \
  "${LOCAL_DIR}/"

echo "[sync-mixed-teacher] complete"
echo "[sync-mixed-teacher] queue=${LOCAL_DIR}/review_queue_1k.csv"
echo "[sync-mixed-teacher] labels=${LOCAL_DIR}/human_labels_vetted.csv"
