#!/usr/bin/env bash
# Rebuild the local human-label audit products for the active S56 1k queue.
#
# This is the post-labeling companion to run_s56_1k_small_pair_vetting_app_local.sh.
# It joins labels to hidden metadata, selects the next useful rows to label, and
# runs the readiness gate. It never edits the queue or label CSV.
set -euo pipefail

cd "${TWIRL_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"

PYTHON_BIN="${PYTHON_BIN:-.venv/bin/python}"
export PYTHONPATH="${PYTHONPATH:-src}"

QUEUE_DIR="${QUEUE_DIR:-reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_pdo}"
QUEUE_CSV="${QUEUE_CSV:-${QUEUE_DIR}/review_queue.csv}"
LABELS_CSV="${LABELS_CSV:-${QUEUE_DIR}/human_labels_vetted.csv}"
LABEL_SUMMARY_DIR="${LABEL_SUMMARY_DIR:-${QUEUE_DIR}/human_label_summary}"
TRAINING_DIR="${TRAINING_DIR:-${QUEUE_DIR}/human_training_table}"
PRIORITY_DIR="${PRIORITY_DIR:-${QUEUE_DIR}/human_label_priority_next}"
READINESS_DIR="${READINESS_DIR:-${QUEUE_DIR}/human_training_readiness}"
NEXT_N_ROWS="${NEXT_N_ROWS:-200}"
TARGET_PER_CELL="${TARGET_PER_CELL:-25}"

if [[ ! -f "${QUEUE_CSV}" ]]; then
  echo "[human-label-audit] missing queue: ${QUEUE_CSV}" >&2
  exit 2
fi
if [[ ! -f "${LABELS_CSV}" ]]; then
  echo "[human-label-audit] missing labels: ${LABELS_CSV}" >&2
  echo "[human-label-audit] run the vetting app and label at least one row first" >&2
  exit 2
fi

echo "[human-label-audit] queue=${QUEUE_CSV}"
echo "[human-label-audit] labels=${LABELS_CSV}"
echo "[human-label-audit] label_summary=${LABEL_SUMMARY_DIR}"
echo "[human-label-audit] training_table=${TRAINING_DIR}"
echo "[human-label-audit] priority=${PRIORITY_DIR}"
echo "[human-label-audit] readiness=${READINESS_DIR}"

"${PYTHON_BIN}" scripts/stage5_validation/summarize_human_vetting_labels.py \
  --queue-csv "${QUEUE_CSV}" \
  --labels-csv "${LABELS_CSV}" \
  --out-dir "${LABEL_SUMMARY_DIR}"

"${PYTHON_BIN}" scripts/stage5_validation/build_human_vetting_training_table.py \
  --queue-csv "${QUEUE_CSV}" \
  --labels-csv "${LABELS_CSV}" \
  --out-dir "${TRAINING_DIR}"

"${PYTHON_BIN}" scripts/stage5_validation/select_next_human_labels.py \
  --training-table "${TRAINING_DIR}/human_vetting_training_table.csv" \
  --out-dir "${PRIORITY_DIR}" \
  --n-rows "${NEXT_N_ROWS}" \
  --target-per-cell "${TARGET_PER_CELL}"

"${PYTHON_BIN}" scripts/stage5_validation/audit_human_training_readiness.py \
  --training-table "${TRAINING_DIR}/human_vetting_training_table.csv" \
  --priority-table "${PRIORITY_DIR}/next_label_priority.csv" \
  --out-dir "${READINESS_DIR}" \
  --target-per-coverage-cell "${TARGET_PER_CELL}"

echo "[human-label-audit] complete"
echo "[human-label-audit] readiness summary: ${READINESS_DIR}/summary.md"
