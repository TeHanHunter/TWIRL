#!/usr/bin/env bash
# Rebuild human-label audit products for the verified all-host real-candidate queue.
#
# Run this on PDO after labeling rows in run_s56_allhost_real_vetting_app_pdo.sh.
# It never edits the review queue or label CSV.
set -euo pipefail

REPO_ROOT="${TWIRL_REPO_ROOT:-/pdo/users/tehan/TWIRL}"
cd "${REPO_ROOT}"

PYTHON_BIN="${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}"
export LD_LIBRARY_PATH="/pdo/app/anaconda/anaconda2-4.4.0/lib:/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${PYTHONPATH:-src}"

QUEUE_DIR="${QUEUE_DIR:-reports/stage5_validation/s56_allhost_ranker_selected_real_leo_queue_pdo}"
QUEUE_CSV="${QUEUE_CSV:-${QUEUE_DIR}/review_queue.csv}"
LABELS_CSV="${LABELS_CSV:-${QUEUE_DIR}/human_labels_vetted.csv}"
VERIFY_JSON="${VERIFY_JSON:-${QUEUE_DIR}/verification.json}"
LABEL_SUMMARY_DIR="${LABEL_SUMMARY_DIR:-${QUEUE_DIR}/human_label_summary}"
TRAINING_DIR="${TRAINING_DIR:-${QUEUE_DIR}/human_training_table}"
PRIORITY_DIR="${PRIORITY_DIR:-${QUEUE_DIR}/human_label_priority_next}"
READINESS_DIR="${READINESS_DIR:-${QUEUE_DIR}/human_training_readiness}"
NEXT_N_ROWS="${NEXT_N_ROWS:-200}"
TARGET_PER_CELL="${TARGET_PER_CELL:-25}"

if [[ ! -s "${VERIFY_JSON}" ]]; then
  echo "[allhost-real-label-audit] missing verifier: ${VERIFY_JSON}" >&2
  exit 2
fi
if ! "${PYTHON_BIN}" - "${VERIFY_JSON}" <<'PY'
import json
import sys
from pathlib import Path

payload = json.loads(Path(sys.argv[1]).read_text())
raise SystemExit(0 if payload.get("passed") else 1)
PY
then
  echo "[allhost-real-label-audit] verifier did not pass: ${VERIFY_JSON}" >&2
  exit 2
fi
if [[ ! -s "${LABELS_CSV}" ]]; then
  echo "[allhost-real-label-audit] missing labels: ${LABELS_CSV}" >&2
  echo "[allhost-real-label-audit] label rows in the vetting app before running this audit" >&2
  exit 2
fi

echo "[allhost-real-label-audit] queue=${QUEUE_CSV}"
echo "[allhost-real-label-audit] labels=${LABELS_CSV}"
echo "[allhost-real-label-audit] label_summary=${LABEL_SUMMARY_DIR}"
echo "[allhost-real-label-audit] training_table=${TRAINING_DIR}"
echo "[allhost-real-label-audit] priority=${PRIORITY_DIR}"
echo "[allhost-real-label-audit] readiness=${READINESS_DIR}"

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

echo "[allhost-real-label-audit] complete"
echo "[allhost-real-label-audit] readiness summary: ${READINESS_DIR}/summary.md"
