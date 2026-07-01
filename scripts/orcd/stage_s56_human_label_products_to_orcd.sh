#!/usr/bin/env bash
# Stage small human-label/readiness products to ORCD for H200 training.
#
# This intentionally does not transfer light curves or injection HDF5 files. It
# copies only the all-host real-candidate queue metadata, labels, teacher rows,
# and readiness products needed by the bounded H200 real-label training smoke.
set -euo pipefail

DRY_RUN=1
ALLOW_INCOMPLETE=0

usage() {
  cat <<'EOF'
Usage: stage_s56_human_label_products_to_orcd.sh [--run] [--allow-incomplete]

Stages small human-label products to ORCD. Default mode is dry-run. By default,
the human-training readiness summary must recommend either:

  ready_for_binary_teacher_smoke
  ready_for_object_teacher_training

Use --allow-incomplete only for manual debugging; do not use it before the real
H200 training smoke.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run)
      DRY_RUN=0
      shift
      ;;
    --allow-incomplete)
      ALLOW_INCOMPLETE=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "[stage-labels-orcd] unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

LOCAL_REPO="${LOCAL_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
ORCD_HOST="${ORCD_HOST:-tehan@orcd-login.mit.edu}"
ORCD_REPO="${ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
ORCD_CONTROL_PATH="${ORCD_CONTROL_PATH:-$HOME/.ssh/cm/%r@%h:%p}"
QUEUE_DIR="${QUEUE_DIR:-reports/stage5_validation/s56_allhost_ranker_selected_real_leo_queue_pdo}"

ORCD_SSH=(
  ssh
  -o BatchMode=yes
  -o PasswordAuthentication=no
  -o KbdInteractiveAuthentication=no
  -o NumberOfPasswordPrompts=0
  -o ControlMaster=auto
  -o ControlPath="${ORCD_CONTROL_PATH}"
)

print_cmd() {
  printf '+'
  printf ' %q' "$@"
  printf '\n'
}

echo "[stage-labels-orcd] local=${LOCAL_REPO}/${QUEUE_DIR}"
echo "[stage-labels-orcd] orcd=${ORCD_HOST}:${ORCD_REPO}/${QUEUE_DIR}"
if [[ "${DRY_RUN}" == "1" ]]; then
  echo "[stage-labels-orcd] dry run. Re-run with --run after ORCD socket is open and readiness passes."
fi

required=(
  "${QUEUE_DIR}/review_queue.csv"
  "${QUEUE_DIR}/verification.json"
  "${QUEUE_DIR}/summary.json"
  "${QUEUE_DIR}/human_labels_vetted.csv"
  "${QUEUE_DIR}/human_label_summary/summary.json"
  "${QUEUE_DIR}/human_training_table/human_vetting_training_table.csv"
  "${QUEUE_DIR}/human_training_table/teacher_labeled_rows.csv"
  "${QUEUE_DIR}/human_training_table/audit_labeled_rows.csv"
  "${QUEUE_DIR}/human_training_table/summary.json"
  "${QUEUE_DIR}/human_label_priority_next/next_label_priority.csv"
  "${QUEUE_DIR}/human_label_priority_next/summary.json"
  "${QUEUE_DIR}/human_training_readiness/summary.json"
  "${QUEUE_DIR}/human_training_readiness/summary.md"
  "${QUEUE_DIR}/human_training_readiness/label_deficits.csv"
)

missing=0
for rel in "${required[@]}"; do
  if [[ ! -s "${LOCAL_REPO}/${rel}" ]]; then
    echo "[stage-labels-orcd] missing required artifact: ${rel}" >&2
    missing=1
  fi
done
if [[ "${missing}" == "1" && "${ALLOW_INCOMPLETE}" != "1" && "${DRY_RUN}" != "1" ]]; then
  echo "[stage-labels-orcd] refusing to stage incomplete label products" >&2
  exit 2
fi

readiness_json="${LOCAL_REPO}/${QUEUE_DIR}/human_training_readiness/summary.json"
recommendation="$(
  python - "${readiness_json}" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
if not path.exists():
    print("")
    raise SystemExit(0)
payload = json.loads(path.read_text())
print(payload.get("recommendation", ""))
PY
)"
echo "[stage-labels-orcd] readiness=${recommendation:-missing}"
case "${recommendation}" in
  ready_for_binary_teacher_smoke|ready_for_object_teacher_training)
    ;;
  *)
    if [[ "${ALLOW_INCOMPLETE}" != "1" && "${DRY_RUN}" != "1" ]]; then
      echo "[stage-labels-orcd] readiness is not sufficient for real H200 training smoke" >&2
      exit 3
    fi
    ;;
esac

if [[ "${DRY_RUN}" == "1" ]]; then
  printf '%s\n' "${required[@]}" | sed 's/^/[stage-labels-orcd] would stage /'
  print_cmd "${ORCD_SSH[@]}" "${ORCD_HOST}" "mkdir -p '${ORCD_REPO}'"
  print_cmd tar -C "${LOCAL_REPO}" -cf - "${required[@]}"
  print_cmd "${ORCD_SSH[@]}" "${ORCD_HOST}" "tar -C '${ORCD_REPO}' -xf -"
else
  "${ORCD_SSH[@]}" "${ORCD_HOST}" "mkdir -p '${ORCD_REPO}'"
  tar -C "${LOCAL_REPO}" -cf - "${required[@]}" |
    "${ORCD_SSH[@]}" "${ORCD_HOST}" "tar -C '${ORCD_REPO}' -xf -"
fi

echo "[stage-labels-orcd] done"
echo "[stage-labels-orcd] next command after this succeeds:"
echo "  scripts/orcd/run_s56_orcd_pilot.sh --run h200-tensor-train-smoke"
