#!/usr/bin/env bash
# Render WD-tuned LEO sheets for the all-host ranker-selected real S56 queue.
#
# This is intentionally a second-stage command. It refuses to render LEO sheets
# until the cheap skip-LEO all-host queue has passed verification.
set -euo pipefail

REPO_ROOT="${TWIRL_REPO_ROOT:-/pdo/users/tehan/TWIRL}"
cd "${REPO_ROOT}"

PYTHON="${PYTHON:-/sw/qlp-environment/.venv/bin/python}"
export LD_LIBRARY_PATH="/pdo/app/anaconda/anaconda2-4.4.0/lib:/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${PYTHONPATH:-src}"

SKIP_QUEUE_VERIFY="${TWIRL_SKIP_QUEUE_VERIFY:-reports/stage5_validation/s56_allhost_ranker_selected_real_review_queue_pdo/verification.json}"
SELECTED_EPHEMERIDES="${TWIRL_SELECTED_EPHEMERIDES:-reports/stage5_validation/s56_allhost_ranker_selected_real_candidates_pdo/selected_ephemerides.csv}"
LEO_OUT_DIR="${TWIRL_RANKER_LEO_OUT_DIR:-reports/stage5_validation/s56_allhost_ranker_selected_real_leo_queue_pdo}"

if [[ ! -s "${SKIP_QUEUE_VERIFY}" ]]; then
  echo "[allhost-ranker-real-leo] missing skip-LEO verifier: ${SKIP_QUEUE_VERIFY}" >&2
  echo "[allhost-ranker-real-leo] wait for run_s56_allhost_peak_ranker_review_pdo.sh to pass first" >&2
  exit 2
fi
if ! "${PYTHON}" - "${SKIP_QUEUE_VERIFY}" <<'PY'
import json
import sys
from pathlib import Path

payload = json.loads(Path(sys.argv[1]).read_text())
raise SystemExit(0 if payload.get("passed") else 1)
PY
then
  echo "[allhost-ranker-real-leo] skip-LEO verifier did not pass: ${SKIP_QUEUE_VERIFY}" >&2
  exit 2
fi
if [[ ! -s "${SELECTED_EPHEMERIDES}" ]]; then
  echo "[allhost-ranker-real-leo] missing selected ephemerides: ${SELECTED_EPHEMERIDES}" >&2
  exit 2
fi

echo "[allhost-ranker-real-leo] selected=${SELECTED_EPHEMERIDES}"
echo "[allhost-ranker-real-leo] skip_queue_verify=${SKIP_QUEUE_VERIFY}"
echo "[allhost-ranker-real-leo] leo_out=${LEO_OUT_DIR}"

TWIRL_SELECTED_EPHEMERIDES="${SELECTED_EPHEMERIDES}" \
TWIRL_RANKER_LEO_OUT_DIR="${LEO_OUT_DIR}" \
  bash scripts/stage5_validation/run_s56_ranker_selected_real_leo_pdo.sh
