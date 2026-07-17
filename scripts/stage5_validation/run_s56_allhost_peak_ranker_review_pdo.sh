#!/usr/bin/env bash
# Run the post-BLS gate/ranker/review handoff for the all-S56-host injection product.
#
# This wrapper intentionally reuses run_s56_peak_ranker_review_pdo.sh and only
# changes paths/minimums to the all-host product. It should be launched after
# run_s56_allhost_peak_training_sharded_pdo.sh has produced and verified the
# merged all-host injected peak table.
set -euo pipefail

REPO_ROOT="${TWIRL_REPO_ROOT:-/pdo/users/tehan/TWIRL}"
cd "${REPO_ROOT}"

ROOT="${TWIRL_ALLHOST_ROOT:-reports/stage5_validation/s56_allhost_predetrend_batman_periodradius_grid_sharded_pdo}"
PEAK_TABLE="${TWIRL_PEAK_TABLE:-${ROOT}/peak_training/s56_allhost_injection_bls_peaks.csv}"
VERIFY_JSON="${PEAK_TABLE%.*}_verification.json"

if [[ ! -s "${PEAK_TABLE}" ]]; then
  echo "[allhost-ranker-review] missing merged peak table: ${PEAK_TABLE}" >&2
  echo "[allhost-ranker-review] monitor active chunks with:" >&2
  echo "  PYTHONPATH=src /sw/qlp-environment/.venv/bin/python scripts/stage5_validation/summarize_injection_peak_chunk_progress.py --chunk-root ${ROOT}/peak_training/chunks --out-dir ${ROOT}/peak_training/chunk_progress" >&2
  exit 2
fi
if [[ ! -s "${VERIFY_JSON}" ]]; then
  echo "[allhost-ranker-review] missing merged peak-table verification: ${VERIFY_JSON}" >&2
  echo "[allhost-ranker-review] rerun the all-host peak builder; it verifies after merge." >&2
  exit 2
fi

export TWIRL_20K_ROOT="${ROOT}"
export TWIRL_PEAK_ROOT="${ROOT}/peak_training"
export TWIRL_PEAK_TABLE="${PEAK_TABLE}"
export TWIRL_RANKER_OUT_DIR="${TWIRL_RANKER_OUT_DIR:-${ROOT}/peak_ranker_pdo}"
export TWIRL_GATE_OUT_DIR="${TWIRL_GATE_OUT_DIR:-${ROOT}/peak_training_gate_pdo}"
export TWIRL_RANKER_SELECTED_OUT_DIR="${TWIRL_RANKER_SELECTED_OUT_DIR:-reports/stage5_validation/s56_allhost_ranker_selected_real_candidates_pdo}"
export TWIRL_RANKER_REVIEW_OUT_DIR="${TWIRL_RANKER_REVIEW_OUT_DIR:-reports/stage5_validation/s56_allhost_ranker_selected_real_review_queue_pdo}"
export TWIRL_INJECTED_PEAK_MIN_INJECTIONS="${TWIRL_INJECTED_PEAK_MIN_INJECTIONS:-19000}"
export TWIRL_INJECTED_PEAK_MIN_CANDIDATE_ROWS="${TWIRL_INJECTED_PEAK_MIN_CANDIDATE_ROWS:-190000}"
export TWIRL_INJECTED_PEAK_MIN_APERTURES="${TWIRL_INJECTED_PEAK_MIN_APERTURES:-2}"
export TWIRL_INJECTED_PEAK_MIN_POSITIVE_RANKS="${TWIRL_INJECTED_PEAK_MIN_POSITIVE_RANKS:-10}"
export TWIRL_RANKER_TOP_N="${TWIRL_RANKER_TOP_N:-3}"
export TWIRL_REVIEW_N_REAL="${TWIRL_REVIEW_N_REAL:-1000}"

echo "[allhost-ranker-review] root=${ROOT}"
echo "[allhost-ranker-review] peak_table=${PEAK_TABLE}"
echo "[allhost-ranker-review] ranker_out=${TWIRL_RANKER_OUT_DIR}"
echo "[allhost-ranker-review] gate_out=${TWIRL_GATE_OUT_DIR}"
echo "[allhost-ranker-review] selected_out=${TWIRL_RANKER_SELECTED_OUT_DIR}"
echo "[allhost-ranker-review] review_out=${TWIRL_RANKER_REVIEW_OUT_DIR}"

bash scripts/stage5_validation/run_s56_peak_ranker_review_pdo.sh
