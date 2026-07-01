#!/usr/bin/env bash
# Render WD-tuned LEO sheets for the ranker-selected real S56 ephemerides.
#
# This is the second stage after run_s56_peak_ranker_review_pdo.sh has produced
# and verified a cheap skip-LEO queue. It reuses the selected ephemerides and
# writes a sibling LEO-backed queue so the preflight artifact stays untouched.
set -euo pipefail

REPO_ROOT="${TWIRL_REPO_ROOT:-/pdo/users/tehan/TWIRL}"
cd "${REPO_ROOT}"

PYTHON="${PYTHON:-/sw/qlp-environment/.venv/bin/python}"
export LD_LIBRARY_PATH="/pdo/app/anaconda/anaconda2-4.4.0/lib:/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="/pdo/users/tehan/LEO-Vetter-twirl:src"
export MPLCONFIGDIR=/pdo/users/tehan/.cache/matplotlib
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

SELECTED_EPHEMERIDES="${TWIRL_SELECTED_EPHEMERIDES:-reports/stage5_validation/s56_ranker_selected_real_candidates_pdo/selected_ephemerides.csv}"
OUT_DIR="${TWIRL_RANKER_LEO_OUT_DIR:-reports/stage5_validation/s56_ranker_selected_real_leo_queue_pdo}"
HLSP_ROOT="${TWIRL_HLSP_ROOT:-/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl_fs_v2_compare}"
INJECTION_H5="${INJECTION_H5:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_20k_predetrend_batman_periodradius_bright_balanced/injected_lightcurves.h5}"
STAR_CATALOG="${TWIRL_STAR_CATALOG:-data_local/catalogs/twirl_master_catalog/twirl_wd_master_catalog_v0_ticmatched.fits}"
STAR_PRIORITY="${TWIRL_STAR_PRIORITY:-H,He,mixed}"
N_REAL="${TWIRL_REVIEW_N_REAL:-1000}"
LEO_WORKERS="${TWIRL_LEO_WORKERS:-8}"
LEO_TIMEOUT_S="${TWIRL_LEO_TIMEOUT_S:-300}"
MAX_LEO_REPORTS="${TWIRL_MAX_LEO_REPORTS:-0}"
RANDOM_STATE="${TWIRL_REVIEW_RANDOM_STATE:-5606}"

if [[ ! -f "${SELECTED_EPHEMERIDES}" ]]; then
  echo "[ranker-real-leo] missing selected ephemerides: ${SELECTED_EPHEMERIDES}" >&2
  echo "[ranker-real-leo] run run_s56_peak_ranker_review_pdo.sh after the peak table is ready first" >&2
  exit 2
fi

echo "[ranker-real-leo] repo=${REPO_ROOT}"
echo "[ranker-real-leo] selected=${SELECTED_EPHEMERIDES}"
echo "[ranker-real-leo] out=${OUT_DIR}"
echo "[ranker-real-leo] n_real=${N_REAL} leo_workers=${LEO_WORKERS}"

REVIEW_ARGS=(
  scripts/stage5_validation/build_s56_pretriage_review_queue.py
  --real-candidates "${SELECTED_EPHEMERIDES}"
  --real-selection ranker_selected
  --n-real "${N_REAL}"
  --n-injections 0
  --hlsp-root "${HLSP_ROOT}"
  --injection-h5 "${INJECTION_H5}"
  --out-dir "${OUT_DIR}"
  --apertures DET_FLUX_ADP_SML DET_FLUX_SML
  --star-atmosphere-priority "${STAR_PRIORITY}"
  --shuffle-review-rows
  --random-state "${RANDOM_STATE}"
  --leo-workers "${LEO_WORKERS}"
  --leo-timeout-s "${LEO_TIMEOUT_S}"
  --max-leo-reports "${MAX_LEO_REPORTS}"
  --overwrite
)
if [[ -f "${STAR_CATALOG}" ]]; then
  REVIEW_ARGS+=(--star-catalog "${STAR_CATALOG}")
fi
"${PYTHON}" "${REVIEW_ARGS[@]}"

"${PYTHON}" scripts/stage5_validation/verify_ranker_selected_real_queue.py \
  --queue "${OUT_DIR}/review_queue.csv" \
  --summary-json "${OUT_DIR}/summary.json" \
  --reports-dir "${OUT_DIR}/vet_reports" \
  --out-json "${OUT_DIR}/verification.json" \
  --min-rows "${N_REAL}" \
  --expect-real "${N_REAL}" \
  --expect-injected 0 \
  --require-reports

"${PYTHON}" scripts/stage5_validation/summarize_ranker_selected_real_queue.py \
  --selected-ephemerides "${SELECTED_EPHEMERIDES}" \
  --review-queue "${OUT_DIR}/review_queue.csv" \
  --verification-json "${OUT_DIR}/verification.json" \
  --out-dir "${OUT_DIR}/ranker_selection_summary"

echo "[ranker-real-leo] complete $(date -Is)"
