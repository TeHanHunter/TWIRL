#!/usr/bin/env bash
# Train/apply the injected-truth peak ranker and build a ranker-selected LEO queue.
#
# This script is intended for PDO after the 20k injected peak table has been
# produced, either by the one-shot runner or by the restartable chunked runner.
set -euo pipefail

REPO_ROOT="${TWIRL_REPO_ROOT:-/pdo/users/tehan/TWIRL}"
cd "${REPO_ROOT}"

PYTHON="${PYTHON:-/sw/qlp-environment/.venv/bin/python}"
export LD_LIBRARY_PATH="/pdo/app/anaconda/anaconda2-4.4.0/lib:/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

ROOT="${TWIRL_20K_ROOT:-reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo}"
PEAK_ROOT="${TWIRL_PEAK_ROOT:-${ROOT}/peak_training}"
DEFAULT_CHUNKED="${PEAK_ROOT}/s56_20k_injection_bls_peaks_chunked.csv"
DEFAULT_ONESHOT="${PEAK_ROOT}/s56_20k_injection_bls_peaks.csv"
PEAK_TABLE="${TWIRL_PEAK_TABLE:-}"
if [[ -z "${PEAK_TABLE}" ]]; then
  if [[ -f "${DEFAULT_CHUNKED}" ]]; then
    PEAK_TABLE="${DEFAULT_CHUNKED}"
  elif [[ -f "${DEFAULT_ONESHOT}" ]]; then
    PEAK_TABLE="${DEFAULT_ONESHOT}"
  else
    echo "[ranker-review] missing peak table. Expected one of:" >&2
    echo "  ${DEFAULT_CHUNKED}" >&2
    echo "  ${DEFAULT_ONESHOT}" >&2
    exit 2
  fi
fi

RANKER_OUT_DIR="${TWIRL_RANKER_OUT_DIR:-${ROOT}/peak_ranker_pdo}"
GATE_OUT_DIR="${TWIRL_GATE_OUT_DIR:-${ROOT}/peak_training_gate_pdo}"
DEFAULT_REAL_PEAK_CSV="data_local/stage2/bls_first_pass_v2/sector_0056/candidates.csv"
DEFAULT_REAL_PEAK_PARQUET="data_local/stage2/bls_first_pass_v2/sector_0056/candidates.parquet"
REAL_PEAK_TABLE="${TWIRL_REAL_PEAK_TABLE:-}"
if [[ -z "${REAL_PEAK_TABLE}" ]]; then
  if [[ -f "${DEFAULT_REAL_PEAK_CSV}" ]]; then
    REAL_PEAK_TABLE="${DEFAULT_REAL_PEAK_CSV}"
  else
    REAL_PEAK_TABLE="${DEFAULT_REAL_PEAK_PARQUET}"
  fi
fi
SELECTED_OUT_DIR="${TWIRL_RANKER_SELECTED_OUT_DIR:-reports/stage5_validation/s56_ranker_selected_real_candidates_pdo}"
REVIEW_OUT_DIR="${TWIRL_RANKER_REVIEW_OUT_DIR:-reports/stage5_validation/s56_ranker_selected_real_review_queue_pdo}"
HLSP_ROOT="${TWIRL_HLSP_ROOT:-/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl_fs_v2_compare}"
INJECTION_H5="${INJECTION_H5:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_20k_predetrend_batman_periodradius_bright_balanced/injected_lightcurves.h5}"
STAR_CATALOG="${TWIRL_STAR_CATALOG:-data_local/catalogs/twirl_master_catalog/twirl_wd_master_catalog_v0_ticmatched.fits}"
TOP_N="${TWIRL_RANKER_TOP_N:-3}"
N_REAL="${TWIRL_REVIEW_N_REAL:-1000}"
LEO_WORKERS="${TWIRL_LEO_WORKERS:-8}"
LEO_TIMEOUT_S="${TWIRL_LEO_TIMEOUT_S:-300}"
MAX_LEO_REPORTS="${TWIRL_MAX_LEO_REPORTS:-0}"
SKIP_LEO="${TWIRL_SKIP_LEO:-0}"
WRITE_SCORED="${TWIRL_WRITE_SCORED:-0}"
INJECTED_MIN_INJECTIONS="${TWIRL_INJECTED_PEAK_MIN_INJECTIONS:-10000}"
INJECTED_MIN_CANDIDATE_ROWS="${TWIRL_INJECTED_PEAK_MIN_CANDIDATE_ROWS:-100000}"
INJECTED_MIN_APERTURES="${TWIRL_INJECTED_PEAK_MIN_APERTURES:-2}"
INJECTED_MIN_POSITIVE_RANKS="${TWIRL_INJECTED_PEAK_MIN_POSITIVE_RANKS:-10}"
REAL_MIN_ROWS="${TWIRL_REAL_PEAK_MIN_ROWS:-1000000}"
REAL_MIN_TICS="${TWIRL_REAL_PEAK_MIN_TICS:-18000}"
REAL_MIN_PEAK_RANKS="${TWIRL_REAL_PEAK_MIN_RANKS:-10}"
REAL_MIN_APERTURES="${TWIRL_REAL_PEAK_MIN_APERTURES:-2}"
REAL_MIN_ROWS_PER_TIC_MAX="${TWIRL_REAL_PEAK_MIN_ROWS_PER_TIC_MAX:-10}"

echo "[ranker-review] repo=${REPO_ROOT}"
echo "[ranker-review] python=${PYTHON}"
echo "[ranker-review] peak_table=${PEAK_TABLE}"
echo "[ranker-review] real_peak_table=${REAL_PEAK_TABLE}"
echo "[ranker-review] ranker_out=${RANKER_OUT_DIR}"
echo "[ranker-review] gate_out=${GATE_OUT_DIR}"
echo "[ranker-review] selected_out=${SELECTED_OUT_DIR}"
echo "[ranker-review] review_out=${REVIEW_OUT_DIR}"
echo "[ranker-review] top_n=${TOP_N} n_real=${N_REAL}"
mkdir -p "${RANKER_OUT_DIR}" "${SELECTED_OUT_DIR}" "${REVIEW_OUT_DIR}"

"${PYTHON}" scripts/stage5_validation/verify_injection_peak_training_table.py \
  --peak-table "${PEAK_TABLE}" \
  --out-json "${RANKER_OUT_DIR}/peak_table_verification.json" \
  --min-injections "${INJECTED_MIN_INJECTIONS}" \
  --min-candidate-rows "${INJECTED_MIN_CANDIDATE_ROWS}" \
  --min-apertures "${INJECTED_MIN_APERTURES}" \
  --min-positive-peak-ranks "${INJECTED_MIN_POSITIVE_RANKS}" \
  --require-cadence-diagnostics

TRAIN_ARGS=(
  scripts/stage5_validation/train_injection_peak_ranker.py
  --peak-table "${PEAK_TABLE}"
  --out-dir "${RANKER_OUT_DIR}"
)
if [[ "${WRITE_SCORED}" != "1" ]]; then
  TRAIN_ARGS+=(--no-scored-table)
fi
"${PYTHON}" "${TRAIN_ARGS[@]}"

"${PYTHON}" scripts/stage5_validation/summarize_injection_peak_gate.py \
  --peak-table "${PEAK_TABLE}" \
  --out-dir "${GATE_OUT_DIR}" \
  --top-k 20 \
  --ranker-summary-json "${RANKER_OUT_DIR}/summary.json"

"${PYTHON}" scripts/stage5_validation/audit_injection_bls_failure_modes.py \
  --peak-table "${PEAK_TABLE}" \
  --out-dir "${GATE_OUT_DIR}/failure_modes" \
  --top-k 20

"${PYTHON}" scripts/stage5_validation/verify_real_bls_peak_table.py \
  --peak-table "${REAL_PEAK_TABLE}" \
  --out-json "${SELECTED_OUT_DIR}/real_peak_table_verification.json" \
  --min-rows "${REAL_MIN_ROWS}" \
  --min-tics "${REAL_MIN_TICS}" \
  --min-peak-ranks "${REAL_MIN_PEAK_RANKS}" \
  --min-apertures "${REAL_MIN_APERTURES}" \
  --min-rows-per-tic-max "${REAL_MIN_ROWS_PER_TIC_MAX}"

"${PYTHON}" scripts/stage5_validation/apply_injection_peak_ranker.py \
  --model "${RANKER_OUT_DIR}/peak_ranker_model.npz" \
  --peak-table "${REAL_PEAK_TABLE}" \
  --out-dir "${SELECTED_OUT_DIR}" \
  --id-column tic \
  --top-n "${TOP_N}" \
  --no-write-scored

REVIEW_ARGS=(
  scripts/stage5_validation/build_s56_pretriage_review_queue.py
  --real-candidates "${SELECTED_OUT_DIR}/selected_ephemerides.csv"
  --real-selection ranker_selected
  --n-real "${N_REAL}"
  --n-injections 0
  --hlsp-root "${HLSP_ROOT}"
  --injection-h5 "${INJECTION_H5}"
  --out-dir "${REVIEW_OUT_DIR}"
  --apertures DET_FLUX_ADP_SML DET_FLUX_SML
  --shuffle-review-rows
  --star-atmosphere-priority H,He,mixed
  --leo-workers "${LEO_WORKERS}"
  --leo-timeout-s "${LEO_TIMEOUT_S}"
  --max-leo-reports "${MAX_LEO_REPORTS}"
  --overwrite
)
if [[ -f "${STAR_CATALOG}" ]]; then
  REVIEW_ARGS+=(--star-catalog "${STAR_CATALOG}")
fi
if [[ "${SKIP_LEO}" == "1" ]]; then
  REVIEW_ARGS+=(--skip-leo)
fi
"${PYTHON}" "${REVIEW_ARGS[@]}"

VERIFY_ARGS=(
  scripts/stage5_validation/verify_ranker_selected_real_queue.py
  --queue "${REVIEW_OUT_DIR}/review_queue.csv"
  --summary-json "${REVIEW_OUT_DIR}/summary.json"
  --out-json "${REVIEW_OUT_DIR}/verification.json"
  --min-rows "${N_REAL}"
  --expect-real "${N_REAL}"
  --expect-injected 0
)
if [[ "${SKIP_LEO}" == "1" ]]; then
  VERIFY_ARGS+=(--expect-skip-leo)
else
  VERIFY_ARGS+=(--require-reports --reports-dir "${REVIEW_OUT_DIR}/vet_reports")
fi
"${PYTHON}" "${VERIFY_ARGS[@]}"

"${PYTHON}" scripts/stage5_validation/summarize_ranker_selected_real_queue.py \
  --selected-ephemerides "${SELECTED_OUT_DIR}/selected_ephemerides.csv" \
  --review-queue "${REVIEW_OUT_DIR}/review_queue.csv" \
  --verification-json "${REVIEW_OUT_DIR}/verification.json" \
  --out-dir "${REVIEW_OUT_DIR}/ranker_selection_summary"

echo "[ranker-review] complete $(date -Is)"
