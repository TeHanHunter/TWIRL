#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-s56-ranker-apply
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
MODEL="${TWIRL_RANKER_MODEL:-reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_ranker_orcd/peak_ranker_model.npz}"
PEAK_TABLE="${TWIRL_REAL_PEAK_TABLE:-data_local/stage2/bls_first_pass_v2/sector_0056/candidates.csv}"
OUT_DIR="${TWIRL_RANKER_SELECTED_OUT_DIR:-reports/stage5_validation/s56_ranker_selected_real_candidates_orcd}"
TOP_N="${TWIRL_RANKER_TOP_N:-3}"
REAL_MIN_ROWS="${TWIRL_REAL_PEAK_MIN_ROWS:-1000000}"
REAL_MIN_TICS="${TWIRL_REAL_PEAK_MIN_TICS:-18000}"
REAL_MIN_PEAK_RANKS="${TWIRL_REAL_PEAK_MIN_RANKS:-10}"
REAL_MIN_APERTURES="${TWIRL_REAL_PEAK_MIN_APERTURES:-2}"
REAL_MIN_ROWS_PER_TIC_MAX="${TWIRL_REAL_PEAK_MIN_ROWS_PER_TIC_MAX:-10}"
REAL_REQUIRE_BLS_PROVENANCE="${TWIRL_REAL_REQUIRE_BLS_PROVENANCE:-0}"
SELECTED_MIN_ROWS="${TWIRL_SELECTED_MIN_ROWS:-1000}"
SELECTED_MIN_TARGETS="${TWIRL_SELECTED_MIN_TARGETS:-300}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs
cd "${REPO}"

export PYTHONPATH=src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[orcd-ranker-apply] host=$(hostname)"
echo "[orcd-ranker-apply] repo=${REPO}"
echo "[orcd-ranker-apply] python=${PYTHON}"
echo "[orcd-ranker-apply] model=${MODEL}"
echo "[orcd-ranker-apply] peak_table=${PEAK_TABLE}"
echo "[orcd-ranker-apply] out_dir=${OUT_DIR}"
echo "[orcd-ranker-apply] top_n=${TOP_N}"

mkdir -p "${OUT_DIR}"

REAL_VERIFY_ARGS=(
  scripts/stage5_validation/verify_real_bls_peak_table.py
  --peak-table "${PEAK_TABLE}"
  --out-json "${OUT_DIR}/real_peak_table_verification.json"
  --min-rows "${REAL_MIN_ROWS}"
  --min-tics "${REAL_MIN_TICS}"
  --min-peak-ranks "${REAL_MIN_PEAK_RANKS}"
  --min-apertures "${REAL_MIN_APERTURES}"
  --min-rows-per-tic-max "${REAL_MIN_ROWS_PER_TIC_MAX}"
)
if [[ "${REAL_REQUIRE_BLS_PROVENANCE}" == "1" ]]; then
  REAL_VERIFY_ARGS+=(--require-bls-provenance)
fi
"${PYTHON}" "${REAL_VERIFY_ARGS[@]}"

"${PYTHON}" scripts/stage5_validation/apply_injection_peak_ranker.py \
  --model "${MODEL}" \
  --peak-table "${PEAK_TABLE}" \
  --out-dir "${OUT_DIR}" \
  --id-column tic \
  --top-n "${TOP_N}" \
  --no-write-scored

"${PYTHON}" scripts/stage5_validation/verify_ranker_selected_ephemerides.py \
  --selected-ephemerides "${OUT_DIR}/selected_ephemerides.csv" \
  --out-json "${OUT_DIR}/selected_ephemerides_verification.json" \
  --id-column tic \
  --top-n "${TOP_N}" \
  --min-rows "${SELECTED_MIN_ROWS}" \
  --min-targets "${SELECTED_MIN_TARGETS}"

"${PYTHON}" scripts/stage5_validation/summarize_ranker_selected_real_queue.py \
  --selected-ephemerides "${OUT_DIR}/selected_ephemerides.csv" \
  --review-queue "${OUT_DIR}/missing_review_queue.csv" \
  --verification-json "${OUT_DIR}/selected_ephemerides_verification.json" \
  --out-dir "${OUT_DIR}/ranker_selection_summary"

echo "[orcd-ranker-apply] complete $(date -Is)"
