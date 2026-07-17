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
PYTHON="${TWIRL_ORCD_PYTHON:-python3}"
MODEL="${TWIRL_RANKER_MODEL:-reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_ranker_orcd/peak_ranker_model.npz}"
PEAK_TABLE="${TWIRL_REAL_PEAK_TABLE:-data_local/stage2/bls_first_pass_v2/sector_0056/candidates.csv}"
OUT_DIR="${TWIRL_RANKER_SELECTED_OUT_DIR:-reports/stage5_validation/s56_ranker_selected_real_candidates_orcd}"
TOP_N="${TWIRL_RANKER_TOP_N:-3}"

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

"${PYTHON}" scripts/stage5_validation/apply_injection_peak_ranker.py \
  --model "${MODEL}" \
  --peak-table "${PEAK_TABLE}" \
  --out-dir "${OUT_DIR}" \
  --id-column tic \
  --top-n "${TOP_N}" \
  --no-write-scored

echo "[orcd-ranker-apply] complete $(date -Is)"
