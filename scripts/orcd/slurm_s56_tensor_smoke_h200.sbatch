#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-s56-tensor-smoke
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 00:20:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
LC_H5="${TWIRL_LC_H5:-data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_lc_export_pdo.h5}"
CANDIDATE_TABLE="${TWIRL_CANDIDATE_TABLE:-reports/stage5_validation/s56_allhost_ranker_selected_real_leo_queue_pdo/review_queue.csv}"
OUT_DIR="${TWIRL_TENSOR_OUT_DIR:-reports/stage5_validation/s56_candidate_tensor_smoke_h200}"
MAX_CANDIDATES="${TWIRL_TENSOR_MAX_CANDIDATES:-128}"
N_POINTS="${TWIRL_TENSOR_N_POINTS:-257}"
WINDOW_DURATIONS="${TWIRL_TENSOR_WINDOW_DURATIONS:-8}"
REQUIRE_LABELS="${TWIRL_REQUIRE_LABELS:-0}"
REQUIRE_TORCH="${TWIRL_REQUIRE_TORCH:-0}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs

echo "[tensor-smoke-h200] host=$(hostname)"
echo "[tensor-smoke-h200] date=$(date -Is)"
echo "[tensor-smoke-h200] repo=${REPO}"
echo "[tensor-smoke-h200] python=${PYTHON}"
echo "[tensor-smoke-h200] lc_h5=${LC_H5}"
echo "[tensor-smoke-h200] candidate_table=${CANDIDATE_TABLE}"
echo "[tensor-smoke-h200] out_dir=${OUT_DIR}"
echo "[tensor-smoke-h200] max_candidates=${MAX_CANDIDATES}"

module load cuda/13.1.0 || module load cuda/12.9.1 || true
nvidia-smi

cd "${REPO}"
ARGS=(
  scripts/stage5_validation/build_s56_candidate_tensor_smoke.py
  --lc-h5 "${LC_H5}"
  --candidate-table "${CANDIDATE_TABLE}"
  --out-dir "${OUT_DIR}"
  --apertures DET_FLUX_SML DET_FLUX DET_FLUX_LAG
  --max-candidates "${MAX_CANDIDATES}"
  --n-points "${N_POINTS}"
  --window-durations "${WINDOW_DURATIONS}"
  --progress-every 25
  --torch-smoke
)
if [[ "${REQUIRE_LABELS}" == "1" ]]; then
  ARGS+=(--require-labels)
fi
if [[ "${REQUIRE_TORCH}" == "1" ]]; then
  ARGS+=(--require-torch)
fi

"${PYTHON}" "${ARGS[@]}"

echo "[tensor-smoke-h200] complete $(date -Is)"
