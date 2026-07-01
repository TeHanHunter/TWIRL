#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-s56-tensor-train
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 00:20:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56-torch/bin/python}"
TENSOR_NPZ="${TWIRL_TENSOR_NPZ:-reports/stage5_validation/s56_candidate_tensor_smoke_h200_torch/candidate_tensor_smoke.npz}"
TENSOR_ROWS="${TWIRL_TENSOR_ROWS:-reports/stage5_validation/s56_candidate_tensor_smoke_h200_torch/candidate_tensor_smoke_rows.csv}"
LABELS="${TWIRL_TENSOR_LABELS:-reports/stage5_validation/s56_allhost_ranker_selected_real_leo_queue_pdo/human_training_table/teacher_labeled_rows.csv}"
OUT_DIR="${TWIRL_TENSOR_TRAIN_OUT_DIR:-reports/stage5_validation/s56_candidate_tensor_train_smoke_h200}"
EPOCHS="${TWIRL_TENSOR_TRAIN_EPOCHS:-12}"
BATCH_SIZE="${TWIRL_TENSOR_TRAIN_BATCH_SIZE:-32}"
MIN_TEACHER_ROWS="${TWIRL_TENSOR_MIN_TEACHER_ROWS:-20}"
SYNTHETIC_LABEL_SMOKE="${TWIRL_TENSOR_SYNTHETIC_LABEL_SMOKE:-0}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs

echo "[tensor-train-h200] host=$(hostname)"
echo "[tensor-train-h200] date=$(date -Is)"
echo "[tensor-train-h200] repo=${REPO}"
echo "[tensor-train-h200] python=${PYTHON}"
echo "[tensor-train-h200] tensor_npz=${TENSOR_NPZ}"
echo "[tensor-train-h200] tensor_rows=${TENSOR_ROWS}"
echo "[tensor-train-h200] labels=${LABELS}"
echo "[tensor-train-h200] out_dir=${OUT_DIR}"
echo "[tensor-train-h200] synthetic_label_smoke=${SYNTHETIC_LABEL_SMOKE}"

module load cuda/13.1.0 || module load cuda/12.9.1 || true
nvidia-smi

cd "${REPO}"
ARGS=(
  scripts/stage5_validation/train_candidate_tensor_classifier_smoke.py
  --tensor-npz "${TENSOR_NPZ}"
  --tensor-rows "${TENSOR_ROWS}"
  --labels "${LABELS}"
  --out-dir "${OUT_DIR}"
  --epochs "${EPOCHS}"
  --batch-size "${BATCH_SIZE}"
  --min-teacher-rows "${MIN_TEACHER_ROWS}"
  --require-cuda
)
if [[ "${SYNTHETIC_LABEL_SMOKE}" == "1" ]]; then
  ARGS+=(--synthetic-label-smoke)
fi

"${PYTHON}" "${ARGS[@]}"

echo "[tensor-train-h200] complete $(date -Is)"
