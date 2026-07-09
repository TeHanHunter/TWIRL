#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-rec50-cnn
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 00:45:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56-torch/bin/python}"
ROOT="${TWIRL_RECOVERY50_ROOT:-reports/stage5_validation/s56_recovery50_teacher_queue_2k}"
TENSOR_DIR="${TWIRL_RECOVERY50_TENSOR_DIR:-${ROOT}/cnn_tensors_adp_only}"
TRAINING_TABLE="${TWIRL_RECOVERY50_TRAINING_TABLE:-${ROOT}/human_training_table_adp_only/human_vetting_training_table.csv}"
CNN_OUT_DIR="${TWIRL_RECOVERY50_CNN_OUT_DIR:-${ROOT}/cnn_teacher_adp_only}"
MIN_CLASS_COUNT="${TWIRL_RECOVERY50_MIN_CLASS_COUNT:-40}"
EPOCHS="${TWIRL_RECOVERY50_CNN_EPOCHS:-80}"
BATCH_SIZE="${TWIRL_RECOVERY50_CNN_BATCH_SIZE:-32}"
ALLOW_CPU="${TWIRL_RECOVERY50_CNN_ALLOW_CPU:-0}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs
cd "${REPO}"

export PYTHONPATH=src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[rec50-cnn-h200] host=$(hostname)"
echo "[rec50-cnn-h200] date=$(date -Is)"
echo "[rec50-cnn-h200] repo=${REPO}"
echo "[rec50-cnn-h200] python=${PYTHON}"
echo "[rec50-cnn-h200] root=${ROOT}"
echo "[rec50-cnn-h200] tensor_dir=${TENSOR_DIR}"
echo "[rec50-cnn-h200] training_table=${TRAINING_TABLE}"
echo "[rec50-cnn-h200] out_dir=${CNN_OUT_DIR}"
echo "[rec50-cnn-h200] epochs=${EPOCHS} batch_size=${BATCH_SIZE}"

module load cuda/13.1.0 || module load cuda/12.9.1 || true
nvidia-smi

train_args=(
  --tensor-npz "${TENSOR_DIR}/recovery50_cnn_tensors.npz"
  --tensor-rows "${TENSOR_DIR}/recovery50_cnn_tensor_rows.csv"
  --training-table "${TRAINING_TABLE}"
  --out-dir "${CNN_OUT_DIR}"
  --min-class-count "${MIN_CLASS_COUNT}"
  --epochs "${EPOCHS}"
  --batch-size "${BATCH_SIZE}"
)
if [[ "${ALLOW_CPU}" == "1" ]]; then
  train_args+=(--allow-cpu)
fi

"${PYTHON}" scripts/stage5_validation/train_s56_recovery50_cnn_teacher.py "${train_args[@]}"

echo "[rec50-cnn-h200] complete $(date -Is)"
