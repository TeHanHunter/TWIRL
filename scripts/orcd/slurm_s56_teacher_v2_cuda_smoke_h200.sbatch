#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=s56-t2-smoke
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --mem=128G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 08:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_teacher_v2}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56-torch/bin/python}"
OUT_ROOT="${TWIRL_TEACHER_V2_ROOT:-${REPO}/reports/stage5_validation/s56_A2v1_teacher_v2}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2 "${OUT_ROOT}/cuda_smoke"
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
module load cuda/13.1.0 || module load cuda/12.9.1 || true
nvidia-smi

start_epoch=$(date +%s)
"${PYTHON}" scripts/stage5_validation/train_s56_teacher_v2.py \
  --training-table "${OUT_ROOT}/roles/training_rows.parquet" \
  --out-dir "${OUT_ROOT}/cuda_smoke" \
  --profiles full_combined \
  --epochs 3 \
  --patience 2 \
  --batch-size 32 \
  --workers 4 \
  --max-development-rows 128
end_epoch=$(date +%s)
printf 'start_epoch=%s\nend_epoch=%s\nelapsed_s=%s\n' \
  "${start_epoch}" "${end_epoch}" "$((end_epoch-start_epoch))" \
  > "${OUT_ROOT}/cuda_smoke/runtime.txt"
echo "[s56-t2-smoke] complete=$(date -Is)"
