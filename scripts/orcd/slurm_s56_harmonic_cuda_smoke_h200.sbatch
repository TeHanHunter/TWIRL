#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-harmonic-smoke
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 01:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56-torch/bin/python}"
OUT="${TWIRL_HARMONIC_SMOKE_ROOT:-reports/stage5_validation/s56_harmonic_cnn_v1/cuda_smoke_${SLURM_JOB_ID}}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs
cd "${REPO}"
export PYTHONPATH="${REPO}:${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

module load cuda/13.1.0 || module load cuda/12.9.1 || true
nvidia-smi

"${PYTHON}" scripts/stage5_validation/smoke_s56_harmonic_cnn_v1.py \
  --out-dir "${OUT}" \
  --profiles metadata_only,full_combined \
  --epochs 1 \
  --pretrain-epochs 1 \
  --batch-size 16 \
  --workers 2

printf '[harmonic-cuda-smoke] complete: %s\n' "${OUT}"
