#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-harmonic-figure
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56-torch/bin/python}"
ROOT="${TWIRL_HARMONIC_TRAINING_ROOT:-reports/stage5_validation/s56_label_adjudication_real343/retrained/s56_harmonic_cnn_v1}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs "/tmp/twirl-matplotlib-${SLURM_JOB_ID}"
cd "${REPO}"
export PYTHONPATH=src
export MPLBACKEND=Agg
export MPLCONFIGDIR="/tmp/twirl-matplotlib-${SLURM_JOB_ID}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

"${PYTHON}" scripts/stage5_validation/plot_s56_harmonic_cnn_performance.py \
  --root "${ROOT}" \
  --out-dir "${ROOT}/performance_figure"

echo "[harmonic-figure] complete $(date -Is): ${ROOT}/performance_figure"
