#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=s56-t2-freeze
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_teacher_v2}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
OUT_ROOT="${TWIRL_TEACHER_V2_ROOT:-${REPO}/reports/stage5_validation/s56_A2v1_teacher_v2}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2 "${OUT_ROOT}/frozen"
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

"${PYTHON}" scripts/stage5_validation/freeze_s56_teacher_v2_selection.py \
  --training-root "${OUT_ROOT}/training" \
  --real-score-dir "${OUT_ROOT}/real_scores" \
  --out-dir "${OUT_ROOT}/frozen" \
  --max-real-tic-fraction 0.05
echo "[s56-t2-freeze] complete=$(date -Is)"
