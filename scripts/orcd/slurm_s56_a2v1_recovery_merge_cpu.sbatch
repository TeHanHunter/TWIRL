#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=s56-a2rec-merge
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=192G
#SBATCH -t 1-00:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/s56_a2v1_recovery/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/s56_a2v1_recovery/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
OUT_ROOT="${TWIRL_RECOVERY_ROOT:-${REPO}/reports/stage5_validation/s56_A2v1_teacher_v1_recovery}"
FULL="${OUT_ROOT}/full"

mkdir -p "${FULL}/teacher_inputs"
cd "${REPO}"
export PYTHONPATH=src HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

"${PYTHON}" scripts/stage5_validation/merge_s56_a2v1_recovery_teacher_inputs.py \
  --candidate-shards "${FULL}"/teacher_input_shards/candidates_*.csv \
  --native-shards "${FULL}"/teacher_input_shards/native_*.h5 \
  --out-candidates "${FULL}/teacher_inputs/candidates_top5.csv" \
  --out-native-h5 "${FULL}/teacher_inputs/native_top5.h5" \
  --maximum-injections 20000

echo "[s56-a2rec-merge] complete $(date -Is)"
