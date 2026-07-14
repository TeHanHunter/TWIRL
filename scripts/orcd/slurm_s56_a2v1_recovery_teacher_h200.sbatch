#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=s56-a2rec-teacher
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=192G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 12:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/s56_a2v1_recovery/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/s56_a2v1_recovery/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56-torch/bin/python}"
SECTOR="${TWIRL_RECOVERY_SECTOR:-56}"
SECTOR_TAG="s${SECTOR}"
LOG_ROOT="${TWIRL_RECOVERY_LOG_ROOT:-/orcd/data/mki_aryeh/001/twirl/logs/${SECTOR_TAG}_a2v1_recovery}"
OUT_ROOT="${TWIRL_RECOVERY_ROOT:-${REPO}/reports/stage5_validation/${SECTOR_TAG}_A2v1_teacher_v1_recovery}"
CHECKPOINT_ROOT="${TWIRL_TEACHER_CHECKPOINT_ROOT:-${REPO}/reports/stage5_validation/s56_label_adjudication_real343/retrained/s56_harmonic_cnn_v1/shape_plus_periodogram_bls}"
MODE="${TWIRL_RECOVERY_MODE:-full}"

if [[ "${MODE}" == "smoke" ]]; then
  INPUT_ROOT="${OUT_ROOT}/smoke100/teacher_inputs"
  CANDIDATES="${INPUT_ROOT}/candidates_00.csv"
  NATIVE_H5="${INPUT_ROOT}/native_00.h5"
  OUT_SCORES="${OUT_ROOT}/smoke100/teacher_scores_128.parquet"
  LIMIT_ARGS=(--limit 128)
else
  INPUT_ROOT="${OUT_ROOT}/full/teacher_inputs"
  CANDIDATES="${INPUT_ROOT}/candidates_top5.csv"
  NATIVE_H5="${INPUT_ROOT}/native_top5.h5"
  OUT_SCORES="${OUT_ROOT}/full/teacher_scores.parquet"
  LIMIT_ARGS=()
fi

mkdir -p "${LOG_ROOT}"
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE

"${PYTHON}" scripts/stage5_validation/score_s56_a2v1_teacher_recovery.py \
  --candidates "${CANDIDATES}" \
  --native-h5 "${NATIVE_H5}" \
  --checkpoints \
    "${CHECKPOINT_ROOT}/fold_0/teacher.pt" \
    "${CHECKPOINT_ROOT}/fold_1/teacher.pt" \
    "${CHECKPOINT_ROOT}/fold_2/teacher.pt" \
    "${CHECKPOINT_ROOT}/fold_3/teacher.pt" \
    "${CHECKPOINT_ROOT}/fold_4/teacher.pt" \
  --out-scores "${OUT_SCORES}" \
  --batch-size 32 \
  --workers 4 \
  "${LIMIT_ARGS[@]}"

echo "[${SECTOR_TAG}-a2rec-teacher] mode=${MODE} complete $(date -Is)"
