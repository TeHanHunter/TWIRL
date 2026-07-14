#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=s56-t2-profile
#SBATCH --array=0-5%1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=192G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 48:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2/%x-%A_%a.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2/%x-%A_%a.err

set -euo pipefail

profiles=(
  metadata_only
  single_period_native_fold
  seven_harmonic_shape
  shape_plus_raw_chronology
  shape_plus_periodogram_bls
  full_combined
)
PROFILE="${profiles[${SLURM_ARRAY_TASK_ID}]}"
REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_teacher_v2}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56-torch/bin/python}"
OUT_ROOT="${TWIRL_TEACHER_V2_ROOT:-${REPO}/reports/stage5_validation/s56_A2v1_teacher_v2}"
TRAIN_ROOT="${OUT_ROOT}/training"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2 \
  "${TRAIN_ROOT}" "${OUT_ROOT}/real_scores" "${OUT_ROOT}/runtime"
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
  --out-dir "${TRAIN_ROOT}" \
  --profiles "${PROFILE}" \
  --epochs 100 \
  --patience 12 \
  --batch-size 32 \
  --workers 8
train_end_epoch=$(date +%s)

checkpoint_args=()
for fold in 0 1 2 3 4; do
  checkpoint_args+=(--checkpoint "${TRAIN_ROOT}/${PROFILE}/fold_${fold}/teacher_v2.pt")
done
"${PYTHON}" scripts/stage5_validation/score_teacher_v2_candidates.py \
  --candidates "${OUT_ROOT}/human_transfer/s56_a2v1_real_candidates_top5.parquet" \
  --native-h5 "${OUT_ROOT}/real_native/s56_a2v1_real_native.h5" \
  "${checkpoint_args[@]}" \
  --profile "${PROFILE}" \
  --out-path "${OUT_ROOT}/real_scores/${PROFILE}.parquet" \
  --batch-size 32 \
  --workers 8
end_epoch=$(date +%s)
printf 'profile=%s\nstart_epoch=%s\ntrain_end_epoch=%s\nend_epoch=%s\ntrain_elapsed_s=%s\ntotal_elapsed_s=%s\n' \
  "${PROFILE}" "${start_epoch}" "${train_end_epoch}" "${end_epoch}" \
  "$((train_end_epoch-start_epoch))" "$((end_epoch-start_epoch))" \
  > "${OUT_ROOT}/runtime/${PROFILE}.txt"
echo "[s56-t2-profile] profile=${PROFILE} complete=$(date -Is)"
