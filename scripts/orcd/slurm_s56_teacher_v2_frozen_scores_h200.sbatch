#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=s56-t2-frozen-score
#SBATCH --array=0-8%1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=192G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 24:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2/%x-%A_%a.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2/%x-%A_%a.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_teacher_v2}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56-torch/bin/python}"
OUT_ROOT="${TWIRL_TEACHER_V2_ROOT:-${REPO}/reports/stage5_validation/s56_A2v1_teacher_v2}"
S56_RECOVERY="${TWIRL_S56_RECOVERY_ROOT:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_s56_a2v1_recovery_2b4f57bb/reports/stage5_validation/s56_A2v1_teacher_v1_recovery}"
S57_RECOVERY="${TWIRL_S57_RECOVERY_ROOT:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_a2v1_recovery_1090f011/reports/stage5_validation/s57_A2v1_teacher_v1_recovery}"
V1_ROOT="${TWIRL_TEACHER_V1_ROOT:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL/reports/stage5_validation/s56_label_adjudication_real343/retrained/s56_harmonic_cnn_v1/shape_plus_periodogram_bls}"
FROZEN="${OUT_ROOT}/frozen/frozen_selection.json"
TASK="${SLURM_ARRAY_TASK_ID}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2 "${OUT_ROOT}/frozen_scores"
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
module load cuda/13.1.0 || module load cuda/12.9.1 || true
nvidia-smi

SELECTED_PROFILE=$("${PYTHON}" -c \
  'import json,sys; print(json.load(open(sys.argv[1]))["selected_profile"])' "${FROZEN}")

v2_score() {
  local profile="$1" candidates="$2" native_h5="$3" manifest="$4" output="$5"
  local checkpoint_args=()
  for fold in 0 1 2 3 4; do
    checkpoint_args+=(--checkpoint "${OUT_ROOT}/training/${profile}/fold_${fold}/teacher_v2.pt")
  done
  local filter_args=()
  if [[ -n "${manifest}" ]]; then
    filter_args+=(--filter-manifest "${manifest}")
  fi
  "${PYTHON}" scripts/stage5_validation/score_teacher_v2_candidates.py \
    --candidates "${candidates}" \
    --native-h5 "${native_h5}" \
    "${checkpoint_args[@]}" \
    --profile "${profile}" \
    --out-path "${output}" \
    --batch-size 32 \
    --workers 8 \
    "${filter_args[@]}"
}

v1_score() {
  local candidates="$1" native_h5="$2" manifest="$3" output="$4"
  local filter_args=()
  if [[ -n "${manifest}" ]]; then
    filter_args+=(--filter-manifest "${manifest}")
  fi
  "${PYTHON}" scripts/stage5_validation/score_s56_a2v1_teacher_recovery.py \
    --candidates "${candidates}" \
    --native-h5 "${native_h5}" \
    --checkpoints \
      "${V1_ROOT}/fold_0/teacher.pt" \
      "${V1_ROOT}/fold_1/teacher.pt" \
      "${V1_ROOT}/fold_2/teacher.pt" \
      "${V1_ROOT}/fold_3/teacher.pt" \
      "${V1_ROOT}/fold_4/teacher.pt" \
    --out-scores "${output}" \
    --batch-size 32 \
    --workers 8 \
    "${filter_args[@]}"
}

S56_CANDIDATES="${S56_RECOVERY}/full/teacher_inputs/candidates_top5.csv"
S56_NATIVE="${S56_RECOVERY}/full/teacher_inputs/native_top5.h5"
S56_HOLDOUT="${OUT_ROOT}/roles/s56_injection_locked_holdout_manifest.parquet"
S57_CANDIDATES="${S57_RECOVERY}/full/teacher_inputs/candidates_top5.csv"
S57_NATIVE="${S57_RECOVERY}/full/teacher_inputs/native_top5.h5"
S57_MANIFEST="${OUT_ROOT}/roles/s57_external_role_manifest.parquet"
HUMAN_CANDIDATES="${OUT_ROOT}/roles/s56_human_role_manifest.parquet"
REAL_NATIVE="${OUT_ROOT}/real_native/s56_a2v1_real_native.h5"

case "${TASK}" in
  0)
    v2_score "${SELECTED_PROFILE}" "${S56_CANDIDATES}" "${S56_NATIVE}" \
      "${S56_HOLDOUT}" "${OUT_ROOT}/frozen_scores/s56_holdout_teacher_v2.parquet"
    ;;
  1)
    v2_score metadata_only "${S56_CANDIDATES}" "${S56_NATIVE}" \
      "${S56_HOLDOUT}" "${OUT_ROOT}/frozen_scores/s56_holdout_metadata_only.parquet"
    ;;
  2)
    v1_score "${S56_CANDIDATES}" "${S56_NATIVE}" "${S56_HOLDOUT}" \
      "${OUT_ROOT}/frozen_scores/s56_holdout_teacher_v1.parquet"
    ;;
  3)
    v2_score "${SELECTED_PROFILE}" "${S57_CANDIDATES}" "${S57_NATIVE}" \
      "${S57_MANIFEST}" "${OUT_ROOT}/frozen_scores/s57_external_teacher_v2.parquet"
    ;;
  4)
    v2_score metadata_only "${S57_CANDIDATES}" "${S57_NATIVE}" \
      "${S57_MANIFEST}" "${OUT_ROOT}/frozen_scores/s57_external_metadata_only.parquet"
    ;;
  5)
    v1_score "${S57_CANDIDATES}" "${S57_NATIVE}" "${S57_MANIFEST}" \
      "${OUT_ROOT}/frozen_scores/s57_external_teacher_v1.parquet"
    ;;
  6)
    v1_score "${OUT_ROOT}/human_transfer/s56_a2v1_real_candidates_top5.parquet" \
      "${REAL_NATIVE}" "" "${OUT_ROOT}/frozen_scores/s56_real_teacher_v1.parquet"
    "${PYTHON}" scripts/stage5_validation/freeze_real_tic_workload_threshold.py \
      --scores "${OUT_ROOT}/frozen_scores/s56_real_teacher_v1.parquet" \
      --score-column p_preserve \
      --max-fraction 0.05 \
      --model-name teacher_v1_shape_plus_periodogram_bls \
      --out-json "${OUT_ROOT}/frozen/teacher_v1_workload_threshold.json"
    ;;
  7)
    v2_score "${SELECTED_PROFILE}" "${HUMAN_CANDIDATES}" "${REAL_NATIVE}" "" \
      "${OUT_ROOT}/frozen_scores/s56_human_teacher_v2.parquet"
    ;;
  8)
    v1_score "${HUMAN_CANDIDATES}" "${REAL_NATIVE}" "" \
      "${OUT_ROOT}/frozen_scores/s56_human_teacher_v1.parquet"
    ;;
  *)
    echo "Unknown frozen-score task ${TASK}" >&2
    exit 64
    ;;
esac
echo "[s56-t2-frozen-score] task=${TASK} complete=$(date -Is)"
