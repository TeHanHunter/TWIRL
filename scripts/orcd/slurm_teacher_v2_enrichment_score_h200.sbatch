#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-enrich-score
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=192G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 08:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/teacher_enrichment/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/teacher_enrichment/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_teacher_v2}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
MODEL_ROOT="${TWIRL_TEACHER_V2_MODEL_ROOT:-${REPO}/reports/stage5_validation/s56_A2v1_teacher_v2/training}"
CANDIDATES="${TWIRL_ENRICH_CANDIDATES:?set TWIRL_ENRICH_CANDIDATES}"
NATIVE_H5="${TWIRL_ENRICH_NATIVE_H5:?set TWIRL_ENRICH_NATIVE_H5}"
SCORE_ROOT="${TWIRL_ENRICH_SCORE_ROOT:?set TWIRL_ENRICH_SCORE_ROOT}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs/teacher_enrichment "${SCORE_ROOT}"
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

score_profile() {
  local profile="$1"
  local checkpoint_args=()
  local fold
  for fold in 0 1 2 3 4; do
    checkpoint_args+=(--checkpoint "${MODEL_ROOT}/${profile}/fold_${fold}/teacher_v2.pt")
  done
  "${PYTHON}" scripts/stage5_validation/score_teacher_v2_candidates.py \
    --candidates "${CANDIDATES}" \
    --native-h5 "${NATIVE_H5}" \
    "${checkpoint_args[@]}" \
    --profile "${profile}" \
    --out-path "${SCORE_ROOT}/${profile}.parquet" \
    --batch-size "${TWIRL_ENRICH_BATCH_SIZE:-32}" \
    --workers "${TWIRL_ENRICH_WORKERS:-8}"
}

score_profile shape_plus_raw_chronology
score_profile shape_plus_periodogram_bls
echo "[enrichment-score] complete=$(date -Is)"
