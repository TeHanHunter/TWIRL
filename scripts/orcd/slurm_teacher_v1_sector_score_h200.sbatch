#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-v1-score
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
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56-torch/bin/python}"
SECTOR="${TWIRL_ENRICH_SECTOR:?set TWIRL_ENRICH_SECTOR}"
WORK_ROOT="${TWIRL_ENRICH_WORK_ROOT:?set TWIRL_ENRICH_WORK_ROOT}"
MODEL_ROOT="${TWIRL_TEACHER_V1_MODEL_ROOT:?set TWIRL_TEACHER_V1_MODEL_ROOT}"
CANDIDATES="${TWIRL_ENRICH_CANDIDATES:-${WORK_ROOT}/candidates/real_candidates_top5.parquet}"
NATIVE_H5="${TWIRL_ENRICH_NATIVE_H5:-${WORK_ROOT}/real_native/s${SECTOR}_a2v1_real_native.h5}"
SCORE_ROOT="${TWIRL_ENRICH_SCORE_ROOT:-${WORK_ROOT}/teacher_v1_scores}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs/teacher_enrichment "${SCORE_ROOT}"
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
module load cuda/13.1.0 || module load cuda/12.9.1 || true
nvidia-smi
"${PYTHON}" -c 'import torch; assert torch.cuda.is_available(), "CUDA is unavailable in the Teacher-v1 scoring environment"'

checkpoints=()
for fold in 0 1 2 3 4; do
  checkpoint="${MODEL_ROOT}/fold_${fold}/teacher.pt"
  [[ -s "${checkpoint}" ]] || {
    echo "Missing Teacher-v1 checkpoint: ${checkpoint}" >&2
    exit 2
  }
  checkpoints+=("${checkpoint}")
done

"${PYTHON}" scripts/stage5_validation/score_s56_teacher_v1_candidates.py \
  --candidates "${CANDIDATES}" \
  --native-h5 "${NATIVE_H5}" \
  --checkpoints "${checkpoints[@]}" \
  --out-dir "${SCORE_ROOT}" \
  --profile shape_plus_periodogram_bls \
  --batch-size "${TWIRL_ENRICH_BATCH_SIZE:-32}" \
  --workers "${TWIRL_ENRICH_WORKERS:-8}"

echo "[teacher-v1-sector-score] sector=${SECTOR} complete=$(date -Is)"
