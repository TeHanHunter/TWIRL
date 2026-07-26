#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-teacher-v3
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=192G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 24:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:-${TWIRL_ROOT}/code/TWIRL}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python}"
CONFIG="${TWIRL_TEACHER_V3_CONFIG:-${REPO}/configs/models/teacher_v3_s56_s62_a2v1_current_adp.yaml}"
OUTPUT_ROOT="${TWIRL_TEACHER_V3_OUTPUT_ROOT:-${TWIRL_ROOT}/results/models/teacher_v3_s56_s62_a2v1_current_adp}"
METADATA_ONLY="${TWIRL_TEACHER_V3_METADATA_ONLY_BOOTSTRAP:-0}"

for input in \
  "${CONFIG}" \
  "${TWIRL_ROOT}/reports/stage5_validation/teacher_v3_s56_s62_a2v1_current_adp/frozen/observation_morphology_corpus.csv" \
  "${TWIRL_ROOT}/reports/stage5_validation/teacher_v3_s56_s62_a2v1_current_adp/frozen/tic_split_registry.csv" \
  "${TWIRL_ROOT}/reports/stage5_validation/teacher_v3_s56_s62_a2v1_current_adp/native_preparation/training_table_with_splits.csv"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing frozen Teacher-v3 input: ${input}" >&2
    exit 3
  fi
done
if [[ "${METADATA_ONLY}" != "0" && "${METADATA_ONLY}" != "1" ]]; then
  echo "TWIRL_TEACHER_V3_METADATA_ONLY_BOOTSTRAP must be 0 or 1." >&2
  exit 2
fi
if [[ "${METADATA_ONLY}" == "0" ]]; then
  for input in \
    "${TWIRL_ROOT}/reports/stage5_validation/teacher_v3_s56_s62_a2v1_current_adp/native/native_input_registry.csv" \
    "${TWIRL_ROOT}/reports/stage5_validation/teacher_v3_s56_s62_a2v1_current_adp/native/native_input_registry.summary.json"; do
    if [[ ! -s "${input}" ]]; then
      echo "Missing full Teacher-v3 native input: ${input}" >&2
      echo "Use TWIRL_TEACHER_V3_METADATA_ONLY_BOOTSTRAP=1 only for the declared real-label baseline." >&2
      exit 3
    fi
  done
fi

mkdir -p "${TWIRL_ROOT}/logs" "${OUTPUT_ROOT}"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

module load cuda/13.1.0 || module load cuda/12.9.1 || true
nvidia-smi

args=(
  --config "${CONFIG}"
  --input-root "${TWIRL_ROOT}"
  --output-root "${OUTPUT_ROOT}"
)
if [[ "${METADATA_ONLY}" == "1" ]]; then
  args+=(--metadata-only-bootstrap)
fi

echo "[teacher-v3] metadata_only_bootstrap=${METADATA_ONLY} start=$(date -Is)"
"${PYTHON}" scripts/stage5_validation/train_teacher_v3.py "${args[@]}"
echo "[teacher-v3] complete=$(date -Is) output=${OUTPUT_ROOT}"
