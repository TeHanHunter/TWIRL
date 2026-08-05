#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-teacher-v4-ssl
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=192G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 08:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:-${TWIRL_ROOT}/code/TWIRL}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python}"
CONFIG="${TWIRL_TEACHER_SSL_V1_CONFIG:-${REPO}/configs/models/teacher_ssl_v1_s56_s62_a2v1_current_adp.yaml}"
REPORT_ROOT="${TWIRL_TEACHER_SSL_V1_REPORT_ROOT:-${TWIRL_ROOT}/reports/stage5_validation/teacher_v3_s56_s62_a2v1_current_adp}"
TRAINING_TABLE="${TWIRL_TEACHER_SSL_V1_TRAINING_TABLE:-${REPORT_ROOT}/native_preparation/training_table_with_splits.csv}"
NATIVE_REGISTRY="${TWIRL_TEACHER_SSL_V1_NATIVE_REGISTRY:-${REPORT_ROOT}/native/native_input_registry.csv}"
NATIVE_REGISTRY_SUMMARY="${TWIRL_TEACHER_SSL_V1_NATIVE_REGISTRY_SUMMARY:-${REPORT_ROOT}/native/native_input_registry.summary.json}"
BASELINE_TEACHER_V3_ROOT="${TWIRL_TEACHER_SSL_V1_BASELINE_ROOT:-${TWIRL_ROOT}/results/models/teacher_v3_s56_s62_a2v1_current_adp/full}"
OUT_DIR="${TWIRL_TEACHER_SSL_V1_OUT_DIR:-${TWIRL_ROOT}/results/models/teacher_ssl_v1_s56_s62_a2v1_current_adp/pilot}"

for input in \
  "${CONFIG}" \
  "${TRAINING_TABLE}" \
  "${NATIVE_REGISTRY}" \
  "${NATIVE_REGISTRY_SUMMARY}" \
  "${BASELINE_TEACHER_V3_ROOT}/summary.json"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing Teacher v4-SSL pilot input: ${input}" >&2
    exit 3
  fi
done
if [[ -e "${OUT_DIR}" ]]; then
  echo "Refusing to reuse non-fresh Teacher v4-SSL output: ${OUT_DIR}" >&2
  exit 4
fi

mkdir -p "${TWIRL_ROOT}/logs" "$(dirname "${OUT_DIR}")"
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

echo "[teacher-v4-ssl] start=$(date -Is) output=${OUT_DIR}"
"${PYTHON}" scripts/stage5_validation/train_teacher_ssl_v1.py \
  --config "${CONFIG}" \
  --input-root "${TWIRL_ROOT}" \
  --training-table "${TRAINING_TABLE}" \
  --native-registry "${NATIVE_REGISTRY}" \
  --native-registry-summary "${NATIVE_REGISTRY_SUMMARY}" \
  --baseline-teacher-v3-root "${BASELINE_TEACHER_V3_ROOT}" \
  --out-dir "${OUT_DIR}"
echo "[teacher-v4-ssl] complete=$(date -Is) output=${OUT_DIR}"
