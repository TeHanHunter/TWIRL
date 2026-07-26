#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-tv3-registry
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:-${TWIRL_ROOT}/code/TWIRL}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python}"
RUN_REL="reports/stage5_validation/teacher_v3_s56_s62_a2v1_current_adp"
SOURCE="${TWIRL_TEACHER_V3_NATIVE_REGISTRY_SOURCE:-${TWIRL_ROOT}/${RUN_REL}/native_preparation/real_native_registry_source.csv}"
REGISTRY="${TWIRL_TEACHER_V3_NATIVE_REGISTRY:-${TWIRL_ROOT}/${RUN_REL}/native/native_input_registry.csv}"
SUMMARY="${TWIRL_TEACHER_V3_NATIVE_REGISTRY_SUMMARY:-${TWIRL_ROOT}/${RUN_REL}/native/native_input_registry.summary.json}"

if [[ ! -s "${SOURCE}" ]]; then
  echo "Missing Teacher-v3 native-registry source: ${SOURCE}" >&2
  exit 3
fi
mkdir -p "${TWIRL_ROOT}/logs" "$(dirname "${REGISTRY}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

"${PYTHON}" scripts/stage5_validation/build_teacher_v1_native_input_registry.py \
  build \
  --source-table "${SOURCE}" \
  --path-base "${TWIRL_ROOT}" \
  --out-registry "${REGISTRY}" \
  --out-summary "${SUMMARY}"

"${PYTHON}" scripts/stage5_validation/build_teacher_v1_native_input_registry.py \
  validate \
  --registry "${REGISTRY}" \
  --summary "${SUMMARY}"

echo "[teacher-v3-native-registry] complete=$(date -Is) registry=${REGISTRY}"
