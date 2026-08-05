#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-ssl-fullpool
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=192G
#SBATCH --gres=gpu:h200:1
#SBATCH --array=0-4%4
#SBATCH -t 2-00:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.err

set -euo pipefail

TWIRL_ROOT="/orcd/data/mki_aryeh/001/twirl"
: "${TWIRL_ORCD_REPO:?controller must export TWIRL_ORCD_REPO}"
: "${TWIRL_EXPECTED_GIT_SHA:?controller must export TWIRL_EXPECTED_GIT_SHA}"
REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}"
PYTHON="${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python"
EXPECTED_SOURCE_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp"
EXPECTED_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v2_s56_s62_a2v1_current_adp_bls_eligible"
SOURCE_RUN_ROOT="$(realpath -m -- "${TWIRL_SSL_FULLPOOL_V1_RUN_ROOT:-${EXPECTED_SOURCE_RUN_ROOT}}")"
RUN_ROOT="$(realpath -m -- "${TWIRL_SSL_FULLPOOL_V2_RUN_ROOT:-${EXPECTED_RUN_ROOT}}")"
ELIGIBILITY_DIR="${RUN_ROOT}/frozen/native_model_eligibility"
ELIGIBILITY_EXCLUSIONS="${ELIGIBILITY_DIR}/native_model_exclusions.csv"
ELIGIBILITY_SUMMARY="${ELIGIBILITY_DIR}/native_model_eligibility.summary.json"
NATIVE_REGISTRY_DIR="${RUN_ROOT}/frozen/native_registry"
NATIVE_REGISTRY="${NATIVE_REGISTRY_DIR}/native_input_registry.parquet"
NATIVE_REGISTRY_SUMMARY="${NATIVE_REGISTRY_DIR}/native_input_registry.summary.json"
NATIVE_RELEASE_SUMMARY="${NATIVE_REGISTRY_DIR}/native_fullpool_release.summary.json"
REGISTRY_DIR="${RUN_ROOT}/frozen/ssl_registry"
REGISTRY="${REGISTRY_DIR}/teacher_ssl_fullpool_registry.parquet"
REGISTRY_SUMMARY="${REGISTRY_DIR}/teacher_ssl_fullpool_registry.summary.json"
SMOKE_SUMMARY="${RUN_ROOT}/smoke/one_epoch/encoder_pretraining/fold_0/summary.json"
OUT_ROOT="${RUN_ROOT}/training/five_fold"
readonly EPOCHS=20
readonly BATCH_SIZE=64
readonly WORKERS=8
readonly SEED=560064

for forbidden_var in \
  TWIRL_SSL_FULLPOOL_EPOCHS \
  TWIRL_SSL_FULLPOOL_BATCH_SIZE \
  TWIRL_SSL_FULLPOOL_WORKERS \
  TWIRL_SSL_FULLPOOL_SEED \
  TWIRL_SSL_FULLPOOL_OUT_ROOT \
  TWIRL_SSL_FULLPOOL_MAX_ROWS \
  TWIRL_SSL_FULLPOOL_SMOKE_OUT_ROOT \
  TWIRL_SSL_FULLPOOL_SMOKE_FOLD \
  TWIRL_SSL_FULLPOOL_SMOKE_EPOCHS \
  TWIRL_SSL_FULLPOOL_SMOKE_MAX_ROWS \
  TWIRL_SSL_FULLPOOL_SMOKE_BATCH_SIZE \
  TWIRL_SSL_FULLPOOL_SMOKE_WORKERS; do
  if [[ -n "${!forbidden_var+x}" ]]; then
    echo "Production fold contract forbids environment override ${forbidden_var}." >&2
    exit 3
  fi
done
if [[ "${SOURCE_RUN_ROOT}" != "${EXPECTED_SOURCE_RUN_ROOT}" || "${RUN_ROOT}" != "${EXPECTED_RUN_ROOT}" ]]; then
  echo "Fold roots differ from the exact Teacher v4-SSL v2 contract." >&2
  exit 3
fi
if [[ "${OUT_ROOT}" != "${EXPECTED_RUN_ROOT}/training/five_fold" ]]; then
  echo "Training output must use the exact v2 five-fold namespace." >&2
  exit 3
fi
if [[ ! "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]]; then
  echo "TWIRL_EXPECTED_GIT_SHA must be a full lowercase Git SHA." >&2
  exit 3
fi
if [[ "$(git -C "${REPO}" rev-parse HEAD)" != "${EXPECTED_SHA}" ]] || [[ -n "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]]; then
  echo "Fold repository is not the clean expected deployed revision." >&2
  exit 3
fi
if [[ "${OUT_ROOT}" == "${SOURCE_RUN_ROOT}" || "${OUT_ROOT}" == "${SOURCE_RUN_ROOT}/"* ]]; then
  echo "Training output may not be published under the v1 run root." >&2
  exit 3
fi

for input in \
  "${ELIGIBILITY_EXCLUSIONS}" \
  "${ELIGIBILITY_SUMMARY}" \
  "${NATIVE_REGISTRY}" \
  "${NATIVE_REGISTRY_SUMMARY}" \
  "${NATIVE_RELEASE_SUMMARY}" \
  "${REGISTRY}" \
  "${REGISTRY_SUMMARY}" \
  "${SMOKE_SUMMARY}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing full-pool SSL input: ${input}" >&2
    exit 3
  fi
done
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "This wrapper must run as the declared fold array." >&2
  exit 4
fi

mkdir -p "${TWIRL_ROOT}/logs" "${OUT_ROOT}"
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

echo "[teacher-ssl-fullpool] start=$(date -Is) fold=${SLURM_ARRAY_TASK_ID} epochs=${EPOCHS} max_rows=full out_root=${OUT_ROOT}"
"${PYTHON}" \
  scripts/stage5_validation/validate_teacher_ssl_fullpool_v2_smoke.py \
  --summary "${SMOKE_SUMMARY}" \
  --expected-code-revision "${EXPECTED_SHA}" \
  --eligibility-exclusions "${ELIGIBILITY_EXCLUSIONS}" \
  --eligibility-summary "${ELIGIBILITY_SUMMARY}" \
  --native-registry "${NATIVE_REGISTRY}" \
  --native-registry-summary "${NATIVE_REGISTRY_SUMMARY}" \
  --native-release-summary "${NATIVE_RELEASE_SUMMARY}" \
  --registry "${REGISTRY}" \
  --registry-summary "${REGISTRY_SUMMARY}" \
  --expected-fold 0
"${PYTHON}" scripts/stage5_validation/train_teacher_ssl_fullpool_fold.py \
  --eligibility-exclusions "${ELIGIBILITY_EXCLUSIONS}" \
  --eligibility-summary "${ELIGIBILITY_SUMMARY}" \
  --native-registry "${NATIVE_REGISTRY}" \
  --native-registry-summary "${NATIVE_REGISTRY_SUMMARY}" \
  --native-release-summary "${NATIVE_RELEASE_SUMMARY}" \
  --registry "${REGISTRY}" \
  --registry-summary "${REGISTRY_SUMMARY}" \
  --out-root "${OUT_ROOT}" \
  --fold "${SLURM_ARRAY_TASK_ID}" \
  --epochs "${EPOCHS}" \
  --batch-size "${BATCH_SIZE}" \
  --workers "${WORKERS}" \
  --seed "${SEED}" \
  --resume
echo "[teacher-ssl-fullpool] complete=$(date -Is) fold=${SLURM_ARRAY_TASK_ID}"
