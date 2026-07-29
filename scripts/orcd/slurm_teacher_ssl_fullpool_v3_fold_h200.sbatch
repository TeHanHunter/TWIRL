#!/bin/bash
# Five independent native-v3r1 encoders, one H200 each, at most four concurrent.
# No prior fold checkpoint or resume path is admitted.
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-ssl-v3r1-fold
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

readonly TWIRL_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly EXPECTED_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v3r1_s56_s62_a2v1_effective_quality_adp_btjd_bls_eligible"
readonly EPOCHS=20
readonly BATCH_SIZE=64
readonly WORKERS=8
readonly SEED=560064
: "${TWIRL_ORCD_REPO:?controller must export TWIRL_ORCD_REPO}"
: "${TWIRL_EXPECTED_GIT_SHA:?controller must export TWIRL_EXPECTED_GIT_SHA}"
: "${TWIRL_SSL_FULLPOOL_V3R1_RUN_ROOT:?controller must export the exact fresh v3r1 run root}"

REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}"
RUN_ROOT="$(realpath -m -- "${TWIRL_SSL_FULLPOOL_V3R1_RUN_ROOT}")"
PYTHON="${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python"
MODEL_RUN_ROOT="${RUN_ROOT}/model_runs/effective_quality_adp_btjd_v2"
SMOKE_SUMMARY="${MODEL_RUN_ROOT}/smoke/one_epoch/encoder_pretraining/fold_2/summary.json"
OUT_ROOT="${MODEL_RUN_ROOT}/training/five_fold"

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
    echo "Native-v3r1 production fold contract forbids environment override ${forbidden_var}." >&2
    exit 2
  fi
done
if [[ "${RUN_ROOT}" != "${EXPECTED_RUN_ROOT}" ]]; then
  echo "Native-v3r1 fold root differs from the exact production contract." >&2
  exit 2
fi
if [[ "${MODEL_RUN_ROOT}" != "${EXPECTED_RUN_ROOT}/model_runs/effective_quality_adp_btjd_v2" ]]; then
  echo "Native-v3r1 fold model namespace differs from the exact contract." >&2
  exit 2
fi
if [[ "${OUT_ROOT}" != "${MODEL_RUN_ROOT}/training/five_fold" ]]; then
  echo "Native-v3r1 folds must use the fresh five-fold namespace." >&2
  exit 2
fi
if [[ ! "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]]; then
  echo "TWIRL_EXPECTED_GIT_SHA must be a full lowercase Git SHA." >&2
  exit 2
fi
if [[ "$(git -C "${REPO}" rev-parse HEAD)" != "${EXPECTED_SHA}" ]] || [[ -n "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]]; then
  echo "Native-v3r1 fold repository is not the exact clean deployed revision." >&2
  exit 2
fi
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]] || (( SLURM_ARRAY_TASK_ID < 0 || SLURM_ARRAY_TASK_ID > 4 )); then
  echo "Native-v3r1 folds must run as the declared exact 0-4%4 array." >&2
  exit 2
fi

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
NUMERIC_GATE_RELEASE="${RUN_ROOT}/frozen/model_input_numeric_gate_v3r1/model_input_numeric_gate.release.json"

for input in \
  "${ELIGIBILITY_EXCLUSIONS}" \
  "${ELIGIBILITY_SUMMARY}" \
  "${NATIVE_REGISTRY}" \
  "${NATIVE_REGISTRY_SUMMARY}" \
  "${NATIVE_RELEASE_SUMMARY}" \
  "${REGISTRY}" \
  "${REGISTRY_SUMMARY}" \
  "${NUMERIC_GATE_RELEASE}" \
  "${NUMERIC_GATE_RELEASE}.sha256" \
  "${SMOKE_SUMMARY}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing native-v3r1 fold authority: ${input}" >&2
    exit 3
  fi
done
PYTHONPATH="${REPO}/src" "${PYTHON}" -c \
  'from pathlib import Path; import sys; from twirl.vetting.ssl_full_pool_numeric import validate_numeric_gate_release; validate_numeric_gate_release(Path(sys.argv[1]), expected_code_revision=sys.argv[2])' \
  "${NUMERIC_GATE_RELEASE}" "${EXPECTED_SHA}"
FOLD_OUT="${OUT_ROOT}/encoder_pretraining/fold_${SLURM_ARRAY_TASK_ID}"
if [[ -e "${FOLD_OUT}" ]]; then
  echo "Refusing to overwrite or resume existing native-v3r1 fold ${SLURM_ARRAY_TASK_ID}." >&2
  exit 4
fi

mkdir -p "${TWIRL_ROOT}/logs" "${OUT_ROOT}"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

module load cuda/13.1.0 || module load cuda/12.9.1 || true
nvidia-smi

"${PYTHON}" scripts/stage5_validation/validate_teacher_ssl_fullpool_v3_smoke.py \
  --summary "${SMOKE_SUMMARY}" \
  --expected-code-revision "${EXPECTED_SHA}" \
  --eligibility-exclusions "${ELIGIBILITY_EXCLUSIONS}" \
  --eligibility-summary "${ELIGIBILITY_SUMMARY}" \
  --native-registry "${NATIVE_REGISTRY}" \
  --native-registry-summary "${NATIVE_REGISTRY_SUMMARY}" \
  --native-release-summary "${NATIVE_RELEASE_SUMMARY}" \
  --registry "${REGISTRY}" \
  --registry-summary "${REGISTRY_SUMMARY}" \
  --numeric-gate-release "${NUMERIC_GATE_RELEASE}" \
  --expected-fold 2

echo "[teacher-ssl-fullpool-v3r1] start=$(date -Is) fold=${SLURM_ARRAY_TASK_ID} epochs=${EPOCHS}"
"${PYTHON}" scripts/stage5_validation/train_teacher_ssl_fullpool_fold.py \
  --eligibility-exclusions "${ELIGIBILITY_EXCLUSIONS}" \
  --eligibility-summary "${ELIGIBILITY_SUMMARY}" \
  --native-registry "${NATIVE_REGISTRY}" \
  --native-registry-summary "${NATIVE_REGISTRY_SUMMARY}" \
  --native-release-summary "${NATIVE_RELEASE_SUMMARY}" \
  --registry "${REGISTRY}" \
  --registry-summary "${REGISTRY_SUMMARY}" \
  --numeric-gate-release "${NUMERIC_GATE_RELEASE}" \
  --out-root "${OUT_ROOT}" \
  --fold "${SLURM_ARRAY_TASK_ID}" \
  --epochs "${EPOCHS}" \
  --batch-size "${BATCH_SIZE}" \
  --workers "${WORKERS}" \
  --seed "${SEED}"
echo "[teacher-ssl-fullpool-v3r1] complete=$(date -Is) fold=${SLURM_ARRAY_TASK_ID}"
