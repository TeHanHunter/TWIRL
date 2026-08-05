#!/bin/bash
# Validate all five completed v4 SSL folds and publish one immutable release.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-ssl-v4-validate-folds
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

readonly TWIRL_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly EXPECTED_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v4r2_s56_s62_a2v1_raw_v1_detector_consistent_adp_bls_eligible"
: "${TWIRL_ORCD_REPO:?controller must export TWIRL_ORCD_REPO}"
: "${TWIRL_VALIDATOR_GIT_SHA:?controller must export TWIRL_VALIDATOR_GIT_SHA}"
: "${TWIRL_EXPECTED_TRAINING_GIT_SHA:?controller must export TWIRL_EXPECTED_TRAINING_GIT_SHA}"
: "${TWIRL_SSL_FULLPOOL_V4_RUN_ROOT:?controller must export the exact fresh v4 run root}"

REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
VALIDATOR_SHA="${TWIRL_VALIDATOR_GIT_SHA}"
TRAINING_SHA="${TWIRL_EXPECTED_TRAINING_GIT_SHA}"
RUN_ROOT="$(realpath -m -- "${TWIRL_SSL_FULLPOOL_V4_RUN_ROOT}")"
PYTHON="${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python"
MODEL_ROOT="${RUN_ROOT}/model_runs/effective_quality_adp_btjd_v2/training/five_fold"
OUTPUT_RELEASE="${RUN_ROOT}/frozen/model/teacher_ssl_fullpool_v4_five_fold.complete.json"

if [[ "${RUN_ROOT}" != "${EXPECTED_RUN_ROOT}" ]]; then
  echo "Native-v4 fold-validation root differs from the exact contract." >&2
  exit 2
fi
if [[ ! "${VALIDATOR_SHA}" =~ ^[0-9a-f]{40}$ ]] || [[ ! "${TRAINING_SHA}" =~ ^[0-9a-f]{40}$ ]] || [[ "$(git -C "${REPO}" rev-parse HEAD)" != "${VALIDATOR_SHA}" ]] || [[ -n "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]]; then
  echo "Native-v4 fold-validation repository is not the exact clean revision." >&2
  exit 2
fi
for fold in {0..4}; do
  for artifact in summary.json run_contract.json history.csv encoder.pt; do
    if [[ ! -s "${MODEL_ROOT}/encoder_pretraining/fold_${fold}/${artifact}" ]]; then
      echo "Missing native-v4 fold-${fold} artifact: ${artifact}" >&2
      exit 3
    fi
  done
done
for output in "${OUTPUT_RELEASE}" "${OUTPUT_RELEASE}.sha256"; do
  if [[ -e "${output}" ]]; then
    echo "Refusing to overwrite an existing native-v4 completion release: ${output}" >&2
    exit 4
  fi
done

mkdir -p "${TWIRL_ROOT}/logs" "$(dirname "${OUTPUT_RELEASE}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

echo "[teacher-ssl-fullpool-v4-fold-validation] start=$(date -Is)"
"${PYTHON}" scripts/stage5_validation/validate_teacher_ssl_fullpool_v4_training.py \
  --model-root "${MODEL_ROOT}" \
  --expected-code-revision "${TRAINING_SHA}" \
  --validator-code-revision "${VALIDATOR_SHA}" \
  --output-release "${OUTPUT_RELEASE}"

test -s "${OUTPUT_RELEASE}"
test -s "${OUTPUT_RELEASE}.sha256"
(
  cd "$(dirname "${OUTPUT_RELEASE}")"
  sha256sum -c "$(basename "${OUTPUT_RELEASE}").sha256"
)
"${PYTHON}" -c 'import json,sys; d=json.load(open(sys.argv[1])); assert d["passed"] is True; assert d["schema_version"]=="twirl_teacher_ssl_fullpool_v4_five_fold_completion_release_v1"; assert d["release_binding"]=="teacher_ssl_fullpool_v4_detector_consistent_raw_v1_effective_quality_adp_btjd_five_fold_complete_v1"; assert d["expected_code_revision"]==sys.argv[2]; assert d["validator_code_revision"]==sys.argv[3]; assert d["training_hyperparameters"]["folds"]==[0,1,2,3,4]; assert d["training_hyperparameters"]["epochs_per_fold"]==20; assert d["counts"]=={"folds":5,"completed_epochs":100}' "${OUTPUT_RELEASE}" "${TRAINING_SHA}" "${VALIDATOR_SHA}"
echo "[teacher-ssl-fullpool-v4-fold-validation] complete=$(date -Is) release=${OUTPUT_RELEASE}"
