#!/bin/bash
# Independently validate all five completed native-v3r1 SSL folds and publish
# one immutable aggregate completion release.  This stage trains no model.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-ssl-v3r1-validate-folds
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

readonly TWIRL_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly EXPECTED_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v3r1_s56_s62_a2v1_effective_quality_adp_btjd_bls_eligible"

: "${TWIRL_ORCD_REPO:?controller must export TWIRL_ORCD_REPO}"
: "${TWIRL_EXPECTED_GIT_SHA:?controller must export TWIRL_EXPECTED_GIT_SHA}"
: "${TWIRL_SSL_FULLPOOL_V3R1_RUN_ROOT:?controller must export the exact fresh v3r1 run root}"

REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}"
RUN_ROOT="$(realpath -m -- "${TWIRL_SSL_FULLPOOL_V3R1_RUN_ROOT}")"
PYTHON="${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python"
MODEL_ROOT="${RUN_ROOT}/model_runs/effective_quality_adp_btjd_v2/training/five_fold"
OUTPUT_RELEASE="${RUN_ROOT}/frozen/model/teacher_ssl_fullpool_v3r1_five_fold.complete.json"

if [[ "${RUN_ROOT}" != "${EXPECTED_RUN_ROOT}" ]]; then
  echo "Native-v3r1 fold-validation root differs from the exact contract." >&2
  exit 2
fi
if [[ "${MODEL_ROOT}" != "${EXPECTED_RUN_ROOT}/model_runs/effective_quality_adp_btjd_v2/training/five_fold" ]]; then
  echo "Native-v3r1 fold-validation model namespace differs." >&2
  exit 2
fi
if [[ ! "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]]; then
  echo "TWIRL_EXPECTED_GIT_SHA must be a full lowercase Git SHA." >&2
  exit 2
fi
if [[ "$(git -C "${REPO}" rev-parse HEAD)" != "${EXPECTED_SHA}" ]] || [[ -n "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]]; then
  echo "Native-v3r1 fold-validation repository is not the exact clean revision." >&2
  exit 2
fi
if [[ ! -d "${MODEL_ROOT}/encoder_pretraining" ]]; then
  echo "Native-v3r1 five-fold model root is missing." >&2
  exit 3
fi
for fold in {0..4}; do
  for artifact in summary.json run_contract.json history.csv encoder.pt; do
    if [[ ! -s "${MODEL_ROOT}/encoder_pretraining/fold_${fold}/${artifact}" ]]; then
      echo "Missing native-v3r1 fold-${fold} artifact: ${artifact}" >&2
      exit 3
    fi
  done
done
for output in "${OUTPUT_RELEASE}" "${OUTPUT_RELEASE}.sha256"; do
  if [[ -e "${output}" ]]; then
    echo "Refusing to overwrite or reuse native-v3r1 completion release: ${output}" >&2
    exit 4
  fi
done

mkdir -p "${TWIRL_ROOT}/logs"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[teacher-ssl-fullpool-v3r1-fold-validation] start=$(date -Is)"
"${PYTHON}" scripts/stage5_validation/validate_teacher_ssl_fullpool_v3_training.py \
  --model-root "${MODEL_ROOT}" \
  --expected-code-revision "${EXPECTED_SHA}" \
  --output-release "${OUTPUT_RELEASE}"

test -s "${OUTPUT_RELEASE}"
test -s "${OUTPUT_RELEASE}.sha256"
(
  cd "$(dirname "${OUTPUT_RELEASE}")"
  sha256sum -c "$(basename "${OUTPUT_RELEASE}").sha256"
)
"${PYTHON}" -c 'import json,sys; d=json.load(open(sys.argv[1])); assert d["passed"] is True; assert d["schema_version"]=="twirl_teacher_ssl_fullpool_v3r1_five_fold_completion_release_v2"; assert d["release_binding"]=="teacher_ssl_fullpool_v3r1_effective_quality_adp_btjd_five_fold_complete_v2"; assert d["expected_code_revision"]==sys.argv[2]; assert d["model_namespace"]=="effective_quality_adp_btjd_v2"; assert d["training_hyperparameters"]["folds"]==[0,1,2,3,4]; assert d["training_hyperparameters"]["epochs_per_fold"]==20; assert d["counts"]=={"folds":5,"completed_epochs":100}; assert [f["fold"] for f in d["folds"]]==[0,1,2,3,4]' "${OUTPUT_RELEASE}" "${EXPECTED_SHA}"
echo "[teacher-ssl-fullpool-v3r1-fold-validation] complete=$(date -Is) release=${OUTPUT_RELEASE}"
