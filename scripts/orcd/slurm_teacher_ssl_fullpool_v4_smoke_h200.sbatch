#!/bin/bash
# Bounded detector-consistent v4 throughput smoke: fold 2, 4,096 rows,
# one epoch, one H200, and the locked diagnostic observation.
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-ssl-v4-smoke
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=192G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 08:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

readonly TWIRL_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly EXPECTED_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v4r2_s56_s62_a2v1_raw_v1_detector_consistent_adp_bls_eligible"
readonly FOLD=2
readonly EPOCHS=1
readonly MAX_ROWS=4096
readonly BATCH_SIZE=64
readonly WORKERS=8
readonly SEED=560064
readonly REQUIRED_OBSERVATION_ID="s0060-tic0000000722078603"
: "${TWIRL_ORCD_REPO:?controller must export TWIRL_ORCD_REPO}"
: "${TWIRL_EXPECTED_GIT_SHA:?controller must export TWIRL_EXPECTED_GIT_SHA}"
: "${TWIRL_SSL_FULLPOOL_V4_RUN_ROOT:?controller must export the exact fresh v4 run root}"

REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}"
RUN_ROOT="$(realpath -m -- "${TWIRL_SSL_FULLPOOL_V4_RUN_ROOT}")"
PYTHON="${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python"
MODEL_RUN_ROOT="${RUN_ROOT}/model_runs/effective_quality_adp_btjd_v2"
OUT_ROOT="${MODEL_RUN_ROOT}/smoke/one_epoch"

if [[ "${RUN_ROOT}" != "${EXPECTED_RUN_ROOT}" ]]; then
  echo "Native-v4 smoke root differs from the exact production contract." >&2
  exit 2
fi
if [[ ! "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]] || [[ "$(git -C "${REPO}" rev-parse HEAD)" != "${EXPECTED_SHA}" ]] || [[ -n "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]]; then
  echo "Native-v4 smoke repository is not the exact clean deployed revision." >&2
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
NUMERIC_GATE_RELEASE="${RUN_ROOT}/frozen/model_input_numeric_gate_v4/model_input_numeric_gate.release.json"

for input in \
  "${ELIGIBILITY_EXCLUSIONS}" \
  "${ELIGIBILITY_SUMMARY}" \
  "${NATIVE_REGISTRY}" \
  "${NATIVE_REGISTRY_SUMMARY}" \
  "${NATIVE_RELEASE_SUMMARY}" \
  "${REGISTRY}" \
  "${REGISTRY_SUMMARY}" \
  "${NUMERIC_GATE_RELEASE}" \
  "${NUMERIC_GATE_RELEASE}.sha256"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing native-v4 smoke authority: ${input}" >&2
    exit 3
  fi
done
PYTHONPATH="${REPO}/src" "${PYTHON}" -c \
  'from pathlib import Path; import sys; from twirl.vetting.ssl_full_pool_numeric import validate_numeric_gate_release; validate_numeric_gate_release(Path(sys.argv[1]), expected_code_revision=sys.argv[2])' \
  "${NUMERIC_GATE_RELEASE}" "${EXPECTED_SHA}"
if [[ -e "${OUT_ROOT}" ]]; then
  echo "Refusing to overwrite an existing native-v4 smoke namespace." >&2
  exit 4
fi

mkdir -p "${TWIRL_ROOT}/logs" "${OUT_ROOT}"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
module load cuda/13.1.0 || module load cuda/12.9.1 || true
nvidia-smi

echo "[teacher-ssl-fullpool-v4-smoke] start=$(date -Is) fold=${FOLD} rows=${MAX_ROWS}"
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
  --fold "${FOLD}" \
  --epochs "${EPOCHS}" \
  --batch-size "${BATCH_SIZE}" \
  --workers "${WORKERS}" \
  --seed "${SEED}" \
  --max-rows "${MAX_ROWS}" \
  --required-observation-id "${REQUIRED_OBSERVATION_ID}"

SMOKE_SUMMARY="${OUT_ROOT}/encoder_pretraining/fold_${FOLD}/summary.json"
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
  --expected-fold "${FOLD}"
echo "[teacher-ssl-fullpool-v4-smoke] complete=$(date -Is) summary=${SMOKE_SUMMARY}"
