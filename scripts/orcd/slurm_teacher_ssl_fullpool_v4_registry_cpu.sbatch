#!/bin/bash
# Build the immutable 175,366-row SSL audit registry from the corrected v4
# BLS eligibility partition and the fresh detector-consistent native release.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-ssl-v4-registry
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

readonly TWIRL_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly EXPECTED_SOURCE_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp"
readonly EXPECTED_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v4r2_s56_s62_a2v1_raw_v1_detector_consistent_adp_bls_eligible"
: "${TWIRL_ORCD_REPO:?controller must export TWIRL_ORCD_REPO}"
: "${TWIRL_EXPECTED_GIT_SHA:?controller must export TWIRL_EXPECTED_GIT_SHA}"
: "${TWIRL_SSL_FULLPOOL_V1_RUN_ROOT:?controller must export the immutable source root}"
: "${TWIRL_SSL_FULLPOOL_V4_RUN_ROOT:?controller must export the exact fresh v4 run root}"

REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}"
SOURCE_RUN_ROOT="$(realpath -m -- "${TWIRL_SSL_FULLPOOL_V1_RUN_ROOT}")"
RUN_ROOT="$(realpath -m -- "${TWIRL_SSL_FULLPOOL_V4_RUN_ROOT}")"
PYTHON="${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python"

if [[ "${SOURCE_RUN_ROOT}" != "${EXPECTED_SOURCE_RUN_ROOT}" || "${RUN_ROOT}" != "${EXPECTED_RUN_ROOT}" ]]; then
  echo "SSL registry-v4 roots differ from the exact production contract." >&2
  exit 2
fi
if [[ ! "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]]; then
  echo "TWIRL_EXPECTED_GIT_SHA must be a full lowercase Git SHA." >&2
  exit 2
fi
if [[ "$(git -C "${REPO}" rev-parse HEAD)" != "${EXPECTED_SHA}" ]] || [[ -n "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]]; then
  echo "SSL registry-v4 repository is not the exact clean deployed revision." >&2
  exit 2
fi

POOL_DIR="${SOURCE_RUN_ROOT}/frozen/full_pool"
POOL_TABLE="${POOL_DIR}/teacher_ssl_full_pool_observations.csv"
POOL_SUMMARY="${POOL_DIR}/teacher_ssl_full_pool_manifest.summary.json"
GLOBAL_BLS_DIR="${RUN_ROOT}/real_bls/global"
BLS_PEAKS="${GLOBAL_BLS_DIR}/teacher_ssl_full_pool_bls.parquet"
BLS_SUMMARY="${GLOBAL_BLS_DIR}/teacher_ssl_full_pool_bls.summary.json"
ELIGIBILITY_DIR="${RUN_ROOT}/frozen/native_model_eligibility"
ELIGIBILITY_EXCLUSIONS="${ELIGIBILITY_DIR}/native_model_exclusions.csv"
ELIGIBILITY_SUMMARY="${ELIGIBILITY_DIR}/native_model_eligibility.summary.json"
NATIVE_REGISTRY_DIR="${RUN_ROOT}/frozen/native_registry"
NATIVE_REGISTRY="${NATIVE_REGISTRY_DIR}/native_input_registry.parquet"
NATIVE_REGISTRY_SUMMARY="${NATIVE_REGISTRY_DIR}/native_input_registry.summary.json"
NATIVE_RELEASE_SUMMARY="${NATIVE_REGISTRY_DIR}/native_fullpool_release.summary.json"
SPLIT_REGISTRY="${TWIRL_ROOT}/reports/stage5_validation/teacher_v3_s56_s62_a2v1_current_adp/frozen/tic_split_registry.csv"
S63_TICS="${SOURCE_RUN_ROOT}/preregistered/s63_reserved_tics.txt"
REGISTRY_DIR="${RUN_ROOT}/frozen/ssl_registry"
REGISTRY_OUT="${REGISTRY_DIR}/teacher_ssl_fullpool_registry.parquet"
SUMMARY_OUT="${REGISTRY_DIR}/teacher_ssl_fullpool_registry.summary.json"

inputs=(
  "${POOL_TABLE}"
  "${POOL_SUMMARY}"
  "${BLS_PEAKS}"
  "${BLS_SUMMARY}"
  "${ELIGIBILITY_EXCLUSIONS}"
  "${ELIGIBILITY_SUMMARY}"
  "${NATIVE_REGISTRY}"
  "${NATIVE_REGISTRY_SUMMARY}"
  "${NATIVE_RELEASE_SUMMARY}"
  "${SPLIT_REGISTRY}"
  "${S63_TICS}"
)
for input in "${inputs[@]}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing SSL registry-v4 authority: ${input}" >&2
    exit 3
  fi
done
for sector in {56..62}; do
  if [[ ! -s "${POOL_DIR}/allowlists/s${sector}_tics.csv" ]]; then
    echo "Missing frozen S${sector} allowlist." >&2
    exit 3
  fi
done
for output in "${REGISTRY_OUT}" "${SUMMARY_OUT}"; do
  if [[ -e "${output}" ]]; then
    echo "Refusing to overwrite an existing SSL registry-v4 artifact: ${output}" >&2
    exit 4
  fi
done

mkdir -p "${TWIRL_ROOT}/logs" "${REGISTRY_DIR}"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[ssl-full-pool-registry-v4] host=$(hostname) start=$(date -Is)"
"${PYTHON}" scripts/stage5_validation/build_teacher_ssl_fullpool_registry.py \
  --frozen-pool "${POOL_TABLE}" \
  --frozen-pool-summary "${POOL_SUMMARY}" \
  --bls-summary "${BLS_SUMMARY}" \
  --bls-peaks "${BLS_PEAKS}" \
  --eligibility-exclusions "${ELIGIBILITY_EXCLUSIONS}" \
  --eligibility-summary "${ELIGIBILITY_SUMMARY}" \
  --native-registry "${NATIVE_REGISTRY}" \
  --native-registry-summary "${NATIVE_REGISTRY_SUMMARY}" \
  --native-release-summary "${NATIVE_RELEASE_SUMMARY}" \
  --frozen-split-registry "${SPLIT_REGISTRY}" \
  --reserved-hosts "${S63_TICS}" \
  --registry-out "${REGISTRY_OUT}" \
  --summary-out "${SUMMARY_OUT}"

for output in "${REGISTRY_OUT}" "${SUMMARY_OUT}"; do
  if [[ ! -s "${output}" ]]; then
    echo "SSL registry-v4 producer did not publish: ${output}" >&2
    exit 5
  fi
done
echo "[ssl-full-pool-registry-v4] complete=$(date -Is) output=${REGISTRY_OUT}"
