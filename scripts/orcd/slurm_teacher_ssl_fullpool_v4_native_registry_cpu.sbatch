#!/bin/bash
# Freeze the exact 112 detector-consistent native-v4 shards.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-ssl-v4-native-registry
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:?TWIRL_ORCD_REPO is required}"
EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA:?TWIRL_EXPECTED_GIT_SHA is required}"
SOURCE_RUN_ROOT="${TWIRL_SSL_FULLPOOL_SOURCE_RUN_ROOT:-${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp}"
RUN_ROOT="${TWIRL_SSL_FULLPOOL_RUN_ROOT:?TWIRL_SSL_FULLPOOL_RUN_ROOT is required}"
EXPECTED_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v4r2_s56_s62_a2v1_raw_v1_detector_consistent_adp_bls_eligible"
PYTHON="${TWIRL_ORCD_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python}"

if [[ "$(realpath -m -- "${RUN_ROOT}")" != "${EXPECTED_RUN_ROOT}" ]]; then
  echo "Native-v4 registry root differs from the exact production contract." >&2
  exit 2
fi
if [[ ! "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]] || [[ "$(git -C "${REPO}" rev-parse HEAD)" != "${EXPECTED_SHA}" ]] || [[ -n "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]]; then
  echo "Native-v4 registry repository is not the exact clean revision." >&2
  exit 2
fi

POOL_DIR="${SOURCE_RUN_ROOT}/frozen/full_pool"
POOL_TABLE="${POOL_DIR}/teacher_ssl_full_pool_observations.csv"
POOL_SUMMARY="${POOL_DIR}/teacher_ssl_full_pool_manifest.summary.json"
ELIGIBILITY_DIR="${RUN_ROOT}/frozen/native_model_eligibility"
ELIGIBILITY_EXCLUSIONS="${ELIGIBILITY_DIR}/native_model_exclusions.csv"
ELIGIBILITY_SUMMARY="${ELIGIBILITY_DIR}/native_model_eligibility.summary.json"
NATIVE_ROOT="${RUN_ROOT}/native"
REGISTRY_DIR="${RUN_ROOT}/frozen/native_registry"
SOURCE_OUT="${REGISTRY_DIR}/native_registry_source.csv"
REGISTRY_OUT="${REGISTRY_DIR}/native_input_registry.parquet"
SUMMARY_OUT="${REGISTRY_DIR}/native_input_registry.summary.json"
RELEASE_SUMMARY_OUT="${REGISTRY_DIR}/native_fullpool_release.summary.json"

for input in "${POOL_TABLE}" "${POOL_SUMMARY}" "${ELIGIBILITY_EXCLUSIONS}" "${ELIGIBILITY_SUMMARY}"; do
  test -s "${input}"
done
native_shards=()
native_summaries=()
for sector in {56..62}; do
  for shard_index in {0..15}; do
    shard="${NATIVE_ROOT}/s${sector}/native_${shard_index}.h5"
    summary="${NATIVE_ROOT}/s${sector}/native_${shard_index}.summary.json"
    sidecar="${shard}.sha256"
    test -s "${shard}" && test -s "${summary}" && test -s "${sidecar}"
    sha256sum -c "${sidecar}"
    native_shards+=("${shard}")
    native_summaries+=("${summary}")
  done
done
if [[ "${#native_shards[@]}" -ne 112 ]]; then
  echo "Native-v4 registry requires exactly 112 shards." >&2
  exit 3
fi
for output in "${SOURCE_OUT}" "${REGISTRY_OUT}" "${SUMMARY_OUT}" "${RELEASE_SUMMARY_OUT}"; do
  if [[ -e "${output}" ]]; then
    echo "Refusing to overwrite native-v4 registry artifact: ${output}" >&2
    exit 4
  fi
done

mkdir -p "${TWIRL_ROOT}/logs" "${REGISTRY_DIR}"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

"${PYTHON}" scripts/stage5_validation/build_teacher_ssl_fullpool_native_registry.py \
  --frozen-pool "${POOL_TABLE}" \
  --frozen-pool-summary "${POOL_SUMMARY}" \
  --eligibility-exclusions "${ELIGIBILITY_EXCLUSIONS}" \
  --eligibility-summary "${ELIGIBILITY_SUMMARY}" \
  --native-shards "${native_shards[@]}" \
  --native-shard-summaries "${native_summaries[@]}" \
  --source-out "${SOURCE_OUT}" \
  --registry-out "${REGISTRY_OUT}" \
  --summary-out "${SUMMARY_OUT}" \
  --release-summary-out "${RELEASE_SUMMARY_OUT}"

for output in "${SOURCE_OUT}" "${REGISTRY_OUT}" "${SUMMARY_OUT}" "${RELEASE_SUMMARY_OUT}"; do
  test -s "${output}"
done
echo "[ssl-full-pool-native-registry-v4] output=${REGISTRY_OUT}"
