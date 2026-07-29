#!/bin/bash
# Freeze the exact 112-shard effective-quality-aware native-v3 inventory.
# The controller submits this only after all seven arrays are terminal
# COMPLETED, every H5/summary/SHA product exists, and all stderr files are empty.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-ssl-v3-native-registry
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

readonly TWIRL_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly EXPECTED_SOURCE_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp"
readonly EXPECTED_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v3_s56_s62_a2v1_effective_quality_adp_bls_eligible"
: "${TWIRL_ORCD_REPO:?controller must export TWIRL_ORCD_REPO}"
: "${TWIRL_EXPECTED_GIT_SHA:?controller must export TWIRL_EXPECTED_GIT_SHA}"
: "${TWIRL_SSL_FULLPOOL_V1_RUN_ROOT:?controller must export the immutable raw-v1 source root}"
: "${TWIRL_SSL_FULLPOOL_V3_RUN_ROOT:?controller must export the exact fresh v3 run root}"
: "${TWIRL_FULLPOOL_NATIVE_SHARDS:?controller must export TWIRL_FULLPOOL_NATIVE_SHARDS}"

REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}"
SOURCE_RUN_ROOT="$(realpath -m -- "${TWIRL_SSL_FULLPOOL_V1_RUN_ROOT}")"
RUN_ROOT="$(realpath -m -- "${TWIRL_SSL_FULLPOOL_V3_RUN_ROOT}")"
N_SHARDS="${TWIRL_FULLPOOL_NATIVE_SHARDS}"
PYTHON="${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python"

if [[ "${SOURCE_RUN_ROOT}" != "${EXPECTED_SOURCE_RUN_ROOT}" || "${RUN_ROOT}" != "${EXPECTED_RUN_ROOT}" ]]; then
  echo "Native-v3 registry roots differ from the exact production contract." >&2
  exit 2
fi
if [[ "${N_SHARDS}" != "16" ]]; then
  echo "Native-v3 registry is locked to 16 shards per sector." >&2
  exit 2
fi
if [[ ! "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]]; then
  echo "TWIRL_EXPECTED_GIT_SHA must be a full lowercase Git SHA." >&2
  exit 2
fi
if [[ "$(git -C "${REPO}" rev-parse HEAD)" != "${EXPECTED_SHA}" ]] || [[ -n "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]]; then
  echo "Native-v3 registry repository is not the exact clean deployed revision." >&2
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
  if [[ ! -s "${input}" ]]; then
    echo "Missing native-v3 registry authority: ${input}" >&2
    exit 3
  fi
done

native_shards=()
native_summaries=()
for sector in {56..62}; do
  for shard_index in {0..15}; do
    shard="${NATIVE_ROOT}/s${sector}/native_${shard_index}.h5"
    shard_summary="${NATIVE_ROOT}/s${sector}/native_${shard_index}.summary.json"
    shard_sha="${shard}.sha256"
    for input in "${shard}" "${shard_summary}" "${shard_sha}"; do
      if [[ ! -s "${input}" ]]; then
        echo "Missing native-v3 shard product: ${input}" >&2
        exit 4
      fi
    done
    sha256sum -c "${shard_sha}"
    native_shards+=("${shard}")
    native_summaries+=("${shard_summary}")
  done
done
if [[ "${#native_shards[@]}" -ne 112 || "${#native_summaries[@]}" -ne 112 ]]; then
  echo "Native-v3 registry requires exactly 112 shards and 112 summaries." >&2
  exit 5
fi

for output in "${SOURCE_OUT}" "${REGISTRY_OUT}" "${SUMMARY_OUT}" "${RELEASE_SUMMARY_OUT}"; do
  if [[ -e "${output}" ]]; then
    echo "Refusing to overwrite or resume an existing native-v3 registry artifact: ${output}" >&2
    exit 6
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

echo "[ssl-full-pool-native-registry-v3] host=$(hostname) shards=${#native_shards[@]} start=$(date -Is)"
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
  if [[ ! -s "${output}" ]]; then
    echo "Native-v3 registry producer did not publish: ${output}" >&2
    exit 7
  fi
done
echo "[ssl-full-pool-native-registry-v3] complete=$(date -Is) output=${REGISTRY_OUT}"
