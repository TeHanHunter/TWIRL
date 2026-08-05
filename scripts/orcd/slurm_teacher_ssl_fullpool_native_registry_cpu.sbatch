#!/bin/bash
# Freeze the exact S56--S62 native-shard inventory into one registry.
# Submit only after every sector's native array has completed successfully.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-ssl-native-registry
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:-${TWIRL_ROOT}/code/TWIRL}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python}"
: "${TWIRL_EXPECTED_GIT_SHA:?controller must export TWIRL_EXPECTED_GIT_SHA}"
EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}"
SOURCE_RUN_ROOT="${TWIRL_SSL_FULLPOOL_V1_RUN_ROOT:-${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp}"
RUN_ROOT="${TWIRL_SSL_FULLPOOL_V2_RUN_ROOT:-${TWIRL_SSL_FULLPOOL_RUN_ROOT:-${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v2_s56_s62_a2v1_current_adp_bls_eligible}}"
POOL_DIR="${SOURCE_RUN_ROOT}/frozen/full_pool"
POOL_TABLE="${POOL_DIR}/teacher_ssl_full_pool_observations.csv"
POOL_SUMMARY="${POOL_DIR}/teacher_ssl_full_pool_manifest.summary.json"
ELIGIBILITY_DIR="${RUN_ROOT}/frozen/native_model_eligibility"
ELIGIBILITY_EXCLUSIONS="${ELIGIBILITY_DIR}/native_model_exclusions.csv"
ELIGIBILITY_SUMMARY="${ELIGIBILITY_DIR}/native_model_eligibility.summary.json"
N_SHARDS="${TWIRL_FULLPOOL_NATIVE_SHARDS:-16}"
NATIVE_ROOT="${RUN_ROOT}/native"
REGISTRY_DIR="${RUN_ROOT}/frozen/native_registry"
SOURCE_OUT="${REGISTRY_DIR}/native_registry_source.csv"
REGISTRY_OUT="${REGISTRY_DIR}/native_input_registry.parquet"
SUMMARY_OUT="${REGISTRY_DIR}/native_input_registry.summary.json"
RELEASE_SUMMARY_OUT="${REGISTRY_DIR}/native_fullpool_release.summary.json"

if [[ "${N_SHARDS}" != "16" ]]; then
  echo "The production native release is locked to 16 shards per sector." >&2
  exit 2
fi
SOURCE_RUN_ROOT="$(realpath -m -- "${SOURCE_RUN_ROOT}")"
RUN_ROOT="$(realpath -m -- "${RUN_ROOT}")"
EXPECTED_SOURCE_RUN_ROOT="$(realpath -m -- "${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp")"
EXPECTED_RUN_ROOT="$(realpath -m -- "${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v2_s56_s62_a2v1_current_adp_bls_eligible")"
if [[ "${SOURCE_RUN_ROOT}" != "${EXPECTED_SOURCE_RUN_ROOT}" || "${RUN_ROOT}" != "${EXPECTED_RUN_ROOT}" ]]; then
  echo "Native registry roots differ from the production contract." >&2
  exit 2
fi
if [[ ! "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]]; then
  echo "TWIRL_EXPECTED_GIT_SHA must be a full lowercase Git SHA." >&2
  exit 2
fi
if [[ "$(git -C "${REPO}" rev-parse HEAD)" != "${EXPECTED_SHA}" ]] || [[ -n "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]]; then
  echo "Native-registry code checkout differs from the exact clean deployed revision." >&2
  exit 2
fi
for input in \
  "${POOL_TABLE}" \
  "${POOL_SUMMARY}" \
  "${ELIGIBILITY_EXCLUSIONS}" \
  "${ELIGIBILITY_SUMMARY}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing frozen full-pool authority: ${input}" >&2
    exit 3
  fi
done

native_shards=()
native_summaries=()
for sector in {56..62}; do
  for ((shard_index = 0; shard_index < 10#${N_SHARDS}; shard_index++)); do
    shard="${NATIVE_ROOT}/s${sector}/native_${shard_index}.h5"
    shard_summary="${NATIVE_ROOT}/s${sector}/native_${shard_index}.summary.json"
    for input in "${shard}" "${shard_summary}"; do
      if [[ ! -s "${input}" ]]; then
        echo "Missing full-pool native shard artifact: ${input}" >&2
        exit 4
      fi
    done
    native_shards+=("${shard}")
    native_summaries+=("${shard_summary}")
  done
done
if (( ${#native_shards[@]} != 7 * 10#${N_SHARDS} )); then
  echo "Native shard inventory count is not exact." >&2
  exit 5
fi
if (( ${#native_summaries[@]} != 7 * 10#${N_SHARDS} )); then
  echo "Native shard-summary inventory count is not exact." >&2
  exit 5
fi
for output in \
  "${SOURCE_OUT}" \
  "${REGISTRY_OUT}" \
  "${SUMMARY_OUT}" \
  "${RELEASE_SUMMARY_OUT}"; do
  if [[ -e "${output}" ]]; then
    echo "Refusing to overwrite an existing v2 native registry artifact: ${output}" >&2
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

echo "[ssl-full-pool-native-registry] host=$(hostname) shards=${#native_shards[@]} start=$(date -Is)"
"${PYTHON}" \
  scripts/stage5_validation/build_teacher_ssl_fullpool_native_registry.py \
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
for output in "${REGISTRY_OUT}" "${SUMMARY_OUT}" "${RELEASE_SUMMARY_OUT}"; do
  if [[ ! -s "${output}" ]]; then
    echo "Native registry producer did not publish: ${output}" >&2
    exit 7
  fi
done
echo "[ssl-full-pool-native-registry] complete=$(date -Is) output=${REGISTRY_OUT}"
