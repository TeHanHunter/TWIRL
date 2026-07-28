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
RUN_ROOT="${TWIRL_SSL_FULLPOOL_RUN_ROOT:-${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp}"
POOL_DIR="${TWIRL_SSL_FULLPOOL_POOL_DIR:-${RUN_ROOT}/frozen/full_pool}"
POOL_TABLE="${TWIRL_SSL_FULLPOOL_POOL_TABLE:-${POOL_DIR}/teacher_ssl_full_pool_observations.parquet}"
POOL_SUMMARY="${POOL_DIR}/teacher_ssl_full_pool_manifest.summary.json"
N_SHARDS="${TWIRL_FULLPOOL_NATIVE_SHARDS:-16}"
NATIVE_ROOT="${TWIRL_FULLPOOL_NATIVE_ROOT:-${RUN_ROOT}/native}"
REGISTRY_DIR="${TWIRL_FULLPOOL_NATIVE_REGISTRY_DIR:-${RUN_ROOT}/frozen/native_registry}"
SOURCE_OUT="${REGISTRY_DIR}/native_registry_source.csv"
REGISTRY_OUT="${REGISTRY_DIR}/native_input_registry.parquet"
SUMMARY_OUT="${REGISTRY_DIR}/native_input_registry.summary.json"

if [[ ! "${N_SHARDS}" =~ ^[0-9]+$ ]] || (( 10#${N_SHARDS} < 1 )); then
  echo "TWIRL_FULLPOOL_NATIVE_SHARDS must be a positive integer." >&2
  exit 2
fi
for input in "${POOL_TABLE}" "${POOL_SUMMARY}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing frozen full-pool authority: ${input}" >&2
    exit 3
  fi
done

native_shards=()
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
  done
done
if (( ${#native_shards[@]} != 7 * 10#${N_SHARDS} )); then
  echo "Native shard inventory count is not exact." >&2
  exit 5
fi

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
  --native-shards "${native_shards[@]}" \
  --source-out "${SOURCE_OUT}" \
  --registry-out "${REGISTRY_OUT}" \
  --summary-out "${SUMMARY_OUT}"
echo "[ssl-full-pool-native-registry] complete=$(date -Is) output=${REGISTRY_OUT}"
