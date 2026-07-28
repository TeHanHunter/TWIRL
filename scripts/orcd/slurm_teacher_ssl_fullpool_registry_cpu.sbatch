#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-ssl-fullpool-registry
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:-${TWIRL_ROOT}/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python}"
RUN_ROOT="${TWIRL_SSL_FULLPOOL_RUN_ROOT:-${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp}"
POOL_DIR="${TWIRL_SSL_FULLPOOL_POOL_DIR:-${RUN_ROOT}/frozen/full_pool}"
POOL_TABLE="${TWIRL_SSL_FULLPOOL_POOL_TABLE:-${POOL_DIR}/teacher_ssl_full_pool_observations.parquet}"
POOL_SUMMARY="${POOL_DIR}/teacher_ssl_full_pool_manifest.summary.json"
BLS_ROOT="${TWIRL_SSL_FULLPOOL_BLS_ROOT:-${RUN_ROOT}/real_bls}"
GLOBAL_BLS_DIR="${TWIRL_SSL_FULLPOOL_GLOBAL_BLS_DIR:-${BLS_ROOT}/global}"
BLS_PEAKS="${TWIRL_SSL_FULLPOOL_GLOBAL_BLS:-${GLOBAL_BLS_DIR}/teacher_ssl_full_pool_bls.parquet}"
BLS_SUMMARY="${TWIRL_SSL_FULLPOOL_GLOBAL_BLS_SUMMARY:-${GLOBAL_BLS_DIR}/teacher_ssl_full_pool_bls.summary.json}"
NATIVE_REGISTRY="${TWIRL_SSL_FULLPOOL_NATIVE_REGISTRY:?set TWIRL_SSL_FULLPOOL_NATIVE_REGISTRY to the external observation-keyed native registry}"
SPLIT_REGISTRY="${TWIRL_SSL_FULLPOOL_SPLIT_REGISTRY:-${TWIRL_ROOT}/reports/stage5_validation/teacher_v3_s56_s62_a2v1_current_adp/frozen/tic_split_registry.csv}"
S63_TICS="${TWIRL_SSL_FULLPOOL_S63_TICS:-${RUN_ROOT}/preregistered/s63_reserved_tics.txt}"
REGISTRY_DIR="${TWIRL_SSL_FULLPOOL_REGISTRY_DIR:-${RUN_ROOT}/frozen/ssl_registry}"
REGISTRY_OUT="${REGISTRY_DIR}/teacher_ssl_fullpool_registry.parquet"
SUMMARY_OUT="${REGISTRY_DIR}/teacher_ssl_fullpool_registry.summary.json"

inputs=(
  "${POOL_TABLE}"
  "${POOL_SUMMARY}"
  "${BLS_PEAKS}"
  "${BLS_SUMMARY}"
  "${NATIVE_REGISTRY}"
  "${SPLIT_REGISTRY}"
  "${S63_TICS}"
)
for input in "${inputs[@]}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing full-pool registry input: ${input}" >&2
    exit 3
  fi
done

mkdir -p "${TWIRL_ROOT}/logs" "${REGISTRY_DIR}"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

echo "[ssl-full-pool-registry] host=$(hostname) start=$(date -Is)"
"${PYTHON}" scripts/stage5_validation/build_teacher_ssl_fullpool_registry.py \
  --frozen-pool "${POOL_TABLE}" \
  --frozen-pool-summary "${POOL_SUMMARY}" \
  --bls-summary "${BLS_SUMMARY}" \
  --bls-peaks "${BLS_PEAKS}" \
  --native-registry "${NATIVE_REGISTRY}" \
  --frozen-split-registry "${SPLIT_REGISTRY}" \
  --reserved-hosts "${S63_TICS}" \
  --registry-out "${REGISTRY_OUT}" \
  --summary-out "${SUMMARY_OUT}"
echo "[ssl-full-pool-registry] complete=$(date -Is) output=${REGISTRY_OUT}"
