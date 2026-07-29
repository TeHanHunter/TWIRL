#!/bin/bash
# Audit every eligible S60/native_4 row through the exact production
# HarmonicNativeDataset/collator/numerical-envelope path before arrays launch.
# This diagnostic publishes no native registry, numeric release, or model.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-ssl-v3-canary-numeric
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

readonly TWIRL_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly EXPECTED_SOURCE_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp"
readonly EXPECTED_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v3_s56_s62_a2v1_effective_quality_adp_bls_eligible"

: "${TWIRL_ORCD_REPO:?controller must export TWIRL_ORCD_REPO}"
: "${TWIRL_EXPECTED_GIT_SHA:?controller must export TWIRL_EXPECTED_GIT_SHA}"
: "${TWIRL_SSL_FULLPOOL_V1_RUN_ROOT:?controller must export the immutable source root}"
: "${TWIRL_SSL_FULLPOOL_V3_RUN_ROOT:?controller must export the exact fresh v3 run root}"

REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}"
SOURCE_RUN_ROOT="$(realpath -m -- "${TWIRL_SSL_FULLPOOL_V1_RUN_ROOT}")"
RUN_ROOT="$(realpath -m -- "${TWIRL_SSL_FULLPOOL_V3_RUN_ROOT}")"
PYTHON="${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python"

if [[ "${SOURCE_RUN_ROOT}" != "${EXPECTED_SOURCE_RUN_ROOT}" || "${RUN_ROOT}" != "${EXPECTED_RUN_ROOT}" ]]; then
  echo "Native-v3 canary-audit roots differ from the exact contract." >&2
  exit 2
fi
if [[ ! "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]]; then
  echo "TWIRL_EXPECTED_GIT_SHA must be a full lowercase Git SHA." >&2
  exit 2
fi
if [[ "$(git -C "${REPO}" rev-parse HEAD)" != "${EXPECTED_SHA}" ]] || [[ -n "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]]; then
  echo "Native-v3 canary-audit repository is not the exact clean revision." >&2
  exit 2
fi

POOL_DIR="${SOURCE_RUN_ROOT}/frozen/full_pool"
POOL="${POOL_DIR}/teacher_ssl_full_pool_observations.csv"
POOL_SUMMARY="${POOL_DIR}/teacher_ssl_full_pool_manifest.summary.json"
GLOBAL_BLS_DIR="${SOURCE_RUN_ROOT}/real_bls/global"
BLS_PEAKS="${GLOBAL_BLS_DIR}/teacher_ssl_full_pool_bls.parquet"
BLS_SUMMARY="${GLOBAL_BLS_DIR}/teacher_ssl_full_pool_bls.summary.json"
CANARY_ROOT="${RUN_ROOT}/native_canary"
NATIVE_H5="${CANARY_ROOT}/s60/native_4.h5"
NATIVE_SUMMARY="${CANARY_ROOT}/s60/native_4.summary.json"
NATIVE_SHA_SIDECAR="${NATIVE_H5}.sha256"
REPORT_OUT="${CANARY_ROOT}/s60/native_4.model_facing_numeric.json"

for input in \
  "${POOL}" \
  "${POOL_SUMMARY}" \
  "${BLS_PEAKS}" \
  "${BLS_SUMMARY}" \
  "${NATIVE_H5}" \
  "${NATIVE_SUMMARY}" \
  "${NATIVE_SHA_SIDECAR}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing native-v3 canary-audit authority: ${input}" >&2
    exit 3
  fi
done
for output in "${REPORT_OUT}" "${REPORT_OUT}.sha256"; do
  if [[ -e "${output}" ]]; then
    echo "Refusing to overwrite or reuse native-v3 canary audit: ${output}" >&2
    exit 4
  fi
done

mkdir -p "${TWIRL_ROOT}/logs" "$(dirname "${REPORT_OUT}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[teacher-ssl-fullpool-v3-canary-numeric] start=$(date -Is) sector=60 shard=4 tic=704538814"
"${PYTHON}" scripts/orcd/audit_teacher_ssl_fullpool_v3_canary.py \
  --frozen-pool "${POOL}" \
  --frozen-pool-summary "${POOL_SUMMARY}" \
  --bls-peaks "${BLS_PEAKS}" \
  --bls-summary "${BLS_SUMMARY}" \
  --native-h5 "${NATIVE_H5}" \
  --native-summary "${NATIVE_SUMMARY}" \
  --native-sha-sidecar "${NATIVE_SHA_SIDECAR}" \
  --report-out "${REPORT_OUT}" \
  --expected-code-revision "${EXPECTED_SHA}"

test -s "${REPORT_OUT}"
test -s "${REPORT_OUT}.sha256"
(
  cd "$(dirname "${REPORT_OUT}")"
  sha256sum -c "$(basename "${REPORT_OUT}").sha256"
)
echo "[teacher-ssl-fullpool-v3-canary-numeric] complete=$(date -Is) report=${REPORT_OUT}"
