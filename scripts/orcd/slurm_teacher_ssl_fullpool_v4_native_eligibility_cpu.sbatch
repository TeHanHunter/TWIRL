#!/bin/bash
# Freeze the exact detector-consistent BLS-derived native/model partition.
# First publish a separate rederived audit authority, then require the corrected
# BLS complement to satisfy the preregistered production partition before the
# final eligibility authority is published.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-ssl-v4-eligibility
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:?TWIRL_ORCD_REPO is required}"
EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA:?TWIRL_EXPECTED_GIT_SHA is required}"
SOURCE_RUN_ROOT="${TWIRL_SSL_FULLPOOL_SOURCE_RUN_ROOT:-${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp}"
RUN_ROOT="${TWIRL_SSL_FULLPOOL_RUN_ROOT:?TWIRL_SSL_FULLPOOL_RUN_ROOT is required}"
EXPECTED_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v4r2_s56_s62_a2v1_raw_v1_detector_consistent_adp_bls_eligible"
PYTHON="${TWIRL_ORCD_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56/bin/python}"

if [[ "$(realpath -m -- "${RUN_ROOT}")" != "${EXPECTED_RUN_ROOT}" ]]; then
  echo "Eligibility output root differs from the exact v4 contract." >&2
  exit 2
fi
if [[ ! "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]] || [[ "$(git -C "${REPO}" rev-parse HEAD)" != "${EXPECTED_SHA}" ]] || [[ -n "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]]; then
  echo "Eligibility repository is not the exact clean deployed revision." >&2
  exit 2
fi

POOL_DIR="${SOURCE_RUN_ROOT}/frozen/full_pool"
POOL_TABLE="${POOL_DIR}/teacher_ssl_full_pool_observations.csv"
POOL_SUMMARY="${POOL_DIR}/teacher_ssl_full_pool_manifest.summary.json"
GLOBAL_BLS_DIR="${RUN_ROOT}/real_bls/global"
BLS_PEAKS="${GLOBAL_BLS_DIR}/teacher_ssl_full_pool_bls.parquet"
BLS_SUMMARY="${GLOBAL_BLS_DIR}/teacher_ssl_full_pool_bls.summary.json"
ELIGIBILITY_DIR="${RUN_ROOT}/frozen/native_model_eligibility"
EXCLUSIONS_OUT="${ELIGIBILITY_DIR}/native_model_exclusions.csv"
SUMMARY_OUT="${ELIGIBILITY_DIR}/native_model_eligibility.summary.json"
REDERIVED_DIR="${RUN_ROOT}/frozen/native_model_eligibility_rederived"
REDERIVED_EXCLUSIONS="${REDERIVED_DIR}/native_model_exclusions.csv"
REDERIVED_SUMMARY="${REDERIVED_DIR}/native_model_eligibility.summary.json"

for input in "${POOL_TABLE}" "${POOL_SUMMARY}" "${BLS_PEAKS}" "${BLS_SUMMARY}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing corrected native/model eligibility authority: ${input}" >&2
    exit 3
  fi
done
for sector in {56..62}; do
  test -s "${POOL_DIR}/allowlists/s${sector}_tics.csv"
done
for output in \
  "${REDERIVED_EXCLUSIONS}" \
  "${REDERIVED_SUMMARY}" \
  "${EXCLUSIONS_OUT}" \
  "${SUMMARY_OUT}"; do
  if [[ -e "${output}" ]]; then
    echo "Refusing to overwrite an existing v4 eligibility artifact: ${output}" >&2
    exit 4
  fi
done

mkdir -p "${TWIRL_ROOT}/logs" "${REDERIVED_DIR}" "${ELIGIBILITY_DIR}"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

"${PYTHON}" scripts/stage5_validation/build_teacher_ssl_fullpool_native_eligibility.py \
  --frozen-pool "${POOL_TABLE}" \
  --frozen-pool-summary "${POOL_SUMMARY}" \
  --bls-peaks "${BLS_PEAKS}" \
  --bls-summary "${BLS_SUMMARY}" \
  --exclusions-out "${REDERIVED_EXCLUSIONS}" \
  --summary-out "${REDERIVED_SUMMARY}"

"${PYTHON}" scripts/stage5_validation/build_teacher_ssl_fullpool_native_eligibility.py \
  --frozen-pool "${POOL_TABLE}" \
  --frozen-pool-summary "${POOL_SUMMARY}" \
  --bls-peaks "${BLS_PEAKS}" \
  --bls-summary "${BLS_SUMMARY}" \
  --exclusions-out "${EXCLUSIONS_OUT}" \
  --summary-out "${SUMMARY_OUT}" \
  --production-lock

test -s "${REDERIVED_EXCLUSIONS}"
test -s "${REDERIVED_SUMMARY}"
test -s "${EXCLUSIONS_OUT}"
test -s "${SUMMARY_OUT}"
echo "[ssl-full-pool-eligibility-v4] summary=${SUMMARY_OUT}"
