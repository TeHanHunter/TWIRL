#!/bin/bash
# Freeze the exact BLS-derived native/model eligibility authority.
#
# This is the first v2-only stage.  It reads the immutable v1 frozen pool and
# global BLS release, but publishes only beneath a wholly fresh v2 run root.
# No native build should start until this job completes successfully and both
# outputs below pass their consumer-side validation.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-ssl-eligibility-v2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:-${TWIRL_ROOT}/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56/bin/python}"
: "${TWIRL_EXPECTED_GIT_SHA:?controller must export TWIRL_EXPECTED_GIT_SHA}"
EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}"
SOURCE_RUN_ROOT="${TWIRL_SSL_FULLPOOL_V1_RUN_ROOT:-${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp}"
RUN_ROOT="${TWIRL_SSL_FULLPOOL_V2_RUN_ROOT:-${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v2_s56_s62_a2v1_current_adp_bls_eligible}"
POOL_DIR="${SOURCE_RUN_ROOT}/frozen/full_pool"
POOL_TABLE="${POOL_DIR}/teacher_ssl_full_pool_observations.csv"
POOL_SUMMARY="${POOL_DIR}/teacher_ssl_full_pool_manifest.summary.json"
GLOBAL_BLS_DIR="${SOURCE_RUN_ROOT}/real_bls/global"
BLS_PEAKS="${GLOBAL_BLS_DIR}/teacher_ssl_full_pool_bls.parquet"
BLS_SUMMARY="${GLOBAL_BLS_DIR}/teacher_ssl_full_pool_bls.summary.json"
ELIGIBILITY_DIR="${RUN_ROOT}/frozen/native_model_eligibility"
EXCLUSIONS_OUT="${ELIGIBILITY_DIR}/native_model_exclusions.csv"
SUMMARY_OUT="${ELIGIBILITY_DIR}/native_model_eligibility.summary.json"

SOURCE_RUN_ROOT="$(realpath -m -- "${SOURCE_RUN_ROOT}")"
RUN_ROOT="$(realpath -m -- "${RUN_ROOT}")"
EXPECTED_SOURCE_RUN_ROOT="$(realpath -m -- "${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp")"
EXPECTED_RUN_ROOT="$(realpath -m -- "${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v2_s56_s62_a2v1_current_adp_bls_eligible")"
if [[ "${SOURCE_RUN_ROOT}" != "${EXPECTED_SOURCE_RUN_ROOT}" || "${RUN_ROOT}" != "${EXPECTED_RUN_ROOT}" ]]; then
  echo "Eligibility source/output roots differ from the production contract." >&2
  exit 2
fi
if [[ ! "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]]; then
  echo "TWIRL_EXPECTED_GIT_SHA must be a full lowercase Git SHA." >&2
  exit 2
fi
if [[ "$(git -C "${REPO}" rev-parse HEAD)" != "${EXPECTED_SHA}" ]] || [[ -n "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]]; then
  echo "Eligibility code checkout differs from the exact clean deployed revision." >&2
  exit 2
fi

for input in "${POOL_TABLE}" "${POOL_SUMMARY}" "${BLS_PEAKS}" "${BLS_SUMMARY}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing native/model eligibility authority: ${input}" >&2
    exit 3
  fi
done

# The pool consumer validates these hashes through POOL_SUMMARY.  Require the
# complete seven-sector inventory here as an additional launch-time fail-safe.
for sector in {56..62}; do
  allowlist="${POOL_DIR}/allowlists/s${sector}_tics.csv"
  if [[ ! -s "${allowlist}" ]]; then
    echo "Missing frozen S${sector} allowlist: ${allowlist}" >&2
    exit 3
  fi
done

for output in "${EXCLUSIONS_OUT}" "${SUMMARY_OUT}"; do
  if [[ -e "${output}" ]]; then
    echo "Refusing to overwrite an existing v2 eligibility artifact: ${output}" >&2
    exit 4
  fi
done

mkdir -p "${TWIRL_ROOT}/logs" "${ELIGIBILITY_DIR}"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[ssl-full-pool-eligibility-v2] host=$(hostname) start=$(date -Is)"
"${PYTHON}" \
  scripts/stage5_validation/build_teacher_ssl_fullpool_native_eligibility.py \
  --frozen-pool "${POOL_TABLE}" \
  --frozen-pool-summary "${POOL_SUMMARY}" \
  --bls-peaks "${BLS_PEAKS}" \
  --bls-summary "${BLS_SUMMARY}" \
  --exclusions-out "${EXCLUSIONS_OUT}" \
  --summary-out "${SUMMARY_OUT}" \
  --production-lock

for output in "${EXCLUSIONS_OUT}" "${SUMMARY_OUT}"; do
  if [[ ! -s "${output}" ]]; then
    echo "Eligibility producer did not publish its required output: ${output}" >&2
    exit 5
  fi
done
echo "[ssl-full-pool-eligibility-v2] complete=$(date -Is) summary=${SUMMARY_OUT}"
