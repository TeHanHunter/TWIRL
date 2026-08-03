#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-ssl-v4-bls-global
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:?TWIRL_ORCD_REPO is required}"
EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA:?TWIRL_EXPECTED_GIT_SHA is required}"
RUN_ROOT="${TWIRL_SSL_FULLPOOL_RUN_ROOT:?TWIRL_SSL_FULLPOOL_RUN_ROOT is required}"
SOURCE_RUN_ROOT="${TWIRL_SSL_FULLPOOL_SOURCE_RUN_ROOT:?TWIRL_SSL_FULLPOOL_SOURCE_RUN_ROOT is required}"
PYTHON="${TWIRL_ORCD_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56/bin/python}"
POOL_SUMMARY="${SOURCE_RUN_ROOT}/frozen/full_pool/teacher_ssl_full_pool_manifest.summary.json"
GLOBAL_DIR="${RUN_ROOT}/real_bls/global"

if [[ "$(git -C "${REPO}" rev-parse HEAD)" != "${EXPECTED_SHA}" ]]; then
  echo "Deployed Git SHA does not match TWIRL_EXPECTED_GIT_SHA" >&2
  exit 3
fi
args=(
  --pool-summary "${POOL_SUMMARY}"
  --out-parquet "${GLOBAL_DIR}/teacher_ssl_full_pool_bls.parquet"
  --out-summary "${GLOBAL_DIR}/teacher_ssl_full_pool_bls.summary.json"
)
for sector in {56..62}; do
  args+=(--sector-bls "${sector}=${RUN_ROOT}/real_bls/s${sector}/real_adp_bls_peaks.parquet")
  args+=(--sector-summary "${sector}=${RUN_ROOT}/real_bls/s${sector}/real_adp_bls_peaks.summary.json")
done
mkdir -p "${TWIRL_ROOT}/logs" "${GLOBAL_DIR}"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
"${PYTHON}" scripts/stage5_validation/consolidate_teacher_ssl_full_pool_bls.py "${args[@]}"
