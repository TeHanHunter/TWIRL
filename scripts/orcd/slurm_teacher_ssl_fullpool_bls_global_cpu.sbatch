#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-ssl-fullpool-bls-global
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
PYTHON="${TWIRL_ORCD_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56/bin/python}"
RUN_ROOT="${TWIRL_SSL_FULLPOOL_RUN_ROOT:-${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp}"
POOL_DIR="${TWIRL_SSL_FULLPOOL_POOL_DIR:-${RUN_ROOT}/frozen/full_pool}"
POOL_SUMMARY="${POOL_DIR}/teacher_ssl_full_pool_manifest.summary.json"
BLS_ROOT="${TWIRL_SSL_FULLPOOL_BLS_ROOT:-${RUN_ROOT}/real_bls}"
GLOBAL_DIR="${TWIRL_SSL_FULLPOOL_GLOBAL_BLS_DIR:-${BLS_ROOT}/global}"
OUT_PARQUET="${GLOBAL_DIR}/teacher_ssl_full_pool_bls.parquet"
OUT_SUMMARY="${GLOBAL_DIR}/teacher_ssl_full_pool_bls.summary.json"

if [[ ! -s "${POOL_SUMMARY}" ]]; then
  echo "Missing frozen full-pool summary: ${POOL_SUMMARY}" >&2
  exit 3
fi

args=(
  --pool-summary "${POOL_SUMMARY}"
  --out-parquet "${OUT_PARQUET}"
  --out-summary "${OUT_SUMMARY}"
)
for sector in {56..62}; do
  sector_table="${BLS_ROOT}/s${sector}/real_adp_bls_peaks.parquet"
  sector_summary="${BLS_ROOT}/s${sector}/real_adp_bls_peaks.summary.json"
  for input in "${sector_table}" "${sector_summary}"; do
    if [[ ! -s "${input}" ]]; then
      echo "Missing merged S${sector} BLS input: ${input}" >&2
      exit 4
    fi
  done
  args+=(--sector-bls "${sector}=${sector_table}")
  args+=(--sector-summary "${sector}=${sector_summary}")
done

mkdir -p "${TWIRL_ROOT}/logs" "${GLOBAL_DIR}"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[ssl-full-pool-bls-global] host=$(hostname) start=$(date -Is)"
"${PYTHON}" \
  scripts/stage5_validation/consolidate_teacher_ssl_full_pool_bls.py \
  "${args[@]}"
echo "[ssl-full-pool-bls-global] complete=$(date -Is) output=${OUT_PARQUET}"
