#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-ssl-fullpool-bls-merge
#SBATCH --array=0-6%7
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.err

set -euo pipefail

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:-${TWIRL_ROOT}/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56/bin/python}"
RUN_ROOT="${TWIRL_SSL_FULLPOOL_RUN_ROOT:-${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp}"
BLS_ROOT="${TWIRL_SSL_FULLPOOL_BLS_ROOT:-${RUN_ROOT}/real_bls}"
N_SHARDS=8

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
if (( TASK_ID < 0 || TASK_ID > 6 )); then
  echo "Merge task ${TASK_ID} is outside the locked 0-6 range" >&2
  exit 3
fi
SECTOR=$((56 + TASK_ID))
SHARD_DIR="${BLS_ROOT}/s${SECTOR}/shards"
OUT_PATH="${BLS_ROOT}/s${SECTOR}/real_adp_bls_peaks.parquet"

mkdir -p "${TWIRL_ROOT}/logs" "$(dirname "${OUT_PATH}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[ssl-full-pool-bls-merge] host=$(hostname) sector=${SECTOR} start=$(date -Is)"
"${PYTHON}" scripts/stage5_validation/merge_adp_real_bls_shards.py \
  --shard-dir "${SHARD_DIR}" \
  --out-path "${OUT_PATH}" \
  --n-shards "${N_SHARDS}"
echo "[ssl-full-pool-bls-merge] sector=${SECTOR} complete=$(date -Is)"
