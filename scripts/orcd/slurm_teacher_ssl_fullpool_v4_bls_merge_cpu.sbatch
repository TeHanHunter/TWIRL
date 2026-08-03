#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-ssl-v4-bls-merge
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
REPO="${TWIRL_ORCD_REPO:?TWIRL_ORCD_REPO is required}"
EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA:?TWIRL_EXPECTED_GIT_SHA is required}"
RUN_ROOT="${TWIRL_SSL_FULLPOOL_RUN_ROOT:?TWIRL_SSL_FULLPOOL_RUN_ROOT is required}"
PYTHON="${TWIRL_ORCD_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56/bin/python}"
N_SHARDS=16
SECTOR=$((56 + ${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}))
SHARD_DIR="${RUN_ROOT}/real_bls/s${SECTOR}/shards"
OUT_PATH="${RUN_ROOT}/real_bls/s${SECTOR}/real_adp_bls_peaks.parquet"

if [[ "$(git -C "${REPO}" rev-parse HEAD)" != "${EXPECTED_SHA}" ]]; then
  echo "Deployed Git SHA does not match TWIRL_EXPECTED_GIT_SHA" >&2
  exit 3
fi
mkdir -p "${TWIRL_ROOT}/logs" "$(dirname "${OUT_PATH}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
"${PYTHON}" scripts/stage5_validation/merge_adp_real_bls_shards.py \
  --shard-dir "${SHARD_DIR}" \
  --out-path "${OUT_PATH}" \
  --n-shards "${N_SHARDS}"
