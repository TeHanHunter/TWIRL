#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-eb-tables
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
OUT_ROOT="${TWIRL_EB_MINER_ROOT:-reports/stage5_validation/s56_eb_miner_adp_only}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs
cd "${REPO}"

export PYTHONPATH=src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[eb-tables] host=$(hostname)"
echo "[eb-tables] date=$(date -Is)"
echo "[eb-tables] repo=${REPO}"
echo "[eb-tables] python=${PYTHON}"
echo "[eb-tables] out_root=${OUT_ROOT}"

"${PYTHON}" scripts/stage5_validation/build_s56_eb_miner_tables.py \
  --out-root "${OUT_ROOT}"

echo "[eb-tables] complete $(date -Is)"
