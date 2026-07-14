#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=s56-t2-real-bls
#SBATCH --array=0-39%16
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2/%x-%A_%a.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2/%x-%A_%a.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_teacher_v2}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
DATA_ROOT="${TWIRL_A2V1_S56_ROOT:-/orcd/data/mki_aryeh/001/twirl/data_local/stage1_lightcurves/s56_A2v1}"
COMPACT="${TWIRL_A2V1_S56_COMPACT:-${DATA_ROOT}/s56_A2v1_adp_pair_rebuilt.h5}"
OUT_ROOT="${TWIRL_TEACHER_V2_ROOT:-${REPO}/reports/stage5_validation/s56_A2v1_teacher_v2}"
N_SHARDS="${TWIRL_TEACHER_V2_BLS_SHARDS:-40}"
SHARD_INDEX="${SLURM_ARRAY_TASK_ID}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2 "${OUT_ROOT}/real_bls/shards"
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

echo "[s56-t2-real-bls] host=$(hostname) shard=${SHARD_INDEX}/${N_SHARDS} start=$(date -Is)"
"${PYTHON}" scripts/stage5_validation/build_s56_adp_real_bls_peaks.py \
  --compact-lc "${COMPACT}" \
  --out-dir "${OUT_ROOT}/real_bls/shards" \
  --workers "${SLURM_CPUS_PER_TASK:-8}" \
  --n-periods 200000 \
  --n-peaks 10 \
  --shard-index "${SHARD_INDEX}" \
  --n-shards "${N_SHARDS}" \
  --progress-every 25 \
  --resume
echo "[s56-t2-real-bls] shard=${SHARD_INDEX} complete=$(date -Is)"
