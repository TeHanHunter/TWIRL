#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=s56-t2-bls-merge
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_teacher_v2}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
OUT_ROOT="${TWIRL_TEACHER_V2_ROOT:-${REPO}/reports/stage5_validation/s56_A2v1_teacher_v2}"
SOURCE_ROOT="${TWIRL_V1_HUMAN_ROOT:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL/reports/stage5_validation/s56_label_adjudication_real343}"
N_SHARDS="${TWIRL_TEACHER_V2_BLS_SHARDS:-40}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2 \
  "${OUT_ROOT}/real_bls" "${OUT_ROOT}/source_snapshots" "${OUT_ROOT}/human_transfer"
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

cp -p \
  "${SOURCE_ROOT}/adjudicated_training_table/human_vetting_training_table_adjudicated.csv" \
  "${OUT_ROOT}/source_snapshots/human_vetting_training_table_adjudicated.csv"
cp -p \
  "${SOURCE_ROOT}/retrained/s56_harmonic_cnn_v1/training_rows_with_fixed_splits.csv" \
  "${OUT_ROOT}/source_snapshots/training_rows_with_fixed_splits.csv"

"${PYTHON}" scripts/stage5_validation/merge_adp_real_bls_shards.py \
  --shard-dir "${OUT_ROOT}/real_bls/shards" \
  --out-path "${OUT_ROOT}/real_bls/real_adp_bls_peaks.parquet" \
  --n-shards "${N_SHARDS}"

"${PYTHON}" scripts/stage5_validation/build_s56_a2v1_human_transfer.py \
  --human-table "${OUT_ROOT}/source_snapshots/human_vetting_training_table_adjudicated.csv" \
  --real-bls-peaks "${OUT_ROOT}/real_bls/real_adp_bls_peaks.parquet" \
  --out-dir "${OUT_ROOT}/human_transfer" \
  --small-peaks 10

sha256sum \
  "${OUT_ROOT}/source_snapshots/human_vetting_training_table_adjudicated.csv" \
  "${OUT_ROOT}/source_snapshots/training_rows_with_fixed_splits.csv" \
  "${OUT_ROOT}/real_bls/real_adp_bls_peaks.parquet" \
  > "${OUT_ROOT}/source_snapshots/preparation_sha256.txt"
echo "[s56-t2-bls-merge] complete=$(date -Is)"
