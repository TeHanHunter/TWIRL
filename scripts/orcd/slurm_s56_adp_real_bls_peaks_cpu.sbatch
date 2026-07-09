#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-adp-real-bls
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH -t 08:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
COMPACT_LC="${TWIRL_ADP_REAL_BLS_COMPACT_LC:-data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp_lc_export_pdo.h5}"
OUT_DIR="${TWIRL_ADP_REAL_BLS_OUT_DIR:-reports/stage5_validation/s56_adp_real_bls_peaks}"
WORKERS="${TWIRL_ADP_REAL_BLS_WORKERS:-${SLURM_CPUS_PER_TASK:-16}}"
N_PERIODS="${TWIRL_ADP_REAL_BLS_N_PERIODS:-50000}"
N_PEAKS="${TWIRL_ADP_REAL_BLS_N_PEAKS:-10}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[adp-real-bls] host=$(hostname)"
echo "[adp-real-bls] date=$(date -Is)"
echo "[adp-real-bls] compact_lc=${COMPACT_LC}"
echo "[adp-real-bls] out_dir=${OUT_DIR}"
echo "[adp-real-bls] workers=${WORKERS} n_periods=${N_PERIODS} n_peaks=${N_PEAKS}"

"${PYTHON}" scripts/stage5_validation/build_s56_adp_real_bls_peaks.py \
  --compact-lc "${COMPACT_LC}" \
  --out-dir "${OUT_DIR}" \
  --workers "${WORKERS}" \
  --n-periods "${N_PERIODS}" \
  --n-peaks "${N_PEAKS}"

echo "[adp-real-bls] complete $(date -Is)"
