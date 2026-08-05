#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=s56-t2-native
#SBATCH --array=0-31%16
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 24:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2/%x-%A_%a.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2/%x-%A_%a.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_teacher_v2}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
DATA_ROOT="${TWIRL_A2V1_S56_ROOT:-/orcd/data/mki_aryeh/001/twirl/data_local/stage1_lightcurves/s56_A2v1}"
OUT_ROOT="${TWIRL_TEACHER_V2_ROOT:-${REPO}/reports/stage5_validation/s56_A2v1_teacher_v2}"
N_SHARDS="${TWIRL_TEACHER_V2_NATIVE_SHARDS:-32}"
SHARD_INDEX="${SLURM_ARRAY_TASK_ID}"
CADENCE_REFERENCE_TABLE="${TWIRL_A2V1_CADENCE_REFERENCE:?set TWIRL_A2V1_CADENCE_REFERENCE}"
CADENCE_REFERENCE_MANIFEST="${TWIRL_A2V1_CADENCE_REFERENCE_MANIFEST:?set TWIRL_A2V1_CADENCE_REFERENCE_MANIFEST}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2 "${OUT_ROOT}/real_native/shards"
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

"${PYTHON}" scripts/stage5_validation/build_s56_harmonic_native_inputs.py \
  --training-table "${OUT_ROOT}/human_transfer/s56_a2v1_real_candidates_top5.parquet" \
  --raw-source-h5 "${DATA_ROOT}/s56_A2v1_tglc_raw_sources.h5" \
  --compact-adp-h5 "${DATA_ROOT}/s56_A2v1_adp_pair_rebuilt.h5" \
  --cadence-reference-table "${CADENCE_REFERENCE_TABLE}" \
  --cadence-reference-manifest "${CADENCE_REFERENCE_MANIFEST}" \
  --out-h5 "${OUT_ROOT}/real_native/shards/native_${SHARD_INDEX}.h5" \
  --repo-root "${REPO}" \
  --n-periods 4096 \
  --shard-index "${SHARD_INDEX}" \
  --n-shards "${N_SHARDS}"
echo "[s56-t2-native] shard=${SHARD_INDEX}/${N_SHARDS} complete=$(date -Is)"
