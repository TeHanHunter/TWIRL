#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=s56-a2rec-native
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=48G
#SBATCH -t 1-00:00:00
#SBATCH --array=0-39%40
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/s56_a2v1_recovery/%x-%A_%a.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/s56_a2v1_recovery/%x-%A_%a.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
SECTOR="${TWIRL_RECOVERY_SECTOR:-56}"
SECTOR_TAG="s${SECTOR}"
LOG_ROOT="${TWIRL_RECOVERY_LOG_ROOT:-/orcd/data/mki_aryeh/001/twirl/logs/${SECTOR_TAG}_a2v1_recovery}"
INPUT_ROOT="${TWIRL_A2V1_INPUT_ROOT:-/orcd/data/mki_aryeh/001/twirl/data_local/stage1_lightcurves/${SECTOR_TAG}_A2v1}"
OUT_ROOT="${TWIRL_RECOVERY_ROOT:-${REPO}/reports/stage5_validation/${SECTOR_TAG}_A2v1_teacher_v1_recovery}"
ADP_H5="${TWIRL_A2V1_ADP_H5:-${INPUT_ROOT}/${SECTOR_TAG}_A2v1_adp_pair_rebuilt.h5}"
SHARD="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
TAG="$(printf '%02d' "${SHARD}")"
FULL="${OUT_ROOT}/full"

mkdir -p "${LOG_ROOT}" "${FULL}/teacher_input_shards"
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE

"${PYTHON}" scripts/stage5_validation/prepare_s56_a2v1_recovery_teacher_shard.py \
  --peak-table "${FULL}/peak_shards/peaks_${TAG}.parquet" \
  --pair-h5 "${FULL}/injection_shards/injections_${TAG}.h5" \
  --raw-h5 "${INPUT_ROOT}/${SECTOR_TAG}_A2v1_tglc_raw_sources.h5" \
  --adp-h5 "${ADP_H5}" \
  --out-candidates "${FULL}/teacher_input_shards/candidates_${TAG}.csv" \
  --out-native-h5 "${FULL}/teacher_input_shards/native_${TAG}.h5" \
  --overwrite

echo "[${SECTOR_TAG}-a2rec-native] shard=${SHARD} complete $(date -Is)"
