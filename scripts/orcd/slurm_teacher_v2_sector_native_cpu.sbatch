#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-sector-native
#SBATCH --array=0-31%32
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 24:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/teacher_enrichment/%x-%A_%a.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/teacher_enrichment/%x-%A_%a.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_teacher_v2}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
SECTOR="${TWIRL_ENRICH_SECTOR:?set TWIRL_ENRICH_SECTOR}"
DATA_ROOT="${TWIRL_ENRICH_DATA_ROOT:-/orcd/data/mki_aryeh/001/twirl/data_local/stage1_lightcurves/s${SECTOR}_A2v1}"
RAW_H5="${TWIRL_ENRICH_RAW_H5:-${DATA_ROOT}/s${SECTOR}_A2v1_tglc_raw_sources.h5}"
COMPACT="${TWIRL_ENRICH_COMPACT_H5:-${DATA_ROOT}/s${SECTOR}_A2v1_adp_pair_rebuilt.h5}"
WORK_ROOT="${TWIRL_ENRICH_WORK_ROOT:?set TWIRL_ENRICH_WORK_ROOT}"
N_SHARDS="${TWIRL_ENRICH_NATIVE_SHARDS:-32}"

mkdir -p "${WORK_ROOT}/real_native/shards"
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

"${PYTHON}" scripts/stage5_validation/build_s56_harmonic_native_inputs.py \
  --training-table "${WORK_ROOT}/candidates/real_candidates_top5.parquet" \
  --raw-source-h5 "${RAW_H5}" \
  --compact-adp-h5 "${COMPACT}" \
  --out-h5 "${WORK_ROOT}/real_native/shards/native_${SLURM_ARRAY_TASK_ID}.h5" \
  --repo-root "${REPO}" \
  --n-periods 4096 \
  --shard-index "${SLURM_ARRAY_TASK_ID}" \
  --n-shards "${N_SHARDS}"
