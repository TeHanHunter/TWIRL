#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-sector-prep
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/teacher_enrichment/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/teacher_enrichment/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_teacher_v2}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
SECTOR="${TWIRL_ENRICH_SECTOR:?set TWIRL_ENRICH_SECTOR}"
DATA_ROOT="${TWIRL_ENRICH_DATA_ROOT:-/orcd/data/mki_aryeh/001/twirl/data_local/stage1_lightcurves/s${SECTOR}_A2v1}"
RAW_H5="${TWIRL_ENRICH_RAW_H5:-${DATA_ROOT}/s${SECTOR}_A2v1_tglc_raw_sources.h5}"
WORK_ROOT="${TWIRL_ENRICH_WORK_ROOT:?set TWIRL_ENRICH_WORK_ROOT}"
N_SHARDS="${TWIRL_ENRICH_BLS_SHARDS:-40}"
SMALL_PEAKS="${TWIRL_ENRICH_SMALL_PEAKS:-3}"

mkdir -p "${WORK_ROOT}/real_bls" "${WORK_ROOT}/candidates"
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"${PYTHON}" scripts/stage5_validation/merge_adp_real_bls_shards.py \
  --shard-dir "${WORK_ROOT}/real_bls/shards" \
  --out-path "${WORK_ROOT}/real_bls/real_adp_bls_peaks.parquet" \
  --n-shards "${N_SHARDS}"
"${PYTHON}" scripts/stage5_validation/prepare_teacher_v2_real_candidates.py \
  --real-bls-peaks "${WORK_ROOT}/real_bls/real_adp_bls_peaks.parquet" \
  --raw-source-h5 "${RAW_H5}" \
  --out-dir "${WORK_ROOT}/candidates" \
  --sector "${SECTOR}" \
  --small-peaks "${SMALL_PEAKS}"
