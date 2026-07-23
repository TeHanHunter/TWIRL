#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-sector-bls
#SBATCH --array=0-39%32
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/teacher_enrichment/%x-%A_%a.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/teacher_enrichment/%x-%A_%a.err

set -euo pipefail

if [[ "${TWIRL_ALLOW_EXPLORATORY_UNGATED_SECTOR_CHAIN:-0}" != "1" ]]; then
  echo "Blocked: generic sector Teacher-v2 preparation lacks the accepted Tier-1 candidate/native provenance chain." >&2
  echo "Use the quality-gated production launcher; set TWIRL_ALLOW_EXPLORATORY_UNGATED_SECTOR_CHAIN=1 only for quarantined diagnostics." >&2
  exit 64
fi

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_teacher_v2}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
SECTOR="${TWIRL_ENRICH_SECTOR:?set TWIRL_ENRICH_SECTOR}"
DATA_ROOT="${TWIRL_ENRICH_DATA_ROOT:-/orcd/data/mki_aryeh/001/twirl/data_local/stage1_lightcurves/s${SECTOR}_A2v1}"
COMPACT="${TWIRL_ENRICH_COMPACT_H5:-${DATA_ROOT}/s${SECTOR}_A2v1_adp_pair_rebuilt.h5}"
CADENCE_REFERENCE="${TWIRL_ENRICH_CADENCE_REFERENCE_TABLE:?set TWIRL_ENRICH_CADENCE_REFERENCE_TABLE}"
CADENCE_REFERENCE_MANIFEST="${TWIRL_ENRICH_CADENCE_REFERENCE_MANIFEST:?set TWIRL_ENRICH_CADENCE_REFERENCE_MANIFEST}"
WORK_ROOT="${TWIRL_ENRICH_WORK_ROOT:?set TWIRL_ENRICH_WORK_ROOT}"
N_SHARDS="${TWIRL_ENRICH_BLS_SHARDS:-40}"
N_PERIODS="${TWIRL_ENRICH_BLS_N_PERIODS:-50000}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs/teacher_enrichment \
  "${WORK_ROOT}/real_bls/shards"
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

"${PYTHON}" scripts/stage5_validation/build_s56_adp_real_bls_peaks.py \
  --compact-lc "${COMPACT}" \
  --sector "${SECTOR}" \
  --cadence-reference "${CADENCE_REFERENCE}" \
  --cadence-reference-manifest "${CADENCE_REFERENCE_MANIFEST}" \
  --out-dir "${WORK_ROOT}/real_bls/shards" \
  --workers "${SLURM_CPUS_PER_TASK:-8}" \
  --n-periods "${N_PERIODS}" \
  --n-peaks 10 \
  --shard-index "${SLURM_ARRAY_TASK_ID}" \
  --n-shards "${N_SHARDS}" \
  --progress-every 25 \
  --source-product-tag A2v1 \
  --resume
