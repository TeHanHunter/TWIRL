#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-native-merge
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=96G
#SBATCH -t 06:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/teacher_enrichment/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/teacher_enrichment/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_teacher_v2}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
SECTOR="${TWIRL_ENRICH_SECTOR:?set TWIRL_ENRICH_SECTOR}"
WORK_ROOT="${TWIRL_ENRICH_WORK_ROOT:?set TWIRL_ENRICH_WORK_ROOT}"
N_SHARDS="${TWIRL_ENRICH_NATIVE_SHARDS:-32}"

cd "${REPO}"
export PYTHONPATH=src
export HDF5_USE_FILE_LOCKING=FALSE
shards=()
for ((index=0; index<N_SHARDS; index++)); do
  shards+=("${WORK_ROOT}/real_native/shards/native_${index}.h5")
done
"${PYTHON}" scripts/stage5_validation/merge_teacher_v2_native_shards.py \
  --shards "${shards[@]}" \
  --out-h5 "${WORK_ROOT}/real_native/s${SECTOR}_a2v1_real_native.h5" \
  --candidate-table "${WORK_ROOT}/candidates/real_candidates_top5.parquet"
