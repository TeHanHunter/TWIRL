#!/bin/bash
# Build one ORCD native shard; submit as an array after raw shards are staged.
# Example: sbatch --array=0-15%8 --export=ALL,TWIRL_FULLPOOL_SECTOR=62,...
#   scripts/orcd/slurm_teacher_ssl_fullpool_native_cpu.sbatch
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-ssl-native
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.err

set -euo pipefail

SECTOR="${TWIRL_FULLPOOL_SECTOR:?set TWIRL_FULLPOOL_SECTOR to 56-62}"
N_SHARDS="${TWIRL_FULLPOOL_NATIVE_SHARDS:?set TWIRL_FULLPOOL_NATIVE_SHARDS}"
SHARD_INDEX="${SLURM_ARRAY_TASK_ID:-${TWIRL_FULLPOOL_SHARD_INDEX:-}}"
if [[ -z "${SHARD_INDEX}" ]]; then
  echo "Submit as a Slurm array or set TWIRL_FULLPOOL_SHARD_INDEX." >&2
  exit 2
fi
if (( SECTOR < 56 || SECTOR > 62 )); then
  echo "Full-pool native preparation is bounded to sectors 56-62." >&2
  exit 2
fi
if (( N_SHARDS < 1 || SHARD_INDEX < 0 || SHARD_INDEX >= N_SHARDS )); then
  echo "Shard index must satisfy 0 <= index < n-shards." >&2
  exit 2
fi

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:-${TWIRL_ROOT}/code/TWIRL}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python}"
RUN_REL="reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp"
RUN="${TWIRL_FULLPOOL_RUN_ROOT:-${TWIRL_ROOT}/${RUN_REL}}"
FROZEN="${RUN}/frozen/full_pool"
POOL="${TWIRL_FULLPOOL_FROZEN_POOL:-${FROZEN}/teacher_ssl_full_pool_observations.parquet}"
POOL_SUMMARY="${TWIRL_FULLPOOL_FROZEN_SUMMARY:-${FROZEN}/teacher_ssl_full_pool_manifest.summary.json}"
ALLOWLIST="${TWIRL_FULLPOOL_SECTOR_ALLOWLIST:-${FROZEN}/allowlists/s${SECTOR}_tics.csv}"
RAW_SOURCE="${TWIRL_FULLPOOL_RAW_SOURCE:-${RUN}/native_preparation/raw_sources/s${SECTOR}/raw_source_${SHARD_INDEX}.h5}"
COMPACT_ADP="${TWIRL_FULLPOOL_COMPACT_ADP_H5:?set the exact S${SECTOR} compact ADP HDF5}"
CADENCE_TABLE="${TWIRL_FULLPOOL_CADENCE_REFERENCE_TABLE:?set the exact S${SECTOR} cadence-reference table}"
CADENCE_MANIFEST="${TWIRL_FULLPOOL_CADENCE_REFERENCE_MANIFEST:?set the exact S${SECTOR} cadence-reference manifest}"
OUT_H5="${TWIRL_FULLPOOL_NATIVE_OUT:-${RUN}/native/s${SECTOR}/native_${SHARD_INDEX}.h5}"
ORBITID_POLICY="strict"
if [[ "${SECTOR}" == "62" ]]; then
  ORBITID_POLICY="reference_by_cadence"
fi

for input in \
  "${POOL}" \
  "${POOL_SUMMARY}" \
  "${ALLOWLIST}" \
  "${RAW_SOURCE}" \
  "${COMPACT_ADP}" \
  "${CADENCE_TABLE}" \
  "${CADENCE_MANIFEST}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing full-pool native input: ${input}" >&2
    exit 3
  fi
done
if [[ -e "${OUT_H5}" || -e "${OUT_H5%.h5}.summary.json" ]]; then
  echo "Refusing to overwrite an existing full-pool native artifact." >&2
  exit 4
fi

mkdir -p "${TWIRL_ROOT}/logs" "$(dirname "${OUT_H5}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

"${PYTHON}" scripts/stage5_validation/build_teacher_ssl_fullpool_native_shard.py \
  --sector "${SECTOR}" \
  --frozen-pool "${POOL}" \
  --frozen-pool-summary "${POOL_SUMMARY}" \
  --sector-allowlist "${ALLOWLIST}" \
  --raw-source-h5 "${RAW_SOURCE}" \
  --compact-adp-h5 "${COMPACT_ADP}" \
  --cadence-reference-table "${CADENCE_TABLE}" \
  --cadence-reference-manifest "${CADENCE_MANIFEST}" \
  --out-h5 "${OUT_H5}" \
  --shard-index "${SHARD_INDEX}" \
  --n-shards "${N_SHARDS}" \
  --n-periods 4096 \
  --orbitid-policy "${ORBITID_POLICY}"

sha256sum "${OUT_H5}" > "${OUT_H5}.sha256"
echo "[fullpool-native] S${SECTOR} shard ${SHARD_INDEX}/${N_SHARDS}: ${OUT_H5}"
