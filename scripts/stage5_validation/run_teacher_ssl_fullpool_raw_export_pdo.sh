#!/usr/bin/env bash
# Export one frozen full-pool sector/shard from the user-owned PDO TGLC tree.
set -euo pipefail

SECTOR="${1:?usage: $0 <sector> <orbit-one> <orbit-two> <shard-index> <n-shards>}"
ORBIT_ONE="${2:?usage: $0 <sector> <orbit-one> <orbit-two> <shard-index> <n-shards>}"
ORBIT_TWO="${3:?usage: $0 <sector> <orbit-one> <orbit-two> <shard-index> <n-shards>}"
SHARD_INDEX="${4:?usage: $0 <sector> <orbit-one> <orbit-two> <shard-index> <n-shards>}"
N_SHARDS="${5:?usage: $0 <sector> <orbit-one> <orbit-two> <shard-index> <n-shards>}"

if (( SECTOR < 56 || SECTOR > 62 )); then
  echo "Full-pool raw export is bounded to sectors 56-62." >&2
  exit 2
fi
if (( ORBIT_ONE <= 0 || ORBIT_TWO <= 0 || ORBIT_ONE == ORBIT_TWO )); then
  echo "Orbit IDs must be distinct positive integers." >&2
  exit 2
fi
if (( N_SHARDS < 1 || SHARD_INDEX < 0 || SHARD_INDEX >= N_SHARDS )); then
  echo "Shard index must satisfy 0 <= index < n-shards." >&2
  exit 2
fi

REPO="${TWIRL_PDO_REPO:-/pdo/users/tehan/TWIRL}"
PYTHON="${TWIRL_PDO_PYTHON:-/sw/qlp-environment/.venv/bin/python}"
STAGE_ROOT="${TWIRL_FULLPOOL_PDO_STAGE_ROOT:-/pdo/users/tehan/twirl-teacher-ssl-fullpool}"
FROZEN_POOL="${TWIRL_FULLPOOL_FROZEN_POOL:-${STAGE_ROOT}/frozen/teacher_ssl_full_pool_observations.csv}"
FROZEN_SUMMARY="${TWIRL_FULLPOOL_FROZEN_SUMMARY:-${STAGE_ROOT}/frozen/teacher_ssl_full_pool_manifest.summary.json}"
ALLOWLIST="${TWIRL_FULLPOOL_SECTOR_ALLOWLIST:-${STAGE_ROOT}/frozen/allowlists/s${SECTOR}_tics.csv}"
RAW_ROOT="${TWIRL_FULLPOOL_RAW_ROOT:-/pdo/users/tehan/tglc-gpu-production-A2v1}"
OUT_H5="${TWIRL_FULLPOOL_RAW_OUT:-${STAGE_ROOT}/raw_sources/s${SECTOR}/raw_source_${SHARD_INDEX}.h5}"

if [[ "${OUT_H5}" != /pdo/users/tehan/* ]]; then
  echo "Refusing a PDO output outside /pdo/users/tehan/: ${OUT_H5}" >&2
  exit 2
fi
for input in "${FROZEN_POOL}" "${FROZEN_SUMMARY}" "${ALLOWLIST}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing staged full-pool authority input: ${input}" >&2
    exit 3
  fi
done
for orbit in "${ORBIT_ONE}" "${ORBIT_TWO}"; do
  if [[ ! -d "${RAW_ROOT}/orbit-${orbit}/ffi" ]]; then
    echo "Missing PDO TGLC orbit root: ${RAW_ROOT}/orbit-${orbit}/ffi" >&2
    exit 3
  fi
done
if [[ -e "${OUT_H5}" || -e "${OUT_H5%.h5}.summary.json" ]]; then
  echo "Refusing to overwrite an existing full-pool raw artifact." >&2
  exit 4
fi

mkdir -p "$(dirname "${OUT_H5}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export HDF5_USE_FILE_LOCKING="${HDF5_USE_FILE_LOCKING:-FALSE}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

"${PYTHON}" scripts/stage5_validation/export_teacher_ssl_fullpool_raw_sources.py \
  --sector "${SECTOR}" \
  --frozen-pool "${FROZEN_POOL}" \
  --frozen-pool-summary "${FROZEN_SUMMARY}" \
  --sector-allowlist "${ALLOWLIST}" \
  --raw-root "${RAW_ROOT}" \
  --orbits "${ORBIT_ONE}" "${ORBIT_TWO}" \
  --out-h5 "${OUT_H5}" \
  --shard-index "${SHARD_INDEX}" \
  --n-shards "${N_SHARDS}"

sha256sum "${OUT_H5}" > "${OUT_H5}.sha256"
echo "[fullpool-raw] S${SECTOR} shard ${SHARD_INDEX}/${N_SHARDS}: ${OUT_H5}"
