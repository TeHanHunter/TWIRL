#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-tv3-s63-native
#SBATCH --array=0-31%16
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=48G
#SBATCH -t 24:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.err

set -euo pipefail
umask 077

AUTHORIZATION_CONTRACT="s63_teacher_v3_prospective_v1"
if [[ "${TWIRL_TEACHER_V3_PROSPECTIVE_CONTRACT:-}" != "${AUTHORIZATION_CONTRACT}" ]]; then
  echo "S63 native preparation requires the prospective authorization contract." >&2
  exit 64
fi
TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:?set the clean pinned TWIRL checkout}"
EXPECTED_GIT_SHA="${TWIRL_EXPECTED_GIT_SHA:?set the full deployed Git SHA}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python}"
CANDIDATES="${TWIRL_S63_CANDIDATES:?set exact S63 rank-one candidates}"
EXPECTED_CANDIDATES_SHA="${TWIRL_S63_CANDIDATES_SHA256:?set candidate SHA-256}"
RAW_SOURCE="${TWIRL_S63_RAW_SOURCE_H5:?set staged S63 raw-source HDF5}"
EXPECTED_RAW_SHA="${TWIRL_S63_RAW_SOURCE_H5_SHA256:?set raw-source SHA-256}"
COMPACT="${TWIRL_S63_COMPACT_H5:?set staged S63 compact ADP HDF5}"
EXPECTED_COMPACT_SHA="${TWIRL_S63_COMPACT_H5_SHA256:?set compact SHA-256}"
CADENCE="${TWIRL_S63_CADENCE_TABLE:?set staged S63 cadence table}"
EXPECTED_CADENCE_SHA="${TWIRL_S63_CADENCE_TABLE_SHA256:?set cadence-table SHA-256}"
CADENCE_MANIFEST="${TWIRL_S63_CADENCE_MANIFEST:?set staged S63 cadence manifest}"
EXPECTED_CADENCE_MANIFEST_SHA="${TWIRL_S63_CADENCE_MANIFEST_SHA256:?set cadence-manifest SHA-256}"
OUT_DIR="${TWIRL_S63_NATIVE_SHARD_DIR:?set fresh S63 native shard directory}"
SHARD_INDEX="${SLURM_ARRAY_TASK_ID:?S63 native must run as the locked 0-31 array}"
N_SHARDS=32

if (( SHARD_INDEX < 0 || SHARD_INDEX >= N_SHARDS )); then
  echo "S63 native array index is outside 0-31." >&2
  exit 2
fi
TWIRL_ROOT="$(realpath -e -- "${TWIRL_ROOT}")"
REPO="$(realpath -e -- "${REPO}")"
OUT_DIR="$(realpath -m -- "${OUT_DIR}")"
TWIRL_EXPECTED_GIT_SHA="${EXPECTED_GIT_SHA}" \
  "${REPO}/scripts/assert_clean_checkout.sh" "${REPO}" >/dev/null
if [[ "${OUT_DIR}" != "${TWIRL_ROOT}"/* \
   || "${OUT_DIR}" == "${REPO}" \
   || "${OUT_DIR}" == "${REPO}"/* ]]; then
  echo "S63 native output must be canonical, under TWIRL_ROOT, and outside the checkout." >&2
  exit 3
fi
sha256_check() {
  local path=$1
  local expected=$2
  local observed remainder
  if [[ ! -s "${path}" ]]; then
    echo "Missing S63 native input: ${path}" >&2
    exit 4
  fi
  read -r observed remainder < <(sha256sum "${path}")
  if [[ "${observed}" != "${expected}" ]]; then
    echo "S63 native input SHA-256 mismatch: ${path}" >&2
    exit 5
  fi
}
sha256_check "${CANDIDATES}" "${EXPECTED_CANDIDATES_SHA}"
sha256_check "${RAW_SOURCE}" "${EXPECTED_RAW_SHA}"
sha256_check "${COMPACT}" "${EXPECTED_COMPACT_SHA}"
sha256_check "${CADENCE}" "${EXPECTED_CADENCE_SHA}"
sha256_check "${CADENCE_MANIFEST}" "${EXPECTED_CADENCE_MANIFEST_SHA}"

printf -v suffix '%03d' "${SHARD_INDEX}"
OUTPUT="${OUT_DIR}/native_${suffix}.h5"
SUMMARY="${OUT_DIR}/native_${suffix}.summary.json"
if [[ -e "${OUTPUT}" || -e "${SUMMARY}" ]]; then
  echo "Refusing to overwrite S63 native shard ${SHARD_INDEX}." >&2
  exit 6
fi

CLAIM_ROOT="${OUT_DIR}/.claims"
PENDING_ROOT="${OUT_DIR}/.pending"
CLAIM="${CLAIM_ROOT}/native_${suffix}.claim"
PENDING="${PENDING_ROOT}/native_${suffix}"
PENDING_OUTPUT="${PENDING}/native_${suffix}.h5"
PENDING_SUMMARY="${PENDING}/native_${suffix}.summary.json"
PENDING_SHA="${PENDING}/native_${suffix}.h5.sha256"
mkdir -p "${TWIRL_ROOT}/logs" "${OUT_DIR}" "${CLAIM_ROOT}" "${PENDING_ROOT}"
if ! mkdir "${CLAIM}" 2>/dev/null; then
  echo "S63 native shard is already claimed/interrupted; audit: ${CLAIM}" >&2
  exit 7
fi
if [[ -e "${PENDING}" ]]; then
  echo "S63 native shard has staged interrupted state; audit: ${PENDING}" >&2
  exit 8
fi
mkdir "${PENDING}"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export PYTHONNOUSERSITE=1
unset LD_LIBRARY_PATH
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
"${PYTHON}" -c 'import astropy, h5py, numpy, pandas, torch; print("S63 ORCD native environment OK")'

echo "[teacher-v3-s63-native] shard=${SHARD_INDEX}/${N_SHARDS} start=$(date -Is)"
"${PYTHON}" scripts/stage5_validation/build_teacher_v3_native_inputs.py \
  --sector 63 \
  --prospective-contract "${AUTHORIZATION_CONTRACT}" \
  --training-table "${CANDIDATES}" \
  --raw-source-h5 "${RAW_SOURCE}" \
  --compact-adp-h5 "${COMPACT}" \
  --cadence-reference-table "${CADENCE}" \
  --cadence-reference-manifest "${CADENCE_MANIFEST}" \
  --out-h5 "${PENDING_OUTPUT}" \
  --repo-root "${REPO}" \
  --n-periods 4096 \
  --shard-index "${SHARD_INDEX}" \
  --n-shards "${N_SHARDS}" \
  --orbitid-policy strict \
  --producer-git-sha "${EXPECTED_GIT_SHA}"
if [[ ! -s "${PENDING_OUTPUT}" || ! -s "${PENDING_SUMMARY}" ]]; then
  echo "S63 native producer did not create a complete staged shard." >&2
  exit 9
fi
sha256sum "${PENDING_OUTPUT}" > "${PENDING_SHA}"
mv -- "${PENDING_OUTPUT}" "${OUTPUT}"
mv -- "${PENDING_SHA}" "${OUTPUT}.sha256"
mv -- "${PENDING_SUMMARY}" "${SUMMARY}"
rmdir "${PENDING}" "${CLAIM}"
echo "[teacher-v3-s63-native] shard=${SHARD_INDEX} complete=$(date -Is)"
