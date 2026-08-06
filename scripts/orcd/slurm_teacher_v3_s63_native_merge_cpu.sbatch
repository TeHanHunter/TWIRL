#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-tv3-s63-native-merge
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 077

AUTHORIZATION_CONTRACT="s63_teacher_v3_prospective_v1"
if [[ "${TWIRL_TEACHER_V3_PROSPECTIVE_CONTRACT:-}" != "${AUTHORIZATION_CONTRACT}" ]]; then
  echo "S63 native merge requires the prospective authorization contract." >&2
  exit 64
fi
TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:?set the clean pinned TWIRL checkout}"
EXPECTED_GIT_SHA="${TWIRL_EXPECTED_GIT_SHA:?set the full deployed Git SHA}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python}"
CANDIDATES="${TWIRL_S63_CANDIDATES:?set exact S63 rank-one candidates}"
EXPECTED_CANDIDATES_SHA="${TWIRL_S63_CANDIDATES_SHA256:?set candidate SHA-256}"
SHARD_DIR="${TWIRL_S63_NATIVE_SHARD_DIR:?set completed S63 native shard directory}"
OUTPUT="${TWIRL_S63_NATIVE_H5:?set fresh merged S63 native HDF5}"
N_SHARDS=32

TWIRL_ROOT="$(realpath -e -- "${TWIRL_ROOT}")"
REPO="$(realpath -e -- "${REPO}")"
SHARD_DIR="$(realpath -e -- "${SHARD_DIR}")"
OUTPUT="$(realpath -m -- "${OUTPUT}")"
TWIRL_EXPECTED_GIT_SHA="${EXPECTED_GIT_SHA}" \
  "${REPO}/scripts/assert_clean_checkout.sh" "${REPO}" >/dev/null
if [[ "${OUTPUT}" != "${TWIRL_ROOT}"/* \
   || "${OUTPUT}" == "${REPO}" \
   || "${OUTPUT}" == "${REPO}"/* ]]; then
  echo "S63 native merged output must be canonical, under TWIRL_ROOT, and outside the checkout." >&2
  exit 3
fi
read -r observed_candidates candidate_hash_remainder < <(sha256sum "${CANDIDATES}")
if [[ "${observed_candidates}" != "${EXPECTED_CANDIDATES_SHA}" ]]; then
  echo "S63 native merge candidate SHA-256 mismatch." >&2
  exit 4
fi
if [[ -e "${OUTPUT}" || -e "${OUTPUT%.h5}.summary.json" || -e "${OUTPUT}.merge_pending" ]]; then
  echo "Refusing to overwrite/presume a merged S63 native output." >&2
  exit 5
fi
shards=()
for ((shard=0; shard<N_SHARDS; shard++)); do
  printf -v suffix '%03d' "${shard}"
  path="${SHARD_DIR}/native_${suffix}.h5"
  summary="${SHARD_DIR}/native_${suffix}.summary.json"
  if [[ ! -s "${path}" || ! -s "${summary}" ]]; then
    echo "Missing S63 native shard/summary for index ${shard}." >&2
    exit 6
  fi
  shards+=("${path}")
done

OUTPUT_PARENT="$(dirname "${OUTPUT}")"
CLAIM_ROOT="${OUTPUT_PARENT}/.claims"
PENDING_ROOT="${OUTPUT_PARENT}/.pending"
CLAIM="${CLAIM_ROOT}/$(basename "${OUTPUT}").claim"
PENDING_DIR="${PENDING_ROOT}/$(basename "${OUTPUT}" .h5)"
PENDING_OUTPUT="${PENDING_DIR}/$(basename "${OUTPUT}")"
PENDING_SUMMARY="${PENDING_OUTPUT%.h5}.summary.json"
PENDING_SHA="${PENDING_OUTPUT}.sha256"
mkdir -p "${TWIRL_ROOT}/logs" "${OUTPUT_PARENT}" "${CLAIM_ROOT}" "${PENDING_ROOT}"
if ! mkdir "${CLAIM}" 2>/dev/null; then
  echo "S63 native merge is already claimed/interrupted; audit: ${CLAIM}" >&2
  exit 7
fi
if [[ -e "${PENDING_DIR}" ]]; then
  echo "S63 native merge has staged interrupted state; audit: ${PENDING_DIR}" >&2
  exit 8
fi
mkdir "${PENDING_DIR}"
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
"${PYTHON}" -c 'import h5py, numpy, pandas, torch; print("S63 ORCD native merge environment OK")'
"${PYTHON}" scripts/stage5_validation/merge_teacher_v3_native_shards.py \
  --sector 63 \
  --prospective-contract "${AUTHORIZATION_CONTRACT}" \
  --training-table "${CANDIDATES}" \
  --shards "${shards[@]}" \
  --out-h5 "${PENDING_OUTPUT}" \
  --producer-git-sha "${EXPECTED_GIT_SHA}"
if [[ ! -s "${PENDING_OUTPUT}" || ! -s "${PENDING_SUMMARY}" ]]; then
  echo "S63 native merge did not create complete staged outputs." >&2
  exit 9
fi
sha256sum "${PENDING_OUTPUT}" > "${PENDING_SHA}"
mv -- "${PENDING_OUTPUT}" "${OUTPUT}"
mv -- "${PENDING_SHA}" "${OUTPUT}.sha256"
mv -- "${PENDING_SUMMARY}" "${OUTPUT%.h5}.summary.json"
rmdir "${PENDING_DIR}" "${CLAIM}"
echo "[teacher-v3-s63-native-merge] complete=$(date -Is) output=${OUTPUT}"
