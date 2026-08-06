#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-tv3-s63-bls-merge
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 077

AUTHORIZATION_CONTRACT="s63_teacher_v3_prospective_v1"
if [[ "${TWIRL_TEACHER_V3_PROSPECTIVE_CONTRACT:-}" != "${AUTHORIZATION_CONTRACT}" ]]; then
  echo "S63 BLS merge requires the prospective authorization contract." >&2
  exit 64
fi
TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:?set the clean pinned TWIRL checkout}"
EXPECTED_GIT_SHA="${TWIRL_EXPECTED_GIT_SHA:?set the full deployed Git SHA}"
PYTHON="${TWIRL_ORCD_STAGE2_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56/bin/python}"
SHARD_DIR="${TWIRL_S63_BLS_SHARD_DIR:?set the completed S63 BLS shard directory}"
OUTPUT="${TWIRL_S63_BLS_MERGED:?set the fresh merged S63 BLS Parquet path}"
N_SHARDS=16

TWIRL_ROOT="$(realpath -e -- "${TWIRL_ROOT}")"
REPO="$(realpath -e -- "${REPO}")"
SHARD_DIR="$(realpath -e -- "${SHARD_DIR}")"
OUTPUT="$(realpath -m -- "${OUTPUT}")"
TWIRL_EXPECTED_GIT_SHA="${EXPECTED_GIT_SHA}" \
  "${REPO}/scripts/assert_clean_checkout.sh" "${REPO}" >/dev/null
if [[ "${OUTPUT}" != "${TWIRL_ROOT}"/* \
   || "${OUTPUT}" == "${REPO}" \
   || "${OUTPUT}" == "${REPO}"/* ]]; then
  echo "S63 BLS merged output must be canonical, under TWIRL_ROOT, and outside the checkout." >&2
  exit 3
fi
if [[ -e "${OUTPUT}" || -e "${OUTPUT%.parquet}.summary.json" ]]; then
  echo "Refusing to overwrite merged S63 BLS output." >&2
  exit 4
fi
for ((shard=0; shard<N_SHARDS; shard++)); do
  printf -v suffix '_%03d' "${shard}"
  for path in \
    "${SHARD_DIR}/real_adp_bls_peaks${suffix}.parquet" \
    "${SHARD_DIR}/summary${suffix}.json"; do
    if [[ ! -s "${path}" ]]; then
      echo "Missing S63 BLS shard input: ${path}" >&2
      exit 5
    fi
  done
done

OUTPUT_PARENT="$(dirname "${OUTPUT}")"
CLAIM_ROOT="${OUTPUT_PARENT}/.claims"
PENDING_ROOT="${OUTPUT_PARENT}/.pending"
CLAIM="${CLAIM_ROOT}/$(basename "${OUTPUT}").claim"
PENDING_DIR="${PENDING_ROOT}/$(basename "${OUTPUT}" .parquet)"
PENDING_OUTPUT="${PENDING_DIR}/$(basename "${OUTPUT}")"
PENDING_SUMMARY="${PENDING_OUTPUT%.parquet}.summary.json"
mkdir -p "${TWIRL_ROOT}/logs" "${OUTPUT_PARENT}" "${CLAIM_ROOT}" "${PENDING_ROOT}"
if ! mkdir "${CLAIM}" 2>/dev/null; then
  echo "S63 BLS merge is already claimed/interrupted; audit: ${CLAIM}" >&2
  exit 7
fi
if [[ -e "${PENDING_DIR}" ]]; then
  echo "S63 BLS merge has staged interrupted state; audit: ${PENDING_DIR}" >&2
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
"${PYTHON}" -c 'import numpy, pandas, pyarrow; print("S63 ORCD BLS merge environment OK")'
"${PYTHON}" scripts/stage5_validation/merge_adp_real_bls_shards.py \
  --shard-dir "${SHARD_DIR}" \
  --out-path "${PENDING_OUTPUT}" \
  --n-shards "${N_SHARDS}"
if [[ ! -s "${PENDING_OUTPUT}" || ! -s "${PENDING_SUMMARY}" ]]; then
  echo "S63 BLS merge did not create complete staged outputs." >&2
  exit 9
fi
read -r summary_pre_sha summary_hash_remainder < <(sha256sum "${PENDING_SUMMARY}")
"${PYTHON}" scripts/stage5_validation/attest_s63_teacher_v3_json.py \
  --json "${PENDING_SUMMARY}" \
  --expected-sha256 "${summary_pre_sha}" \
  --producer-git-sha "${EXPECTED_GIT_SHA}" \
  --bls-published-peak-table "${OUTPUT}" \
  --bls-published-summary "${OUTPUT%.parquet}.summary.json" >/dev/null
mv -- "${PENDING_OUTPUT}" "${OUTPUT}"
mv -- "${PENDING_SUMMARY}" "${OUTPUT%.parquet}.summary.json"
rmdir "${PENDING_DIR}" "${CLAIM}"
echo "[teacher-v3-s63-bls-merge] complete=$(date -Is) output=${OUTPUT}"
