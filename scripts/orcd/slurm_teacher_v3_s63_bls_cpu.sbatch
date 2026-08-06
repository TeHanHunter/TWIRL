#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-tv3-s63-bls
#SBATCH --array=0-15%8
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.err

set -euo pipefail
umask 077

AUTHORIZATION_CONTRACT="s63_teacher_v3_prospective_v1"
if [[ "${TWIRL_TEACHER_V3_PROSPECTIVE_CONTRACT:-}" != "${AUTHORIZATION_CONTRACT}" ]]; then
  echo "S63 BLS requires TWIRL_TEACHER_V3_PROSPECTIVE_CONTRACT=${AUTHORIZATION_CONTRACT}." >&2
  exit 64
fi

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:?set the clean pinned TWIRL checkout}"
EXPECTED_GIT_SHA="${TWIRL_EXPECTED_GIT_SHA:?set the full deployed Git SHA}"
PYTHON="${TWIRL_ORCD_STAGE2_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56/bin/python}"
VALIDATION="${TWIRL_S63_ACCEPTED_VALIDATION:?set staged accepted S63 validation JSON}"
EXPECTED_VALIDATION_SHA="${TWIRL_S63_ACCEPTED_VALIDATION_SHA256:?set accepted validation SHA-256}"
COMPACT="${TWIRL_S63_COMPACT_H5:?set staged S63 compact ADP HDF5}"
EXPECTED_COMPACT_SHA="${TWIRL_S63_COMPACT_H5_SHA256:?set compact HDF5 SHA-256}"
CADENCE="${TWIRL_S63_CADENCE_TABLE:?set staged S63 cadence table}"
EXPECTED_CADENCE_SHA="${TWIRL_S63_CADENCE_TABLE_SHA256:?set cadence-table SHA-256}"
CADENCE_MANIFEST="${TWIRL_S63_CADENCE_MANIFEST:?set staged S63 cadence manifest}"
EXPECTED_CADENCE_MANIFEST_SHA="${TWIRL_S63_CADENCE_MANIFEST_SHA256:?set cadence-manifest SHA-256}"
ALLOWLIST="${TWIRL_S63_MODEL_READY_ALLOWLIST:?set exact model-ready TIC allowlist}"
EXPECTED_ALLOWLIST_SHA="${TWIRL_S63_MODEL_READY_ALLOWLIST_SHA256:?set allowlist SHA-256}"
OUT_DIR="${TWIRL_S63_BLS_SHARD_DIR:?set fresh S63 BLS shard directory}"
SHARD_INDEX="${SLURM_ARRAY_TASK_ID:?S63 BLS must run as the locked 0-15 array}"
N_SHARDS=16

if (( SHARD_INDEX < 0 || SHARD_INDEX >= N_SHARDS )); then
  echo "S63 BLS array index is outside 0-15." >&2
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
  echo "S63 BLS output must be canonical, under TWIRL_ROOT, and outside the checkout." >&2
  exit 3
fi
for input in "${VALIDATION}" "${COMPACT}" "${CADENCE}" "${CADENCE_MANIFEST}" "${ALLOWLIST}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing S63 BLS input: ${input}" >&2
    exit 4
  fi
done
sha256_check() {
  local path=$1
  local expected=$2
  local observed remainder
  read -r observed remainder < <(sha256sum "${path}")
  if [[ "${observed}" != "${expected}" ]]; then
    echo "S63 BLS input SHA-256 mismatch: ${path}" >&2
    exit 5
  fi
}
sha256_check "${VALIDATION}" "${EXPECTED_VALIDATION_SHA}"
sha256_check "${COMPACT}" "${EXPECTED_COMPACT_SHA}"
sha256_check "${CADENCE}" "${EXPECTED_CADENCE_SHA}"
sha256_check "${CADENCE_MANIFEST}" "${EXPECTED_CADENCE_MANIFEST_SHA}"
sha256_check "${ALLOWLIST}" "${EXPECTED_ALLOWLIST_SHA}"

printf -v suffix '_%03d' "${SHARD_INDEX}"
PEAKS="${OUT_DIR}/real_adp_bls_peaks${suffix}.parquet"
SUMMARY="${OUT_DIR}/summary${suffix}.json"
if [[ -e "${PEAKS}" || -e "${SUMMARY}" ]]; then
  echo "Refusing to overwrite S63 BLS shard ${SHARD_INDEX}." >&2
  exit 6
fi

CLAIM_ROOT="${OUT_DIR}/.claims"
PENDING_ROOT="${OUT_DIR}/.pending"
CLAIM="${CLAIM_ROOT}/bls${suffix}.claim"
PENDING="${PENDING_ROOT}/bls${suffix}"
mkdir -p "${TWIRL_ROOT}/logs" "${OUT_DIR}" "${CLAIM_ROOT}" "${PENDING_ROOT}"
if ! mkdir "${CLAIM}" 2>/dev/null; then
  echo "S63 BLS shard is already claimed/interrupted; audit before retry: ${CLAIM}" >&2
  exit 7
fi
if [[ -e "${PENDING}" ]]; then
  echo "S63 BLS shard has interrupted staged output; audit before retry: ${PENDING}" >&2
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
"${PYTHON}" -c 'import astropy, h5py, numpy, pandas, pyarrow; print("S63 ORCD BLS environment OK")'
"${PYTHON}" scripts/stage5_validation/build_s63_teacher_v3_launch_manifest.py \
  --accepted-validation "${VALIDATION}" \
  --expected-git-sha "${EXPECTED_GIT_SHA}" \
  --preflight-only >/dev/null

echo "[teacher-v3-s63-bls] shard=${SHARD_INDEX}/${N_SHARDS} start=$(date -Is)"
"${PYTHON}" scripts/stage5_validation/build_s56_adp_real_bls_peaks.py \
  --sector 63 \
  --compact-lc "${COMPACT}" \
  --expected-compact-lc-sha256 "${EXPECTED_COMPACT_SHA}" \
  --cadence-reference "${CADENCE}" \
  --cadence-reference-manifest "${CADENCE_MANIFEST}" \
  --target-allowlist "${ALLOWLIST}" \
  --expected-target-allowlist-sha256 "${EXPECTED_ALLOWLIST_SHA}" \
  --orbitid-policy strict \
  --out-dir "${PENDING}" \
  --workers "${SLURM_CPUS_PER_TASK:-8}" \
  --n-periods 50000 \
  --n-peaks 10 \
  --shard-index "${SHARD_INDEX}" \
  --n-shards "${N_SHARDS}" \
  --progress-every 25 \
  --source-product-tag A2v1
PENDING_PEAKS="${PENDING}/real_adp_bls_peaks${suffix}.parquet"
PENDING_SUMMARY="${PENDING}/summary${suffix}.json"
if [[ ! -s "${PENDING_PEAKS}" || ! -s "${PENDING_SUMMARY}" ]]; then
  echo "S63 BLS producer did not create a complete staged shard." >&2
  exit 9
fi
read -r summary_pre_sha summary_hash_remainder < <(sha256sum "${PENDING_SUMMARY}")
"${PYTHON}" scripts/stage5_validation/attest_s63_teacher_v3_json.py \
  --json "${PENDING_SUMMARY}" \
  --expected-sha256 "${summary_pre_sha}" \
  --producer-git-sha "${EXPECTED_GIT_SHA}" \
  --bls-published-peak-table "${PEAKS}" \
  --bls-published-summary "${SUMMARY}" >/dev/null
mv -- "${PENDING_PEAKS}" "${PEAKS}"
mv -- "${PENDING_SUMMARY}" "${SUMMARY}"
rmdir "${PENDING}" "${CLAIM}"
echo "[teacher-v3-s63-bls] shard=${SHARD_INDEX} complete=$(date -Is)"
