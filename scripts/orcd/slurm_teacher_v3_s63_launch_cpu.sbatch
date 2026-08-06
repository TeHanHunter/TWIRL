#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-tv3-s63-launch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 077

AUTHORIZATION_CONTRACT="s63_teacher_v3_prospective_v1"
if [[ "${TWIRL_TEACHER_V3_PROSPECTIVE_CONTRACT:-}" != "${AUTHORIZATION_CONTRACT}" ]]; then
  echo "S63 launch freeze requires the prospective authorization contract." >&2
  exit 64
fi

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:?set the clean pinned TWIRL checkout}"
EXPECTED_GIT_SHA="${TWIRL_EXPECTED_GIT_SHA:?set the full deployed Git SHA}"
PYTHON="${TWIRL_ORCD_STAGE2_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56/bin/python}"

PLAN="${TWIRL_S63_PROSPECTIVE_PLAN:?set prospective plan}"
PLAN_SHA="${TWIRL_S63_PROSPECTIVE_PLAN_SHA256:?set prospective-plan SHA-256}"
POLICY="${TWIRL_S63_SELECTION_POLICY:?set selection policy}"
POLICY_SHA="${TWIRL_S63_SELECTION_POLICY_SHA256:?set selection-policy SHA-256}"
VALIDATION="${TWIRL_S63_ACCEPTED_VALIDATION:?set attested Stage-1 receipt copy}"
VALIDATION_SHA="${TWIRL_S63_ACCEPTED_VALIDATION_SHA256:?set receipt SHA-256}"
CADENCE="${TWIRL_S63_CADENCE_TABLE:?set cadence table}"
CADENCE_SHA="${TWIRL_S63_CADENCE_TABLE_SHA256:?set cadence-table SHA-256}"
CADENCE_MANIFEST="${TWIRL_S63_CADENCE_MANIFEST:?set cadence manifest}"
CADENCE_MANIFEST_SHA="${TWIRL_S63_CADENCE_MANIFEST_SHA256:?set cadence-manifest SHA-256}"
COMPACT="${TWIRL_S63_COMPACT_H5:?set compact ADP HDF5}"
COMPACT_SHA="${TWIRL_S63_COMPACT_H5_SHA256:?set compact-HDF5 SHA-256}"
COMPACT_MANIFEST="${TWIRL_S63_COMPACT_MANIFEST:?set compact manifest}"
COMPACT_MANIFEST_SHA="${TWIRL_S63_COMPACT_MANIFEST_SHA256:?set compact-manifest SHA-256}"
RESERVED="${TWIRL_S63_RESERVED_TICS:?set reserved TIC inventory}"
RESERVED_SHA="${TWIRL_S63_RESERVED_TICS_SHA256:?set reserved-TIC SHA-256}"
CORPUS="${TWIRL_TEACHER_V3_CORPUS:?set frozen Teacher-v3 corpus}"
CORPUS_SHA="${TWIRL_TEACHER_V3_CORPUS_SHA256:?set corpus SHA-256}"
ALLOWLIST="${TWIRL_S63_MODEL_READY_ALLOWLIST:?set model-ready allowlist}"
ALLOWLIST_SHA="${TWIRL_S63_MODEL_READY_ALLOWLIST_SHA256:?set allowlist SHA-256}"
ALLOWLIST_SUMMARY="${TWIRL_S63_MODEL_READY_SUMMARY:?set model-ready summary}"
ALLOWLIST_SUMMARY_SHA="${TWIRL_S63_MODEL_READY_SUMMARY_SHA256:?set model-ready-summary SHA-256}"
PRIMARY="${TWIRL_S63_PRIMARY_COHORT:?set primary cohort}"
PRIMARY_SHA="${TWIRL_S63_PRIMARY_COHORT_SHA256:?set primary-cohort SHA-256}"
REPEATED="${TWIRL_S63_REPEATED_HOST_COHORT:?set repeated-host cohort}"
REPEATED_SHA="${TWIRL_S63_REPEATED_HOST_COHORT_SHA256:?set repeated-host SHA-256}"
NOT_READY="${TWIRL_S63_RESERVED_NOT_MODEL_READY:?set excluded inventory}"
NOT_READY_SHA="${TWIRL_S63_RESERVED_NOT_MODEL_READY_SHA256:?set excluded-inventory SHA-256}"
COHORT_SUMMARY="${TWIRL_S63_COHORT_SUMMARY:?set cohort summary}"
COHORT_SUMMARY_SHA="${TWIRL_S63_COHORT_SUMMARY_SHA256:?set cohort-summary SHA-256}"
BLS="${TWIRL_S63_BLS_MERGED:?set merged BLS table}"
BLS_SHA="${TWIRL_S63_BLS_MERGED_SHA256:?set merged-BLS SHA-256}"
BLS_SUMMARY="${TWIRL_S63_BLS_SUMMARY:?set merged BLS summary}"
BLS_SUMMARY_SHA="${TWIRL_S63_BLS_SUMMARY_SHA256:?set BLS-summary SHA-256}"
CANDIDATES="${TWIRL_S63_CANDIDATES:?set candidate table}"
CANDIDATES_SHA="${TWIRL_S63_CANDIDATES_SHA256:?set candidate SHA-256}"
CANDIDATE_SUMMARY="${TWIRL_S63_CANDIDATE_SUMMARY:?set candidate summary}"
CANDIDATE_SUMMARY_SHA="${TWIRL_S63_CANDIDATE_SUMMARY_SHA256:?set candidate-summary SHA-256}"
RAW="${TWIRL_S63_RAW_SOURCE_H5:?set raw-source HDF5}"
RAW_SHA="${TWIRL_S63_RAW_SOURCE_H5_SHA256:?set raw-source SHA-256}"
RAW_SUMMARY="${TWIRL_S63_RAW_SOURCE_SUMMARY:?set raw-source summary}"
RAW_SUMMARY_SHA="${TWIRL_S63_RAW_SOURCE_SUMMARY_SHA256:?set raw-summary SHA-256}"
NATIVE="${TWIRL_S63_NATIVE_H5:?set native HDF5}"
NATIVE_SHA="${TWIRL_S63_NATIVE_H5_SHA256:?set native SHA-256}"
NATIVE_SUMMARY="${TWIRL_S63_NATIVE_SUMMARY:?set native summary}"
NATIVE_SUMMARY_SHA="${TWIRL_S63_NATIVE_SUMMARY_SHA256:?set native-summary SHA-256}"
OUTPUT="${TWIRL_S63_LAUNCH_MANIFEST:?set fresh launch-manifest path}"

TWIRL_ROOT="$(realpath -e -- "${TWIRL_ROOT}")"
REPO="$(realpath -e -- "${REPO}")"
OUTPUT="$(realpath -m -- "${OUTPUT}")"
TWIRL_EXPECTED_GIT_SHA="${EXPECTED_GIT_SHA}" \
  "${REPO}/scripts/assert_clean_checkout.sh" "${REPO}" >/dev/null
if [[ "${OUTPUT}" != "${TWIRL_ROOT}"/* \
   || "${OUTPUT}" == "${REPO}" \
   || "${OUTPUT}" == "${REPO}"/* ]]; then
  echo "S63 launch output must be canonical, under TWIRL_ROOT, and outside the checkout." >&2
  exit 3
fi
if [[ -e "${OUTPUT}" ]]; then
  echo "Refusing to overwrite S63 launch manifest: ${OUTPUT}" >&2
  exit 4
fi

sha256_check() {
  local path=$1
  local expected=$2
  local observed remainder
  if [[ ! -s "${path}" ]]; then
    echo "Missing S63 launch input: ${path}" >&2
    exit 5
  fi
  read -r observed remainder < <(sha256sum "${path}")
  if [[ "${observed}" != "${expected}" ]]; then
    echo "S63 launch input SHA-256 mismatch: ${path}" >&2
    exit 6
  fi
}
for binding in \
  "${PLAN}|${PLAN_SHA}" "${POLICY}|${POLICY_SHA}" \
  "${VALIDATION}|${VALIDATION_SHA}" "${CADENCE}|${CADENCE_SHA}" \
  "${CADENCE_MANIFEST}|${CADENCE_MANIFEST_SHA}" "${COMPACT}|${COMPACT_SHA}" \
  "${COMPACT_MANIFEST}|${COMPACT_MANIFEST_SHA}" "${RESERVED}|${RESERVED_SHA}" \
  "${CORPUS}|${CORPUS_SHA}" "${ALLOWLIST}|${ALLOWLIST_SHA}" \
  "${ALLOWLIST_SUMMARY}|${ALLOWLIST_SUMMARY_SHA}" "${PRIMARY}|${PRIMARY_SHA}" \
  "${REPEATED}|${REPEATED_SHA}" "${NOT_READY}|${NOT_READY_SHA}" \
  "${COHORT_SUMMARY}|${COHORT_SUMMARY_SHA}" "${BLS}|${BLS_SHA}" \
  "${BLS_SUMMARY}|${BLS_SUMMARY_SHA}" "${CANDIDATES}|${CANDIDATES_SHA}" \
  "${CANDIDATE_SUMMARY}|${CANDIDATE_SUMMARY_SHA}" "${RAW}|${RAW_SHA}" \
  "${RAW_SUMMARY}|${RAW_SUMMARY_SHA}" "${NATIVE}|${NATIVE_SHA}" \
  "${NATIVE_SUMMARY}|${NATIVE_SUMMARY_SHA}"; do
  sha256_check "${binding%%|*}" "${binding#*|}"
done

OUTPUT_PARENT="$(dirname "${OUTPUT}")"
CLAIM_ROOT="${OUTPUT_PARENT}/.claims"
PENDING_ROOT="${OUTPUT_PARENT}/.pending"
CLAIM="${CLAIM_ROOT}/$(basename "${OUTPUT}").claim"
PENDING_DIR="${PENDING_ROOT}/$(basename "${OUTPUT}" .json)"
PENDING_OUTPUT="${PENDING_DIR}/$(basename "${OUTPUT}")"
mkdir -p "${TWIRL_ROOT}/logs" "${OUTPUT_PARENT}" "${CLAIM_ROOT}" "${PENDING_ROOT}"
if ! mkdir "${CLAIM}" 2>/dev/null; then
  echo "S63 launch freeze is already claimed/interrupted; audit: ${CLAIM}" >&2
  exit 7
fi
if [[ -e "${PENDING_DIR}" ]]; then
  echo "S63 launch freeze has staged interrupted state; audit: ${PENDING_DIR}" >&2
  exit 8
fi
mkdir "${PENDING_DIR}"

cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export PYTHONNOUSERSITE=1
unset LD_LIBRARY_PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
"${PYTHON}" -c 'import h5py, pandas, pyarrow; print("S63 launch environment OK")'
"${PYTHON}" scripts/stage5_validation/build_s63_teacher_v3_launch_manifest.py \
  --preregistered-contract "${PLAN}" --selection-policy "${POLICY}" \
  --accepted-validation "${VALIDATION}" \
  --cadence-table "${CADENCE}" --cadence-manifest "${CADENCE_MANIFEST}" \
  --compact-h5 "${COMPACT}" --compact-manifest "${COMPACT_MANIFEST}" \
  --reserved-tics "${RESERVED}" --teacher-v3-corpus "${CORPUS}" \
  --model-ready-allowlist "${ALLOWLIST}" --model-ready-summary "${ALLOWLIST_SUMMARY}" \
  --primary-cohort "${PRIMARY}" --repeated-host-cohort "${REPEATED}" \
  --reserved-not-model-ready "${NOT_READY}" --cohort-summary "${COHORT_SUMMARY}" \
  --bls-peaks "${BLS}" --bls-summary "${BLS_SUMMARY}" \
  --candidates "${CANDIDATES}" --candidate-summary "${CANDIDATE_SUMMARY}" \
  --raw-source-h5 "${RAW}" --raw-source-summary "${RAW_SUMMARY}" \
  --native-h5 "${NATIVE}" --native-summary "${NATIVE_SUMMARY}" \
  --repo "${REPO}" --expected-git-sha "${EXPECTED_GIT_SHA}" \
  --out "${PENDING_OUTPUT}" >/dev/null
if [[ ! -s "${PENDING_OUTPUT}" ]]; then
  echo "S63 launch builder did not create its staged manifest." >&2
  exit 9
fi
mv -- "${PENDING_OUTPUT}" "${OUTPUT}"
rmdir "${PENDING_DIR}" "${CLAIM}"
echo "[teacher-v3-s63-launch] complete=$(date -Is) output=${OUTPUT}"
