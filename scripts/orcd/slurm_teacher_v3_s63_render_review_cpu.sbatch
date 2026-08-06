#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-tv3-s63-render
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 077

AUTHORIZATION_CONTRACT="s63_teacher_v3_prospective_v1"
if [[ "${TWIRL_TEACHER_V3_PROSPECTIVE_CONTRACT:-}" != "${AUTHORIZATION_CONTRACT}" ]]; then
  echo "S63 review rendering requires the prospective authorization contract." >&2
  exit 64
fi
TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:?set the clean pinned TWIRL checkout}"
EXPECTED_GIT_SHA="${TWIRL_EXPECTED_GIT_SHA:?set the full deployed Git SHA}"
PYTHON="${TWIRL_ORCD_STAGE2_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56/bin/python}"
QUEUE="${TWIRL_S63_REVIEW_QUEUE:?set public review_queue_1100.csv}"
QUEUE_SHA="${TWIRL_S63_REVIEW_QUEUE_SHA256:?set queue SHA-256}"
PUBLIC_SUMMARY="${TWIRL_S63_REVIEW_PUBLIC_SUMMARY:?set public_summary.json}"
PUBLIC_SUMMARY_SHA="${TWIRL_S63_REVIEW_PUBLIC_SUMMARY_SHA256:?set public-summary SHA-256}"
HIDDEN_QUEUE_DIR="${TWIRL_S63_QUEUE_HIDDEN_DIR:?set the paired hidden queue directory}"
BUNDLE_MARKER="${TWIRL_S63_QUEUE_BUNDLE_MARKER:-$(dirname "${QUEUE}")/bundle_complete.json}"
LAUNCH="${TWIRL_S63_LAUNCH_MANIFEST:?set launch manifest}"
LAUNCH_SHA="${TWIRL_S63_LAUNCH_MANIFEST_SHA256:?set launch-manifest SHA-256}"
COMPACT="${TWIRL_S63_COMPACT_H5:?set compact ADP HDF5}"
COMPACT_SHA="${TWIRL_S63_COMPACT_H5_SHA256:?set compact-HDF5 SHA-256}"
CADENCE="${TWIRL_S63_CADENCE_TABLE:?set cadence-reference table}"
CADENCE_SHA="${TWIRL_S63_CADENCE_TABLE_SHA256:?set cadence-table SHA-256}"
CADENCE_MANIFEST="${TWIRL_S63_CADENCE_MANIFEST:?set cadence manifest}"
CADENCE_MANIFEST_SHA="${TWIRL_S63_CADENCE_MANIFEST_SHA256:?set cadence-manifest SHA-256}"
OUTPUT="${TWIRL_S63_REVIEW_RENDER_DIR:?set fresh review-render directory}"

TWIRL_ROOT="$(realpath -e -- "${TWIRL_ROOT}")"
REPO="$(realpath -e -- "${REPO}")"
OUTPUT="$(realpath -m -- "${OUTPUT}")"
TWIRL_EXPECTED_GIT_SHA="${EXPECTED_GIT_SHA}" \
  "${REPO}/scripts/assert_clean_checkout.sh" "${REPO}" >/dev/null
if [[ "${OUTPUT}" != "${TWIRL_ROOT}"/* \
   || "${OUTPUT}" == "${REPO}" \
   || "${OUTPUT}" == "${REPO}"/* ]]; then
  echo "S63 render output must be canonical, under TWIRL_ROOT, and outside checkout." >&2
  exit 3
fi
if [[ -e "${OUTPUT}" ]]; then
  echo "Refusing to overwrite S63 review-render output: ${OUTPUT}" >&2
  exit 4
fi
sha256_check() {
  local path=$1 expected=$2 observed remainder
  if [[ ! -s "${path}" ]]; then
    echo "Missing S63 render input: ${path}" >&2
    exit 5
  fi
  read -r observed remainder < <(sha256sum "${path}")
  if [[ "${observed}" != "${expected}" ]]; then
    echo "S63 render input SHA-256 mismatch: ${path}" >&2
    exit 6
  fi
}
sha256_check "${QUEUE}" "${QUEUE_SHA}"
sha256_check "${PUBLIC_SUMMARY}" "${PUBLIC_SUMMARY_SHA}"
sha256_check "${LAUNCH}" "${LAUNCH_SHA}"
sha256_check "${COMPACT}" "${COMPACT_SHA}"
sha256_check "${CADENCE}" "${CADENCE_SHA}"
sha256_check "${CADENCE_MANIFEST}" "${CADENCE_MANIFEST_SHA}"

OUTPUT_PARENT="$(dirname "${OUTPUT}")"
CLAIM_ROOT="${OUTPUT_PARENT}/.claims"
PENDING_ROOT="${OUTPUT_PARENT}/.pending"
CLAIM="${CLAIM_ROOT}/$(basename "${OUTPUT}").claim"
PENDING="${PENDING_ROOT}/$(basename "${OUTPUT}")"
SHEETS="${PENDING}/sheets"
METRICS="${PENDING}/twirl_vet_metrics.csv"
RENDER_SUMMARY="${PENDING}/twirl_two_aperture_vet_summary.json"
VERIFICATION="${PENDING}/review_sheet_coverage.json"
mkdir -p "${TWIRL_ROOT}/logs" "${OUTPUT_PARENT}" "${CLAIM_ROOT}" "${PENDING_ROOT}"
if ! mkdir "${CLAIM}" 2>/dev/null; then
  echo "S63 render is already claimed/interrupted; audit: ${CLAIM}" >&2
  exit 7
fi
if [[ -e "${PENDING}" ]]; then
  echo "S63 render has staged interrupted state; audit: ${PENDING}" >&2
  exit 8
fi
mkdir "${PENDING}"

cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export PYTHONNOUSERSITE=1
unset LD_LIBRARY_PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg
"${PYTHON}" -c 'import astropy, h5py, matplotlib, numpy, pandas; print("S63 render environment OK")'
"${PYTHON}" -c \
  'from pathlib import Path; import sys; from twirl.vetting.teacher_v3_prospective import verify_s63_prospective_review_bundle; verify_s63_prospective_review_bundle(public_dir=Path(sys.argv[1]), hidden_dir=Path(sys.argv[2]))' \
  "$(dirname "${QUEUE}")" "${HIDDEN_QUEUE_DIR}"
"${PYTHON}" scripts/stage5_validation/render_two_aperture_vet_sheets.py \
  --queue-csv "${QUEUE}" \
  --lc-export-h5 "${COMPACT}" \
  --out-dir "${SHEETS}" \
  --metrics-csv "${METRICS}" \
  --summary-json "${RENDER_SUMMARY}" \
  --branch-name current_adp \
  --apertures DET_FLUX_ADP_SML DET_FLUX_ADP \
  --use-row-ephemeris \
  --no-pdf \
  --cadence-reference "${CADENCE}" \
  --cadence-reference-manifest "${CADENCE_MANIFEST}" \
  --orbitid-policy strict \
  --producer-git-sha "${EXPECTED_GIT_SHA}" \
  --launch-manifest-sha256 "${LAUNCH_SHA}" \
  --workers "${SLURM_CPUS_PER_TASK:-8}" >/dev/null
"${PYTHON}" scripts/stage5_validation/verify_teacher_v3_s63_review_sheets.py \
  --queue-csv "${QUEUE}" \
  --public-summary "${PUBLIC_SUMMARY}" \
  --bundle-completion-marker "${BUNDLE_MARKER}" \
  --metrics-csv "${METRICS}" \
  --render-summary "${RENDER_SUMMARY}" \
  --sheet-dir "${SHEETS}" \
  --compact-h5 "${COMPACT}" \
  --cadence-reference "${CADENCE}" \
  --cadence-reference-manifest "${CADENCE_MANIFEST}" \
  --expected-rows 1100 \
  --producer-git-sha "${EXPECTED_GIT_SHA}" \
  --launch-manifest-sha256 "${LAUNCH_SHA}" \
  --out "${VERIFICATION}" >/dev/null
mv -- "${PENDING}" "${OUTPUT}"
rmdir "${CLAIM}"
echo "[teacher-v3-s63-render] complete=$(date -Is) output=${OUTPUT}"
