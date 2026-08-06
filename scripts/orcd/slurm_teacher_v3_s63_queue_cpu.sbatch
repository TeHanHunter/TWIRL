#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-tv3-s63-queue
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH --exclude=node4900
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 077

AUTHORIZATION_CONTRACT="s63_teacher_v3_prospective_v1"
if [[ "${TWIRL_TEACHER_V3_PROSPECTIVE_CONTRACT:-}" != "${AUTHORIZATION_CONTRACT}" ]]; then
  echo "S63 queue construction requires TWIRL_TEACHER_V3_PROSPECTIVE_CONTRACT=${AUTHORIZATION_CONTRACT}." >&2
  exit 64
fi

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:?set the clean pinned TWIRL checkout}"
EXPECTED_GIT_SHA="${TWIRL_EXPECTED_GIT_SHA:?set the full deployed Git SHA}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python}"

PLAN="${TWIRL_S63_PROSPECTIVE_PLAN:?set the frozen S63 prospective plan}"
EXPECTED_PLAN_SHA="${TWIRL_S63_PROSPECTIVE_PLAN_SHA256:?set prospective-plan SHA-256}"
SELECTION_POLICY="${TWIRL_S63_SELECTION_POLICY:?set the pre-score selection policy}"
EXPECTED_SELECTION_SHA="${TWIRL_S63_SELECTION_POLICY_SHA256:?set selection-policy SHA-256}"
LAUNCH_MANIFEST="${TWIRL_S63_LAUNCH_MANIFEST:?set the frozen S63 launch manifest}"
EXPECTED_LAUNCH_SHA="${TWIRL_S63_LAUNCH_MANIFEST_SHA256:?set launch-manifest SHA-256}"

SCORING_SUMMARY="${TWIRL_S63_SCORING_SUMMARY:?set the exact S63 Teacher-v3 scoring summary}"
EXPECTED_SCORING_SUMMARY_SHA="${TWIRL_S63_SCORING_SUMMARY_SHA256:?set scoring-summary SHA-256}"
SCORES="${TWIRL_S63_SCORE_TABLE:?set the exact S63 Teacher-v3 score table}"
EXPECTED_SCORES_SHA="${TWIRL_S63_SCORE_TABLE_SHA256:?set score-table SHA-256}"
CANDIDATES="${TWIRL_S63_CANDIDATES:?set the exact S63 rank-one candidate table}"
EXPECTED_CANDIDATES_SHA="${TWIRL_S63_CANDIDATES_SHA256:?set candidate-table SHA-256}"
PRIMARY_COHORT="${TWIRL_S63_PRIMARY_COHORT:?set the frozen primary-cohort TIC inventory}"
EXPECTED_PRIMARY_COHORT_SHA="${TWIRL_S63_PRIMARY_COHORT_SHA256:?set primary-cohort SHA-256}"
REPEATED_HOST_COHORT="${TWIRL_S63_REPEATED_HOST_COHORT:?set the frozen repeated-host TIC inventory}"
EXPECTED_REPEATED_HOST_COHORT_SHA="${TWIRL_S63_REPEATED_HOST_COHORT_SHA256:?set repeated-host cohort SHA-256}"
COHORT_SUMMARY="${TWIRL_S63_COHORT_SUMMARY:?set the frozen S63 cohort summary}"
EXPECTED_COHORT_SUMMARY_SHA="${TWIRL_S63_COHORT_SUMMARY_SHA256:?set cohort-summary SHA-256}"

PUBLIC_OUT_DIR="${TWIRL_S63_QUEUE_PUBLIC_OUT_DIR:?set a fresh public queue output directory}"
HIDDEN_OUT_DIR="${TWIRL_S63_QUEUE_HIDDEN_OUT_DIR:?set a fresh hidden provenance output directory}"

TWIRL_ROOT="$(realpath -e -- "${TWIRL_ROOT}")"
REPO="$(realpath -e -- "${REPO}")"
PUBLIC_OUT_DIR="$(realpath -m -- "${PUBLIC_OUT_DIR}")"
HIDDEN_OUT_DIR="$(realpath -m -- "${HIDDEN_OUT_DIR}")"
TWIRL_EXPECTED_GIT_SHA="${EXPECTED_GIT_SHA}" \
  "${REPO}/scripts/assert_clean_checkout.sh" "${REPO}" >/dev/null

for output in "${PUBLIC_OUT_DIR}" "${HIDDEN_OUT_DIR}"; do
  if [[ "${output}" != /* || "${output}" == "${REPO}" || "${output}" == "${REPO}"/* ]]; then
    echo "S63 queue outputs must be absolute and outside the checkout: ${output}" >&2
    exit 3
  fi
  if [[ -e "${output}" ]]; then
    echo "Refusing to overwrite S63 queue output: ${output}" >&2
    exit 4
  fi
  if [[ "${output}" != "${TWIRL_ROOT}"/* ]]; then
    echo "S63 queue outputs must remain under the TWIRL data root." >&2
    exit 5
  fi
done
if [[ "${PUBLIC_OUT_DIR}" == "${HIDDEN_OUT_DIR}" \
   || "${PUBLIC_OUT_DIR}" == "${HIDDEN_OUT_DIR}"/* \
   || "${HIDDEN_OUT_DIR}" == "${PUBLIC_OUT_DIR}"/* ]]; then
  echo "Public and hidden S63 queue outputs must be distinct and non-nested." >&2
  exit 6
fi

sha256_check() {
  local path=$1
  local expected=$2
  local observed remainder
  if [[ ! -s "${path}" ]]; then
    echo "Missing S63 queue input: ${path}" >&2
    exit 7
  fi
  read -r observed remainder < <(sha256sum "${path}")
  if [[ "${observed}" != "${expected}" ]]; then
    echo "S63 queue input SHA-256 mismatch: ${path}" >&2
    exit 8
  fi
}
sha256_check "${PLAN}" "${EXPECTED_PLAN_SHA}"
sha256_check "${SELECTION_POLICY}" "${EXPECTED_SELECTION_SHA}"
sha256_check "${LAUNCH_MANIFEST}" "${EXPECTED_LAUNCH_SHA}"
sha256_check "${SCORING_SUMMARY}" "${EXPECTED_SCORING_SUMMARY_SHA}"
sha256_check "${SCORES}" "${EXPECTED_SCORES_SHA}"
sha256_check "${CANDIDATES}" "${EXPECTED_CANDIDATES_SHA}"
sha256_check "${PRIMARY_COHORT}" "${EXPECTED_PRIMARY_COHORT_SHA}"
sha256_check "${REPEATED_HOST_COHORT}" "${EXPECTED_REPEATED_HOST_COHORT_SHA}"
sha256_check "${COHORT_SUMMARY}" "${EXPECTED_COHORT_SUMMARY_SHA}"

mkdir -p "${TWIRL_ROOT}/logs" "$(dirname "${PUBLIC_OUT_DIR}")" "$(dirname "${HIDDEN_OUT_DIR}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export PYTHONNOUSERSITE=1
unset LD_LIBRARY_PATH
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

"${PYTHON}" -c 'import numpy, pandas, pyarrow, torch; print("S63 queue CPU environment ready")'

echo "[teacher-v3-s63-queue] start=$(date -Is)"
"${PYTHON}" scripts/stage5_validation/build_teacher_v3_s63_review_queue.py \
  --repo "${REPO}" \
  --expected-git-sha "${EXPECTED_GIT_SHA}" \
  --producer-git-sha "${EXPECTED_GIT_SHA}" \
  --plan "${PLAN}" \
  --expected-plan-sha256 "${EXPECTED_PLAN_SHA}" \
  --selection-policy "${SELECTION_POLICY}" \
  --expected-selection-policy-sha256 "${EXPECTED_SELECTION_SHA}" \
  --launch-manifest "${LAUNCH_MANIFEST}" \
  --expected-launch-manifest-sha256 "${EXPECTED_LAUNCH_SHA}" \
  --scoring-summary "${SCORING_SUMMARY}" \
  --expected-scoring-summary-sha256 "${EXPECTED_SCORING_SUMMARY_SHA}" \
  --scores "${SCORES}" \
  --expected-scores-sha256 "${EXPECTED_SCORES_SHA}" \
  --candidates "${CANDIDATES}" \
  --expected-candidates-sha256 "${EXPECTED_CANDIDATES_SHA}" \
  --primary-cohort "${PRIMARY_COHORT}" \
  --expected-primary-cohort-sha256 "${EXPECTED_PRIMARY_COHORT_SHA}" \
  --repeated-host-cohort "${REPEATED_HOST_COHORT}" \
  --expected-repeated-host-cohort-sha256 "${EXPECTED_REPEATED_HOST_COHORT_SHA}" \
  --cohort-summary "${COHORT_SUMMARY}" \
  --expected-cohort-summary-sha256 "${EXPECTED_COHORT_SUMMARY_SHA}" \
  --public-out-dir "${PUBLIC_OUT_DIR}" \
  --hidden-out-dir "${HIDDEN_OUT_DIR}" >/dev/null
echo "[teacher-v3-s63-queue] complete=$(date -Is) public=${PUBLIC_OUT_DIR}"
