#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-tv3-s63-score
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=192G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 08:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 077

AUTHORIZATION_CONTRACT="s63_teacher_v3_prospective_v1"
if [[ "${TWIRL_TEACHER_V3_PROSPECTIVE_CONTRACT:-}" != "${AUTHORIZATION_CONTRACT}" ]]; then
  echo "S63 scoring requires TWIRL_TEACHER_V3_PROSPECTIVE_CONTRACT=${AUTHORIZATION_CONTRACT}." >&2
  exit 64
fi

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:?set the clean pinned TWIRL checkout}"
EXPECTED_GIT_SHA="${TWIRL_EXPECTED_GIT_SHA:?set the full deployed Git SHA}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python}"

LAUNCH_MANIFEST="${TWIRL_S63_LAUNCH_MANIFEST:?set the frozen S63 launch manifest}"
EXPECTED_LAUNCH_SHA="${TWIRL_S63_LAUNCH_MANIFEST_SHA256:?set launch-manifest SHA-256}"
PLAN="${TWIRL_S63_PROSPECTIVE_PLAN:?set the frozen S63 prospective plan}"
EXPECTED_PLAN_SHA="${TWIRL_S63_PROSPECTIVE_PLAN_SHA256:?set prospective-plan SHA-256}"
SELECTION_POLICY="${TWIRL_S63_SELECTION_POLICY:?set the pre-score selection policy}"
EXPECTED_SELECTION_SHA="${TWIRL_S63_SELECTION_POLICY_SHA256:?set selection-policy SHA-256}"

RELEASE_SUMMARY="${TWIRL_TEACHER_V3_RELEASE_SUMMARY:?set frozen Teacher-v3 release summary}"
EXPECTED_RELEASE_SUMMARY_SHA="${TWIRL_TEACHER_V3_RELEASE_SUMMARY_SHA256:?set release-summary SHA-256}"
PRETEST_FREEZE="${TWIRL_TEACHER_V3_PRETEST_FREEZE:?set frozen Teacher-v3 pretest freeze}"
EXPECTED_PRETEST_FREEZE_SHA="${TWIRL_TEACHER_V3_PRETEST_FREEZE_SHA256:?set pretest-freeze SHA-256}"
CHECKPOINT_MANIFEST="${TWIRL_TEACHER_V3_CHECKPOINT_MANIFEST:?set selected Teacher-v3 checkpoint manifest}"
EXPECTED_CHECKPOINT_MANIFEST_SHA="${TWIRL_TEACHER_V3_CHECKPOINT_MANIFEST_SHA256:?set checkpoint-manifest SHA-256}"
CALIBRATION="${TWIRL_TEACHER_V3_CALIBRATION:?set pooled Teacher-v3 calibration}"
EXPECTED_CALIBRATION_SHA="${TWIRL_TEACHER_V3_CALIBRATION_SHA256:?set calibration SHA-256}"

CANDIDATES="${TWIRL_S63_CANDIDATES:?set exact S63 rank-one candidate table}"
EXPECTED_CANDIDATES_SHA="${TWIRL_S63_CANDIDATES_SHA256:?set candidate-table SHA-256}"
RESERVED_TICS="${TWIRL_S63_RESERVED_TICS:?set sealed S63 TIC inventory}"
EXPECTED_RESERVED_TICS_SHA="${TWIRL_S63_RESERVED_TICS_SHA256:?set reserved-TIC SHA-256}"
NATIVE_H5="${TWIRL_S63_NATIVE_H5:?set merged S63 native HDF5}"
EXPECTED_NATIVE_H5_SHA="${TWIRL_S63_NATIVE_H5_SHA256:?set native-HDF5 SHA-256}"
OUT_DIR="${TWIRL_S63_SCORE_OUT_DIR:?set a fresh S63 score output directory}"

TWIRL_ROOT="$(realpath -e -- "${TWIRL_ROOT}")"
REPO="$(realpath -e -- "${REPO}")"
OUT_DIR="$(realpath -m -- "${OUT_DIR}")"
TWIRL_EXPECTED_GIT_SHA="${EXPECTED_GIT_SHA}" \
  "${REPO}/scripts/assert_clean_checkout.sh" "${REPO}" >/dev/null
if [[ "${OUT_DIR}" != "${TWIRL_ROOT}"/* \
   || "${OUT_DIR}" == "${REPO}" \
   || "${OUT_DIR}" == "${REPO}"/* ]]; then
  echo "S63 score output must be canonical, under TWIRL_ROOT, and outside the checkout." >&2
  exit 3
fi
if [[ -e "${OUT_DIR}" ]]; then
  echo "Refusing to overwrite S63 score output: ${OUT_DIR}" >&2
  exit 4
fi

sha256_check() {
  local path=$1
  local expected=$2
  local observed remainder
  if [[ ! -s "${path}" ]]; then
    echo "Missing S63 scoring input: ${path}" >&2
    exit 5
  fi
  read -r observed remainder < <(sha256sum "${path}")
  if [[ "${observed}" != "${expected}" ]]; then
    echo "S63 scoring input SHA-256 mismatch: ${path}" >&2
    exit 6
  fi
}
sha256_check "${LAUNCH_MANIFEST}" "${EXPECTED_LAUNCH_SHA}"
sha256_check "${PLAN}" "${EXPECTED_PLAN_SHA}"
sha256_check "${SELECTION_POLICY}" "${EXPECTED_SELECTION_SHA}"
sha256_check "${RELEASE_SUMMARY}" "${EXPECTED_RELEASE_SUMMARY_SHA}"
sha256_check "${PRETEST_FREEZE}" "${EXPECTED_PRETEST_FREEZE_SHA}"
sha256_check "${CHECKPOINT_MANIFEST}" "${EXPECTED_CHECKPOINT_MANIFEST_SHA}"
sha256_check "${CALIBRATION}" "${EXPECTED_CALIBRATION_SHA}"
sha256_check "${CANDIDATES}" "${EXPECTED_CANDIDATES_SHA}"
sha256_check "${RESERVED_TICS}" "${EXPECTED_RESERVED_TICS_SHA}"
sha256_check "${NATIVE_H5}" "${EXPECTED_NATIVE_H5_SHA}"

mkdir -p "${TWIRL_ROOT}/logs" "$(dirname "${OUT_DIR}")"
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

module load cuda/13.1.0 || module load cuda/12.9.1 || true
nvidia-smi
"${PYTHON}" -c 'import h5py, pandas, pyarrow, torch; assert torch.cuda.is_available(), "CUDA is unavailable for S63 Teacher-v3 scoring"; print(torch.cuda.get_device_name(0))'

echo "[teacher-v3-s63-score] start=$(date -Is)"
"${PYTHON}" scripts/stage5_validation/score_teacher_v3_s63_candidates.py \
  --repo "${REPO}" \
  --expected-git-sha "${EXPECTED_GIT_SHA}" \
  --launch-manifest "${LAUNCH_MANIFEST}" \
  --launch-manifest-sha256 "${EXPECTED_LAUNCH_SHA}" \
  --preregistered-contract "${PLAN}" \
  --selection-policy "${SELECTION_POLICY}" \
  --release-summary "${RELEASE_SUMMARY}" \
  --release-summary-sha256 "${EXPECTED_RELEASE_SUMMARY_SHA}" \
  --pretest-freeze "${PRETEST_FREEZE}" \
  --pretest-freeze-sha256 "${EXPECTED_PRETEST_FREEZE_SHA}" \
  --checkpoint-manifest "${CHECKPOINT_MANIFEST}" \
  --checkpoint-manifest-sha256 "${EXPECTED_CHECKPOINT_MANIFEST_SHA}" \
  --calibration "${CALIBRATION}" \
  --calibration-sha256 "${EXPECTED_CALIBRATION_SHA}" \
  --candidates "${CANDIDATES}" \
  --candidates-sha256 "${EXPECTED_CANDIDATES_SHA}" \
  --s63-reserved-tics "${RESERVED_TICS}" \
  --s63-reserved-tics-sha256 "${EXPECTED_RESERVED_TICS_SHA}" \
  --native-h5 "${NATIVE_H5}" \
  --native-h5-sha256 "${EXPECTED_NATIVE_H5_SHA}" \
  --out-dir "${OUT_DIR}" \
  --batch-size "${TWIRL_S63_SCORE_BATCH_SIZE:-32}" \
  --workers "${TWIRL_S63_SCORE_WORKERS:-4}"
echo "[teacher-v3-s63-score] complete=$(date -Is) output=${OUT_DIR}"
