#!/bin/bash
# Development-only representation export and matched five-fold transfer test.
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-ssl-v4-eval
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=192G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 12:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

readonly TWIRL_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v4r2_s56_s62_a2v1_raw_v1_detector_consistent_adp_bls_eligible"
readonly TEACHER_V3_REPORT_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_v3_s56_s62_a2v1_current_adp"
readonly TEACHER_V3_MODEL_ROOT="${TWIRL_ROOT}/results/models/teacher_v3_s56_s62_a2v1_current_adp/full"
readonly PYTHON="${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python"
readonly TRAINING_TABLE="${TEACHER_V3_REPORT_ROOT}/native_preparation/training_table_with_splits.csv"
readonly REGISTRY="${RUN_ROOT}/frozen/ssl_registry/teacher_ssl_fullpool_registry.parquet"
readonly REGISTRY_SUMMARY="${RUN_ROOT}/frozen/ssl_registry/teacher_ssl_fullpool_registry.summary.json"
readonly COMPLETION_RELEASE="${RUN_ROOT}/frozen/model/teacher_ssl_fullpool_v4_five_fold.complete.json"
readonly EVALUATION_ROOT="${RUN_ROOT}/evaluation/teacher_ssl_fullpool_v4_development_transfer_v1"
readonly FINE_TUNE_EPOCHS=100
readonly BATCH_SIZE=32
readonly WORKERS=8
readonly SEED=560064
readonly LINEAR_PROBE_STEPS=500
readonly BOOTSTRAP_DRAWS=2000

: "${TWIRL_ORCD_REPO:?controller must export TWIRL_ORCD_REPO}"
: "${TWIRL_EXPECTED_GIT_SHA:?controller must export TWIRL_EXPECTED_GIT_SHA}"
readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}"

if [[ ! "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]] || \
   [[ "$(git -C "${REPO}" rev-parse HEAD)" != "${EXPECTED_SHA}" ]] || \
   [[ -n "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]]; then
  echo "Evaluation repository is not the exact clean deployed revision." >&2
  exit 2
fi
for input in \
  "${TRAINING_TABLE}" \
  "${REGISTRY}" \
  "${REGISTRY_SUMMARY}" \
  "${COMPLETION_RELEASE}" \
  "${COMPLETION_RELEASE}.sha256" \
  "${TEACHER_V3_MODEL_ROOT}/summary.json"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing full-pool evaluation authority: ${input}" >&2
    exit 3
  fi
done
if [[ -e "${EVALUATION_ROOT}" ]]; then
  echo "Refusing to overwrite full-pool evaluation: ${EVALUATION_ROOT}" >&2
  exit 4
fi

mkdir -p "${TWIRL_ROOT}/logs" "$(dirname "${EVALUATION_ROOT}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
module load cuda/13.1.0 || module load cuda/12.9.1 || true
nvidia-smi

echo "[teacher-ssl-fullpool-v4-eval] start=$(date -Is) output=${EVALUATION_ROOT}"
"${PYTHON}" scripts/stage5_validation/evaluate_teacher_ssl_fullpool.py \
  --training-table "${TRAINING_TABLE}" \
  --registry "${REGISTRY}" \
  --registry-summary "${REGISTRY_SUMMARY}" \
  --completion-release "${COMPLETION_RELEASE}" \
  --baseline-teacher-v3-root "${TEACHER_V3_MODEL_ROOT}" \
  --out-dir "${EVALUATION_ROOT}" \
  --fine-tune-epochs "${FINE_TUNE_EPOCHS}" \
  --batch-size "${BATCH_SIZE}" \
  --workers "${WORKERS}" \
  --seed "${SEED}" \
  --linear-probe-steps "${LINEAR_PROBE_STEPS}" \
  --bootstrap-draws "${BOOTSTRAP_DRAWS}"
echo "[teacher-ssl-fullpool-v4-eval] complete=$(date -Is)"
