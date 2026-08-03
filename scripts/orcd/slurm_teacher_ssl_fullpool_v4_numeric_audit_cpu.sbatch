#!/bin/bash
# Audit every model-facing sample in all 112 detector-consistent native shards.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-ssl-v4-numeric-audit
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH --array=0-111%4
#SBATCH -t 08:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.err

set -euo pipefail

readonly TWIRL_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly EXPECTED_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v4r2_s56_s62_a2v1_raw_v1_detector_consistent_adp_bls_eligible"
readonly N_SHARDS=16
readonly N_TASKS=112
: "${TWIRL_ORCD_REPO:?controller must export TWIRL_ORCD_REPO}"
: "${TWIRL_EXPECTED_GIT_SHA:?controller must export TWIRL_EXPECTED_GIT_SHA}"
: "${TWIRL_SSL_FULLPOOL_V4_RUN_ROOT:?controller must export the exact fresh v4 run root}"

REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}"
RUN_ROOT="$(realpath -m -- "${TWIRL_SSL_FULLPOOL_V4_RUN_ROOT}")"
PYTHON="${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python"

if [[ "${RUN_ROOT}" != "${EXPECTED_RUN_ROOT}" ]]; then
  echo "Numeric audit-v4 root differs from the exact production contract." >&2
  exit 2
fi
if [[ ! "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]] || [[ "$(git -C "${REPO}" rev-parse HEAD)" != "${EXPECTED_SHA}" ]] || [[ -n "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]]; then
  echo "Numeric audit-v4 repository is not the exact clean deployed revision." >&2
  exit 2
fi
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]] || (( SLURM_ARRAY_TASK_ID < 0 || SLURM_ARRAY_TASK_ID >= N_TASKS )); then
  echo "Numeric audit-v4 must run as the exact 0-111 array." >&2
  exit 2
fi

SECTOR=$((56 + SLURM_ARRAY_TASK_ID / N_SHARDS))
SHARD_INDEX=$((SLURM_ARRAY_TASK_ID % N_SHARDS))
NATIVE_REGISTRY_DIR="${RUN_ROOT}/frozen/native_registry"
NATIVE_REGISTRY="${NATIVE_REGISTRY_DIR}/native_input_registry.parquet"
NATIVE_REGISTRY_SUMMARY="${NATIVE_REGISTRY_DIR}/native_input_registry.summary.json"
NATIVE_RELEASE_SUMMARY="${NATIVE_REGISTRY_DIR}/native_fullpool_release.summary.json"
REGISTRY_DIR="${RUN_ROOT}/frozen/ssl_registry"
REGISTRY="${REGISTRY_DIR}/teacher_ssl_fullpool_registry.parquet"
REGISTRY_SUMMARY="${REGISTRY_DIR}/teacher_ssl_fullpool_registry.summary.json"
NATIVE_H5="${RUN_ROOT}/native/s${SECTOR}/native_${SHARD_INDEX}.h5"
NATIVE_SUMMARY="${RUN_ROOT}/native/s${SECTOR}/native_${SHARD_INDEX}.summary.json"
NATIVE_SHA="${NATIVE_H5}.sha256"
NUMERIC_GATE_DIR="${RUN_ROOT}/frozen/model_input_numeric_gate_v4"
NUMERIC_RELEASE="${NUMERIC_GATE_DIR}/model_input_numeric_gate.release.json"
REPORT_DIR="${NUMERIC_GATE_DIR}/shard_reports/s${SECTOR}"
REPORT_OUT="${REPORT_DIR}/native_${SHARD_INDEX}.numeric.json"

for input in \
  "${NATIVE_REGISTRY}" \
  "${NATIVE_REGISTRY_SUMMARY}" \
  "${NATIVE_RELEASE_SUMMARY}" \
  "${REGISTRY}" \
  "${REGISTRY_SUMMARY}" \
  "${NATIVE_H5}" \
  "${NATIVE_SUMMARY}" \
  "${NATIVE_SHA}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing numeric audit-v4 authority: ${input}" >&2
    exit 3
  fi
done
sha256sum -c "${NATIVE_SHA}"
for output in "${REPORT_OUT}" "${REPORT_OUT}.sha256"; do
  if [[ -e "${output}" ]]; then
    echo "Refusing to overwrite an existing numeric audit-v4 artifact: ${output}" >&2
    exit 4
  fi
done
if [[ -e "${NUMERIC_RELEASE}" || -e "${NUMERIC_RELEASE}.sha256" ]]; then
  echo "Refusing to audit after the numeric-v4 release was published." >&2
  exit 4
fi

mkdir -p "${TWIRL_ROOT}/logs" "${REPORT_DIR}"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

echo "[teacher-ssl-numeric-audit-v4] start=$(date -Is) sector=${SECTOR} shard=${SHARD_INDEX}"
"${PYTHON}" scripts/stage5_validation/audit_teacher_ssl_fullpool_numeric.py \
  --native-registry "${NATIVE_REGISTRY}" \
  --native-registry-summary "${NATIVE_REGISTRY_SUMMARY}" \
  --native-release-summary "${NATIVE_RELEASE_SUMMARY}" \
  --registry "${REGISTRY}" \
  --registry-summary "${REGISTRY_SUMMARY}" \
  --sector "${SECTOR}" \
  --shard-index "${SHARD_INDEX}" \
  --expected-code-revision "${EXPECTED_SHA}" \
  --report-out "${REPORT_OUT}"
test -s "${REPORT_OUT}"
test -s "${REPORT_OUT}.sha256"
echo "[teacher-ssl-numeric-audit-v4] complete=$(date -Is) report=${REPORT_OUT}"
