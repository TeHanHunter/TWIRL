#!/bin/bash
# Merge the complete 112-shard numeric audit into one immutable release.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-ssl-numeric-release
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

TWIRL_ROOT="/orcd/data/mki_aryeh/001/twirl"
: "${TWIRL_ORCD_REPO:?controller must export TWIRL_ORCD_REPO}"
: "${TWIRL_EXPECTED_GIT_SHA:?controller must export TWIRL_EXPECTED_GIT_SHA}"
REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}"
PYTHON="${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python"
EXPECTED_SOURCE_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp"
EXPECTED_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v2_s56_s62_a2v1_current_adp_bls_eligible"
SOURCE_RUN_ROOT="$(realpath -m -- "${TWIRL_SSL_FULLPOOL_V1_RUN_ROOT:-${EXPECTED_SOURCE_RUN_ROOT}}")"
RUN_ROOT="$(realpath -m -- "${TWIRL_SSL_FULLPOOL_V2_RUN_ROOT:-${EXPECTED_RUN_ROOT}}")"
NATIVE_REGISTRY_DIR="${RUN_ROOT}/frozen/native_registry"
NATIVE_REGISTRY="${NATIVE_REGISTRY_DIR}/native_input_registry.parquet"
NATIVE_REGISTRY_SUMMARY="${NATIVE_REGISTRY_DIR}/native_input_registry.summary.json"
NATIVE_RELEASE_SUMMARY="${NATIVE_REGISTRY_DIR}/native_fullpool_release.summary.json"
REGISTRY_DIR="${RUN_ROOT}/frozen/ssl_registry"
REGISTRY="${REGISTRY_DIR}/teacher_ssl_fullpool_registry.parquet"
REGISTRY_SUMMARY="${REGISTRY_DIR}/teacher_ssl_fullpool_registry.summary.json"
NUMERIC_GATE_DIR="${RUN_ROOT}/frozen/model_input_numeric_gate_v1"
NUMERIC_SHARD_DIR="${NUMERIC_GATE_DIR}/shard_reports"
NUMERIC_AUDIT="${NUMERIC_GATE_DIR}/model_input_numeric_audit.parquet"
NUMERIC_RELEASE="${NUMERIC_GATE_DIR}/model_input_numeric_gate.release.json"

if [[ "${SOURCE_RUN_ROOT}" != "${EXPECTED_SOURCE_RUN_ROOT}" || "${RUN_ROOT}" != "${EXPECTED_RUN_ROOT}" ]]; then
  echo "Numeric-release roots differ from the exact Teacher v4-SSL v2 contract." >&2
  exit 2
fi
if [[ ! "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]]; then
  echo "TWIRL_EXPECTED_GIT_SHA must be a full lowercase Git SHA." >&2
  exit 2
fi
if [[ "$(git -C "${REPO}" rev-parse HEAD)" != "${EXPECTED_SHA}" ]] || [[ -n "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]]; then
  echo "Numeric-release repository is not the clean expected deployed revision." >&2
  exit 2
fi

for input in \
  "${NATIVE_REGISTRY}" \
  "${NATIVE_REGISTRY_SUMMARY}" \
  "${NATIVE_RELEASE_SUMMARY}" \
  "${REGISTRY}" \
  "${REGISTRY_SUMMARY}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing numeric-release authority: ${input}" >&2
    exit 3
  fi
done
if [[ -e "${NUMERIC_RELEASE}.sha256" ]]; then
  for completed in \
    "${NUMERIC_AUDIT}" \
    "${NUMERIC_AUDIT}.sha256" \
    "${NUMERIC_RELEASE}"; do
    if [[ ! -s "${completed}" ]]; then
      echo "Numeric release completion sidecar has an incomplete artifact set." >&2
      exit 4
    fi
  done
  export PYTHONPATH="${REPO}/src"
  "${PYTHON}" -c \
    'from pathlib import Path; import sys; from twirl.vetting.ssl_full_pool_numeric import validate_numeric_gate_release; validate_numeric_gate_release(Path(sys.argv[1]), expected_code_revision=sys.argv[2])' \
    "${NUMERIC_RELEASE}" "${EXPECTED_SHA}"
  echo "[teacher-ssl-numeric-release] already complete release=${NUMERIC_RELEASE}"
  exit 0
fi
if [[
  -e "${NUMERIC_AUDIT}" ||
  -e "${NUMERIC_AUDIT}.sha256" ||
  -e "${NUMERIC_RELEASE}"
]]; then
  echo "[teacher-ssl-numeric-release] recovering incomplete publication"
fi

reports=()
for sector in {56..62}; do
  for shard_index in {0..15}; do
    report="${NUMERIC_SHARD_DIR}/s${sector}/native_${shard_index}.numeric.json"
    if [[ ! -s "${report}" || ! -s "${report}.sha256" ]]; then
      echo "Missing numeric shard report: ${report}" >&2
      exit 3
    fi
    reports+=("${report}")
  done
done
if [[ "${#reports[@]}" -ne 112 ]]; then
  echo "Numeric release requires exactly 112 shard reports." >&2
  exit 3
fi

mkdir -p "${TWIRL_ROOT}/logs" "${NUMERIC_GATE_DIR}"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

echo "[teacher-ssl-numeric-release] start=$(date -Is)"
"${PYTHON}" scripts/stage5_validation/merge_teacher_ssl_fullpool_numeric.py \
  --native-registry "${NATIVE_REGISTRY}" \
  --native-registry-summary "${NATIVE_REGISTRY_SUMMARY}" \
  --native-release-summary "${NATIVE_RELEASE_SUMMARY}" \
  --registry "${REGISTRY}" \
  --registry-summary "${REGISTRY_SUMMARY}" \
  --shard-reports "${reports[@]}" \
  --expected-code-revision "${EXPECTED_SHA}" \
  --release-out "${NUMERIC_RELEASE}"
test -s "${NUMERIC_RELEASE}"
"${PYTHON}" -c \
  'from pathlib import Path; import sys; from twirl.vetting.ssl_full_pool_numeric import validate_numeric_gate_release; validate_numeric_gate_release(Path(sys.argv[1]), expected_code_revision=sys.argv[2])' \
  "${NUMERIC_RELEASE}" "${EXPECTED_SHA}"
echo "[teacher-ssl-numeric-release] complete=$(date -Is) release=${NUMERIC_RELEASE}"
