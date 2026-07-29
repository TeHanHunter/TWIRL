#!/bin/bash
# Audit every eligible model-facing sample in one native-v2 shard.
# The 112-task array maps exactly to sectors 56--62 x shards 0--15.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-ssl-numeric-audit
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH --array=0-111%4
#SBATCH -t 08:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.err

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
NUMERIC_RELEASE="${NUMERIC_GATE_DIR}/model_input_numeric_gate.release.json"
readonly N_SHARDS=16
readonly N_TASKS=112

if [[ "${SOURCE_RUN_ROOT}" != "${EXPECTED_SOURCE_RUN_ROOT}" || "${RUN_ROOT}" != "${EXPECTED_RUN_ROOT}" ]]; then
  echo "Numeric-audit roots differ from the exact Teacher v4-SSL v2 contract." >&2
  exit 2
fi
if [[ ! "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]]; then
  echo "TWIRL_EXPECTED_GIT_SHA must be a full lowercase Git SHA." >&2
  exit 2
fi
if [[ "$(git -C "${REPO}" rev-parse HEAD)" != "${EXPECTED_SHA}" ]] || [[ -n "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]]; then
  echo "Numeric-audit repository is not the clean expected deployed revision." >&2
  exit 2
fi
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]] || (( SLURM_ARRAY_TASK_ID < 0 || SLURM_ARRAY_TASK_ID >= N_TASKS )); then
  echo "Numeric audit must run as the exact 0-111 array." >&2
  exit 2
fi

SECTOR=$((56 + SLURM_ARRAY_TASK_ID / N_SHARDS))
SHARD_INDEX=$((SLURM_ARRAY_TASK_ID % N_SHARDS))
NATIVE_H5="${RUN_ROOT}/native/s${SECTOR}/native_${SHARD_INDEX}.h5"
NATIVE_SUMMARY="${RUN_ROOT}/native/s${SECTOR}/native_${SHARD_INDEX}.summary.json"
REPORT_DIR="${NUMERIC_SHARD_DIR}/s${SECTOR}"
REPORT_OUT="${REPORT_DIR}/native_${SHARD_INDEX}.numeric.json"

for input in \
  "${NATIVE_REGISTRY}" \
  "${NATIVE_REGISTRY_SUMMARY}" \
  "${NATIVE_RELEASE_SUMMARY}" \
  "${REGISTRY}" \
  "${REGISTRY_SUMMARY}" \
  "${NATIVE_H5}" \
  "${NATIVE_SUMMARY}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing numeric-audit authority: ${input}" >&2
    exit 3
  fi
done
if [[ -e "${NUMERIC_RELEASE}" || -e "${NUMERIC_RELEASE}.sha256" ]]; then
  echo "Refusing to audit after the numeric release was published." >&2
  exit 4
fi
if [[ -e "${REPORT_OUT}.sha256" ]]; then
  if [[ ! -s "${REPORT_OUT}" ]]; then
    echo "Numeric shard sidecar exists without its report: ${REPORT_OUT}" >&2
    exit 4
  fi
  "${PYTHON}" -c \
    'import hashlib,json,pathlib,sys; report=pathlib.Path(sys.argv[1]); sidecar=pathlib.Path(sys.argv[2]); payload=json.loads(report.read_text(encoding="utf-8")); digest=hashlib.sha256(report.read_bytes()).hexdigest(); expected=f"{digest}  {report.name}\n"; assert sidecar.read_text(encoding="ascii")==expected; assert payload.get("passed") is True; assert payload.get("code_revision")==sys.argv[3]; assert payload.get("sector")==int(sys.argv[4]); assert payload.get("shard_index")==int(sys.argv[5])' \
    "${REPORT_OUT}" "${REPORT_OUT}.sha256" "${EXPECTED_SHA}" \
    "${SECTOR}" "${SHARD_INDEX}"
  echo "[teacher-ssl-numeric-audit] already complete report=${REPORT_OUT}"
  exit 0
fi
if [[ -e "${REPORT_OUT}" ]]; then
  echo "[teacher-ssl-numeric-audit] recovering incomplete report=${REPORT_OUT}"
fi

mkdir -p "${TWIRL_ROOT}/logs" "${REPORT_DIR}"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

echo "[teacher-ssl-numeric-audit] start=$(date -Is) sector=${SECTOR} shard=${SHARD_INDEX}"
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
echo "[teacher-ssl-numeric-audit] complete=$(date -Is) report=${REPORT_OUT}"
