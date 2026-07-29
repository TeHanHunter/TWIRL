#!/bin/bash
# Build one fresh effective-quality-aware ADP native-v3 shard.
# The immutable raw-v1 release is the only photometric source.  The compact
# export supplies cadence/time/geometry alignment only; native-v3 refuses its
# stored ADP photometry.
# Example production submit:
#   sbatch --array=0-15%4 --export=TWIRL_ORCD_REPO=...,TWIRL_EXPECTED_GIT_SHA=...,\
# TWIRL_SSL_FULLPOOL_V1_RUN_ROOT=...,TWIRL_SSL_FULLPOOL_V3_RUN_ROOT=...,\
# TWIRL_FULLPOOL_SECTOR=60,TWIRL_FULLPOOL_NATIVE_SHARDS=16 \
# scripts/orcd/slurm_teacher_ssl_fullpool_v3_native_cpu.sbatch
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-ssl-v3-native
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.err

set -euo pipefail

readonly TWIRL_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly EXPECTED_SOURCE_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp"
readonly EXPECTED_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v3_s56_s62_a2v1_effective_quality_adp_bls_eligible"
: "${TWIRL_ORCD_REPO:?controller must export TWIRL_ORCD_REPO}"
: "${TWIRL_EXPECTED_GIT_SHA:?controller must export TWIRL_EXPECTED_GIT_SHA}"
: "${TWIRL_SSL_FULLPOOL_V1_RUN_ROOT:?controller must export the immutable raw-v1 source root}"
: "${TWIRL_SSL_FULLPOOL_V3_RUN_ROOT:?controller must export the exact fresh v3 run root}"

REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}"
SOURCE_RUN_ROOT="$(realpath -m -- "${TWIRL_SSL_FULLPOOL_V1_RUN_ROOT}")"
RUN_ROOT="$(realpath -m -- "${TWIRL_SSL_FULLPOOL_V3_RUN_ROOT}")"
PYTHON="${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python"
SECTOR="${TWIRL_FULLPOOL_SECTOR:?set TWIRL_FULLPOOL_SECTOR to 56-62}"
N_SHARDS="${TWIRL_FULLPOOL_NATIVE_SHARDS:?set TWIRL_FULLPOOL_NATIVE_SHARDS to 16}"
SHARD_INDEX="${SLURM_ARRAY_TASK_ID:-${TWIRL_FULLPOOL_SHARD_INDEX:-}}"

if [[ -z "${SHARD_INDEX}" ]]; then
  echo "Submit as a Slurm array or set TWIRL_FULLPOOL_SHARD_INDEX for a canary." >&2
  exit 2
fi
if (( SECTOR < 56 || SECTOR > 62 )); then
  echo "Native-v3 preparation is bounded to sectors 56-62." >&2
  exit 2
fi
if [[ "${N_SHARDS}" != "16" ]] || (( SHARD_INDEX < 0 || SHARD_INDEX >= 16 )); then
  echo "Native-v3 is locked to exactly 16 shards per sector." >&2
  exit 2
fi
if [[ "${SOURCE_RUN_ROOT}" != "${EXPECTED_SOURCE_RUN_ROOT}" || "${RUN_ROOT}" != "${EXPECTED_RUN_ROOT}" ]]; then
  echo "Native-v3 source/output roots differ from the exact production contract." >&2
  exit 2
fi
if [[ ! "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]]; then
  echo "TWIRL_EXPECTED_GIT_SHA must be a full lowercase Git SHA." >&2
  exit 2
fi
if [[ "$(git -C "${REPO}" rev-parse HEAD)" != "${EXPECTED_SHA}" ]] || [[ -n "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]]; then
  echo "Native-v3 repository is not the exact clean deployed revision." >&2
  exit 2
fi

POOL_DIR="${SOURCE_RUN_ROOT}/frozen/full_pool"
POOL="${POOL_DIR}/teacher_ssl_full_pool_observations.csv"
POOL_SUMMARY="${POOL_DIR}/teacher_ssl_full_pool_manifest.summary.json"
ALLOWLIST="${POOL_DIR}/allowlists/s${SECTOR}_tics.csv"
RAW_DIR="${SOURCE_RUN_ROOT}/native_preparation/raw_sources"
RAW_SOURCE="${RAW_DIR}/s${SECTOR}/raw_source_${SHARD_INDEX}.h5"
RAW_SOURCE_SUMMARY="${RAW_DIR}/s${SECTOR}/raw_source_${SHARD_INDEX}.summary.json"
RAW_EXPORT_COMPLETE="${RAW_DIR}/raw_export.complete.json"
RAW_TRANSFER_VALIDATION="${RAW_DIR}/s${SECTOR}/transfer_validation.json"
ELIGIBILITY_DIR="${RUN_ROOT}/frozen/native_model_eligibility"
ELIGIBILITY_EXCLUSIONS="${ELIGIBILITY_DIR}/native_model_exclusions.csv"
ELIGIBILITY_SUMMARY="${ELIGIBILITY_DIR}/native_model_eligibility.summary.json"
CADENCE_ROOT="${TWIRL_ROOT}/staging/teacher_v3_s56_s62_pdo_inputs_33f5855e_extract_orcdlane_20260726T0143/a2v1_cadence_reference"
CADENCE_TABLE="${CADENCE_ROOT}/s${SECTOR}/cadence_reference.csv"
CADENCE_MANIFEST="${CADENCE_ROOT}/s${SECTOR}/cadence_reference.manifest.json"
OUT_H5="${TWIRL_FULLPOOL_NATIVE_OUT:-${RUN_ROOT}/native/s${SECTOR}/native_${SHARD_INDEX}.h5}"
OUT_H5="$(realpath -m -- "${OUT_H5}")"
SUMMARY_OUT="${OUT_H5%.h5}.summary.json"
SHA_OUT="${OUT_H5}.sha256"

if [[ "${OUT_H5}" != "${RUN_ROOT}/"* ]]; then
  echo "Native-v3 output must remain beneath the exact fresh v3 run root." >&2
  exit 2
fi

ORBITID_POLICY="strict"
if [[ "${SECTOR}" == "62" ]]; then
  ORBITID_POLICY="reference_by_cadence"
fi

for authority in \
  "${POOL}" \
  "${POOL_SUMMARY}" \
  "${ALLOWLIST}" \
  "${ELIGIBILITY_EXCLUSIONS}" \
  "${ELIGIBILITY_SUMMARY}" \
  "${RAW_EXPORT_COMPLETE}" \
  "${RAW_TRANSFER_VALIDATION}" \
  "${RAW_SOURCE_SUMMARY}"; do
  if [[ ! -s "${authority}" ]]; then
    echo "Missing native-v3 authority: ${authority}" >&2
    exit 3
  fi
done

# This compact path is checksum-bound by POOL_SUMMARY and is admitted only for
# cadence/time/geometry/internal-quality alignment.  The builder's v3 contract
# requires raw_photometry_only=1 and compact_adp_photometry_reused=0.
COMPACT_ADP="$(
  "${PYTHON}" -c \
    'import json,sys; d=json.load(open(sys.argv[1])); s=int(sys.argv[2]); rows=[v for v in d["inputs"]["compact_exports"] if int(v["sector"]) == s]; assert len(rows) == 1; print(rows[0]["compact_h5"]["path"])' \
    "${POOL_SUMMARY}" "${SECTOR}"
)"
for input in "${RAW_SOURCE}" "${COMPACT_ADP}" "${CADENCE_TABLE}" "${CADENCE_MANIFEST}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing native-v3 input: ${input}" >&2
    exit 3
  fi
done
for output in "${OUT_H5}" "${SUMMARY_OUT}" "${SHA_OUT}"; do
  if [[ -e "${output}" ]]; then
    echo "Refusing to overwrite or resume an existing native-v3 artifact: ${output}" >&2
    exit 4
  fi
done

mkdir -p "${TWIRL_ROOT}/logs" "$(dirname "${OUT_H5}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

"${PYTHON}" scripts/stage5_validation/build_teacher_ssl_fullpool_native_shard.py \
  --sector "${SECTOR}" \
  --frozen-pool "${POOL}" \
  --frozen-pool-summary "${POOL_SUMMARY}" \
  --sector-allowlist "${ALLOWLIST}" \
  --eligibility-exclusions "${ELIGIBILITY_EXCLUSIONS}" \
  --eligibility-summary "${ELIGIBILITY_SUMMARY}" \
  --raw-source-h5 "${RAW_SOURCE}" \
  --raw-source-summary "${RAW_SOURCE_SUMMARY}" \
  --raw-export-complete "${RAW_EXPORT_COMPLETE}" \
  --raw-transfer-validation "${RAW_TRANSFER_VALIDATION}" \
  --compact-adp-h5 "${COMPACT_ADP}" \
  --cadence-reference-table "${CADENCE_TABLE}" \
  --cadence-reference-manifest "${CADENCE_MANIFEST}" \
  --out-h5 "${OUT_H5}" \
  --shard-index "${SHARD_INDEX}" \
  --n-shards "${N_SHARDS}" \
  --n-periods 4096 \
  --orbitid-policy "${ORBITID_POLICY}" \
  --expected-code-revision "${EXPECTED_SHA}"

test -s "${OUT_H5}"
test -s "${SUMMARY_OUT}"
SHA_TMP="${SHA_OUT}.tmp.$$"
sha256sum "${OUT_H5}" > "${SHA_TMP}"
mv "${SHA_TMP}" "${SHA_OUT}"
echo "[fullpool-native-v3] S${SECTOR} shard ${SHARD_INDEX}/${N_SHARDS}: ${OUT_H5}"
