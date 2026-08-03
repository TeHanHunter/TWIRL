#!/bin/bash
# Build one detector-consistent raw-v1 native-v4 shard.  The compact export is
# admitted only as the frozen cadence/time inventory reference.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-ssl-v4-native
#SBATCH --array=0-111%28
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.err

set -euo pipefail

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:?TWIRL_ORCD_REPO is required}"
EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA:?TWIRL_EXPECTED_GIT_SHA is required}"
SOURCE_RUN_ROOT="${TWIRL_SSL_FULLPOOL_SOURCE_RUN_ROOT:-${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp}"
RUN_ROOT="${TWIRL_SSL_FULLPOOL_RUN_ROOT:?TWIRL_SSL_FULLPOOL_RUN_ROOT is required}"
EXPECTED_SOURCE_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp"
EXPECTED_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v4r2_s56_s62_a2v1_raw_v1_detector_consistent_adp_bls_eligible"
PYTHON="${TWIRL_ORCD_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python}"
N_SHARDS=16
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

if (( TASK_ID < 0 || TASK_ID >= 7 * N_SHARDS )); then
  echo "Array task ${TASK_ID} is outside the locked 0-$((7 * N_SHARDS - 1)) range" >&2
  exit 2
fi
if [[ "$(realpath -m -- "${SOURCE_RUN_ROOT}")" != "${EXPECTED_SOURCE_RUN_ROOT}" || "$(realpath -m -- "${RUN_ROOT}")" != "${EXPECTED_RUN_ROOT}" ]]; then
  echo "Native-v4 source/output roots differ from the production contract." >&2
  exit 2
fi
if [[ ! "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]] || [[ "$(git -C "${REPO}" rev-parse HEAD)" != "${EXPECTED_SHA}" ]] || [[ -n "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]]; then
  echo "Native-v4 repository is not the exact clean deployed revision." >&2
  exit 2
fi

SECTOR=$((56 + TASK_ID / N_SHARDS))
SHARD_INDEX=$((TASK_ID % N_SHARDS))
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
NATIVE_ROOT="${TWIRL_SSL_FULLPOOL_NATIVE_OUT_ROOT:-${RUN_ROOT}/native}"
OUT_H5="${NATIVE_ROOT}/s${SECTOR}/native_${SHARD_INDEX}.h5"
SUMMARY_OUT="${OUT_H5%.h5}.summary.json"
SHA_OUT="${OUT_H5}.sha256"

if [[ "$(realpath -m -- "${NATIVE_ROOT}")" != "${RUN_ROOT}/"* ]]; then
  echo "Native-v4 output must remain beneath the exact fresh v4 run root." >&2
  exit 2
fi

mapfile -t compact_fields < <(
  "${PYTHON}" -c '
import json,sys
d=json.load(open(sys.argv[1])); s=int(sys.argv[2])
rows=[v for v in d["inputs"]["compact_exports"] if int(v["sector"]) == s]
assert len(rows) == 1
print(rows[0]["compact_h5"]["path"])
' "${POOL_SUMMARY}" "${SECTOR}"
)
if (( ${#compact_fields[@]} != 1 )); then
  echo "Could not resolve the frozen compact cadence/time reference." >&2
  exit 3
fi
COMPACT_ADP="${compact_fields[0]}"
ORBITID_POLICY="strict"
if [[ "${SECTOR}" == "62" ]]; then
  ORBITID_POLICY="reference_by_cadence"
fi

for input in "${POOL}" "${POOL_SUMMARY}" "${ALLOWLIST}" \
  "${ELIGIBILITY_EXCLUSIONS}" "${ELIGIBILITY_SUMMARY}" \
  "${RAW_SOURCE}" "${RAW_SOURCE_SUMMARY}" "${RAW_EXPORT_COMPLETE}" \
  "${RAW_TRANSFER_VALIDATION}" "${COMPACT_ADP}" "${CADENCE_TABLE}" \
  "${CADENCE_MANIFEST}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing native-v4 input: ${input}" >&2
    exit 3
  fi
done
for output in "${OUT_H5}" "${SUMMARY_OUT}" "${SHA_OUT}"; do
  if [[ -e "${output}" ]]; then
    echo "Refusing to overwrite or reuse native-v4 artifact: ${output}" >&2
    exit 4
  fi
done

mkdir -p "${TWIRL_ROOT}/logs" "$(dirname "${OUT_H5}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

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

test -s "${OUT_H5}" && test -s "${SUMMARY_OUT}"
SHA_TMP="${SHA_OUT}.tmp.$$"
sha256sum "${OUT_H5}" > "${SHA_TMP}"
mv "${SHA_TMP}" "${SHA_OUT}"
echo "[fullpool-native-v4] S${SECTOR} shard ${SHARD_INDEX}/${N_SHARDS}: ${OUT_H5}"
