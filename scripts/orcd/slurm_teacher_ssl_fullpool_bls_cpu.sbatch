#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-ssl-fullpool-bls
#SBATCH --array=0-55%28
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.err

set -euo pipefail

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:-${TWIRL_ROOT}/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56/bin/python}"
RUN_ROOT="${TWIRL_SSL_FULLPOOL_RUN_ROOT:-${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp}"
POOL_DIR="${TWIRL_SSL_FULLPOOL_POOL_DIR:-${RUN_ROOT}/frozen/full_pool}"
POOL_SUMMARY="${POOL_DIR}/teacher_ssl_full_pool_manifest.summary.json"
CADENCE_ROOT="${TWIRL_SSL_FULLPOOL_CADENCE_ROOT:-${TWIRL_ROOT}/staging/teacher_v3_s56_s62_pdo_inputs_33f5855e_extract_orcdlane_20260726T0143/a2v1_cadence_reference}"
OUT_ROOT="${TWIRL_SSL_FULLPOOL_BLS_ROOT:-${RUN_ROOT}/real_bls}"
N_SHARDS=8

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
if (( TASK_ID < 0 || TASK_ID >= 7 * N_SHARDS )); then
  echo "Array task ${TASK_ID} is outside the locked 0-$((7 * N_SHARDS - 1)) range" >&2
  exit 3
fi
SECTOR=$((56 + TASK_ID / N_SHARDS))
SHARD_INDEX=$((TASK_ID % N_SHARDS))
CADENCE_REFERENCE="${CADENCE_ROOT}/s${SECTOR}/cadence_reference.csv"
CADENCE_MANIFEST="${CADENCE_ROOT}/s${SECTOR}/cadence_reference.manifest.json"
OUT_DIR="${OUT_ROOT}/s${SECTOR}/shards"

for input in "${POOL_SUMMARY}" "${CADENCE_REFERENCE}" "${CADENCE_MANIFEST}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing full-pool BLS input: ${input}" >&2
    exit 4
  fi
done

mapfile -t pool_fields < <(
  "${PYTHON}" -c '
import json
import sys
from pathlib import Path

summary = json.loads(Path(sys.argv[1]).read_text())
sector = int(sys.argv[2])
if summary.get("passed") is not True:
    raise SystemExit("full-pool manifest did not pass")
if summary.get("schema_version") != "twirl_teacher_ssl_full_pool_summary_v1":
    raise SystemExit("full-pool manifest schema mismatch")
exports = {
    int(item["sector"]): item
    for item in summary["inputs"]["compact_exports"]
}
item = exports[sector]
allowlist = summary["outputs"]["sector_allowlists"][str(sector)]
print(item["compact_h5"]["path"])
print(item["compact_h5"]["sha256"])
print(allowlist["path"])
print(allowlist["sha256"])
' "${POOL_SUMMARY}" "${SECTOR}"
)
if (( ${#pool_fields[@]} != 4 )); then
  echo "Could not resolve the frozen S${SECTOR} pool inputs" >&2
  exit 5
fi
COMPACT="${pool_fields[0]}"
COMPACT_SHA256="${pool_fields[1]}"
ALLOWLIST="${pool_fields[2]}"
ALLOWLIST_SHA256="${pool_fields[3]}"
for input in "${COMPACT}" "${ALLOWLIST}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing frozen S${SECTOR} pool input: ${input}" >&2
    exit 6
  fi
done

ORBITID_POLICY="strict"
if [[ "${SECTOR}" == "62" ]]; then
  ORBITID_POLICY="reference_by_cadence"
fi

mkdir -p "${TWIRL_ROOT}/logs" "${OUT_DIR}"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

echo "[ssl-full-pool-bls] host=$(hostname) sector=${SECTOR} shard=${SHARD_INDEX}/${N_SHARDS} policy=${ORBITID_POLICY} start=$(date -Is)"
"${PYTHON}" scripts/stage5_validation/build_s56_adp_real_bls_peaks.py \
  --compact-lc "${COMPACT}" \
  --expected-compact-lc-sha256 "${COMPACT_SHA256}" \
  --sector "${SECTOR}" \
  --cadence-reference "${CADENCE_REFERENCE}" \
  --cadence-reference-manifest "${CADENCE_MANIFEST}" \
  --target-allowlist "${ALLOWLIST}" \
  --expected-target-allowlist-sha256 "${ALLOWLIST_SHA256}" \
  --orbitid-policy "${ORBITID_POLICY}" \
  --out-dir "${OUT_DIR}" \
  --workers "${SLURM_CPUS_PER_TASK:-8}" \
  --n-periods 50000 \
  --n-peaks 10 \
  --shard-index "${SHARD_INDEX}" \
  --n-shards "${N_SHARDS}" \
  --progress-every 25 \
  --source-product-tag A2v1-fullpool-v1 \
  --resume
echo "[ssl-full-pool-bls] sector=${SECTOR} shard=${SHARD_INDEX} complete=$(date -Is)"
