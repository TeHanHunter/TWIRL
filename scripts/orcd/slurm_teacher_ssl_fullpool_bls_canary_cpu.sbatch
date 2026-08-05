#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-ssl-bls-canary
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:-${TWIRL_ROOT}/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56/bin/python}"
RUN_ROOT="${TWIRL_SSL_FULLPOOL_RUN_ROOT:-${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp}"
POOL_DIR="${TWIRL_SSL_FULLPOOL_POOL_DIR:-${RUN_ROOT}/frozen/full_pool}"
POOL_SUMMARY="${POOL_DIR}/teacher_ssl_full_pool_manifest.summary.json"
CADENCE_ROOT="${TWIRL_SSL_FULLPOOL_CADENCE_ROOT:-${TWIRL_ROOT}/staging/teacher_v3_s56_s62_pdo_inputs_33f5855e_extract_orcdlane_20260726T0143/a2v1_cadence_reference}"
PRODUCTION_ROOT="${TWIRL_SSL_FULLPOOL_BLS_ROOT:-${RUN_ROOT}/real_bls}"
CANARY_ROOT="${TWIRL_SSL_FULLPOOL_BLS_CANARY_ROOT:-${RUN_ROOT}/canary/real_bls_first100}"
SECTOR="${TWIRL_SSL_FULLPOOL_BLS_CANARY_SECTOR:-62}"
N_TARGETS="${TWIRL_SSL_FULLPOOL_BLS_CANARY_TARGETS:-100}"
PRODUCTION_ROOT="$(realpath -m -- "${PRODUCTION_ROOT}")"
CANARY_ROOT="$(realpath -m -- "${CANARY_ROOT}")"

if [[ ! "${SECTOR}" =~ ^[0-9]+$ ]] || (( 10#${SECTOR} < 56 || 10#${SECTOR} > 62 )); then
  echo "Canary sector must be an integer in [56,62]." >&2
  exit 3
fi
if [[ ! "${N_TARGETS}" =~ ^[0-9]+$ ]] || (( 10#${N_TARGETS} < 2 || 10#${N_TARGETS} > 100 )); then
  echo "Canary target count must be an integer in [2,100]." >&2
  exit 3
fi
if [[ "${CANARY_ROOT}" == "${PRODUCTION_ROOT}" ]]; then
  echo "Canary and production BLS roots must be distinct." >&2
  exit 3
fi

CADENCE_REFERENCE="${CADENCE_ROOT}/s${SECTOR}/cadence_reference.csv"
CADENCE_MANIFEST="${CADENCE_ROOT}/s${SECTOR}/cadence_reference.manifest.json"
INPUT_DIR="${CANARY_ROOT}/inputs"
OUT_DIR="${CANARY_ROOT}/s${SECTOR}"
for input in "${POOL_SUMMARY}" "${CADENCE_REFERENCE}" "${CADENCE_MANIFEST}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing full-pool BLS canary input: ${input}" >&2
    exit 4
  fi
done

mkdir -p "${TWIRL_ROOT}/logs" "${INPUT_DIR}" "${OUT_DIR}"
mapfile -t canary_fields < <(
  "${PYTHON}" -c '
import csv
import hashlib
import json
import sys
from pathlib import Path

summary_path = Path(sys.argv[1]).resolve()
sector = int(sys.argv[2])
n_targets = int(sys.argv[3])
input_dir = Path(sys.argv[4]).resolve()
def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
summary = json.loads(summary_path.read_text())
if summary.get("passed") is not True:
    raise SystemExit("full-pool manifest did not pass")
if summary.get("schema_version") != "twirl_teacher_ssl_full_pool_summary_v1":
    raise SystemExit("full-pool manifest schema mismatch")
exports = {int(item["sector"]): item for item in summary["inputs"]["compact_exports"]}
compact = exports[sector]["compact_h5"]
full_allowlist = summary["outputs"]["sector_allowlists"][str(sector)]
compact_path = Path(compact["path"]).resolve()
allowlist_path = Path(full_allowlist["path"]).resolve()
if not compact_path.is_file():
    raise SystemExit(f"missing compact HDF5: {compact_path}")
if not allowlist_path.is_file():
    raise SystemExit(f"missing full allowlist: {allowlist_path}")
if file_sha256(allowlist_path) != full_allowlist["sha256"]:
    raise SystemExit(f"full allowlist hash mismatch: {allowlist_path}")
with allowlist_path.open(newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle)
    if reader.fieldnames != ["tic"]:
        raise SystemExit("full allowlist must contain exactly one tic column")
    tics = [int(row["tic"]) for row in reader]
if tics != sorted(set(tics)):
    raise SystemExit("full allowlist TICs are not sorted and unique")
if len(tics) < n_targets:
    raise SystemExit(
        f"full allowlist has {len(tics)} TICs, fewer than requested {n_targets}"
    )
selected = tics[:n_targets]
payload = ("tic\n" + "".join(f"{tic}\n" for tic in selected)).encode("ascii")
subset_path = input_dir / f"s{sector}_first{n_targets}_tics.csv"
temporary = subset_path.with_suffix(subset_path.suffix + ".tmp")
temporary.write_bytes(payload)
temporary.replace(subset_path)
subset_sha256 = hashlib.sha256(payload).hexdigest()
manifest = {
    "schema_version": "twirl_teacher_ssl_fullpool_bls_canary_allowlist_v1",
    "selection": "first_n_sorted_tics",
    "sector": sector,
    "n_targets": n_targets,
    "parent_allowlist": {
        "path": str(allowlist_path),
        "sha256": full_allowlist["sha256"],
        "n_tics": int(full_allowlist["n_tics"]),
    },
    "subset_allowlist": {
        "path": str(subset_path),
        "sha256": subset_sha256,
        "n_tics": len(selected),
    },
}
manifest_path = subset_path.with_suffix(".manifest.json")
manifest_tmp = manifest_path.with_suffix(manifest_path.suffix + ".tmp")
manifest_tmp.write_text(
    json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
    encoding="utf-8",
)
manifest_tmp.replace(manifest_path)
print(compact_path)
print(compact["sha256"])
print(subset_path)
print(subset_sha256)
' "${POOL_SUMMARY}" "${SECTOR}" "${N_TARGETS}" "${INPUT_DIR}"
)
if (( ${#canary_fields[@]} != 4 )); then
  echo "Could not resolve the frozen S${SECTOR} canary inputs." >&2
  exit 5
fi
COMPACT="${canary_fields[0]}"
COMPACT_SHA256="${canary_fields[1]}"
ALLOWLIST="${canary_fields[2]}"
ALLOWLIST_SHA256="${canary_fields[3]}"

ORBITID_POLICY="strict"
if [[ "${SECTOR}" == "62" ]]; then
  ORBITID_POLICY="reference_by_cadence"
fi

cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

echo "[ssl-full-pool-bls-canary] host=$(hostname) sector=${SECTOR} targets=${N_TARGETS} policy=${ORBITID_POLICY} out=${OUT_DIR} start=$(date -Is)"
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
  --n-shards 1 \
  --progress-every 10 \
  --source-product-tag A2v1-fullpool-v1-canary \
  --resume
echo "[ssl-full-pool-bls-canary] complete=$(date -Is)"
