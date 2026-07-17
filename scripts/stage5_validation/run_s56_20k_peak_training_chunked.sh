#!/usr/bin/env bash
# Restartable peak-level BLS truth table build for the S56 20k injection set.
set -euo pipefail

REPO_ROOT="${TWIRL_REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
cd "${REPO_ROOT}"

if [[ -z "${PYTHON:-}" ]]; then
  if [[ -x /sw/qlp-environment/.venv/bin/python ]]; then
    PYTHON=/sw/qlp-environment/.venv/bin/python
  else
    PYTHON=python3
  fi
fi

if [[ -d /pdo/app/anaconda/anaconda2-4.4.0/lib && -d /sw/python-versions/python-3.11.9/lib ]]; then
  export LD_LIBRARY_PATH="/pdo/app/anaconda/anaconda2-4.4.0/lib:/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
fi
export PYTHONPATH=src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

INJECTION_H5="${INJECTION_H5:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_20k_predetrend_batman_periodradius_bright_balanced/injected_lightcurves.h5}"
OUT_ROOT="${OUT_ROOT:-reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training}"
CHUNK_DIR="${CHUNK_DIR:-${OUT_ROOT}/chunk_ids}"
RUN_CHUNK_DIR="${RUN_CHUNK_DIR:-${OUT_ROOT}/chunks}"
FINAL_TABLE="${FINAL_TABLE:-${OUT_ROOT}/s56_20k_injection_bls_peaks_chunked.csv}"
VERIFY_JSON="${VERIFY_JSON:-${FINAL_TABLE%.*}_verification.json}"
CHUNK_SIZE="${TWIRL_PEAK_CHUNK_SIZE:-250}"
CHUNK_JOBS="${TWIRL_PEAK_CHUNK_JOBS:-1}"
WORKERS="${TWIRL_PEAK_WORKERS:-8}"
N_INJECTIONS="${TWIRL_N_INJECTIONS:-0}"
N_PEAKS="${TWIRL_BLS_N_PEAKS:-20}"
N_PERIODS="${TWIRL_BLS_N_PERIODS:-200000}"
MIN_APERTURES="${TWIRL_PEAK_MIN_APERTURES:-2}"
MIN_POSITIVE_PEAK_RANKS="${TWIRL_PEAK_MIN_POSITIVE_RANKS:-10}"
REQUIRE_SIGNAL_PEAKS="${TWIRL_PEAK_REQUIRE_SIGNAL_PEAKS:-1}"
TABLE_NAME="injection_bls_peaks.csv"

log() { echo "[$(date -Is)] [peak-chunked] $*"; }

mkdir -p "${CHUNK_DIR}" "${RUN_CHUNK_DIR}" logs

log "repo=${REPO_ROOT}"
log "python=${PYTHON}"
log "injection_h5=${INJECTION_H5}"
log "chunk_size=${CHUNK_SIZE} chunk_jobs=${CHUNK_JOBS} workers=${WORKERS}"
log "n_injections=${N_INJECTIONS} n_periods=${N_PERIODS} n_peaks=${N_PEAKS}"

"${PYTHON}" - "${INJECTION_H5}" "${CHUNK_DIR}" "${CHUNK_SIZE}" "${N_INJECTIONS}" <<'PY'
from pathlib import Path
import sys

import h5py

injection_h5 = Path(sys.argv[1])
chunk_dir = Path(sys.argv[2])
chunk_size = int(sys.argv[3])
n_injections = int(sys.argv[4])

with h5py.File(injection_h5, "r") as h5:
    keys = sorted(h5["injections"].keys())
if n_injections > 0:
    keys = keys[:n_injections]

for old in chunk_dir.glob("chunk_*.txt"):
    old.unlink()

for start in range(0, len(keys), chunk_size):
    chunk = keys[start:start + chunk_size]
    path = chunk_dir / f"chunk_{start // chunk_size:03d}.txt"
    path.write_text("\n".join(chunk) + "\n", encoding="utf-8")

(chunk_dir / "manifest.txt").write_text(
    f"n_injections={len(keys)}\nchunk_size={chunk_size}\nn_chunks={(len(keys) + chunk_size - 1) // chunk_size}\n",
    encoding="utf-8",
)
print((chunk_dir / "manifest.txt").read_text(), end="")
PY

run_one_chunk() {
  local id_file="$1"
  local chunk_name chunk_out summary_path run_log
  chunk_name="$(basename "${id_file}" .txt)"
  chunk_out="${RUN_CHUNK_DIR}/${chunk_name}"
  summary_path="${chunk_out}/injection_bls_peaks_summary.json"
  run_log="${chunk_out}/run.log"
  if [[ -f "${summary_path}" ]]; then
    log "skip ${chunk_name}: ${summary_path} exists"
    return 0
  fi
  mkdir -p "${chunk_out}"
  log "run ${chunk_name} with ${WORKERS} workers; log=${run_log}"
  if ! "${PYTHON}" scripts/stage5_validation/build_injection_peak_training_table.py \
    --injection-h5 "${INJECTION_H5}" \
    --out-table "${chunk_out}/${TABLE_NAME}" \
    --apertures DET_FLUX_ADP_SML DET_FLUX_SML \
    --limit-keys-file "${id_file}" \
    --n-injections 0 \
    --n-peaks "${N_PEAKS}" \
    --n-periods "${N_PERIODS}" \
    --workers "${WORKERS}" \
    --overwrite >"${run_log}" 2>&1; then
    log "failed ${chunk_name}; tail ${run_log}"
    tail -n 40 "${run_log}" || true
    return 1
  fi
  log "done ${chunk_name}"
}
export -f run_one_chunk log
export PYTHON INJECTION_H5 RUN_CHUNK_DIR TABLE_NAME WORKERS N_PEAKS N_PERIODS

if [[ "${CHUNK_JOBS}" -le 1 ]]; then
  for id_file in "${CHUNK_DIR}"/chunk_*.txt; do
    run_one_chunk "${id_file}"
  done
else
  find "${CHUNK_DIR}" -maxdepth 1 -name 'chunk_*.txt' -print0 |
    sort -z |
    xargs -0 -n 1 -P "${CHUNK_JOBS}" bash -c 'run_one_chunk "$1"' _
fi

EXPECT_N="$(awk -F= '$1 == "n_injections" {print $2}' "${CHUNK_DIR}/manifest.txt")"
log "merging chunks expect_n=${EXPECT_N}"
"${PYTHON}" scripts/stage5_validation/merge_injection_peak_training_chunks.py \
  --chunk-root "${RUN_CHUNK_DIR}" \
  --out-table "${FINAL_TABLE}" \
  --table-name "${TABLE_NAME}" \
  --expect-injections "${EXPECT_N}"

MIN_INJECTIONS="${TWIRL_PEAK_MIN_INJECTIONS:-${EXPECT_N}}"
MIN_CANDIDATE_ROWS="${TWIRL_PEAK_MIN_CANDIDATE_ROWS:-$(( EXPECT_N * 5 ))}"
VERIFY_ARGS=(
  scripts/stage5_validation/verify_injection_peak_training_table.py
  --peak-table "${FINAL_TABLE}"
  --out-json "${VERIFY_JSON}"
  --min-injections "${MIN_INJECTIONS}"
  --min-candidate-rows "${MIN_CANDIDATE_ROWS}"
  --min-apertures "${MIN_APERTURES}"
  --min-positive-peak-ranks "${MIN_POSITIVE_PEAK_RANKS}"
  --require-cadence-diagnostics
)
if [[ "${REQUIRE_SIGNAL_PEAKS}" != "1" ]]; then
  VERIFY_ARGS+=(--no-require-signal-peaks)
fi
log "verifying merged table"
"${PYTHON}" "${VERIFY_ARGS[@]}"

log "complete"
log "final_table=${FINAL_TABLE}"
log "verification_json=${VERIFY_JSON}"
