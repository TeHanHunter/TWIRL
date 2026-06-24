#!/usr/bin/env bash
# Build a dense S56 pre-detrend BATMAN injection sample and run the canonical
# small-aperture BLS recovery pass. This is for visualizing the BLS sensitivity
# boundary in period/depth space; it intentionally skips LEO and human-queue
# rendering.
set -euo pipefail

cd /pdo/users/tehan/TWIRL

export LD_LIBRARY_PATH="/pdo/app/anaconda/anaconda2-4.4.0/lib:/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH=src
export MPLCONFIGDIR=/pdo/users/tehan/.cache/matplotlib
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

PYTHON="${PYTHON:-/sw/qlp-environment/.venv/bin/python}"
ORBIT_ROOTS=(
  "${ORBIT_119_ROOT:-/pdo/users/tehan/tglc-gpu-production/orbit-119/ffi}"
  "${ORBIT_120_ROOT:-/pdo/users/tehan/tglc-gpu-production/orbit-120/ffi}"
)
N_INJECTIONS="${N_INJECTIONS:-10000}"
RANDOM_STATE="${RANDOM_STATE:-5623}"
GRID_PERIOD_BINS="${GRID_PERIOD_BINS:-50}"
GRID_DEPTH_BINS="${GRID_DEPTH_BINS:-50}"
GRID_DEPTH_SPACING="${GRID_DEPTH_SPACING:-linear}"
BASELINE_SOURCE="${BASELINE_SOURCE:-raw_median}"
INJECTION_DIR="${INJECTION_DIR:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_10k_predetrend_batman_depthgrid_dense_bls_map}"
OUT_DIR="${OUT_DIR:-reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo}"
CHUNK_DIR="${OUT_DIR}/chunk_ids"
RUN_CHUNK_DIR="${OUT_DIR}/chunks"
CHUNK_SIZE="${TWIRL_CHUNK_SIZE:-100}"
WORKERS="${TWIRL_CHUNK_WORKERS:-16}"
SWEEP_NAME="${SWEEP_NAME:-small_pair_200k}"
SWEEP_SPEC="${SWEEP_NAME}:DET_FLUX_ADP_SML+DET_FLUX_SML:200000"
REBUILD_INJECTIONS="${REBUILD_INJECTIONS:-0}"

log() { echo "[$(date -Is)] [dense-bls-map] $*"; }

mkdir -p "${CHUNK_DIR}" "${RUN_CHUNK_DIR}" logs

log "python=${PYTHON}"
log "n_injections=${N_INJECTIONS} grid=${GRID_PERIOD_BINS}x${GRID_DEPTH_BINS} depth_spacing=${GRID_DEPTH_SPACING}"
log "injection_dir=${INJECTION_DIR}"
log "out_dir=${OUT_DIR}"

if [[ "${REBUILD_INJECTIONS}" == "1" || ! -s "${INJECTION_DIR}/injected_lightcurves.h5" ]]; then
  log "building dense raw-flux pre-detrend BATMAN injections"
  "${PYTHON}" scripts/stage3_injections/make_s56_predetrend_injection_set.py \
    --orbit-roots "${ORBIT_ROOTS[@]}" \
    --out-dir "${INJECTION_DIR}" \
    --n-injections "${N_INJECTIONS}" \
    --apertures DET_FLUX_ADP_SML DET_FLUX_SML \
    --sampling-mode period_depth_grid \
    --grid-period-range-d 0.12,13.0 \
    --grid-depth-range 0.003,0.995 \
    --grid-period-bins "${GRID_PERIOD_BINS}" \
    --grid-depth-bins "${GRID_DEPTH_BINS}" \
    --grid-depth-spacing "${GRID_DEPTH_SPACING}" \
    --baseline-source "${BASELINE_SOURCE}" \
    --random-state "${RANDOM_STATE}" \
    --min-in-transit 2 \
    --max-attempts-per-injection 100 \
    --progress-every 100 \
    --overwrite
else
  log "reusing injections: ${INJECTION_DIR}/injected_lightcurves.h5"
fi

log "writing restartable BLS chunk id files"
"${PYTHON}" - "${INJECTION_DIR}/injected_lightcurves.h5" "${CHUNK_DIR}" "${CHUNK_SIZE}" <<'PY'
from pathlib import Path
import sys

import h5py

injection_h5 = Path(sys.argv[1])
chunk_dir = Path(sys.argv[2])
chunk_size = int(sys.argv[3])

with h5py.File(injection_h5, "r") as h5:
    keys = sorted(h5["injections"].keys())

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

for id_file in "${CHUNK_DIR}"/chunk_*.txt; do
  chunk_name="$(basename "${id_file}" .txt)"
  chunk_out="${RUN_CHUNK_DIR}/${chunk_name}"
  summary_path="${chunk_out}/${SWEEP_NAME}/recovery_mode_summary/summary.json"
  if [[ -f "${summary_path}" ]]; then
    log "skip ${chunk_name}: ${summary_path} exists"
    continue
  fi
  log "run ${chunk_name} with ${WORKERS} workers"
  "${PYTHON}" scripts/stage5_validation/run_injection_recovery_mode_sweep.py \
    --injection-h5 "${INJECTION_DIR}/injected_lightcurves.h5" \
    --out-dir "${chunk_out}" \
    --n-injections 0 \
    --workers "${WORKERS}" \
    --injection-id-file "${id_file}" \
    --sweep "${SWEEP_SPEC}"
done

log "merging BLS chunks"
"${PYTHON}" scripts/stage5_validation/merge_recovery_sweep_chunks.py \
  --chunk-root "${RUN_CHUNK_DIR}" \
  --out-dir "${OUT_DIR}" \
  --sweep-name "${SWEEP_NAME}" \
  --apertures DET_FLUX_ADP_SML DET_FLUX_SML \
  --n-periods 200000 \
  --summary-aperture DET_FLUX_ADP_SML \
  --expect-n "${N_INJECTIONS}"

log "complete"
log "recovery_csv=${OUT_DIR}/${SWEEP_NAME}/injection_bls_recoveries.csv"
