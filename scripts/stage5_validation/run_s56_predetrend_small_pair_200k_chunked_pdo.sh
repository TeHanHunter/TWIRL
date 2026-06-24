#!/usr/bin/env bash
set -euo pipefail

cd /pdo/users/tehan/TWIRL

export LD_LIBRARY_PATH="/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH=src
export MPLCONFIGDIR=/pdo/users/tehan/.cache/matplotlib
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

PYTHON=/sw/qlp-environment/.venv/bin/python
INJECTION_H5=data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_1k_predetrend_batman_depthgrid_adp_compare/injected_lightcurves.h5
SURVIVAL_CSV=reports/stage5_validation/s56_1k_predetrend_batman_adp_only_review_queue_pdo/predet_signal_survival_snr_diagnostics.csv
OUT_DIR=reports/stage5_validation/s56_1k_predetrend_small_pair_200k_chunked_pdo
CHUNK_DIR="${OUT_DIR}/chunk_ids"
RUN_CHUNK_DIR="${OUT_DIR}/chunks"
CHUNK_SIZE="${TWIRL_CHUNK_SIZE:-100}"
WORKERS="${TWIRL_CHUNK_WORKERS:-8}"
SWEEP_NAME=small_pair_200k
SWEEP_SPEC="${SWEEP_NAME}:DET_FLUX_ADP_SML+DET_FLUX_SML:200000"

mkdir -p "${CHUNK_DIR}" "${RUN_CHUNK_DIR}" logs

"${PYTHON}" - "${INJECTION_H5}" "${CHUNK_DIR}" "${CHUNK_SIZE}" <<'PY'
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
    echo "[chunked-200k] skip ${chunk_name}: ${summary_path} exists"
    continue
  fi
  echo "[chunked-200k] run ${chunk_name} with ${WORKERS} workers"
  "${PYTHON}" scripts/stage5_validation/run_injection_recovery_mode_sweep.py \
    --injection-h5 "${INJECTION_H5}" \
    --survival-csv "${SURVIVAL_CSV}" \
    --out-dir "${chunk_out}" \
    --n-injections 0 \
    --workers "${WORKERS}" \
    --injection-id-file "${id_file}" \
    --sweep "${SWEEP_SPEC}"
done

"${PYTHON}" scripts/stage5_validation/merge_recovery_sweep_chunks.py \
  --chunk-root "${RUN_CHUNK_DIR}" \
  --out-dir "${OUT_DIR}" \
  --sweep-name "${SWEEP_NAME}" \
  --apertures DET_FLUX_ADP_SML DET_FLUX_SML \
  --n-periods 200000 \
  --survival-csv "${SURVIVAL_CSV}" \
  --summary-aperture DET_FLUX_ADP_SML \
  --expect-n 1000
