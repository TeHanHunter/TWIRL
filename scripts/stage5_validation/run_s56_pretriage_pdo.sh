#!/usr/bin/env bash
set -euo pipefail

cd /pdo/users/tehan/TWIRL

export LD_LIBRARY_PATH="/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="/pdo/users/tehan/LEO-Vetter-twirl:src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

PYTHON_BIN="${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}"
HLSP_ROOT="${HLSP_ROOT:-/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl_fs_v2_compare}"
EXPORT_H5="${EXPORT_H5:-data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_lc_export_pdo.h5}"
INJECTION_DIR="${INJECTION_DIR:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_batman_depthgrid_r1p5j_900}"
REVIEW_DIR="${REVIEW_DIR:-reports/stage5_validation/s56_pretriage_review_queue_pdo}"
CANDIDATES="${CANDIDATES:-data_local/stage2/bls_first_pass_v2/sector_0056/vetted_per_tic_centroid.csv}"
WORKERS="${WORKERS:-8}"
N_PERIODS="${N_PERIODS:-200000}"
LEO_TIMEOUT_S="${LEO_TIMEOUT_S:-300}"

echo "[pdo-pretriage] start $(date -Is)"
echo "[pdo-pretriage] python=${PYTHON_BIN}"
echo "[pdo-pretriage] hlsp_root=${HLSP_ROOT}"

"${PYTHON_BIN}" scripts/stage3_injections/export_s56_lc_training_set.py \
  --hlsp-root "${HLSP_ROOT}" \
  --out-h5 "${EXPORT_H5}" \
  --progress-every 1000 \
  --overwrite

"${PYTHON_BIN}" scripts/stage3_injections/make_s56_lc_injection_training_set.py \
  --export-h5 "${EXPORT_H5}" \
  --out-dir "${INJECTION_DIR}" \
  --apertures DET_FLUX_SML DET_FLUX DET_FLUX_LAG \
  --n-injections 900 \
  --sampling-mode period_depth_grid \
  --grid-period-bins 10 \
  --grid-depth-bins 10 \
  --store-original \
  --overwrite

"${PYTHON_BIN}" scripts/stage5_validation/build_s56_pretriage_review_queue.py \
  --real-candidates "${CANDIDATES}" \
  --hlsp-root "${HLSP_ROOT}" \
  --injection-h5 "${INJECTION_DIR}/injected_lightcurves.h5" \
  --out-dir "${REVIEW_DIR}" \
  --n-real 100 \
  --n-injections 900 \
  --n-periods "${N_PERIODS}" \
  --workers "${WORKERS}" \
  --leo-timeout-s "${LEO_TIMEOUT_S}" \
  --overwrite

echo "[pdo-pretriage] complete $(date -Is)"
echo "[pdo-pretriage] review queue: ${REVIEW_DIR}/review_queue.csv"
