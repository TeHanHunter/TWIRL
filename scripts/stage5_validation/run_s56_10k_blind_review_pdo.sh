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
INJECTION_DIR="${INJECTION_DIR:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_10k_batman_depthgrid_r1p5j_1000}"
REVIEW_DIR="${REVIEW_DIR:-reports/stage5_validation/s56_10k_blind_review_queue_r1p5j_pdo}"
CANDIDATES="${CANDIDATES:-data_local/stage2/bls_first_pass_v2/sector_0056/vetted_per_tic_centroid.csv}"
N_REAL="${N_REAL:-9000}"
N_INJECTIONS="${N_INJECTIONS:-1000}"
WORKERS="${WORKERS:-16}"
N_PERIODS="${N_PERIODS:-200000}"
LEO_WORKERS="${LEO_WORKERS:-8}"
LEO_TIMEOUT_S="${LEO_TIMEOUT_S:-300}"
RANDOM_STATE="${RANDOM_STATE:-5610}"
REBUILD_EXPORT="${REBUILD_EXPORT:-0}"
REBUILD_INJECTIONS="${REBUILD_INJECTIONS:-0}"
REUSE_INJECTION_BLS="${REUSE_INJECTION_BLS:-0}"
SKIP_LEO="${SKIP_LEO:-0}"
MAX_LEO_REPORTS="${MAX_LEO_REPORTS:-0}"

echo "[pdo-10k-review] start $(date -Is)"
echo "[pdo-10k-review] python=${PYTHON_BIN}"
echo "[pdo-10k-review] hlsp_root=${HLSP_ROOT}"
echo "[pdo-10k-review] review_dir=${REVIEW_DIR}"

if [[ "${REBUILD_EXPORT}" == "1" || ! -s "${EXPORT_H5}" ]]; then
  "${PYTHON_BIN}" scripts/stage3_injections/export_s56_lc_training_set.py \
    --hlsp-root "${HLSP_ROOT}" \
    --out-h5 "${EXPORT_H5}" \
    --progress-every 1000 \
    --overwrite
else
  echo "[pdo-10k-review] reusing compact export: ${EXPORT_H5}"
fi

if [[ "${REBUILD_INJECTIONS}" == "1" || ! -s "${INJECTION_DIR}/injected_lightcurves.h5" ]]; then
  "${PYTHON_BIN}" scripts/stage3_injections/make_s56_lc_injection_training_set.py \
    --export-h5 "${EXPORT_H5}" \
    --out-dir "${INJECTION_DIR}" \
    --apertures DET_FLUX_SML DET_FLUX DET_FLUX_LAG \
    --n-injections "${N_INJECTIONS}" \
    --sampling-mode period_depth_grid \
    --grid-period-bins 10 \
    --grid-depth-bins 10 \
    --families wd_small_body wd_earth_size wd_giant_or_bd wd_long_period wd_roche_edge wd1856_like short_deep roche_boundary \
    --random-state "${RANDOM_STATE}" \
    --min-in-transit 2 \
    --max-attempts-per-injection 50 \
    --compression lzf \
    --store-original \
    --overwrite
else
  echo "[pdo-10k-review] reusing injections: ${INJECTION_DIR}/injected_lightcurves.h5"
fi

queue_args=(
  --real-candidates "${CANDIDATES}"
  --hlsp-root "${HLSP_ROOT}"
  --injection-h5 "${INJECTION_DIR}/injected_lightcurves.h5"
  --out-dir "${REVIEW_DIR}"
  --n-real "${N_REAL}"
  --n-injections "${N_INJECTIONS}"
  --real-selection stratified_blind
  --blind-review-metadata
  --shuffle-review-rows
  --apertures DET_FLUX_SML DET_FLUX DET_FLUX_LAG
  --n-periods "${N_PERIODS}"
  --workers "${WORKERS}"
  --leo-workers "${LEO_WORKERS}"
  --leo-timeout-s "${LEO_TIMEOUT_S}"
  --max-leo-reports "${MAX_LEO_REPORTS}"
  --random-state "${RANDOM_STATE}"
  --overwrite
)
if [[ "${REUSE_INJECTION_BLS}" == "1" ]]; then
  queue_args+=(--reuse-injection-bls)
fi
if [[ "${SKIP_LEO}" == "1" ]]; then
  queue_args+=(--skip-leo)
fi

"${PYTHON_BIN}" scripts/stage5_validation/build_s56_pretriage_review_queue.py "${queue_args[@]}"

echo "[pdo-10k-review] complete $(date -Is)"
echo "[pdo-10k-review] review queue: ${REVIEW_DIR}/review_queue.csv"
