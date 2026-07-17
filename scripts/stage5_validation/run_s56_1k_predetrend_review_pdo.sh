#!/usr/bin/env bash
# Build the scientifically useful S56 pre-detrend injection review sample on PDO.
#
# This injects BATMAN signals into raw TGLC RawFlux before TWIRL-FS detrending,
# reruns canonical + ADP detrending, then sends the detrended injected products
# through BLS and WD-tuned LEO-Vetter.
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}"
export LD_LIBRARY_PATH="/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${LEO_VETTER_ROOT:-/pdo/users/tehan/LEO-Vetter-twirl}:src"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

ORBIT_ROOTS=(
  "${ORBIT_119_ROOT:-/pdo/users/tehan/tglc-gpu-production/orbit-119/ffi}"
  "${ORBIT_120_ROOT:-/pdo/users/tehan/tglc-gpu-production/orbit-120/ffi}"
)
N_INJECTIONS="${N_INJECTIONS:-1000}"
RANDOM_STATE="${RANDOM_STATE:-5606}"
INJECTION_DIR="${INJECTION_DIR:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_1k_predetrend_batman_depthgrid_adp_compare}"
REVIEW_DIR="${REVIEW_DIR:-reports/stage5_validation/s56_1k_predetrend_batman_adp_compare_review_queue_pdo}"
REBUILD_INJECTIONS="${REBUILD_INJECTIONS:-1}"
REUSE_INJECTION_BLS="${REUSE_INJECTION_BLS:-0}"
QUEUE_APERTURES="${QUEUE_APERTURES:-DET_FLUX_ADP DET_FLUX}"
SHUFFLE_REVIEW_ROWS="${SHUFFLE_REVIEW_ROWS:-0}"
LEO_WORKERS="${LEO_WORKERS:-8}"
BLS_WORKERS="${BLS_WORKERS:-8}"
N_PERIODS="${N_PERIODS:-5000}"
MAX_LEO_REPORTS="${MAX_LEO_REPORTS:-0}"
read -r -a QUEUE_APERTURE_ARGS <<< "${QUEUE_APERTURES}"

if [[ "${REBUILD_INJECTIONS}" == "1" || ! -s "${INJECTION_DIR}/injected_lightcurves.h5" ]]; then
  echo "[pdo-predetrend-1k] building raw-flux pre-detrend injections"
  "${PYTHON_BIN}" scripts/stage3_injections/make_s56_predetrend_injection_set.py \
    --orbit-roots "${ORBIT_ROOTS[@]}" \
    --out-dir "${INJECTION_DIR}" \
    --n-injections "${N_INJECTIONS}" \
    --apertures DET_FLUX_ADP DET_FLUX DET_FLUX_ADP_SML DET_FLUX_ADP_LAG DET_FLUX_SML DET_FLUX_LAG \
    --sampling-mode period_depth_grid \
    --grid-period-range-d 0.12,13.0 \
    --grid-depth-range 0.003,0.995 \
    --grid-period-bins 10 \
    --grid-depth-bins 10 \
    --grid-depth-spacing linear \
    --random-state "${RANDOM_STATE}" \
    --min-in-transit 2 \
    --max-attempts-per-injection 80 \
    --progress-every 50 \
    --overwrite
else
  echo "[pdo-predetrend-1k] reusing injections: ${INJECTION_DIR}/injected_lightcurves.h5"
fi

queue_args=(
  --injection-h5 "${INJECTION_DIR}/injected_lightcurves.h5"
  --out-dir "${REVIEW_DIR}"
  --n-real 0
  --n-injections "${N_INJECTIONS}"
  --apertures "${QUEUE_APERTURE_ARGS[@]}"
  --workers "${BLS_WORKERS}"
  --n-periods "${N_PERIODS}"
  --leo-workers "${LEO_WORKERS}"
  --max-leo-reports "${MAX_LEO_REPORTS}"
  --overwrite
)
if [[ "${REUSE_INJECTION_BLS}" == "1" ]]; then
  queue_args+=(--reuse-injection-bls)
fi
if [[ "${SHUFFLE_REVIEW_ROWS}" == "1" ]]; then
  queue_args+=(--shuffle-review-rows)
fi

echo "[pdo-predetrend-1k] building BLS + LEO review queue"
"${PYTHON_BIN}" scripts/stage5_validation/build_s56_pretriage_review_queue.py "${queue_args[@]}"

echo "[pdo-predetrend-1k] complete"
echo "  injection_dir=${INJECTION_DIR}"
echo "  review_dir=${REVIEW_DIR}"
