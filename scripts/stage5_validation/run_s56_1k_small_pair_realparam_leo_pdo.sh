#!/usr/bin/env bash
# Re-render the active S56 1k small-pair LEO queue with catalog-backed WD hosts.
#
# This writes a sibling queue directory so the in-progress human labels on the
# canonical-host queue are not touched. The BLS/injection rows are reused; only
# the LEO star parameters and reports are regenerated.
set -euo pipefail

cd /pdo/users/tehan/TWIRL

export LD_LIBRARY_PATH="/pdo/app/anaconda/anaconda2-4.4.0/lib:/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="/pdo/users/tehan/LEO-Vetter-twirl:src"
export MPLCONFIGDIR=/pdo/users/tehan/.cache/matplotlib
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

PYTHON="${PYTHON:-/sw/qlp-environment/.venv/bin/python}"
INJECTION_H5="${INJECTION_H5:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_1k_predetrend_batman_depthgrid_adp_compare/injected_lightcurves.h5}"
SOURCE_BLS="${SOURCE_BLS:-reports/stage5_validation/s56_1k_predetrend_small_pair_200k_chunked_pdo/small_pair_200k/injection_bls_recoveries.csv}"
OUT_DIR="${OUT_DIR:-reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_realparam_pdo}"
HLSP_ROOT="${HLSP_ROOT:-/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl_fs_v2_compare}"
STAR_CATALOG="${STAR_CATALOG:-data_local/catalogs/twirl_master_catalog/twirl_wd_master_catalog_v0_ticmatched.fits}"
STAR_PRIORITY="${STAR_PRIORITY:-H,He,mixed}"

mkdir -p "${OUT_DIR}" logs
cp "${SOURCE_BLS}" "${OUT_DIR}/injection_bls_recoveries.csv"

"${PYTHON}" scripts/stage5_validation/build_s56_pretriage_review_queue.py \
  --injection-h5 "${INJECTION_H5}" \
  --hlsp-root "${HLSP_ROOT}" \
  --star-catalog "${STAR_CATALOG}" \
  --star-atmosphere-priority "${STAR_PRIORITY}" \
  --out-dir "${OUT_DIR}" \
  --n-real 0 \
  --n-injections 1000 \
  --apertures DET_FLUX_ADP_SML DET_FLUX_SML \
  --workers 1 \
  --n-periods 200000 \
  --reuse-injection-bls \
  --shuffle-review-rows \
  --random-state 5606 \
  --blind-review-metadata \
  --leo-workers "${TWIRL_LEO_WORKERS:-8}" \
  --leo-timeout-s 300 \
  --max-leo-reports 0 \
  --overwrite

"${PYTHON}" scripts/stage5_validation/verify_s56_pretriage_queue.py \
  --queue "${OUT_DIR}/review_queue.csv" \
  --reports-dir "${OUT_DIR}/vet_reports" \
  --summary-json "${OUT_DIR}/summary.json" \
  --min-rows 1000 \
  --expect-real 0 \
  --expect-injected 1000 \
  --min-reports 1000 \
  --apertures DET_FLUX_ADP_SML DET_FLUX_SML \
  --skip-wd1856-check \
  --out-json "${OUT_DIR}/verification.json"
