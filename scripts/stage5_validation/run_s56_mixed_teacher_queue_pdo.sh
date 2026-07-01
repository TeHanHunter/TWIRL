#!/usr/bin/env bash
# Build the S56 10k mixed teacher pool and the blinded 1k first-pass queue.
#
# This deliberately does not use the injected-truth peak ranker.
set -euo pipefail

REPO_ROOT="${TWIRL_REPO_ROOT:-/pdo/users/tehan/TWIRL}"
cd "${REPO_ROOT}"

PYTHON_BIN="${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}"
export LD_LIBRARY_PATH="/pdo/app/anaconda/anaconda2-4.4.0/lib:/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="/pdo/users/tehan/LEO-Vetter-twirl:src:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

OUT_DIR="${OUT_DIR:-reports/stage5_validation/s56_mixed_teacher_queue_pdo}"
REAL_CANDIDATES="${REAL_CANDIDATES:-data_local/stage2/bls_first_pass_v2/sector_0056/vetted_per_tic_centroid.csv}"
INJECTED_CANDIDATES="${INJECTED_CANDIDATES:-reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/small_pair_200k/review_queue.csv}"
INJECTION_H5="${INJECTION_H5:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_20k_predetrend_batman_periodradius_bright_balanced/injected_lightcurves.h5}"
HLSP_ROOT="${HLSP_ROOT:-/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl_fs_v2_compare}"
STAR_CATALOG="${STAR_CATALOG:-data_local/catalogs/twirl_master_catalog/twirl_wd_master_catalog_v0_ticmatched.fits}"
N_REAL_POOL="${N_REAL_POOL:-9000}"
N_INJECTED_POOL="${N_INJECTED_POOL:-1000}"
N_REVIEW="${N_REVIEW:-1000}"
REVIEW_REAL="${REVIEW_REAL:-900}"
REVIEW_INJECTED="${REVIEW_INJECTED:-100}"
RANDOM_STATE="${RANDOM_STATE:-5608}"
APERTURE="${APERTURE:-DET_FLUX_ADP_SML}"
MAX_LEO_REPORTS="${MAX_LEO_REPORTS:-1000}"
LEO_WORKERS="${LEO_WORKERS:-8}"
LEO_TIMEOUT_S="${LEO_TIMEOUT_S:-300}"

args=(
  --real-candidates "${REAL_CANDIDATES}"
  --injected-candidates "${INJECTED_CANDIDATES}"
  --injection-h5 "${INJECTION_H5}"
  --hlsp-root "${HLSP_ROOT}"
  --star-catalog "${STAR_CATALOG}"
  --out-dir "${OUT_DIR}"
  --n-real-pool "${N_REAL_POOL}"
  --n-injected-pool "${N_INJECTED_POOL}"
  --n-review "${N_REVIEW}"
  --review-real "${REVIEW_REAL}"
  --review-injected "${REVIEW_INJECTED}"
  --random-state "${RANDOM_STATE}"
  --aperture "${APERTURE}"
  --max-leo-reports "${MAX_LEO_REPORTS}"
  --leo-workers "${LEO_WORKERS}"
  --leo-timeout-s "${LEO_TIMEOUT_S}"
)

if [[ "${SKIP_LEO:-0}" == "1" ]]; then
  args+=(--skip-leo)
fi
if [[ "${OVERWRITE:-0}" == "1" ]]; then
  args+=(--overwrite)
fi

echo "[mixed-teacher-pdo] out_dir=${OUT_DIR}"
echo "[mixed-teacher-pdo] real=${REAL_CANDIDATES}"
echo "[mixed-teacher-pdo] injected=${INJECTED_CANDIDATES}"
echo "[mixed-teacher-pdo] leo_reports=${MAX_LEO_REPORTS} workers=${LEO_WORKERS}"
"${PYTHON_BIN}" scripts/stage5_validation/build_s56_mixed_teacher_queue.py "${args[@]}"

echo "[mixed-teacher-pdo] complete"
echo "[mixed-teacher-pdo] queue=${OUT_DIR}/review_queue_1k.csv"
echo "[mixed-teacher-pdo] verifier=${OUT_DIR}/verification.json"
