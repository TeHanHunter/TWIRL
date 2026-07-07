#!/usr/bin/env bash
# Run the S56 ADP+ detrending/BLS preservation audit on PDO.
set -euo pipefail

cd "${TWIRL_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"

PYTHON_BIN="${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}"
export PYTHONPATH="${PYTHONPATH:-src}"
export LD_LIBRARY_PATH="/pdo/app/anaconda/anaconda2-4.4.0/lib:/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

OUT_DIR="${OUT_DIR:-reports/stage5_validation/s56_adpplus_bls_audit_pdo}"
INJECTION_CSV="${INJECTION_CSV:-reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/small_pair_200k/review_queue.csv}"
INJECTION_H5="${INJECTION_H5:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_20k_predetrend_batman_periodradius_bright_balanced/injected_lightcurves.h5}"
REAL_CANDIDATES="${REAL_CANDIDATES:-data_local/stage2/bls_first_pass_v2/sector_0056/vetted_per_tic_centroid.csv}"
HLSP_ROOT="${HLSP_ROOT:-/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl_fs_v2_compare}"
N_INJECTED="${N_INJECTED:-3000}"
N_REAL="${N_REAL:-2000}"
N_PERIODS="${N_PERIODS:-50000}"
N_PEAKS="${N_PEAKS:-20}"
WORKERS="${WORKERS:-8}"

echo "[adpplus-audit-pdo] out_dir=${OUT_DIR}"
echo "[adpplus-audit-pdo] injected=${N_INJECTED} real=${N_REAL} n_periods=${N_PERIODS} workers=${WORKERS}"

"${PYTHON_BIN}" scripts/stage5_validation/run_s56_adpplus_bls_audit.py \
  --injection-csv "${INJECTION_CSV}" \
  --injection-h5 "${INJECTION_H5}" \
  --real-candidates "${REAL_CANDIDATES}" \
  --hlsp-root "${HLSP_ROOT}" \
  --out-dir "${OUT_DIR}" \
  --n-injected "${N_INJECTED}" \
  --n-real "${N_REAL}" \
  --n-periods "${N_PERIODS}" \
  --n-peaks "${N_PEAKS}" \
  --workers "${WORKERS}"
