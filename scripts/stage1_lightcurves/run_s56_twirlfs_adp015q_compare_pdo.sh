#!/usr/bin/env bash
# Build an S56 TWIRL-FS compare tree with the ADP015Q candidate columns.
#
# The canonical DET_FLUX* columns remain twirl-fs-v2.  The extra compare
# columns are:
#   DET_FLUX_ADP015, DET_FLUX_ADP015_ERR, DET_FLUX_ADP015_SML,
#   DET_FLUX_ADP015_LAG
set -euo pipefail

cd "${TWIRL_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"

PYTHON_BIN="${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}"
ROOT="${ROOT:-/pdo/users/tehan/tglc-gpu-production}"
OUT_ROOT="${OUT_ROOT:-/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl_fs_v2_adp015q_compare}"
if [[ -n "${ORBIT_ROOTS:-}" ]]; then
  read -r -a ORBIT_ROOTS_ARRAY <<<"${ORBIT_ROOTS}"
else
  ORBIT_ROOTS_ARRAY=("${ROOT}/orbit-119/ffi" "${ROOT}/orbit-120/ffi")
fi
WORKERS="${WORKERS:-32}"
SECTOR="${SECTOR:-56}"

export PYTHONPATH="${PYTHONPATH:-src}"
export LD_LIBRARY_PATH="/pdo/app/anaconda/anaconda2-4.4.0/lib:/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export HDF5_USE_FILE_LOCKING="${HDF5_USE_FILE_LOCKING:-FALSE}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

args=(
  --sector "${SECTOR}"
  --orbit-roots "${ORBIT_ROOTS_ARRAY[@]}"
  --out-root "${OUT_ROOT}"
  --workers "${WORKERS}"
  --max-failures "${MAX_FAILURES:-250}"
  --min-success-frac "${MIN_SUCCESS_FRAC:-0.99}"
  --include-adaptive-columns
  --adaptive-bkspace 0.15
  --adaptive-gap-split 0.2
  --adaptive-method-version twirl-fs-v2-adp015q
  --adaptive-column-tag ADP015
)

if [[ -n "${LIMIT:-}" ]]; then
  args+=(--limit "${LIMIT}")
fi
if [[ -n "${TIC:-}" ]]; then
  args+=(--tic "${TIC}")
fi

echo "[s56-adp015q-compare] python=${PYTHON_BIN}"
echo "[s56-adp015q-compare] out_root=${OUT_ROOT}"
echo "[s56-adp015q-compare] orbit_roots=${ORBIT_ROOTS_ARRAY[*]}"
echo "[s56-adp015q-compare] workers=${WORKERS}"
"${PYTHON_BIN}" scripts/stage1_lightcurves/build_twirl_hlsp.py "${args[@]}"
