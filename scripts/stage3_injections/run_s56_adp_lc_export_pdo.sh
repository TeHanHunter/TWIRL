#!/usr/bin/env bash
# Build a compact S56 adaptive light-curve export for ORCD-side vetting/search tests.
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

APERTURE_TAG="${APERTURE_TAG:-ADP015}"
HLSP_ROOT="${HLSP_ROOT:-/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl_fs_v2_adp015q_compare}"
OUT_H5="${OUT_H5:-data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp015_lc_export_pdo.h5}"
SECTOR="${SECTOR:-56}"
LIMIT="${LIMIT:-0}"
OVERWRITE="${OVERWRITE:-0}"

args=(
  --hlsp-root "${HLSP_ROOT}"
  --out-h5 "${OUT_H5}"
  --sector "${SECTOR}"
  --columns "DET_FLUX_${APERTURE_TAG}_SML" "DET_FLUX_${APERTURE_TAG}"
  --compression lzf
  --progress-every 1000
)
if [[ "${LIMIT}" != "0" ]]; then
  args+=(--limit "${LIMIT}")
fi
if [[ "${OVERWRITE}" == "1" ]]; then
  args+=(--overwrite)
fi

echo "[adp-export-pdo] hlsp_root=${HLSP_ROOT}"
echo "[adp-export-pdo] out_h5=${OUT_H5}"
echo "[adp-export-pdo] aperture_tag=${APERTURE_TAG}"
"${PYTHON_BIN}" scripts/stage3_injections/export_s56_lc_training_set.py "${args[@]}"
