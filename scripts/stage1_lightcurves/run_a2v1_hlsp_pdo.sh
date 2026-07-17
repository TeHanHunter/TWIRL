#!/bin/bash
# Build A2v1 ADP/ADP015-only TWIRL-FS FITS products from completed TGLC HDF5s.
#
# Usage:
#   run_a2v1_hlsp_pdo.sh <sector> <orbit> [<orbit> ...]
#
# Example:
#   run_a2v1_hlsp_pdo.sh 56 119 120
#
# Output root defaults to:
#   /pdo/users/tehan/tglc-gpu-production-A2v1/hlsp_s<sector>_A2v1

set -uo pipefail

SECTOR="${1:?usage: $0 <sector> <orbit> [<orbit> ...]}"
shift || true
if [ "$#" -lt 1 ]; then
  echo "usage: $0 <sector> <orbit> [<orbit> ...]" >&2
  exit 1
fi

REPO=${TWIRL_REPO:-/pdo/users/tehan/TWIRL}
export TWIRL_REPO="$REPO"
A2V1_ROOT=${TWIRL_A2V1_ROOT:-/pdo/users/tehan/tglc-gpu-production-A2v1}
SCRIPT_DIR=${TWIRL_SCRIPT_DIR:-$REPO/scripts/stage1_lightcurves}
PY=${TWIRL_QLP_PY:-/sw/qlp-environment/.venv/bin/python}
SECTOR_TAG="s$(printf '%04d' "$SECTOR")"
OUT_ROOT=${TWIRL_A2V1_HLSP_ROOT:-$A2V1_ROOT/hlsp_${SECTOR_TAG}_A2v1}
WORKERS=${TWIRL_A2V1_HLSP_WORKERS:-8}
MAX_FAILURES=${TWIRL_A2V1_HLSP_MAX_FAILURES:-0}
MIN_SUCCESS_FRAC=${TWIRL_A2V1_HLSP_MIN_SUCCESS_FRAC:-1.0}

export HDF5_USE_FILE_LOCKING=${HDF5_USE_FILE_LOCKING:-FALSE}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export VECLIB_MAXIMUM_THREADS=${VECLIB_MAXIMUM_THREADS:-1}
export NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS:-1}
export LD_LIBRARY_PATH=/sw/python-versions/python-3.11.9/lib:/pdo/app/anaconda/anaconda2-4.4.0/lib:${LD_LIBRARY_PATH:-}
export PYTHONPATH=$REPO/src:/sw/qlp-environment/.venv/lib/python3.11/site-packages:${PYTHONPATH:-}

"$REPO/scripts/assert_clean_checkout.sh" "$REPO"

ORBIT_ROOTS=()
for orbit in "$@"; do
  ORBIT_ROOTS+=("$A2V1_ROOT/orbit-$orbit/ffi")
done

echo "[a2v1-hlsp] sector=$SECTOR out_root=$OUT_ROOT workers=$WORKERS"
for root in "${ORBIT_ROOTS[@]}"; do
  if [ ! -d "$root" ]; then
    echo "[a2v1-hlsp] missing orbit root: $root" >&2
    exit 1
  fi
  count=$(find "$root" -path "*/LC/*.h5" 2>/dev/null | wc -l | tr -d " ")
  echo "[a2v1-hlsp] $root h5=$count"
done

"$PY" "$SCRIPT_DIR/build_twirl_hlsp.py" \
  --sector "$SECTOR" \
  --orbit-roots "${ORBIT_ROOTS[@]}" \
  --out-root "$OUT_ROOT" \
  --workers "$WORKERS" \
  --a2v1-only \
  --max-failures "$MAX_FAILURES" \
  --min-success-frac "$MIN_SUCCESS_FRAC"
