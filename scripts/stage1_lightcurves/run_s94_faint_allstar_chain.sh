#!/bin/bash
# Low-priority S94 faint all-star HDF5 + precision-plot chain.
#
# Outputs are confined to /pdo/users/tehan/tglc-s94-faint-allstar.
# Shared /pdo/qlp-data products are read-only inputs.

set -uo pipefail

REPO=${TWIRL_REPO:-/pdo/users/tehan/TWIRL}
ROOT=${TWIRL_S94_ROOT:-/pdo/users/tehan/tglc-s94-faint-allstar}
INPUT_ROOT=${TWIRL_S94_INPUT_ROOT:-/pdo/qlp-data}
SECTOR=94
ORBIT_1=195
ORBIT_2=196
CCD_ARG=${TWIRL_S94_CCD:-all}
WORKERS=${TWIRL_S94_WORKERS:-4}
CATALOG_NPROCS=${TWIRL_S94_CATALOG_NPROCS:-16}
TMAG_MAX=${TWIRL_S94_TMAG_MAX:-18.5}
LIMIT_SOURCES=${TWIRL_S94_LIMIT_SOURCES:-}

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
export PYTHONPATH=/pdo/users/tehan/tess-gaia-light-curve-twirl:$REPO/src:/sw/qlp-environment/.venv/lib/python3.11/site-packages:${PYTHONPATH:-}
export LD_LIBRARY_PATH=/sw/python-versions/python-3.11.9/lib:/pdo/app/anaconda/anaconda2-4.4.0/lib:${LD_LIBRARY_PATH:-}

PY=/sw/qlp-environment/.venv/bin/python
if [ -d "$REPO/scripts/stage1_lightcurves" ]; then
  SCRIPT_DIR="$REPO/scripts/stage1_lightcurves"
else
  SCRIPT_DIR="$REPO/scripts/stage1_lcs"
fi
LOG="$ROOT/logs"
PRECISION_DIR="$ROOT/qc_precision"

mkdir -p "$ROOT" "$LOG" "$PRECISION_DIR"
cd "$REPO" || exit 1

stamp() { date '+%Y-%m-%d %H:%M:%S %Z'; }
log() { echo "[$(stamp)] $*" | tee -a "$LOG/s94_faint_allstar_chain.log"; }
abort() { log "ABORT: $1"; exit 1; }

LC_LIMIT_ARGS=()
if [ -n "$LIMIT_SOURCES" ]; then
  LC_LIMIT_ARGS=(--limit-sources "$LIMIT_SOURCES")
fi

log "step A: faint TIC-backed TGLC HDF5 light curves, ccd=$CCD_ARG, Tmag<=$TMAG_MAX"
nice -n 19 "$PY" "$SCRIPT_DIR/run_faint_allstar_lightcurves.py" \
  --sector "$SECTOR" \
  --orbits "$ORBIT_1" "$ORBIT_2" \
  --input-root "$INPUT_ROOT" \
  --out-root "$ROOT" \
  --ccd "$CCD_ARG" \
  --stages catalogs,lightcurves \
  --tmag-max "$TMAG_MAX" \
  --catalog-nprocs "$CATALOG_NPROCS" \
  --workers "$WORKERS" \
  "${LC_LIMIT_ARGS[@]}" \
  >> "$LOG/faint_allstar_lcs.log" 2>&1 || abort "faint all-star light curves failed"

log "step B: TGLC-style all-star MAD precision plot from HDF5 RawFlux"
nice -n 19 "$PY" "$SCRIPT_DIR/plot_hlsp_precision.py" \
  --hlsp-root "$ROOT" \
  --sector "$SECTOR" \
  --output "$PRECISION_DIR/s94_mad_all_tmag18p5.png" \
  --input-kind tglc-h5 \
  --all \
  --style tglc \
  --label "TGLC Aperture" \
  --h5-aperture Primary \
  --tmag-max "$TMAG_MAX" \
  --noisemodel "$REPO/data_local/refs/noisemodel.dat" \
  >> "$LOG/s94_mad_all_tmag18p5.log" 2>&1 || abort "precision plot failed"

log "done"
