#!/bin/bash
# PDO worker: run the per-sector BLS first pass on one sector and gate it
# through the WD 1856 audit. Mirrors the finalize_worker.sh idiom so a daemon
# loop / tmux session can call this once per ready sector.
#
# Usage:
#   stage2_bls_worker.sh <sector> [hlsp_root]
#
# Defaults:
#   hlsp_root  = /pdo/users/tehan/tglc-gpu-production/hlsp_s{NN:04d}
#   workers    = $TWIRL_BLS_WORKERS or 32
#   out_dir    = $TWIRL_BLS_OUT_DIR or $REPO/data_local/stage2/bls_first_pass
#
# BLAS thread caps are mandatory on PDO multi-core nodes (CLAUDE.md): each
# pool worker would otherwise spawn 128 BLAS threads and the run goes 10-100x
# slower. Keep these exports.

set -euo pipefail

if [ "${1-}" = "" ]; then
  echo "usage: stage2_bls_worker.sh <sector> [hlsp_root]" >&2
  exit 64
fi
SECTOR="$1"
HLSP_OVERRIDE="${2-}"

REPO="${TWIRL_REPO:-/pdo/users/tehan/TWIRL}"
WORKERS="${TWIRL_BLS_WORKERS:-32}"
OUT_DIR="${TWIRL_BLS_OUT_DIR:-$REPO/data_local/stage2/bls_first_pass}"
PYTHON="${TWIRL_PYTHON:-python}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONPATH="${PYTHONPATH:-}:$REPO/src"

cd "$REPO"

log() { echo "[$(date '+%F %H:%M:%S')] [stage2-bls] $*"; }

log "sector=$SECTOR workers=$WORKERS out_dir=$OUT_DIR"

ARGS=(--sector "$SECTOR" --workers "$WORKERS" --out-dir "$OUT_DIR")
if [ -n "$HLSP_OVERRIDE" ]; then
  ARGS+=(--hlsp-root "$HLSP_OVERRIDE")
fi

"$PYTHON" -m twirl.search.sector_run "${ARGS[@]}"

# Audit gate: only relevant if WD 1856 is observed in this sector.
# WD 1856 is in S56 and S59 in the 200 s era; tess-coverage drives the
# auto-decision in audit_wd1856_recovery.py (which simply checks for a row).
if "$PYTHON" "$REPO/scripts/stage2_search/audit_wd1856_recovery.py" \
      --sector "$SECTOR" --out-dir "$OUT_DIR"; then
  log "sector=$SECTOR audit PASS"
else
  rc=$?
  if [ "$rc" -eq 2 ]; then
    log "sector=$SECTOR audit skipped (WD 1856 not in this sector)"
  else
    log "sector=$SECTOR audit FAILED (rc=$rc) — block scale-out and investigate"
    exit "$rc"
  fi
fi

log "sector=$SECTOR DONE"
