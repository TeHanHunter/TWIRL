#!/bin/bash
# Wait for S57 TWIRL-FS v2 to finish, then plot a random dim-target
# detrending sample. Designed to run in tmux on pdogpu6.

set -uo pipefail

REPO=/pdo/users/tehan/TWIRL
ROOT=/pdo/users/tehan/tglc-gpu-production
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SECTOR=57
TAG=57
O1=121
O2=122
PRODUCT_LABEL=twirl_fs_v2
HLSP_ROOT=$ROOT/hlsp_s0057_$PRODUCT_LABEL
OUT_ROOT=$ROOT/reports/twirl_fs_dim_qa
OUT_DIR=$OUT_ROOT/s57_random_dim_${PRODUCT_LABEL}_$(date +%Y%m%dT%H%M%S)
CSV_NAME=random_${PRODUCT_LABEL}_sample.csv
LOG=$ROOT/twirl_logs/s57_dim_twirl_fs_qa.log
PY=/pdo/users/tehan/twirl-gpu-venv/bin/python

export LD_LIBRARY_PATH=/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}

stamp() { date '+%Y-%m-%d %H:%M:%S %Z'; }
log() { echo "[$(stamp)] [s57-dim-qa] $*" | tee -a "$LOG"; }
abort() { log "ABORT: $1"; exit 1; }

mkdir -p "$OUT_ROOT" "$(dirname "$LOG")"

log "waiting for S57 TWIRL-FS done marker and nonempty $HLSP_ROOT"
while true; do
  if [ -f "$ROOT/markers/s${TAG}_done.flag" ] && [ -d "$HLSP_ROOT" ]; then
    if find "$HLSP_ROOT" -name '*.fits' -type f -print -quit 2>/dev/null | grep -q .; then
      break
    fi
  fi
  sleep 300
done

log "S57 product is ready; writing QA to $OUT_DIR"
mkdir -p "$OUT_DIR"

"$PY" "$SCRIPT_DIR/plot_twirl_fs_random_sample.py" \
  --root "$HLSP_ROOT" \
  --sector "$SECTOR" \
  --orbit-roots "$ROOT/orbit-$O1/ffi" "$ROOT/orbit-$O2/ffi" \
  --out-dir "$OUT_DIR" \
  --n 8 \
  --seed 5701 \
  --min-tmag 18.5 \
  --min-neg-q0 100 \
  --product-label "$PRODUCT_LABEL" \
  >> "$LOG" 2>&1

if [ ! -s "$OUT_DIR/$CSV_NAME" ] || [ "$(wc -l < "$OUT_DIR/$CSV_NAME")" -le 1 ]; then
  log "strict sample was empty; retrying with min-neg-q0=25"
  "$PY" "$SCRIPT_DIR/plot_twirl_fs_random_sample.py" \
    --root "$HLSP_ROOT" \
    --sector "$SECTOR" \
    --orbit-roots "$ROOT/orbit-$O1/ffi" "$ROOT/orbit-$O2/ffi" \
    --out-dir "$OUT_DIR" \
    --n 8 \
    --seed 5701 \
    --min-tmag 18.5 \
    --min-neg-q0 25 \
    --product-label "$PRODUCT_LABEL" \
    >> "$LOG" 2>&1 || abort "fallback QA plot failed"
fi

if [ ! -s "$OUT_DIR/$CSV_NAME" ] || [ "$(wc -l < "$OUT_DIR/$CSV_NAME")" -le 1 ]; then
  abort "QA sample is still empty"
fi

printf "%s\n" "$OUT_DIR" > "$OUT_ROOT/latest_s57_dim_${PRODUCT_LABEL}_qa_path.txt"
date -Iseconds > "$ROOT/markers/s57_dim_${PRODUCT_LABEL}_qa_done.flag"
log "DONE: $OUT_DIR"
