#!/bin/bash
# Refresh previously-produced sectors with the patched TGLC lightcurves stage,
# then build TWIRL-FS v2 HLSPs. Designed for pdogpu6.
#
# Existing LC directories are moved aside once per sector under
# LC_legacy_qlp_<stamp>/ so the patched TGLC fork can regenerate HDF5 files
# with RawFlux/RawFluxError. ePSFs and cutouts are reused.
#
# Usage:
#   relight_twirl_fs_sectors.sh 57 58 59

set -uo pipefail

REPO=/pdo/users/tehan/TWIRL
ROOT=/pdo/users/tehan/tglc-gpu-production
LOG=$ROOT/twirl_logs
MARKERS=$ROOT/markers
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY=/sw/qlp-environment/.venv/bin/python
LD_PREFIX=/pdo/app/anaconda/anaconda2-4.4.0/lib:/sw/python-versions/python-3.11.9/lib
STAMP=${TWIRL_RELIGHT_STAMP:-$(date +%Y%m%dT%H%M%S)}

export LD_LIBRARY_PATH=$LD_PREFIX:${LD_LIBRARY_PATH:-}
export PYTHONPATH=/pdo/users/tehan/tess-gaia-light-curve-twirl:/sw/qlp-environment/.venv/lib/python3.11/site-packages

mkdir -p "$LOG" "$MARKERS"

stamp() { date '+%Y-%m-%d %H:%M:%S %Z'; }
log() { echo "[$(stamp)] [twirlfs-relight] $*"; }
abort() { log "ABORT: $1"; exit 1; }

if [ "$#" -eq 0 ]; then
  set -- 57 58 59 60 61 62 63
fi

backup_lightcurves_once() {
  local sector=$1
  local o1=$((2 * sector + 7))
  local o2=$((o1 + 1))
  local tag=$(printf '%02d' "$sector")
  local backup_flag="$MARKERS/s${tag}_twirlfs_relight_lc_backup.flag"

  if [ -f "$backup_flag" ]; then
    log "S${tag}: LC backup already recorded at $backup_flag; leaving current LC dirs in place"
    return 0
  fi

  log "S${tag}: moving existing LC dirs aside with suffix LC_legacy_qlp_$STAMP"
  for orbit in "$o1" "$o2"; do
    for cam in 1 2 3 4; do
      for ccd in 1 2 3 4; do
        local ccd_root="$ROOT/orbit-$orbit/ffi/cam$cam/ccd$ccd"
        local lc_dir="$ccd_root/LC"
        local backup_dir="$ccd_root/LC_legacy_qlp_$STAMP"
        if [ -d "$lc_dir" ]; then
          if [ -e "$backup_dir" ]; then
            abort "backup dir already exists: $backup_dir"
          fi
          mv "$lc_dir" "$backup_dir" || abort "failed to move $lc_dir"
        fi
        mkdir -p "$lc_dir" || abort "failed to recreate $lc_dir"
      done
    done
  done

  date -Iseconds > "$backup_flag"
}

preserve_done_marker_once() {
  local sector=$1
  local tag=$(printf '%02d' "$sector")
  local old_done="$MARKERS/s${tag}_done.flag"
  if [ -f "$old_done" ]; then
    local legacy_done="$MARKERS/s${tag}_done_legacy_qlp_$STAMP.flag"
    log "S${tag}: preserving old done marker as $(basename "$legacy_done")"
    mv "$old_done" "$legacy_done" || abort "failed to preserve $old_done"
  fi
}

refresh_orbit_lightcurves() {
  local sector=$1
  local orbit=$2
  local orbit_tag=$3
  local tag=$(printf '%02d' "$sector")
  local sector_log="$LOG/twirlfs_relight_s${tag}"
  mkdir -p "$sector_log"

  log "S${tag}: refresh lightcurves orbit-$orbit ($orbit_tag)"
  HDF5_USE_FILE_LOCKING=FALSE \
  "$PY" "$SCRIPT_DIR/run_tglc_orbit_pipeline.py" \
    --orbit "$orbit" --sector "$sector" --orbit-tag "$orbit_tag" \
    --tglc-data-dir "$ROOT" \
    --log-dir "$LOG" \
    --python-bin "$PY" \
    --ld-library-prefix "$LD_PREFIX" \
    --catalogs-nprocs 16 \
    --cutouts-nprocs 16 \
    --lightcurves-nprocs 16 \
    --max-parallel-ccd-jobs 4 \
    --max-magnitude 20 \
    --run-label "s${tag}-twirlfs-relight" \
    --stages lightcurves \
    --no-auto-resume \
    >> "$sector_log/relight_orbit_${orbit}.log" 2>&1 \
    || abort "lightcurve refresh failed for S${tag} orbit-$orbit"
}

finalize_twirl_fs() {
  local sector=$1
  local o1=$((2 * sector + 7))
  local o2=$((o1 + 1))
  local tag=$(printf '%02d' "$sector")

  log "S${tag}: build TWIRL-FS product"
  TWIRL_SKIP_STEP_A=1 bash "$SCRIPT_DIR/finalize_sector_gpu.sh" "$sector" "$o1" "$o2" \
    || abort "TWIRL-FS finalize failed for S${tag}"
}

for sector in "$@"; do
  tag=$(printf '%02d' "$sector")
  o1=$((2 * sector + 7))
  o2=$((o1 + 1))
  log "S${tag}: start orbit-$o1 + orbit-$o2"
  preserve_done_marker_once "$sector"
  backup_lightcurves_once "$sector"
  refresh_orbit_lightcurves "$sector" "$o1" o1
  refresh_orbit_lightcurves "$sector" "$o2" o2
  finalize_twirl_fs "$sector"
  log "S${tag}: complete"
done

log "all requested sectors complete"
