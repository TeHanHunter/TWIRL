#!/bin/bash
# GPU-side finalize worker: poll for the next sector whose prep is done
# (cutouts_done.flag present, _done.flag absent) and run
# finalize_sector_gpu.sh on it. ONE sector at a time so QC and review stay
# sequential (per operator request).
#
# Optional QC-pause: if /pdo/users/tehan/tglc-gpu-production/markers/qc_pause.flag
# exists, the worker waits before pulling the next sector. Touch that file
# to halt cleanly between sectors; rm it to resume.
#
# Usage:
#   finalize_worker.sh

set -uo pipefail

REPO=/pdo/users/tehan/TWIRL
ROOT=/pdo/users/tehan/tglc-gpu-production
QUEUE=$REPO/scripts/stage1_lightcurves/sector_queue.txt
MARKERS=$ROOT/markers
PAUSE_FLAG=$MARKERS/qc_pause.flag

log() { echo "[$(date '+%H:%M:%S')] [finalize] $*"; }

next_ready_sector() {
  # Lowest-numbered sector with cutouts_done.flag but not _done.flag.
  while IFS= read -r line; do
    line="${line%%#*}"; line=$(echo "$line" | xargs)
    [ -z "$line" ] && continue
    set -- $line
    local sector=$1 o1=$2 o2=$3
    local tag=$(printf '%02d' "$sector")
    [ -f "$MARKERS/s${tag}_cutouts_done.flag" ] || continue
    [ -f "$MARKERS/s${tag}_done.flag" ] && continue
    echo "$sector $o1 $o2"
    return 0
  done < "$QUEUE"
  return 1
}

while true; do
  if [ -f "$PAUSE_FLAG" ]; then
    log "qc_pause.flag present; sleeping 5 min before re-checking"
    sleep 300
    continue
  fi

  ready=$(next_ready_sector)
  if [ -z "$ready" ]; then
    log "no sector ready; polling every 5 min"
    sleep 300
    continue
  fi
  set -- $ready
  sector=$1 o1=$2 o2=$3
  log "finalize sector=$sector orbit_1=$o1 orbit_2=$o2"

  if bash $REPO/scripts/stage1_lightcurves/finalize_sector_gpu.sh "$sector" "$o1" "$o2"; then
    log "sector=$sector finalize DONE"
  else
    log "sector=$sector finalize FAILED; check $ROOT/twirl_logs/post_lc_chain_s$(printf '%02d' "$sector")/"
    log "qc_pause.flag will be set to give operator time to investigate"
    date -Iseconds > "$PAUSE_FLAG"
  fi
done
