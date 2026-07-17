#!/bin/bash
# CPU-side prep worker: pop sectors off the queue and run prep_sector_cpu.sh.
# Designed to run as 1 or 2 concurrent tmux sessions on pdogpu1.
#
# Each worker takes a "lease" by atomically renaming the next available queue
# line into a per-worker lease file, runs prep, then deletes the lease.
#
# Usage:
#   prep_worker.sh <worker_id>      # e.g. prep_worker.sh A   (or B for the 2nd one)
#
# Stop conditions: queue exhausted, or backlog of unfinalized sectors >= 4
# (don't run too far ahead of pdogpu6).

set -uo pipefail

WORKER_ID="${1:?usage: $0 <worker_id e.g. A>}"

REPO=/pdo/users/tehan/TWIRL
ROOT=/pdo/users/tehan/tglc-gpu-production
QUEUE=$REPO/scripts/stage1_lcs/sector_queue.txt
MARKERS=$ROOT/markers
LEASES=$ROOT/markers/leases
mkdir -p "$MARKERS" "$LEASES"

log()   { echo "[$(date '+%H:%M:%S')] [prep-$WORKER_ID] $*"; }

next_sector() {
  # Atomically pick the next sector from the queue that:
  #   - has no s<NN>_cutouts_done.flag
  #   - has no active lease
  # Print "sector orbit_1 orbit_2" or empty if nothing available.
  while IFS= read -r line; do
    line="${line%%#*}"; line=$(echo "$line" | xargs)
    [ -z "$line" ] && continue
    set -- $line
    local sector=$1 o1=$2 o2=$3
    local tag=$(printf '%02d' "$sector")
    local marker="$MARKERS/s${tag}_cutouts_done.flag"
    local lease="$LEASES/s${tag}.lease"
    [ -f "$marker" ] && continue
    # Try to take the lease atomically.
    if ( set -o noclobber; > "$lease" ) 2>/dev/null; then
      echo "$WORKER_ID:$(date -Iseconds)" > "$lease"
      echo "$sector $o1 $o2"
      return 0
    fi
  done < "$QUEUE"
  return 1
}

backlog_size() {
  local n=0
  for f in "$MARKERS"/s??_cutouts_done.flag; do
    [ -f "$f" ] || continue
    local tag=$(basename "$f" | sed 's/^s\(..\)_cutouts_done.flag/\1/')
    [ -f "$MARKERS/s${tag}_done.flag" ] && continue
    n=$(( n + 1 ))
  done
  echo "$n"
}

while true; do
  bl=$(backlog_size)
  if (( bl >= 4 )); then
    log "backlog of $bl unfinalized sectors; pausing 10 min so pdogpu6 catches up"
    sleep 600
    continue
  fi

  next=$(next_sector)
  if [ -z "$next" ]; then
    log "queue empty or all leased; pausing 10 min"
    sleep 600
    continue
  fi
  set -- $next
  sector=$1 o1=$2 o2=$3
  tag=$(printf '%02d' "$sector")
  log "prep sector=$sector orbit_1=$o1 orbit_2=$o2"

  if bash $REPO/scripts/stage1_lcs/prep_sector_cpu.sh "$sector" "$o1" "$o2"; then
    log "sector=$sector prep DONE"
  else
    log "sector=$sector prep FAILED; check $ROOT/twirl_logs/s${tag}-prep/"
  fi
  rm -f "$LEASES/s${tag}.lease"
done
