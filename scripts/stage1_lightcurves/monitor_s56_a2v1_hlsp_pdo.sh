#!/bin/bash
# Wait for the S56 A2v1 HDF5 reproduction to finish, then build A2v1 FITS.
#
# Intended to run in tmux on PDO. It never cleans outputs. It only launches the
# ADP/ADP015-only FITS build after the reproduction log contains the completion
# marker written by run_a2v1_reproduction_pdo.sh.

set -uo pipefail

REPO=${TWIRL_REPO:-/pdo/users/tehan/TWIRL}
A2V1_ROOT=${TWIRL_A2V1_ROOT:-/pdo/users/tehan/tglc-gpu-production-A2v1}
LOG_DIR=${TWIRL_A2V1_LOG_DIR:-$A2V1_ROOT/twirl_logs}
REPRO_LOG=${TWIRL_A2V1_REPRO_LOG:-$LOG_DIR/s56_a2v1_reproduction.log}
MONITOR_LOG=${TWIRL_A2V1_HLSP_MONITOR_LOG:-$LOG_DIR/s56_a2v1_hlsp_monitor.log}
POLL_SECONDS=${TWIRL_A2V1_HLSP_MONITOR_POLL_SECONDS:-300}

stamp() { date '+%Y-%m-%d %H:%M:%S %Z'; }
log() { echo "[$(stamp)] $*" | tee -a "$MONITOR_LOG"; }

mkdir -p "$LOG_DIR"
log "monitor start; repro_log=$REPRO_LOG"

while true; do
  if grep -q "s56 A2v1 HDF5 reproduction complete" "$REPRO_LOG" 2>/dev/null; then
    log "detected HDF5 completion marker; starting A2v1 FITS build"
    cd "$REPO" || exit 1
    bash scripts/stage1_lightcurves/run_a2v1_hlsp_pdo.sh 56 119 120 \
      2>&1 | tee -a "$MONITOR_LOG"
    rc=${PIPESTATUS[0]}
    log "A2v1 FITS build finished rc=$rc"
    exit "$rc"
  fi

  if grep -q "ABORT:" "$REPRO_LOG" 2>/dev/null; then
    log "detected reproduction ABORT marker; not launching FITS build"
    exit 1
  fi

  if ! pgrep -af "run_a2v1_reproduction_pdo.sh 56" >/dev/null; then
    log "reproduction process is not running and completion marker is absent"
    exit 1
  fi

  o119_h5=$(find "$A2V1_ROOT/orbit-119/ffi" -path "*/LC/*.h5" 2>/dev/null | wc -l | tr -d " ")
  o120_h5=$(find "$A2V1_ROOT/orbit-120/ffi" -path "*/LC/*.h5" 2>/dev/null | wc -l | tr -d " ")
  o119_epsf=$(find "$A2V1_ROOT/orbit-119/ffi" -path "*/epsf/epsf_*.npy" 2>/dev/null | wc -l | tr -d " ")
  o120_epsf=$(find "$A2V1_ROOT/orbit-120/ffi" -path "*/epsf/epsf_*.npy" 2>/dev/null | wc -l | tr -d " ")
  log "waiting; orbit119 epsf=$o119_epsf h5=$o119_h5; orbit120 epsf=$o120_epsf h5=$o120_h5"
  sleep "$POLL_SECONDS"
done
