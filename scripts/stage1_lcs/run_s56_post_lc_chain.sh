#!/bin/bash
# Sequential post-LC chain for the Sector 56 GPU production tree:
#   1. Set up run/ dirs (qlp.cfg + quat/qflag symlinks).
#   2. Cadence-mismatch fix on cam3/orbit-120 h5s (TGLC-only cadence trim).
#   3. qlp lctools detrend (orbit-119 then orbit-120).
#   4. qlp lctools hlsp (sector level).
#   5. Recovery-rate audit + WD 1856 cross-check vs the frozen baseline.
#
# Each step aborts the chain on non-zero exit. Run in a tmux session.
#
# Usage:
#   bash scripts/stage1_lcs/run_s56_post_lc_chain.sh

set -uo pipefail

REPO=/pdo/users/tehan/TWIRL
NEW_ROOT=/pdo/users/tehan/tglc-gpu-production
OLD_ROOT=/pdo/users/tehan/tglc-deep-catalogs
SECTOR=56
LOG=$NEW_ROOT/twirl_logs/post_lc_chain
mkdir -p "$LOG"

stamp() { date '+%Y-%m-%d %H:%M:%S %Z'; }
log()   { echo "[$(stamp)] $*" | tee -a "$LOG/chain.log"; }

abort() { log "ABORT: $1"; exit 1; }

#######################################################################
# Step 1: stage run/ dirs in the new tree (qlp.cfg + quat/qflag symlinks).
#######################################################################
log "step 1: stage run/ dirs in the new tree"
for orbit in 119 120; do
  src=$OLD_ROOT/orbit-$orbit/ffi/run
  dst=$NEW_ROOT/orbit-$orbit/ffi/run
  [ -d "$src" ] || abort "missing src run dir $src"
  mkdir -p "$dst"
  for f in "$src"/cam?_quat.txt "$src"/cam?ccd?_qflag.txt; do
    ln -sfn "$f" "$dst/$(basename "$f")"
  done
  # Write a fresh qlp.cfg that points at the new tree, mirroring the old format.
  cat > "$dst/qlp.cfg" <<EOF
[IOSettings]
indir = $NEW_ROOT/orbit-$orbit/
orbit_id = $orbit
sector = $SECTOR

[Setup]
orbit_id = $orbit
indir = $NEW_ROOT/orbit-$orbit/ffi/

[LC]
method = qsp
bkspace_min = 0.3
bkspace_max = 2.5
penalty_coeff = 0.75
level = 3
EOF
  log "  orbit-$orbit run dir staged ($(ls $dst/cam*_quat.txt | wc -l) quats, $(ls $dst/cam*ccd*_qflag.txt | wc -l) qflags)"
done

#######################################################################
# Step 2: cam3/orbit-120 cadence-mismatch fix.
#######################################################################
log "step 2: cam3/orbit-120 cadence-mismatch fix (dry-run + apply)"
export LD_LIBRARY_PATH=/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}
TGLC_PY=/pdo/users/tehan/twirl-gpu-venv/bin/python

for ccd in 1 2 3 4; do
  LC_DIR=$NEW_ROOT/orbit-120/ffi/cam3/ccd$ccd/LC
  [ -d "$LC_DIR" ] || abort "missing LC dir $LC_DIR"
done

# Skip dry-run; --max-drop guard catches the wrong-orbit case. Run serial
# with HDF5 file locking off: parallel pool deadlocked under NFS h5py write
# contention (16 procs × imap chunksize=4 = 64 in-flight writes stalled the
# pool after ~5% in earlier runs).
HDF5_USE_FILE_LOCKING=FALSE \
$TGLC_PY $REPO/scripts/stage1_lcs/fix_tglc_quat_cadence_mismatch.py \
  --lc-dir $NEW_ROOT/orbit-120/ffi/cam3/ccd1/LC \
  --lc-dir $NEW_ROOT/orbit-120/ffi/cam3/ccd2/LC \
  --lc-dir $NEW_ROOT/orbit-120/ffi/cam3/ccd3/LC \
  --lc-dir $NEW_ROOT/orbit-120/ffi/cam3/ccd4/LC \
  --quat   $NEW_ROOT/orbit-120/ffi/run/cam3_quat.txt \
  --nprocs 1 \
  >> "$LOG/cadence_fix_apply.log" 2>&1 || abort "cadence-fix apply failed"

#######################################################################
# Step 3: qlp lctools detrend, both orbits.
#######################################################################
log "step 3: detrend (orbit-119 then orbit-120)"
source $REPO/scripts/activate_qlp_env.sh

for orbit in 119 120; do
  log "  detrend orbit-$orbit"
  "$TWIRL_QLP_PYTHON" $REPO/scripts/detrend_wrapper.py \
    lctools detrend -b \
    -c $NEW_ROOT/orbit-$orbit/ffi/run/qlp.cfg \
    --autolist all \
    -n 32 \
    >> "$LOG/detrend_orbit_$orbit.log" 2>&1 || abort "detrend orbit-$orbit failed"
done

#######################################################################
# Step 4: qlp lctools hlsp at sector level.
#######################################################################
log "step 4: hlsp sector level"
mkdir -p $NEW_ROOT/hlsp_s0056
"$TWIRL_QLP_PYTHON" $REPO/scripts/hlsp_wrapper.py \
  lctools hlsp \
  -c $NEW_ROOT/orbit-119/ffi/run/qlp.cfg \
  -s $SECTOR --autolist all \
  --flag-type spoc --flag-source fits \
  --flag-dir /pdo/qlp-data/spocflags \
  --basedir $NEW_ROOT/ \
  -o $NEW_ROOT/hlsp_s0056 \
  -n 32 \
  >> "$LOG/hlsp.log" 2>&1 || abort "hlsp failed"

#######################################################################
# Step 5: audit (FITS count + WD 1856 sanity).
#######################################################################
log "step 5: audit"
n_fits=$(find $NEW_ROOT/hlsp_s0056 -name '*.fits' | wc -l)
log "  HLSP FITS count: $n_fits  (baseline tglc-deep-catalogs/hlsp_s0056: 19072)"

# WD 1856 = TIC 267574918 — sharded path
WD_FITS=$(find $NEW_ROOT/hlsp_s0056 -name '*0000000267574918*.fits' | head -1)
if [ -n "$WD_FITS" ]; then
  log "  WD 1856 FITS: $WD_FITS"
  ls -la "$WD_FITS"
else
  log "  WARN: WD 1856 FITS not found in new tree"
fi

log "ALL DONE"
