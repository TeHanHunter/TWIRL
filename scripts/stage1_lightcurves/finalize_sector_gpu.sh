#!/bin/bash
# GPU finalize step: take a sector whose catalogs+cutouts are done (marker
# present) and run ePSFs+lightcurves on both orbits, then post-LC chain
# (cadence-fix + detrend + HLSP + audit + QC). Designed to run on pdogpu6.
# After all stages succeed, writes
#   /pdo/users/tehan/tglc-gpu-production/markers/s<NN>_done.flag
#
# The post-LC chain assumes:
#   - sector_orbit JSON entry exists for $SECTOR (built by build_sector_orbit_json.py)
#   - quat + qflag files exist at /pdo/qlp-data/orbit-${ORBIT_*}/ffi/run/
#
# Usage:
#   finalize_sector_gpu.sh <sector> <orbit_1> <orbit_2>
# Example:
#   finalize_sector_gpu.sh 58 123 124

set -uo pipefail

SECTOR="${1:?usage: $0 <sector> <orbit_1> <orbit_2>}"
ORBIT_1="${2:?usage: $0 <sector> <orbit_1> <orbit_2>}"
ORBIT_2="${3:?usage: $0 <sector> <orbit_1> <orbit_2>}"
SECTOR_TAG=$(printf '%02d' "$SECTOR")

REPO=/pdo/users/tehan/TWIRL
ROOT=/pdo/users/tehan/tglc-gpu-production
LOG=$ROOT/twirl_logs
MARKERS=$ROOT/markers
SECTOR_LOG=$LOG/post_lc_chain_s${SECTOR_TAG}
mkdir -p "$ROOT" "$LOG" "$MARKERS" "$SECTOR_LOG"

export PYTHONPATH=/sw/qlp-environment/.venv/lib/python3.11/site-packages
export LD_LIBRARY_PATH=/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}
export TWIRL_SECTOR_ORBIT_JSON=$ROOT/sector_orbit.json
TGLC_PY=/pdo/users/tehan/twirl-gpu-venv/bin/python

stamp() { date '+%Y-%m-%d %H:%M:%S %Z'; }
log()   { echo "[$(stamp)] $*" | tee -a "$SECTOR_LOG/chain.log"; }
abort() { log "ABORT: $1"; exit 1; }

cd "$REPO"

#######################################################################
# Step A: ePSFs + lightcurves on both orbits via the orbit driver
#         (--stages epsfs,lightcurves; auto-skip on summary.json).
#######################################################################

run_finalize_orbit() {
  local ORBIT=$1
  local ORBIT_TAG=$2
  log "  finalize orbit-$ORBIT (epsfs+lightcurves)"
  GPU_ARGS=$(bash "$REPO/scripts/stage1_lightcurves/pick_gpu_config.sh" 8 4 4)
  log "    gpu_args: $GPU_ARGS"
  $TGLC_PY scripts/stage1_lightcurves/run_tglc_orbit_pipeline.py \
    --orbit "$ORBIT" --sector "$SECTOR" --orbit-tag "$ORBIT_TAG" \
    --tglc-data-dir "$ROOT" \
    --log-dir "$LOG" \
    --python-bin /pdo/users/tehan/twirl-gpu-venv/bin/python \
    --ld-library-prefix /sw/python-versions/python-3.11.9/lib \
    --catalogs-nprocs 16 \
    --cutouts-nprocs 16 \
    --lightcurves-nprocs 16 \
    --max-magnitude 20 \
    --run-label "s${SECTOR_TAG}-gpu" \
    --stages epsfs,lightcurves \
    $GPU_ARGS \
    >> "$SECTOR_LOG/finalize_orbit_$ORBIT.log" 2>&1
}

log "step A: ePSFs+LCs (orbit-$ORBIT_1, orbit-$ORBIT_2)"
MAX_RETRIES=10
attempt=0
while (( attempt < MAX_RETRIES )); do
  attempt=$(( attempt + 1 ))
  if run_finalize_orbit "$ORBIT_1" o1 && run_finalize_orbit "$ORBIT_2" o2; then
    break
  fi
  log "  finalize attempt $attempt/$MAX_RETRIES failed; sleeping 60s"
  sleep 60
done
(( attempt <= MAX_RETRIES )) || abort "finalize stage exhausted retries"

#######################################################################
# Step B: stage run/ dirs + cadence-mismatch fix.
#######################################################################
log "step B: stage run/ dirs"
for orbit in "$ORBIT_1" "$ORBIT_2"; do
  src=/pdo/qlp-data/orbit-$orbit/ffi/run
  dst=$ROOT/orbit-$orbit/ffi/run
  [ -d "$src" ] || abort "missing src run dir $src"
  mkdir -p "$dst"
  for f in "$src"/cam?_quat.txt "$src"/cam?ccd?_qflag.txt; do
    ln -sfn "$f" "$dst/$(basename "$f")"
  done
  cat > "$dst/qlp.cfg" <<EOF
[IOSettings]
indir = $ROOT/orbit-$orbit/
orbit_id = $orbit
sector = $SECTOR

[Setup]
orbit_id = $orbit
indir = $ROOT/orbit-$orbit/ffi/

[LC]
method = qsp
bkspace_min = 0.3
bkspace_max = 2.5
penalty_coeff = 0.75
level = 3
EOF
  log "  orbit-$orbit run dir staged"
done

# Cadence-mismatch trim per orbit per cam, with a high --max-drop to handle
# TICA cross-orbit re-deliveries (S57+ ship o1/o2 dirs that include cadences
# from the other orbit; the trim drops them down to the quat-file boundary).
log "step C: cadence-mismatch fix per orbit per cam"
for orbit in "$ORBIT_1" "$ORBIT_2"; do
  for cam in 1 2 3 4; do
    quat=$ROOT/orbit-$orbit/ffi/run/cam${cam}_quat.txt
    [ -f "$quat" ] || abort "missing quat $quat"
    LCDIRS=()
    for ccd in 1 2 3 4; do
      LCDIRS+=("--lc-dir" "$ROOT/orbit-$orbit/ffi/cam$cam/ccd$ccd/LC")
    done
    HDF5_USE_FILE_LOCKING=FALSE \
    $TGLC_PY $REPO/scripts/stage1_lightcurves/fix_tglc_quat_cadence_mismatch.py \
      "${LCDIRS[@]}" --quat "$quat" --nprocs 1 --max-drop 20000 \
      >> "$SECTOR_LOG/cadence_fix_orbit-${orbit}_cam${cam}.log" 2>&1 \
      || abort "cadence-fix orbit-$orbit cam$cam failed"
  done
  log "  orbit-$orbit cadence-fix done"
done

#######################################################################
# Step D: qlp lctools detrend per orbit.
#######################################################################
log "step D: detrend"
source $REPO/scripts/activate_qlp_env.sh
for orbit in "$ORBIT_1" "$ORBIT_2"; do
  log "  detrend orbit-$orbit"
  "$TWIRL_QLP_PYTHON" $REPO/scripts/detrend_wrapper.py \
    lctools detrend -b \
    -c $ROOT/orbit-$orbit/ffi/run/qlp.cfg \
    --autolist all -n 32 \
    >> "$SECTOR_LOG/detrend_orbit_$orbit.log" 2>&1 \
    || abort "detrend orbit-$orbit failed"
done

#######################################################################
# Step E: qlp lctools hlsp at sector level.
#######################################################################
log "step E: hlsp"
mkdir -p $ROOT/hlsp_s00${SECTOR_TAG}
"$TWIRL_QLP_PYTHON" $REPO/scripts/hlsp_wrapper.py \
  lctools hlsp \
  -c $ROOT/orbit-$ORBIT_1/ffi/run/qlp.cfg \
  -s $SECTOR --autolist all \
  --flag-type spoc --flag-source fits \
  --flag-dir /pdo/qlp-data/spocflags \
  --basedir $ROOT/ \
  -o $ROOT/hlsp_s00${SECTOR_TAG} \
  -n 32 \
  >> "$SECTOR_LOG/hlsp.log" 2>&1 \
  || abort "hlsp failed"

#######################################################################
# Step F: audit + QC plots.
#######################################################################
log "step F: audit"
n_fits=$(find $ROOT/hlsp_s00${SECTOR_TAG} -name '*.fits' | wc -l)
log "  HLSP FITS count: $n_fits"
if (( n_fits == 0 )); then
  touch "$MARKERS/qc_pause.flag"
  abort "HLSP produced 0 FITS — likely detrend/cadence-fix upstream failure (see hlsp.log + detrend_orbit_*.log). qc_pause.flag set."
fi

log "step G: QC plot sample"
mkdir -p $ROOT/qc_s${SECTOR_TAG}_gpu
$TGLC_PY $REPO/scripts/stage1_lightcurves/qc_plot_hlsp_sample.py \
  --hlsp-root $ROOT/hlsp_s00${SECTOR_TAG} \
  --output-dir $ROOT/qc_s${SECTOR_TAG}_gpu \
  --n-per-bin 64 --seed 0 \
  >> "$SECTOR_LOG/qc_plot.log" 2>&1 \
  || log "  WARN: QC plot failed (non-fatal)"

# Step H: per-sector QC PDF (precision vs Tmag, sky map, example LCs) + raw
# data NPZ for re-plotting. Non-fatal — production can ship without it.
log "step H: QC sector PDF"
mkdir -p $ROOT/qc_pdf
$TGLC_PY $REPO/scripts/stage1_lightcurves/qc_sector_pdf.py \
  --hlsp-root $ROOT/hlsp_s00${SECTOR_TAG} \
  --sector $SECTOR \
  --output $ROOT/qc_pdf/s${SECTOR_TAG}_qc.pdf \
  --n-precision 6000 --n-per-bin 6 --seed 0 \
  >> "$SECTOR_LOG/qc_sector_pdf.log" 2>&1 \
  || log "  WARN: QC sector PDF failed (non-fatal)"

marker=$MARKERS/s${SECTOR_TAG}_done.flag
date -Iseconds > "$marker"
log "ALL DONE; wrote $marker"
