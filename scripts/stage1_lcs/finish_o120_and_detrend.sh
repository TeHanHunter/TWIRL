#!/usr/bin/env bash
# Finish orbit 120 (cam3_ccd2 + cam4_ccd4) then auto-launch detrend for orbit 120.
# Run inside a tmux session on pdogpu1.
set -euo pipefail

REPO=/pdo/users/tehan/TWIRL
TGLC_DATA=/pdo/users/tehan/tglc-deep-catalogs
PYTHON=/sw/qlp-environment/.venv/bin/python
QLP_PYTHON=/pdo/app/qlp-environment/.venv/bin/python
FORK=/pdo/users/tehan/tess-gaia-light-curve-twirl
LD_PREFIX="/pdo/app/anaconda/anaconda2-4.4.0/lib:/pdo/app/python-versions/python-3.11.9/lib"
LOG_DIR=${REPO}/logs/o120_finish

mkdir -p "${LOG_DIR}"

export LD_LIBRARY_PATH="${LD_PREFIX}:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${FORK}:${PYTHONPATH:-}"

echo "[finish-o120] $(date): Starting pipeline for cam3_ccd2 + cam4_ccd4"

${PYTHON} ${REPO}/scripts/stage1_lcs/run_tglc_orbit_pipeline.py \
    --orbit 120 --sector 56 --orbit-tag o2 \
    --ccd 3,2 --ccd 4,4 \
    --max-parallel-ccd-jobs 2 \
    --epsfs-nprocs 16 --lightcurves-nprocs 16 \
    --log-dir "${LOG_DIR}" \
    --run-label twirl-s56-o120-finish \
    2>&1 | tee "${LOG_DIR}/twirl-s56-o120-finish.log"

echo "[finish-o120] $(date): Pipeline complete. Verifying h5 counts..."

TOTAL_H5=0
for cam in 1 2 3 4; do
    for ccd in 1 2 3 4; do
        count=$(ls "${TGLC_DATA}/orbit-120/ffi/cam${cam}/ccd${ccd}/LC/"*.h5 2>/dev/null | wc -l)
        echo "  cam${cam}_ccd${ccd}: ${count} h5"
        TOTAL_H5=$((TOTAL_H5 + count))
    done
done
echo "[finish-o120] Total orbit 120 h5 files: ${TOTAL_H5}"

if [ "${TOTAL_H5}" -lt 1000 ]; then
    echo "[finish-o120] ERROR: Total h5 count unexpectedly low (${TOTAL_H5}). Aborting detrend."
    exit 1
fi

echo "[finish-o120] $(date): Launching orbit 120 detrend (nprocs=1)..."

export PYTHONPATH="${FORK}"
${QLP_PYTHON} ${REPO}/scripts/detrend_wrapper.py \
    lctools detrend \
    -c "${TGLC_DATA}/orbit-120/ffi/run/qlp.cfg" \
    --autolist 1:1,2,3,4 2:1,2,3,4 3:1,2,3,4 4:1,2,3,4 \
    --best -n 1 \
    --logfile "${REPO}/logs/detrend_orbit120.log" \
    2>&1 | tee "${REPO}/logs/detrend_orbit120_stdout.log"

echo "[finish-o120] $(date): Orbit 120 detrend complete."
echo "DONE"
