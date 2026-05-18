#!/usr/bin/env bash
# Convenience launcher for the S96 vanilla-TGLC vs MAST-QLP QA benchmark on PDO.
#
# This runs the upstream TGLC code at TeHanHunter/TESS_Gaia_Light_Curve
# (mitpdo branch) — NOT the MIT QLP-integrated fork. It ingests the TICA FFIs
# at /pdo/qlp-data/tica-delivery/s0096/camC-ccdD/ via tglc.ffi.ffi_tica,
# builds 196 source_XX_YY.pkl patches per CCD, and runs tglc.run.lc_per_ccd
# for the specific CCDs that host the selected transit targets.
#
#   tmux new -s twirl-s96-qa
#   cd /pdo/users/tehan/TWIRL
#   bash scripts/stage1_lightcurves/run_s96_qa.sh 2>&1 | tee logs/s96_qa.log
#
# After it finishes, rsync the TGLC FITS back locally:
#   rsync -av pdogpu1:/pdo/users/tehan/tglc-vanilla-s96qa/lc/ \
#       ~/PycharmProjects/TWIRL/reports/stage1_lightcurves/s96_qa/tglc_vanilla_lc/
# then locally (plotter consumes vanilla TGLC FITS by Gaia DR3 id):
#   python reports/exploratory/_scripts/plot_tglc_tica_vs_qlp_s96.py \
#       --tglc-vanilla-lc-dir reports/stage1_lightcurves/s96_qa/tglc_vanilla_lc

set -euo pipefail

TWIRL_ROOT=${TWIRL_ROOT:-/pdo/users/tehan/TWIRL}
TARGETS_JSON=${TARGETS_JSON:-$TWIRL_ROOT/reports/stage1_lightcurves/s96_qa/targets.json}
TGLC_DATA_DIR=${TGLC_DATA_DIR:-/pdo/users/tehan/tglc-vanilla-s96qa}
TGLC_FORK_PATH=${TGLC_FORK_PATH:-/pdo/users/tehan/TESS_Gaia_Light_Curve}

# Match the shell env the PDO python needs; the vanilla tglc package imports
# astroquery etc. from /sw/qlp-environment/.venv and libffi.so.6 from the
# anaconda2 tree (CLAUDE.md PDO-Specific Notes).
export LD_LIBRARY_PATH=/pdo/app/anaconda/anaconda2-4.4.0/lib:/pdo/app/python-versions/python-3.11.9/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export PYTHONPATH=$TGLC_FORK_PATH${PYTHONPATH:+:$PYTHONPATH}
TGLC_PYTHON=${TGLC_PYTHON:-/sw/qlp-environment/.venv/bin/python}

# Cap BLAS threads before the Pool spawns (CLAUDE.md: multiprocessing on PDO).
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# PDO injects a self-signed cert chain; certifi's bundle rejects it. Point
# astroquery / requests at the system CA bundle so Gaia TAP + MAST queries
# made per Source() work. Without this, Source(x=0,y=0,...) crashes with
# "TypeError: Cannot convert None to table column" after an SSL retry.
export REQUESTS_CA_BUNDLE=/etc/pki/tls/certs/ca-bundle.crt
export SSL_CERT_FILE=/etc/pki/tls/certs/ca-bundle.crt
export CURL_CA_BUNDLE=/etc/pki/tls/certs/ca-bundle.crt

mkdir -p "$TGLC_DATA_DIR" "$TWIRL_ROOT/logs"

"$TGLC_PYTHON" "$TWIRL_ROOT/scripts/stage1_lightcurves/run_s96_qa_vanilla.py" \
    --targets-json "$TARGETS_JSON" \
    --local-directory "$TGLC_DATA_DIR" \
    --sector 96 \
    --size 150 \
    --cores 16
