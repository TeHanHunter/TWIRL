#!/bin/bash
set -e
cd /pdo/users/tehan/tglc-deep-catalogs
source /pdo/users/tehan/TWIRL/scripts/activate_qlp_env.sh
mkdir -p /pdo/users/tehan/tglc-deep-catalogs/twirl_logs/hlsp_s0056_r2
exec "$TWIRL_QLP_PYTHON" /pdo/users/tehan/TWIRL/scripts/hlsp_wrapper.py \
    lctools hlsp \
    -c /pdo/users/tehan/tglc-deep-catalogs/orbit-119/ffi/run/qlp.cfg \
    -s 56 --autolist all \
    --flag-type spoc --flag-source fits \
    --basedir /pdo/users/tehan/tglc-deep-catalogs/ \
    --flag-dir /pdo/qlp-data/spocflags/ \
    -o /pdo/users/tehan/tglc-deep-catalogs/hlsp_s0056 \
    -n 32
