#!/usr/bin/env bash
set -euo pipefail

cd /pdo/users/tehan/TWIRL

export LD_LIBRARY_PATH="/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH=src
export MPLCONFIGDIR=/pdo/users/tehan/.cache/matplotlib
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

/sw/qlp-environment/.venv/bin/python \
  scripts/stage5_validation/sweep_predetrend_detrending_methods.py \
  --injection-h5 data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_1k_predetrend_batman_depthgrid_adp_compare/injected_lightcurves.h5 \
  --out-dir reports/stage5_validation/s56_1k_predetrend_detrending_method_sweep_extended_small_pdo \
  --raw-aperture Small \
  --method constant \
  --method poly2_gap05 \
  --method poly3_gap05 \
  --method median06_gap05 \
  --method median10_gap05 \
  --method median20_gap05 \
  --method savgol06_p2_gap05 \
  --method savgol10_p2_gap05 \
  --method savgol20_p2_gap05 \
  --method current_adp03q \
  --method spline_q05_gap02 \
  --method current_uniform08_gap05 \
  --method oracle_adp03q \
  --method oracle_q05_gap02
