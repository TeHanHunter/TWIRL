#!/usr/bin/env bash
# Build a publication-oriented S56 pre-detrend BATMAN injection-recovery map:
#
# - period-radius grid, not period-depth grid
# - transit-conditioned random impact parameters / inclinations
# - equal target sampling in Tmag < 17, 17-18, 18-19, and > 19
#
# This keeps the y-axis physical (companion radius) while marginalizing the
# recovery fraction over the physically allowed chord-length/depth scatter at
# fixed period and radius.
set -euo pipefail

export N_INJECTIONS="${N_INJECTIONS:-20000}"
export RANDOM_STATE="${RANDOM_STATE:-5624}"
export SAMPLING_MODE="${SAMPLING_MODE:-period_radius_grid}"
export GRID_PERIOD_RANGE_D="${GRID_PERIOD_RANGE_D:-0.12,13.0}"
export GRID_RADIUS_RANGE_REARTH="${GRID_RADIUS_RANGE_REARTH:-0.18,18.0}"
export GRID_PERIOD_BINS="${GRID_PERIOD_BINS:-50}"
export GRID_RADIUS_BINS="${GRID_RADIUS_BINS:-50}"
export TARGET_TMAG_BIN_EDGES="${TARGET_TMAG_BIN_EDGES:-0,17,18,19,inf}"
export TARGET_TMAG_BIN_WEIGHTS="${TARGET_TMAG_BIN_WEIGHTS:-1,1,1,1}"
export INJECTION_DIR="${INJECTION_DIR:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_20k_predetrend_batman_periodradius_bright_balanced}"
export OUT_DIR="${OUT_DIR:-reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo}"
export SWEEP_NAME="${SWEEP_NAME:-small_pair_200k}"

exec scripts/stage5_validation/run_s56_predetrend_dense_bls_map_pdo.sh
