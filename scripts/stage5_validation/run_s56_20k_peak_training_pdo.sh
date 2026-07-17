#!/usr/bin/env bash
# Build a peak-level BLS truth table for the S56 20k pre-detrend injections.
#
# This is a detector/ranker training product: one row per BLS peak per aperture,
# labeled against the known injected ephemeris.
set -euo pipefail

cd /pdo/users/tehan/TWIRL

export LD_LIBRARY_PATH="/pdo/app/anaconda/anaconda2-4.4.0/lib:/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

PYTHON="${PYTHON:-/sw/qlp-environment/.venv/bin/python}"
INJECTION_H5="${INJECTION_H5:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_20k_predetrend_batman_periodradius_bright_balanced/injected_lightcurves.h5}"
OUT_TABLE="${OUT_TABLE:-reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training/s56_20k_injection_bls_peaks.csv}"
VERIFY_JSON="${VERIFY_JSON:-${OUT_TABLE%.*}_verification.json}"
WORKERS="${TWIRL_PEAK_WORKERS:-16}"
N_PEAKS="${TWIRL_BLS_N_PEAKS:-20}"
N_PERIODS="${TWIRL_BLS_N_PERIODS:-200000}"
MIN_INJECTIONS="${TWIRL_PEAK_MIN_INJECTIONS:-10000}"
MIN_CANDIDATE_ROWS="${TWIRL_PEAK_MIN_CANDIDATE_ROWS:-100000}"
MIN_APERTURES="${TWIRL_PEAK_MIN_APERTURES:-2}"
MIN_POSITIVE_PEAK_RANKS="${TWIRL_PEAK_MIN_POSITIVE_RANKS:-10}"

"${PYTHON}" scripts/stage5_validation/build_injection_peak_training_table.py \
  --injection-h5 "${INJECTION_H5}" \
  --out-table "${OUT_TABLE}" \
  --apertures DET_FLUX_ADP_SML DET_FLUX_SML \
  --n-peaks "${N_PEAKS}" \
  --n-periods "${N_PERIODS}" \
  --workers "${WORKERS}" \
  --overwrite

"${PYTHON}" scripts/stage5_validation/verify_injection_peak_training_table.py \
  --peak-table "${OUT_TABLE}" \
  --out-json "${VERIFY_JSON}" \
  --min-injections "${MIN_INJECTIONS}" \
  --min-candidate-rows "${MIN_CANDIDATE_ROWS}" \
  --min-apertures "${MIN_APERTURES}" \
  --min-positive-peak-ranks "${MIN_POSITIVE_PEAK_RANKS}" \
  --require-cadence-diagnostics

echo "[peak-training] verification passed: ${VERIFY_JSON}"
