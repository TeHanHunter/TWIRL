#!/usr/bin/env bash
# Run recovery-mode aperture/grid sweeps on the S56 pre-detrend injection set.
#
# This does not render LEO reports. It reruns BLS on the injected HDF5 and
# summarizes strict top-1, exact top-N, harmonic top-N, and unmatched recovery
# modes by empirical signal SNR.
set -euo pipefail

cd "${TWIRL_ROOT:-/pdo/users/tehan/TWIRL}"

PYTHON_BIN="${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}"
export LD_LIBRARY_PATH="/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="src"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/pdo/users/tehan/.cache/matplotlib}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

mkdir -p "${MPLCONFIGDIR}"

INJECTION_H5="${INJECTION_H5:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_1k_predetrend_batman_depthgrid_adp_compare/injected_lightcurves.h5}"
SURVIVAL_CSV="${SURVIVAL_CSV:-reports/stage5_validation/s56_1k_predetrend_batman_adp_only_review_queue_pdo/predet_signal_survival_snr_diagnostics.csv}"
OUT_DIR="${OUT_DIR:-reports/stage5_validation/s56_1k_predetrend_recovery_mode_sweep_pdo}"
N_INJECTIONS="${N_INJECTIONS:-1000}"
WORKERS="${WORKERS:-16}"
REUSE="${REUSE:-0}"

sweep_args=(
  --sweep "adp_sml_5k:DET_FLUX_ADP_SML:5000"
  --sweep "adp_5k:DET_FLUX_ADP:5000"
  --sweep "adp_lag_5k:DET_FLUX_ADP_LAG:5000"
  --sweep "det_sml_5k:DET_FLUX_SML:5000"
  --sweep "det_5k:DET_FLUX:5000"
  --sweep "det_lag_5k:DET_FLUX_LAG:5000"
  --sweep "adp_sml_200k:DET_FLUX_ADP_SML:200000"
  --sweep "adp_200k:DET_FLUX_ADP:200000"
)

cmd=(
  "${PYTHON_BIN}" scripts/stage5_validation/run_injection_recovery_mode_sweep.py
  --injection-h5 "${INJECTION_H5}"
  --survival-csv "${SURVIVAL_CSV}"
  --out-dir "${OUT_DIR}"
  --n-injections "${N_INJECTIONS}"
  --workers "${WORKERS}"
  "${sweep_args[@]}"
)
if [[ "${REUSE}" == "1" ]]; then
  cmd+=(--reuse)
fi

echo "[pdo-recovery-sweep] start $(date -Is)"
echo "[pdo-recovery-sweep] python=${PYTHON_BIN}"
echo "[pdo-recovery-sweep] injection_h5=${INJECTION_H5}"
echo "[pdo-recovery-sweep] survival_csv=${SURVIVAL_CSV}"
echo "[pdo-recovery-sweep] out_dir=${OUT_DIR}"
"${cmd[@]}"
echo "[pdo-recovery-sweep] complete $(date -Is)"
echo "[pdo-recovery-sweep] overview: ${OUT_DIR}/sweep_overview.csv"
