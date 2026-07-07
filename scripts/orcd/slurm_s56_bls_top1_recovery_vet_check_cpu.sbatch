#!/bin/bash
#SBATCH -p mit_normal
#SBATCH --job-name=twirl-top1-vet
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
WORKERS="${TWIRL_TOP1_VET_WORKERS:-${SLURM_CPUS_PER_TASK:-8}}"
QUEUE_CSV="${TWIRL_TOP1_VET_QUEUE_CSV:-reports/stage5_validation/s56_bls_top1_recovery_vet_check/review_queue_100.csv}"
INJECTION_H5="${TWIRL_TOP1_VET_INJECTION_H5:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_20k_predetrend_batman_periodradius_bright_balanced/injected_lightcurves.h5}"
HLSP_ROOT="${TWIRL_TOP1_VET_HLSP_ROOT:-data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare}"
OUT_DIR="${TWIRL_TOP1_VET_OUT_DIR:-reports/stage5_validation/s56_bls_top1_recovery_vet_check/twirl_vet_sheets}"
METRICS_CSV="${TWIRL_TOP1_VET_METRICS_CSV:-reports/stage5_validation/s56_bls_top1_recovery_vet_check/twirl_vet_metrics.csv}"
SUMMARY_JSON="${TWIRL_TOP1_VET_SUMMARY_JSON:-reports/stage5_validation/s56_bls_top1_recovery_vet_check/twirl_two_aperture_vet_summary.json}"
N_PERIODS="${TWIRL_TOP1_VET_N_PERIODS:-20000}"
N_PEAKS="${TWIRL_TOP1_VET_N_PEAKS:-10}"
OVERWRITE="${TWIRL_TOP1_VET_OVERWRITE:-1}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs
mkdir -p /orcd/data/mki_aryeh/001/twirl/mpl-cache
cd "${REPO}"

export PYTHONPATH=src
export MPLCONFIGDIR=/orcd/data/mki_aryeh/001/twirl/mpl-cache
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[orcd-top1-vet] host=$(hostname)"
echo "[orcd-top1-vet] repo=${REPO}"
echo "[orcd-top1-vet] python=${PYTHON}"
echo "[orcd-top1-vet] queue=${QUEUE_CSV}"
echo "[orcd-top1-vet] injection_h5=${INJECTION_H5}"
echo "[orcd-top1-vet] out_dir=${OUT_DIR}"
echo "[orcd-top1-vet] workers=${WORKERS} n_periods=${N_PERIODS}"

"${PYTHON}" - <<'PY'
import importlib
import sys

missing = []
for name in ("numpy", "pandas", "h5py", "astropy", "scipy", "matplotlib"):
    try:
        importlib.import_module(name)
    except Exception as exc:  # pragma: no cover - cluster preflight
        missing.append(f"{name}: {type(exc).__name__}: {exc}")
if missing:
    print("[orcd-top1-vet] missing required Python packages:", file=sys.stderr)
    for item in missing:
        print(f"  {item}", file=sys.stderr)
    sys.exit(3)
PY

args=(
  scripts/stage5_validation/render_two_aperture_vet_sheets.py
  --queue-csv "${QUEUE_CSV}"
  --hlsp-root "${HLSP_ROOT}"
  --injection-h5 "${INJECTION_H5}"
  --out-dir "${OUT_DIR}"
  --metrics-csv "${METRICS_CSV}"
  --summary-json "${SUMMARY_JSON}"
  --branch-name current_adp
  --apertures DET_FLUX_ADP_SML DET_FLUX_SML
  --workers "${WORKERS}"
  --n-periods "${N_PERIODS}"
  --n-peaks "${N_PEAKS}"
)

if [[ "${OVERWRITE}" == "1" ]]; then
  args+=(--overwrite)
fi

"${PYTHON}" "${args[@]}"

echo "[orcd-top1-vet] complete $(date -Is)"
