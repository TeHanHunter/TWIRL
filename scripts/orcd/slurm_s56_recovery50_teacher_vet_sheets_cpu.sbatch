#!/bin/bash
#SBATCH -p mit_normal
#SBATCH --job-name=twirl-rec50-vet
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=96G
#SBATCH -t 04:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
WORKERS="${TWIRL_REC50_VET_WORKERS:-${SLURM_CPUS_PER_TASK:-8}}"
QUEUE_DIR="${TWIRL_REC50_QUEUE_DIR:-reports/stage5_validation/s56_recovery50_teacher_queue}"
REAL_QUEUE_CSV="${TWIRL_REC50_REAL_QUEUE_CSV:-${QUEUE_DIR}/review_queue_1k_real.csv}"
INJ_QUEUE_CSV="${TWIRL_REC50_INJ_QUEUE_CSV:-${QUEUE_DIR}/review_queue_1k_injected.csv}"
REAL_LC_EXPORT_H5="${TWIRL_REC50_REAL_LC_EXPORT_H5:-data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp_lc_export_pdo.h5}"
INJECTION_H5="${TWIRL_REC50_INJECTION_H5:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_20k_predetrend_batman_periodradius_bright_balanced/injected_lightcurves.h5}"
HLSP_ROOT="${TWIRL_REC50_HLSP_ROOT:-data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare}"
OUT_DIR="${TWIRL_REC50_OUT_DIR:-${QUEUE_DIR}/twirl_vet_sheets}"
N_PERIODS="${TWIRL_REC50_N_PERIODS:-20000}"
N_PEAKS="${TWIRL_REC50_N_PEAKS:-10}"
OVERWRITE="${TWIRL_REC50_OVERWRITE:-1}"
LIMIT="${TWIRL_REC50_LIMIT:-0}"

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

echo "[orcd-rec50-vet] host=$(hostname)"
echo "[orcd-rec50-vet] repo=${REPO}"
echo "[orcd-rec50-vet] python=${PYTHON}"
echo "[orcd-rec50-vet] real_queue=${REAL_QUEUE_CSV}"
echo "[orcd-rec50-vet] inj_queue=${INJ_QUEUE_CSV}"
echo "[orcd-rec50-vet] real_lc_export_h5=${REAL_LC_EXPORT_H5}"
echo "[orcd-rec50-vet] injection_h5=${INJECTION_H5}"
echo "[orcd-rec50-vet] out_dir=${OUT_DIR}"
echo "[orcd-rec50-vet] workers=${WORKERS} n_periods=${N_PERIODS} limit=${LIMIT}"

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
    print("[orcd-rec50-vet] missing required Python packages:", file=sys.stderr)
    for item in missing:
        print(f"  {item}", file=sys.stderr)
    sys.exit(3)
PY

common_args=(
  --hlsp-root "${HLSP_ROOT}"
  --out-dir "${OUT_DIR}"
  --branch-name current_adp
  --workers "${WORKERS}"
  --n-periods "${N_PERIODS}"
  --n-peaks "${N_PEAKS}"
)

if [[ "${OVERWRITE}" == "1" ]]; then
  common_args+=(--overwrite)
fi
if [[ "${LIMIT}" != "0" ]]; then
  common_args+=(--limit "${LIMIT}")
fi

echo "[orcd-rec50-vet] rendering real rows $(date -Is)"
"${PYTHON}" scripts/stage5_validation/render_two_aperture_vet_sheets.py \
  --queue-csv "${REAL_QUEUE_CSV}" \
  --lc-export-h5 "${REAL_LC_EXPORT_H5}" \
  --metrics-csv "${QUEUE_DIR}/twirl_vet_metrics_real.csv" \
  --summary-json "${QUEUE_DIR}/twirl_two_aperture_vet_summary_real.json" \
  --apertures DET_FLUX_ADP_SML DET_FLUX_ADP \
  "${common_args[@]}"

echo "[orcd-rec50-vet] rendering injected rows $(date -Is)"
"${PYTHON}" scripts/stage5_validation/render_two_aperture_vet_sheets.py \
  --queue-csv "${INJ_QUEUE_CSV}" \
  --injection-h5 "${INJECTION_H5}" \
  --metrics-csv "${QUEUE_DIR}/twirl_vet_metrics_injected.csv" \
  --summary-json "${QUEUE_DIR}/twirl_two_aperture_vet_summary_injected.json" \
  --apertures DET_FLUX_ADP_SML DET_FLUX_SML \
  "${common_args[@]}"

echo "[orcd-rec50-vet] complete $(date -Is)"
