#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-twoap-vet
#SBATCH -N 1
#SBATCH -c 48
#SBATCH --mem=192G
#SBATCH -t 1-00:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
WORKERS="${TWIRL_TWOAP_WORKERS:-${SLURM_CPUS_PER_TASK:-16}}"
QUEUE_CSV="${TWIRL_TWOAP_QUEUE_CSV:-reports/stage5_validation/s56_mixed_teacher_queue_pdo/review_queue_1k.csv}"
LC_EXPORT_H5="${TWIRL_TWOAP_LC_EXPORT_H5:-data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp015_lc_export_pdo.h5}"
HLSP_ROOT="${TWIRL_TWOAP_HLSP_ROOT:-data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_adp015q_compare}"
OUT_DIR="${TWIRL_TWOAP_OUT_DIR:-reports/stage5_validation/s56_mixed_teacher_queue_pdo/twirl_vet_sheets_adp015q_orcd}"
METRICS_CSV="${TWIRL_TWOAP_METRICS_CSV:-reports/stage5_validation/s56_mixed_teacher_queue_pdo/twirl_vet_metrics_adp015q_orcd.csv}"
SUMMARY_JSON="${TWIRL_TWOAP_SUMMARY_JSON:-reports/stage5_validation/s56_mixed_teacher_queue_pdo/twirl_two_aperture_vet_summary_adp015q_orcd.json}"
BRANCH_NAME="${TWIRL_TWOAP_BRANCH_NAME:-twirl_fs_v2_adp015q}"
APERTURE_TAG="${TWIRL_TWOAP_APERTURE_TAG:-ADP015}"
LIMIT="${TWIRL_TWOAP_LIMIT:-0}"
N_PERIODS="${TWIRL_TWOAP_N_PERIODS:-20000}"
N_PEAKS="${TWIRL_TWOAP_N_PEAKS:-10}"
OVERWRITE="${TWIRL_TWOAP_OVERWRITE:-0}"

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

echo "[orcd-twoap] host=$(hostname)"
echo "[orcd-twoap] repo=${REPO}"
echo "[orcd-twoap] python=${PYTHON}"
echo "[orcd-twoap] queue=${QUEUE_CSV}"
echo "[orcd-twoap] lc_export_h5=${LC_EXPORT_H5}"
echo "[orcd-twoap] out_dir=${OUT_DIR}"
echo "[orcd-twoap] summary_json=${SUMMARY_JSON}"
echo "[orcd-twoap] workers=${WORKERS} limit=${LIMIT} n_periods=${N_PERIODS}"

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
    print("[orcd-twoap] missing required Python packages:", file=sys.stderr)
    for item in missing:
        print(f"  {item}", file=sys.stderr)
    sys.exit(3)
PY

args=(
  scripts/stage5_validation/render_two_aperture_vet_sheets.py
  --queue-csv "${QUEUE_CSV}"
  --hlsp-root "${HLSP_ROOT}"
  --lc-export-h5 "${LC_EXPORT_H5}"
  --out-dir "${OUT_DIR}"
  --metrics-csv "${METRICS_CSV}"
  --summary-json "${SUMMARY_JSON}"
  --branch-name "${BRANCH_NAME}"
  --aperture-tag "${APERTURE_TAG}"
  --limit "${LIMIT}"
  --workers "${WORKERS}"
  --n-periods "${N_PERIODS}"
  --n-peaks "${N_PEAKS}"
)

if [[ "${OVERWRITE}" == "1" ]]; then
  args+=(--overwrite)
fi

"${PYTHON}" "${args[@]}"

echo "[orcd-twoap] complete $(date -Is)"
