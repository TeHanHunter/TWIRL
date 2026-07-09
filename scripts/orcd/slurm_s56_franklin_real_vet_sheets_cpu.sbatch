#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-franklin-real5k
#SBATCH -N 1
#SBATCH -c 48
#SBATCH --mem=160G
#SBATCH -t 10:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
QUEUE_DIR="${TWIRL_FRANKLIN_QUEUE_DIR:-reports/stage5_validation/s56_franklin_real5k_handoff}"
REAL_LC_EXPORT_H5="${TWIRL_FRANKLIN_REAL_LC_EXPORT_H5:-data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp_lc_export_pdo.h5}"
HLSP_ROOT="${TWIRL_FRANKLIN_HLSP_ROOT:-data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare}"
WORKERS="${TWIRL_FRANKLIN_WORKERS:-${SLURM_CPUS_PER_TASK:-8}}"
N_PERIODS="${TWIRL_FRANKLIN_N_PERIODS:-20000}"
N_PEAKS="${TWIRL_FRANKLIN_N_PEAKS:-10}"
N_REVIEW="${TWIRL_FRANKLIN_N_REVIEW:-5000}"
N_SOURCE_POOL="${TWIRL_FRANKLIN_N_SOURCE_POOL:-9000}"
RANDOM_STATE="${TWIRL_FRANKLIN_RANDOM_STATE:-561208}"
OVERWRITE="${TWIRL_FRANKLIN_OVERWRITE:-0}"
WRITE_PDF="${TWIRL_FRANKLIN_WRITE_PDF:-0}"

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

echo "[franklin-real5k] host=$(hostname)"
echo "[franklin-real5k] repo=${REPO}"
echo "[franklin-real5k] python=${PYTHON}"
echo "[franklin-real5k] queue_dir=${QUEUE_DIR}"
echo "[franklin-real5k] real_lc_export_h5=${REAL_LC_EXPORT_H5}"
echo "[franklin-real5k] workers=${WORKERS} n_periods=${N_PERIODS} n_review=${N_REVIEW}"

"${PYTHON}" - <<'PY'
import importlib
import sys

missing = []
for name in ("numpy", "pandas", "h5py", "astropy", "scipy", "matplotlib"):
    try:
        importlib.import_module(name)
    except Exception as exc:
        missing.append(f"{name}: {type(exc).__name__}: {exc}")
if missing:
    print("[franklin-real5k] missing required Python packages:", file=sys.stderr)
    for item in missing:
        print(f"  {item}", file=sys.stderr)
    sys.exit(3)
PY

echo "[franklin-real5k] building fresh real-only handoff $(date -Is)"
"${PYTHON}" scripts/stage5_validation/build_s56_franklin_real_handoff.py \
  --out-dir "${QUEUE_DIR}" \
  --n-review "${N_REVIEW}" \
  --n-source-pool "${N_SOURCE_POOL}" \
  --random-state "${RANDOM_STATE}"

mkdir -p "${QUEUE_DIR}/vet_sheets"
echo "[franklin-real5k] removing redundant PDF sheets $(date -Is)"
find "${QUEUE_DIR}/vet_sheets" -maxdepth 1 -type f -name "*.pdf" -delete

render_args=(
  --queue-csv "${QUEUE_DIR}/franklin_review_queue_5k_real.csv"
  --hlsp-root "${HLSP_ROOT}"
  --lc-export-h5 "${REAL_LC_EXPORT_H5}"
  --out-dir "${QUEUE_DIR}/vet_sheets"
  --metrics-csv "${QUEUE_DIR}/twirl_vet_metrics_current_adp.csv"
  --summary-json "${QUEUE_DIR}/twirl_two_aperture_vet_summary_current_adp.json"
  --branch-name current_adp
  --apertures DET_FLUX_ADP_SML DET_FLUX_ADP
  --workers "${WORKERS}"
  --n-periods "${N_PERIODS}"
  --n-peaks "${N_PEAKS}"
  --use-row-ephemeris
)
if [[ "${OVERWRITE}" == "1" ]]; then
  render_args+=(--overwrite)
fi
if [[ "${WRITE_PDF}" != "1" ]]; then
  render_args+=(--no-pdf)
fi

echo "[franklin-real5k] rendering newest-style sheets $(date -Is)"
"${PYTHON}" scripts/stage5_validation/render_two_aperture_vet_sheets.py "${render_args[@]}"

echo "[franklin-real5k] validating standalone browser package $(date -Is)"
"${PYTHON}" scripts/stage5_validation/franklin_vetting_app.py \
  --queue "${QUEUE_DIR}/franklin_review_queue_5k_real.csv" \
  --labels-out "${QUEUE_DIR}/franklin_labels_vetted.csv" \
  --sheet-root "${QUEUE_DIR}/vet_sheets" \
  --check-only

echo "[franklin-real5k] complete $(date -Is)"
