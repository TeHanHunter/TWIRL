#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-s56-allhost-peak
#SBATCH -N 1
#SBATCH -c 48
#SBATCH --mem=192G
#SBATCH -t 2-00:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
INJECTION_ROOT="${TWIRL_ALLHOST_INJECTION_ROOT:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_allhost_predetrend_batman_periodradius_grid_sharded}"
OUT_ROOT="${TWIRL_ALLHOST_PEAK_OUT_ROOT:-reports/stage5_validation/s56_allhost_predetrend_batman_periodradius_grid_sharded_pdo/peak_training_orcd}"
FINAL_TABLE="${TWIRL_ALLHOST_FINAL_TABLE:-${OUT_ROOT}/s56_allhost_injection_bls_peaks_orcd.csv}"
CHUNK_JOBS="${TWIRL_ALLHOST_PEAK_CHUNK_JOBS:-4}"
WORKERS="${TWIRL_PEAK_WORKERS:-8}"
N_PEAKS="${TWIRL_BLS_N_PEAKS:-20}"
N_PERIODS="${TWIRL_BLS_N_PERIODS:-200000}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs
cd "${REPO}"

export PYTHONPATH=src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[orcd-allhost-peak] host=$(hostname)"
echo "[orcd-allhost-peak] repo=${REPO}"
echo "[orcd-allhost-peak] python=${PYTHON}"
echo "[orcd-allhost-peak] injection_root=${INJECTION_ROOT}"
echo "[orcd-allhost-peak] out_root=${OUT_ROOT}"
echo "[orcd-allhost-peak] final_table=${FINAL_TABLE}"
echo "[orcd-allhost-peak] chunk_jobs=${CHUNK_JOBS} workers=${WORKERS} n_periods=${N_PERIODS}"

"${PYTHON}" - <<'PY'
import importlib
import sys

required = ("numpy", "pandas", "h5py", "astropy")
missing = []
for name in required:
    try:
        importlib.import_module(name)
    except Exception as exc:  # pragma: no cover - cluster preflight
        missing.append(f"{name}: {type(exc).__name__}: {exc}")
if missing:
    print("[orcd-allhost-peak] missing required Python packages:", file=sys.stderr)
    for item in missing:
        print(f"  {item}", file=sys.stderr)
    sys.exit(3)
PY

PYTHON_BIN="${PYTHON}" \
INJECTION_ROOT="${INJECTION_ROOT}" \
OUT_ROOT="${OUT_ROOT}" \
FINAL_TABLE="${FINAL_TABLE}" \
TWIRL_ALLHOST_PEAK_CHUNK_JOBS="${CHUNK_JOBS}" \
TWIRL_PEAK_WORKERS="${WORKERS}" \
TWIRL_BLS_N_PEAKS="${N_PEAKS}" \
TWIRL_BLS_N_PERIODS="${N_PERIODS}" \
  bash scripts/stage5_validation/run_s56_allhost_peak_training_sharded_pdo.sh

echo "[orcd-allhost-peak] complete $(date -Is)"
