#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-s56-peak
#SBATCH -N 1
#SBATCH -c 48
#SBATCH --mem=192G
#SBATCH -t 2-00:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
WORKERS="${TWIRL_PEAK_WORKERS:-${SLURM_CPUS_PER_TASK:-16}}"
N_INJECTIONS="${TWIRL_N_INJECTIONS:-0}"
N_PEAKS="${TWIRL_BLS_N_PEAKS:-20}"
N_PERIODS="${TWIRL_BLS_N_PERIODS:-200000}"
INJECTION_H5="${TWIRL_INJECTION_H5:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_20k_predetrend_batman_periodradius_bright_balanced/injected_lightcurves.h5}"
OUT_ROOT="${TWIRL_PEAK_OUT_ROOT:-reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training_orcd_full}"
FINAL_TABLE="${TWIRL_FINAL_TABLE:-${OUT_ROOT}/s56_20k_injection_bls_peaks_orcd_full.csv}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs
cd "${REPO}"

export PYTHONPATH=src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[orcd-peak] host=$(hostname)"
echo "[orcd-peak] repo=${REPO}"
echo "[orcd-peak] python=${PYTHON}"
echo "[orcd-peak] workers=${WORKERS} n_injections=${N_INJECTIONS} n_periods=${N_PERIODS}"
echo "[orcd-peak] injection_h5=${INJECTION_H5}"
echo "[orcd-peak] out_root=${OUT_ROOT}"
echo "[orcd-peak] final_table=${FINAL_TABLE}"

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
    print("[orcd-peak] missing required Python packages:", file=sys.stderr)
    for item in missing:
        print(f"  {item}", file=sys.stderr)
    sys.exit(3)
PY

export TWIRL_REPO_ROOT="${REPO}"
export PYTHON
export INJECTION_H5
export OUT_ROOT
export FINAL_TABLE
export TWIRL_N_INJECTIONS="${N_INJECTIONS}"
export TWIRL_PEAK_WORKERS="${WORKERS}"
export TWIRL_BLS_N_PEAKS="${N_PEAKS}"
export TWIRL_BLS_N_PERIODS="${N_PERIODS}"

scripts/stage5_validation/run_s56_20k_peak_training_chunked.sh

echo "[orcd-peak] complete $(date -Is)"
