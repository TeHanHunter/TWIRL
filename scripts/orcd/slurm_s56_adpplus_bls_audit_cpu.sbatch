#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-adpplus-audit
#SBATCH -N 1
#SBATCH -c 48
#SBATCH --mem=192G
#SBATCH -t 2-00:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
WORKERS="${TWIRL_ADPPLUS_WORKERS:-${SLURM_CPUS_PER_TASK:-16}}"
N_INJECTED="${TWIRL_ADPPLUS_N_INJECTED:-3000}"
N_REAL="${TWIRL_ADPPLUS_N_REAL:-2000}"
N_PERIODS="${TWIRL_ADPPLUS_N_PERIODS:-50000}"
N_PEAKS="${TWIRL_ADPPLUS_N_PEAKS:-20}"
OUT_DIR="${TWIRL_ADPPLUS_OUT_DIR:-reports/stage5_validation/s56_adpplus_bls_audit_orcd}"
INJECTION_CSV="${TWIRL_ADPPLUS_INJECTION_CSV:-reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/small_pair_200k/review_queue.csv}"
INJECTION_H5="${TWIRL_ADPPLUS_INJECTION_H5:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_20k_predetrend_batman_periodradius_bright_balanced/injected_lightcurves.h5}"
REAL_CANDIDATES="${TWIRL_ADPPLUS_REAL_CANDIDATES:-data_local/stage2/bls_first_pass_v2/sector_0056/vetted_per_tic_centroid.csv}"
EXPORT_H5="${TWIRL_ADPPLUS_EXPORT_H5:-data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp_lc_export_pdo.h5}"
HLSP_ROOT="${TWIRL_ADPPLUS_HLSP_ROOT:-data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs
cd "${REPO}"

export PYTHONPATH=src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[orcd-adpplus] host=$(hostname)"
echo "[orcd-adpplus] repo=${REPO}"
echo "[orcd-adpplus] python=${PYTHON}"
echo "[orcd-adpplus] workers=${WORKERS} injected=${N_INJECTED} real=${N_REAL} n_periods=${N_PERIODS}"
echo "[orcd-adpplus] export_h5=${EXPORT_H5}"
echo "[orcd-adpplus] out_dir=${OUT_DIR}"

"${PYTHON}" - <<'PY'
import importlib
import sys

missing = []
for name in ("numpy", "pandas", "h5py", "astropy", "scipy"):
    try:
        importlib.import_module(name)
    except Exception as exc:
        missing.append(f"{name}: {type(exc).__name__}: {exc}")
if missing:
    print("[orcd-adpplus] missing required Python packages:", file=sys.stderr)
    for item in missing:
        print(f"  {item}", file=sys.stderr)
    sys.exit(3)
PY

"${PYTHON}" scripts/stage5_validation/run_s56_adpplus_bls_audit.py \
  --injection-csv "${INJECTION_CSV}" \
  --injection-h5 "${INJECTION_H5}" \
  --real-candidates "${REAL_CANDIDATES}" \
  --hlsp-root "${HLSP_ROOT}" \
  --export-h5 "${EXPORT_H5}" \
  --out-dir "${OUT_DIR}" \
  --n-injected "${N_INJECTED}" \
  --n-real "${N_REAL}" \
  --n-periods "${N_PERIODS}" \
  --n-peaks "${N_PEAKS}" \
  --workers "${WORKERS}"

"${PYTHON}" scripts/stage5_validation/summarize_s56_adpplus_bls_audit.py \
  --audit-dir "${OUT_DIR}"

echo "[orcd-adpplus] complete $(date -Is)"
