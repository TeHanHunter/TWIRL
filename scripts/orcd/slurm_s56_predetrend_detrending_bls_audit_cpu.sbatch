#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-predetrend-bls
#SBATCH -N 1
#SBATCH -c 48
#SBATCH --mem=192G
#SBATCH -t 1-00:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
WORKERS="${TWIRL_PREDET_BLS_WORKERS:-${SLURM_CPUS_PER_TASK:-16}}"
N_INJECTIONS="${TWIRL_PREDET_BLS_N_INJECTIONS:-3000}"
N_PERIODS="${TWIRL_PREDET_BLS_N_PERIODS:-50000}"
N_PEAKS="${TWIRL_PREDET_BLS_N_PEAKS:-20}"
INJECTION_H5="${TWIRL_PREDET_BLS_INJECTION_H5:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_20k_predetrend_batman_periodradius_bright_balanced/injected_lightcurves.h5}"
INJECTION_MANIFEST="${TWIRL_PREDET_BLS_INJECTION_MANIFEST:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_20k_predetrend_batman_periodradius_bright_balanced/injection_manifest.csv}"
OUT_DIR="${TWIRL_PREDET_BLS_OUT_DIR:-reports/stage5_validation/s56_predetrend_detrending_bls_audit_orcd}"
METHODS="${TWIRL_PREDET_BLS_METHODS:-current_adp03q,twirl_fs_v2_adp015q,spline_q02_gap02,spline_q025_gap02,spline_q05_gap02,spline_q08_gap02,current_uniform08_gap05}"
RAW_APERTURES="${TWIRL_PREDET_BLS_RAW_APERTURES:-Small,Primary}"

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

echo "[orcd-predetrend-bls] host=$(hostname)"
echo "[orcd-predetrend-bls] repo=${REPO}"
echo "[orcd-predetrend-bls] python=${PYTHON}"
echo "[orcd-predetrend-bls] workers=${WORKERS} n_injections=${N_INJECTIONS} n_periods=${N_PERIODS}"
echo "[orcd-predetrend-bls] injection_h5=${INJECTION_H5}"
echo "[orcd-predetrend-bls] injection_manifest=${INJECTION_MANIFEST}"
echo "[orcd-predetrend-bls] out_dir=${OUT_DIR}"
echo "[orcd-predetrend-bls] methods=${METHODS}"
echo "[orcd-predetrend-bls] raw_apertures=${RAW_APERTURES}"

"${PYTHON}" - <<'PY'
import importlib
import sys

required = ("numpy", "pandas", "h5py", "astropy", "scipy")
missing = []
for name in required:
    try:
        importlib.import_module(name)
    except Exception as exc:  # pragma: no cover - cluster preflight
        missing.append(f"{name}: {type(exc).__name__}: {exc}")
if missing:
    print("[orcd-predetrend-bls] missing required Python packages:", file=sys.stderr)
    for item in missing:
        print(f"  {item}", file=sys.stderr)
    sys.exit(3)
PY

args=(
  scripts/stage5_validation/run_s56_predetrend_detrending_bls_audit.py
  --injection-h5 "${INJECTION_H5}"
  --injection-manifest "${INJECTION_MANIFEST}"
  --out-dir "${OUT_DIR}"
  --n-injections "${N_INJECTIONS}"
  --workers "${WORKERS}"
  --n-periods "${N_PERIODS}"
  --n-peaks "${N_PEAKS}"
)

IFS=',' read -r -a method_array <<<"${METHODS}"
for method in "${method_array[@]}"; do
  [[ -n "${method}" ]] && args+=(--method "${method}")
done

IFS=',' read -r -a aperture_array <<<"${RAW_APERTURES}"
for aperture in "${aperture_array[@]}"; do
  [[ -n "${aperture}" ]] && args+=(--raw-aperture "${aperture}")
done

"${PYTHON}" "${args[@]}"

echo "[orcd-predetrend-bls] complete $(date -Is)"
