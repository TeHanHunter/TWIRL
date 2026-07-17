#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-s56-ranker
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-python3}"
PEAK_TABLE="${TWIRL_PEAK_TABLE:-reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training/s56_20k_injection_bls_peaks_orcd_chunked.csv}"
OUT_DIR="${TWIRL_RANKER_OUT_DIR:-reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_ranker_orcd}"
GATE_OUT_DIR="${TWIRL_GATE_OUT_DIR:-reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training_gate_orcd}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs
cd "${REPO}"

export PYTHONPATH=src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[orcd-ranker] host=$(hostname)"
echo "[orcd-ranker] repo=${REPO}"
echo "[orcd-ranker] python=${PYTHON}"
echo "[orcd-ranker] peak_table=${PEAK_TABLE}"
echo "[orcd-ranker] out_dir=${OUT_DIR}"
echo "[orcd-ranker] gate_out_dir=${GATE_OUT_DIR}"

"${PYTHON}" - <<'PY'
import importlib
import sys

required = ("numpy", "pandas", "scipy")
missing = []
for name in required:
    try:
        importlib.import_module(name)
    except Exception as exc:  # pragma: no cover - cluster preflight
        missing.append(f"{name}: {type(exc).__name__}: {exc}")
if missing:
    print("[orcd-ranker] missing required Python packages:", file=sys.stderr)
    for item in missing:
        print(f"  {item}", file=sys.stderr)
    sys.exit(3)
PY

mkdir -p "${OUT_DIR}"
"${PYTHON}" scripts/stage5_validation/verify_injection_peak_training_table.py \
  --peak-table "${PEAK_TABLE}" \
  --out-json "${OUT_DIR}/peak_table_verification.json" \
  --min-injections 10000 \
  --min-candidate-rows 100000 \
  --min-apertures 2 \
  --min-positive-peak-ranks 10 \
  --require-cadence-diagnostics

"${PYTHON}" scripts/stage5_validation/train_injection_peak_ranker.py \
  --peak-table "${PEAK_TABLE}" \
  --out-dir "${OUT_DIR}"

"${PYTHON}" scripts/stage5_validation/summarize_injection_peak_gate.py \
  --peak-table "${PEAK_TABLE}" \
  --out-dir "${GATE_OUT_DIR}" \
  --top-k 20 \
  --ranker-summary-json "${OUT_DIR}/summary.json"

echo "[orcd-ranker] complete $(date -Is)"
