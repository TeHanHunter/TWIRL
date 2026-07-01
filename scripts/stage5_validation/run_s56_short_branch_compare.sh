#!/usr/bin/env bash
# Build and compare a short-period BLS branch against the standard S56 peak table.
#
# Intended use after the standard 20k injected peak table has finished and
# verified. This script does not modify the standard table; it audits its
# high-observability misses, reruns BLS only on those injection IDs with a
# short-period cap, then compares standard vs short-branch recovery.
set -euo pipefail

REPO_ROOT="${TWIRL_REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
cd "${REPO_ROOT}"

if [[ -z "${PYTHON:-}" ]]; then
  if [[ -x /sw/qlp-environment/.venv/bin/python ]]; then
    PYTHON=/sw/qlp-environment/.venv/bin/python
  else
    PYTHON=python3
  fi
fi

if [[ -d /pdo/app/anaconda/anaconda2-4.4.0/lib && -d /sw/python-versions/python-3.11.9/lib ]]; then
  export LD_LIBRARY_PATH="/pdo/app/anaconda/anaconda2-4.4.0/lib:/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
fi
export PYTHONPATH=src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

ROOT="${TWIRL_20K_ROOT:-reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo}"
PEAK_ROOT="${TWIRL_PEAK_ROOT:-${ROOT}/peak_training}"
DEFAULT_CHUNKED="${PEAK_ROOT}/s56_20k_injection_bls_peaks_chunked.csv"
DEFAULT_ONESHOT="${PEAK_ROOT}/s56_20k_injection_bls_peaks.csv"
STANDARD_PEAK_TABLE="${TWIRL_STANDARD_PEAK_TABLE:-}"
if [[ -z "${STANDARD_PEAK_TABLE}" ]]; then
  if [[ -f "${DEFAULT_CHUNKED}" ]]; then
    STANDARD_PEAK_TABLE="${DEFAULT_CHUNKED}"
  elif [[ -f "${DEFAULT_ONESHOT}" ]]; then
    STANDARD_PEAK_TABLE="${DEFAULT_ONESHOT}"
  else
    echo "[short-branch] missing standard peak table. Expected one of:" >&2
    echo "  ${DEFAULT_CHUNKED}" >&2
    echo "  ${DEFAULT_ONESHOT}" >&2
    exit 2
  fi
fi

INJECTION_H5="${INJECTION_H5:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_20k_predetrend_batman_periodradius_bright_balanced/injected_lightcurves.h5}"
OUT_ROOT="${TWIRL_SHORT_BRANCH_OUT_ROOT:-${ROOT}/peak_training_branch_compare/short_pmax2}"
STANDARD_AUDIT_DIR="${TWIRL_SHORT_BRANCH_STANDARD_AUDIT_DIR:-${OUT_ROOT}/standard_failure_modes}"
KEY_FILE="${TWIRL_SHORT_BRANCH_KEY_FILE:-${OUT_ROOT}/high_observability_miss_ids.txt}"
BRANCH_TABLE="${TWIRL_SHORT_BRANCH_PEAK_TABLE:-${OUT_ROOT}/short_pmax2_peaks.csv}"
BRANCH_SUMMARY="${TWIRL_SHORT_BRANCH_SUMMARY_JSON:-${BRANCH_TABLE%.*}_summary.json}"
BRANCH_GATE_DIR="${TWIRL_SHORT_BRANCH_GATE_DIR:-${OUT_ROOT}/short_pmax2_gate}"
COMPARE_OUT_DIR="${TWIRL_SHORT_BRANCH_COMPARE_OUT_DIR:-${OUT_ROOT}/compare_standard_vs_short_pmax2}"

TOP_K="${TWIRL_SHORT_BRANCH_TOP_K:-20}"
MAX_IDS="${TWIRL_SHORT_BRANCH_MAX_IDS:-0}"
MIN_STANDARD_INJECTIONS="${TWIRL_SHORT_BRANCH_MIN_STANDARD_INJECTIONS:-10000}"
MIN_STANDARD_ROWS="${TWIRL_SHORT_BRANCH_MIN_STANDARD_ROWS:-100000}"
MIN_APERTURES="${TWIRL_SHORT_BRANCH_MIN_APERTURES:-2}"
MIN_POSITIVE_PEAK_RANKS="${TWIRL_SHORT_BRANCH_MIN_POSITIVE_PEAK_RANKS:-10}"
WORKERS="${TWIRL_SHORT_BRANCH_WORKERS:-8}"
N_PERIODS="${TWIRL_SHORT_BRANCH_N_PERIODS:-100000}"
N_PEAKS="${TWIRL_SHORT_BRANCH_N_PEAKS:-20}"
P_MAX_CAP_D="${TWIRL_SHORT_BRANCH_P_MAX_CAP_D:-2.0}"
SEARCH_BRANCH="${TWIRL_SHORT_BRANCH_NAME:-short_pmax2}"
PERIOD_BIN_EDGES="${TWIRL_SHORT_BRANCH_PERIOD_BIN_EDGES:-}"
MAX_PEAKS_PER_PERIOD_BIN="${TWIRL_SHORT_BRANCH_MAX_PEAKS_PER_PERIOD_BIN:-0}"
OVERWRITE="${TWIRL_SHORT_BRANCH_OVERWRITE:-0}"
SKIP_BUILD="${TWIRL_SHORT_BRANCH_SKIP_BUILD:-0}"

log() { echo "[$(date -Is)] [short-branch] $*"; }

mkdir -p "${OUT_ROOT}" "${STANDARD_AUDIT_DIR}" "${BRANCH_GATE_DIR}" "${COMPARE_OUT_DIR}"

log "repo=${REPO_ROOT}"
log "python=${PYTHON}"
log "standard_peak_table=${STANDARD_PEAK_TABLE}"
log "injection_h5=${INJECTION_H5}"
log "out_root=${OUT_ROOT}"
log "branch=${SEARCH_BRANCH} p_max_cap_d=${P_MAX_CAP_D} n_periods=${N_PERIODS} n_peaks=${N_PEAKS}"
log "period_bin_edges=${PERIOD_BIN_EDGES:-none} max_peaks_per_period_bin=${MAX_PEAKS_PER_PERIOD_BIN}"
log "branch_table=${BRANCH_TABLE} skip_build=${SKIP_BUILD} overwrite=${OVERWRITE}"

"${PYTHON}" scripts/stage5_validation/verify_injection_peak_training_table.py \
  --peak-table "${STANDARD_PEAK_TABLE}" \
  --out-json "${OUT_ROOT}/standard_peak_table_verification.json" \
  --min-injections "${MIN_STANDARD_INJECTIONS}" \
  --min-candidate-rows "${MIN_STANDARD_ROWS}" \
  --min-apertures "${MIN_APERTURES}" \
  --min-positive-peak-ranks "${MIN_POSITIVE_PEAK_RANKS}" \
  --require-cadence-diagnostics

"${PYTHON}" scripts/stage5_validation/audit_injection_bls_failure_modes.py \
  --peak-table "${STANDARD_PEAK_TABLE}" \
  --out-dir "${STANDARD_AUDIT_DIR}" \
  --top-k "${TOP_K}"

"${PYTHON}" - "${STANDARD_AUDIT_DIR}/high_observability_misses.csv" "${KEY_FILE}" "${MAX_IDS}" <<'PY'
from pathlib import Path
import sys

import pandas as pd

source = Path(sys.argv[1])
out = Path(sys.argv[2])
max_ids = int(sys.argv[3])
if not source.exists():
    raise SystemExit(f"missing high-observability miss table: {source}")
df = pd.read_csv(source)
ids = df.get("injection_id", pd.Series(dtype=object)).dropna().astype(str).drop_duplicates().tolist()
if max_ids > 0:
    ids = ids[:max_ids]
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("\n".join(ids) + ("\n" if ids else ""), encoding="utf-8")
print(f"[short-branch] wrote {len(ids)} high-observability miss IDs to {out}")
PY

N_IDS="$(wc -l < "${KEY_FILE}" | tr -d ' ')"
if [[ "${N_IDS}" -eq 0 ]]; then
  log "no high-observability misses; skipping branch build"
  exit 0
fi

if [[ "${SKIP_BUILD}" == "1" ]]; then
  if [[ ! -f "${BRANCH_TABLE}" ]]; then
    echo "[short-branch] TWIRL_SHORT_BRANCH_SKIP_BUILD=1 but branch table is missing: ${BRANCH_TABLE}" >&2
    exit 3
  fi
  log "skip branch build; using existing ${BRANCH_TABLE}"
elif [[ -f "${BRANCH_TABLE}" && "${OVERWRITE}" != "1" ]]; then
  log "branch table exists and TWIRL_SHORT_BRANCH_OVERWRITE!=1; reusing ${BRANCH_TABLE}"
else
  BUILD_ARGS=(
    --injection-h5 "${INJECTION_H5}"
    --out-table "${BRANCH_TABLE}"
    --summary-json "${BRANCH_SUMMARY}"
    --apertures DET_FLUX_ADP_SML DET_FLUX_SML
    --limit-keys-file "${KEY_FILE}"
    --n-injections 0
    --n-peaks "${N_PEAKS}"
    --n-periods "${N_PERIODS}"
    --p-max-cap-d "${P_MAX_CAP_D}"
    --workers "${WORKERS}"
    --search-branch "${SEARCH_BRANCH}"
    --max-peaks-per-period-bin "${MAX_PEAKS_PER_PERIOD_BIN}"
    --overwrite
  )
  if [[ -n "${PERIOD_BIN_EDGES}" ]]; then
    BUILD_ARGS+=(--period-bin-edges "${PERIOD_BIN_EDGES}")
  fi
  "${PYTHON}" scripts/stage5_validation/build_injection_peak_training_table.py \
    "${BUILD_ARGS[@]}"
fi

"${PYTHON}" scripts/stage5_validation/summarize_injection_peak_gate.py \
  --peak-table "${BRANCH_TABLE}" \
  --out-dir "${BRANCH_GATE_DIR}" \
  --top-k "${TOP_K}" \
  --min-cell-count 1

"${PYTHON}" scripts/stage5_validation/compare_injection_bls_branches.py \
  --baseline-peak-table "${STANDARD_PEAK_TABLE}" \
  --branch-peak-table "${BRANCH_TABLE}" \
  --out-dir "${COMPARE_OUT_DIR}" \
  --top-k "${TOP_K}" \
  --baseline-label standard \
  --branch-label "${SEARCH_BRANCH}" \
  --branch-search-branch "${SEARCH_BRANCH}" \
  --comparison-scope intersection

log "complete"
log "branch_table=${BRANCH_TABLE}"
log "branch_gate=${BRANCH_GATE_DIR}"
log "comparison=${COMPARE_OUT_DIR}"
