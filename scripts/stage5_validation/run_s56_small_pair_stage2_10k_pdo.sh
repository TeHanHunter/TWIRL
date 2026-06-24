#!/usr/bin/env bash
# Build the S56 small-aperture real-candidate table and mixed 10k review queue.
#
# This is the current science path after the 2026-06-23 detrending audit:
#   1. search the S56 TWIRL-FS v2 compare product with DET_FLUX_ADP_SML+DET_FLUX_SML,
#   2. run the heuristic WD vetter and centroid enrichment on that matching search table,
#   3. merge 9k real candidates with the verified 1k raw-flux pre-detrend BATMAN injections,
#   4. render WD-tuned LEO reports for the mixed blinded queue.
set -euo pipefail

cd /pdo/users/tehan/TWIRL

export LD_LIBRARY_PATH="/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
QLP_LD_LIBRARY_PATH="/pdo/app/anaconda/anaconda2-4.4.0/lib:/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="/pdo/users/tehan/LEO-Vetter-twirl:src"
export MPLCONFIGDIR=/pdo/users/tehan/.cache/matplotlib
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

PY_STAGE2="${PY_STAGE2:-/pdo/users/tehan/twirl-stage2-venv/bin/python}"
PY_QLP="${PY_QLP:-/sw/qlp-environment/.venv/bin/python}"
HLSP_ROOT="${HLSP_ROOT:-/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl_fs_v2_compare}"
STAGE2_OUT="${STAGE2_OUT:-data_local/stage2/bls_s56_twirlfs_v2_small_pair_200k}"
SECTOR_DIR="${STAGE2_OUT}/sector_0056"
TARGET_LIST="${TARGET_LIST:-${STAGE2_OUT}/s56_twirlfs_v2_compare_target_paths.txt}"
TIC_LIST="${TIC_LIST:-${HLSP_ROOT}/tic_ids_s0056_lightcurves_orbits119_120.txt}"
INJECTION_H5="${INJECTION_H5:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_1k_predetrend_batman_depthgrid_adp_compare/injected_lightcurves.h5}"
SOURCE_BLS="${SOURCE_BLS:-reports/stage5_validation/s56_1k_predetrend_small_pair_200k_chunked_pdo/small_pair_200k/injection_bls_recoveries.csv}"
REVIEW_DIR="${REVIEW_DIR:-reports/stage5_validation/s56_10k_predetrend_small_pair_200k_review_queue_pdo}"
N_REAL="${N_REAL:-9000}"
N_INJECTIONS="${N_INJECTIONS:-1000}"
BLS_WORKERS="${BLS_WORKERS:-64}"
CENTROID_WORKERS="${CENTROID_WORKERS:-32}"
LEO_WORKERS="${LEO_WORKERS:-8}"
LEO_TIMEOUT_S="${LEO_TIMEOUT_S:-300}"
RANDOM_STATE="${RANDOM_STATE:-5610}"
MAX_LEO_REPORTS="${MAX_LEO_REPORTS:-0}"
REBUILD_STAGE2="${REBUILD_STAGE2:-0}"
REBUILD_VET="${REBUILD_VET:-0}"
REBUILD_CENTROID="${REBUILD_CENTROID:-0}"
REUSE_INJECTION_BLS="${REUSE_INJECTION_BLS:-1}"
SKIP_LEO="${SKIP_LEO:-0}"
SKIP_QUEUE="${SKIP_QUEUE:-0}"

log() { echo "[$(date -Is)] [s56-smallpair-10k] $*"; }

mkdir -p logs "${REVIEW_DIR}"

log "stage2_python=${PY_STAGE2}"
log "qlp_python=${PY_QLP}"
log "hlsp_root=${HLSP_ROOT}"
log "stage2_out=${STAGE2_OUT}"
log "review_dir=${REVIEW_DIR}"

if [[ "${REBUILD_STAGE2}" == "1" || ! -s "${SECTOR_DIR}/consolidated.parquet" ]]; then
  if [[ "${REBUILD_STAGE2}" == "1" || ! -s "${TARGET_LIST}" ]]; then
    log "building explicit HLSP target list"
    mkdir -p "$(dirname "${TARGET_LIST}")"
    "${PY_STAGE2}" - "${HLSP_ROOT}" "${TIC_LIST}" "${TARGET_LIST}" <<'PY'
from pathlib import Path
import sys

hlsp_root = Path(sys.argv[1])
tic_list = Path(sys.argv[2])
target_list = Path(sys.argv[3])

if not tic_list.exists():
    raise SystemExit(f"missing TIC list: {tic_list}")

paths = []
for raw in tic_list.read_text().splitlines():
    value = raw.strip()
    if not value or value == "tic_id":
        continue
    tic = int(value.split(",", 1)[0])
    tic16 = f"{tic:016d}"
    paths.append(
        hlsp_root
        / tic16[0:4]
        / tic16[4:8]
        / tic16[8:12]
        / tic16[12:16]
        / f"hlsp_twirlfs_tess_ffi_s0056-{tic16}_tess_v01_llc.fits"
    )

target_list.write_text("\n".join(str(path) for path in paths) + "\n")
print(f"wrote {len(paths)} paths to {target_list}")
PY
    log "target_list=${TARGET_LIST} n=$(wc -l < "${TARGET_LIST}")"
  else
    log "reusing target_list=${TARGET_LIST} n=$(wc -l < "${TARGET_LIST}")"
  fi
  log "running small-pair S56 BLS"
  "${PY_STAGE2}" -m twirl.search.sector_run \
    --sector 56 \
    --hlsp-root "${HLSP_ROOT}" \
    --target-list "${TARGET_LIST}" \
    --out-dir "${STAGE2_OUT}" \
    --workers "${BLS_WORKERS}" \
    --apertures DET_FLUX_ADP_SML,DET_FLUX_SML \
    --n-periods 200000 \
    --n-peaks 10
else
  log "reusing ${SECTOR_DIR}/consolidated.parquet"
fi

if [[ "${REBUILD_VET}" == "1" || ! -s "${SECTOR_DIR}/vetted_per_tic.parquet" ]]; then
  log "running heuristic vetter"
  "${PY_STAGE2}" scripts/stage5_validation/vet_s56.py \
    --sector 56 \
    --out-dir "${STAGE2_OUT}" \
    --grid-max-duration-min 30
else
  log "reusing ${SECTOR_DIR}/vetted_per_tic.parquet"
fi

if [[ "${REBUILD_CENTROID}" == "1" || ! -s "${SECTOR_DIR}/vetted_per_tic_centroid.csv" ]]; then
  log "running centroid enrichment"
  TWIRL_HLSP_ROOT="${HLSP_ROOT}" \
  TWIRL_INPUT_PARQUET="${SECTOR_DIR}/vetted_per_tic.parquet" \
  TWIRL_OUTPUT_PARQUET="${SECTOR_DIR}/vetted_per_tic_centroid.parquet" \
  TWIRL_WORKERS="${CENTROID_WORKERS}" \
    "${PY_STAGE2}" scripts/stage5_validation/centroid_check.py
else
  log "reusing ${SECTOR_DIR}/vetted_per_tic_centroid.csv"
fi

if [[ "${SKIP_QUEUE}" == "1" ]]; then
  log "SKIP_QUEUE=1; stopping after Stage 2/vetter/centroid products"
  exit 0
fi

if [[ ! -s "${INJECTION_H5}" ]]; then
  echo "[s56-smallpair-10k] missing injection HDF5: ${INJECTION_H5}" >&2
  exit 2
fi
if [[ ! -s "${SOURCE_BLS}" ]]; then
  echo "[s56-smallpair-10k] missing source injected BLS table: ${SOURCE_BLS}" >&2
  exit 2
fi

if [[ "${REUSE_INJECTION_BLS}" == "1" ]]; then
  cp "${SOURCE_BLS}" "${REVIEW_DIR}/injection_bls_recoveries.csv"
fi

queue_args=(
  --real-candidates "${SECTOR_DIR}/vetted_per_tic_centroid.csv"
  --hlsp-root "${HLSP_ROOT}"
  --injection-h5 "${INJECTION_H5}"
  --out-dir "${REVIEW_DIR}"
  --n-real "${N_REAL}"
  --n-injections "${N_INJECTIONS}"
  --real-selection stratified_blind
  --blind-review-metadata
  --shuffle-review-rows
  --apertures DET_FLUX_ADP_SML DET_FLUX_SML
  --workers 1
  --n-periods 200000
  --leo-workers "${LEO_WORKERS}"
  --leo-timeout-s "${LEO_TIMEOUT_S}"
  --max-leo-reports "${MAX_LEO_REPORTS}"
  --random-state "${RANDOM_STATE}"
  --overwrite
)
if [[ "${REUSE_INJECTION_BLS}" == "1" ]]; then
  queue_args+=(--reuse-injection-bls)
fi
if [[ "${SKIP_LEO}" == "1" ]]; then
  queue_args+=(--skip-leo)
fi

log "building mixed review queue"
LD_LIBRARY_PATH="${QLP_LD_LIBRARY_PATH}" "${PY_QLP}" scripts/stage5_validation/build_s56_pretriage_review_queue.py "${queue_args[@]}"

if [[ "${SKIP_LEO}" != "1" && "${MAX_LEO_REPORTS}" == "0" ]]; then
  log "verifying mixed queue"
  LD_LIBRARY_PATH="${QLP_LD_LIBRARY_PATH}" "${PY_QLP}" scripts/stage5_validation/verify_s56_pretriage_queue.py \
    --queue "${REVIEW_DIR}/review_queue.csv" \
    --reports-dir "${REVIEW_DIR}/vet_reports" \
    --summary-json "${REVIEW_DIR}/summary.json" \
    --min-rows "$((N_REAL + N_INJECTIONS))" \
    --expect-real "${N_REAL}" \
    --expect-injected "${N_INJECTIONS}" \
    --min-reports "$((N_REAL + N_INJECTIONS))" \
    --apertures DET_FLUX_ADP_SML DET_FLUX_SML \
    --out-json "${REVIEW_DIR}/verification.json"
else
  log "skipping strict verifier because SKIP_LEO=${SKIP_LEO} MAX_LEO_REPORTS=${MAX_LEO_REPORTS}"
fi

log "complete"
log "stage2_dir=${SECTOR_DIR}"
log "review_queue=${REVIEW_DIR}/review_queue.csv"
