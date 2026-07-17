#!/usr/bin/env bash
# Run robust BLS peak-table training over the sharded all-host S56 injections.
set -euo pipefail

log() {
  printf '[s56-allhost-peak] %s\n' "$*" >&2
}

REPO_ROOT="${TWIRL_REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
cd "${REPO_ROOT}"

PYTHON_BIN="${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}"
export LD_LIBRARY_PATH="/sw/python-versions/python-3.11.9/lib:/pdo/app/anaconda/anaconda2-4.4.0/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${PYTHONPATH:-src}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

INJECTION_ROOT="${INJECTION_ROOT:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_allhost_predetrend_batman_periodradius_grid_sharded}"
OUT_ROOT="${OUT_ROOT:-reports/stage5_validation/s56_allhost_predetrend_batman_periodradius_grid_sharded_pdo/peak_training}"
PEAK_CHUNK_ROOT="${PEAK_CHUNK_ROOT:-${OUT_ROOT}/chunks}"
FINAL_TABLE="${FINAL_TABLE:-${OUT_ROOT}/s56_allhost_injection_bls_peaks.csv}"
VERIFY_JSON="${VERIFY_JSON:-${FINAL_TABLE%.*}_verification.json}"
TABLE_NAME="${TABLE_NAME:-injection_bls_peaks.csv}"

CHUNK_JOBS="${TWIRL_ALLHOST_PEAK_CHUNK_JOBS:-1}"
WORKERS="${TWIRL_PEAK_WORKERS:-8}"
N_PERIODS="${TWIRL_BLS_N_PERIODS:-200000}"
N_PEAKS="${TWIRL_BLS_N_PEAKS:-20}"
N_INJECTIONS_PER_SHARD="${TWIRL_N_INJECTIONS_PER_SHARD:-0}"
MAX_SHARDS="${TWIRL_ALLHOST_PEAK_MAX_SHARDS:-0}"
MERGE_WHEN_COMPLETE="${TWIRL_ALLHOST_PEAK_MERGE:-1}"
REQUIRE_ALL_SHARDS_FOR_MERGE="${TWIRL_ALLHOST_REQUIRE_ALL_SHARDS_FOR_MERGE:-1}"
MIN_APERTURES="${TWIRL_PEAK_MIN_APERTURES:-2}"
MIN_POSITIVE_PEAK_RANKS="${TWIRL_PEAK_MIN_POSITIVE_RANKS:-10}"
REQUIRE_SIGNAL_PEAKS="${TWIRL_PEAK_REQUIRE_SIGNAL_PEAKS:-1}"

mkdir -p "${PEAK_CHUNK_ROOT}" "${OUT_ROOT}"

log "repo=${REPO_ROOT}"
log "python=${PYTHON_BIN}"
log "injection_root=${INJECTION_ROOT}"
log "out_root=${OUT_ROOT}"
log "chunk_jobs=${CHUNK_JOBS} workers=${WORKERS} n_periods=${N_PERIODS} n_peaks=${N_PEAKS}"

chunk_list="${OUT_ROOT}/ready_injection_shards.txt"
find "${INJECTION_ROOT}/chunks" -mindepth 2 -maxdepth 2 -name injected_lightcurves.h5 2>/dev/null |
  sort > "${chunk_list}"
if [[ "${MAX_SHARDS}" -gt 0 ]]; then
  tmp_list="${OUT_ROOT}/ready_injection_shards.limited.txt"
  head -n "${MAX_SHARDS}" "${chunk_list}" > "${tmp_list}"
  mv "${tmp_list}" "${chunk_list}"
fi

n_ready="$(wc -l < "${chunk_list}" | tr -d ' ')"
log "ready injection shards=${n_ready}"
if [[ "${n_ready}" == "0" ]]; then
  log "no completed injection shards found yet"
  exit 0
fi
mkdir -p "${OUT_ROOT}/chunk_ids"
{
  echo "n_injections=0"
  echo "chunk_size=sharded"
  echo "n_chunks=${n_ready}"
  echo "source_list=${chunk_list}"
} > "${OUT_ROOT}/chunk_ids/manifest.txt"

run_one_peak_chunk() {
  local h5_path="$1"
  local chunk_name chunk_out summary_path run_log
  chunk_name="$(basename "$(dirname "${h5_path}")")"
  chunk_out="${PEAK_CHUNK_ROOT}/${chunk_name}"
  summary_path="${chunk_out}/injection_bls_peaks_summary.json"
  run_log="${chunk_out}/run.log"
  if [[ -s "${summary_path}" && -s "${chunk_out}/${TABLE_NAME}" ]]; then
    log "skip ${chunk_name}: completed"
    return 0
  fi
  mkdir -p "${chunk_out}"
  log "run ${chunk_name}; h5=${h5_path}; log=${run_log}"
  if ! "${PYTHON_BIN}" scripts/stage5_validation/build_injection_peak_training_table.py \
    --injection-h5 "${h5_path}" \
    --out-table "${chunk_out}/${TABLE_NAME}" \
    --apertures DET_FLUX_ADP_SML DET_FLUX_SML \
    --n-injections "${N_INJECTIONS_PER_SHARD}" \
    --n-peaks "${N_PEAKS}" \
    --n-periods "${N_PERIODS}" \
    --workers "${WORKERS}" \
    --overwrite >"${run_log}" 2>&1; then
    log "failed ${chunk_name}; tail ${run_log}"
    tail -n 40 "${run_log}" || true
    return 1
  fi
  log "done ${chunk_name}"
}
export -f run_one_peak_chunk log
export PYTHON_BIN PEAK_CHUNK_ROOT TABLE_NAME N_INJECTIONS_PER_SHARD N_PEAKS N_PERIODS WORKERS

if [[ "${CHUNK_JOBS}" -le 1 ]]; then
  while IFS= read -r h5_path; do
    run_one_peak_chunk "${h5_path}"
  done < "${chunk_list}"
else
  xargs -n 1 -P "${CHUNK_JOBS}" bash -c 'run_one_peak_chunk "$1"' _ < "${chunk_list}"
fi

n_peak_chunks="$(find "${PEAK_CHUNK_ROOT}" -mindepth 2 -maxdepth 2 -name "${TABLE_NAME}" 2>/dev/null | wc -l | tr -d ' ')"
log "peak chunks completed=${n_peak_chunks}/${n_ready}"

if [[ "${MERGE_WHEN_COMPLETE}" != "1" ]]; then
  log "merge disabled"
  exit 0
fi
if [[ "${REQUIRE_ALL_SHARDS_FOR_MERGE}" == "1" && "${n_peak_chunks}" != "${n_ready}" ]]; then
  log "merge skipped until all ready injection shards have peak tables"
  exit 0
fi

EXPECT_N=0
if [[ -s "${INJECTION_ROOT}/summary.json" ]]; then
  EXPECT_N="$("${PYTHON_BIN}" -c 'import json,sys; print(json.load(open(sys.argv[1])).get("n_injections", 0))' "${INJECTION_ROOT}/summary.json")"
elif [[ "${N_INJECTIONS_PER_SHARD}" -gt 0 ]]; then
  EXPECT_N=$(( n_peak_chunks * N_INJECTIONS_PER_SHARD ))
fi

log "merging peak chunks expect_n=${EXPECT_N}"
"${PYTHON_BIN}" scripts/stage5_validation/merge_injection_peak_training_chunks.py \
  --chunk-root "${PEAK_CHUNK_ROOT}" \
  --out-table "${FINAL_TABLE}" \
  --table-name "${TABLE_NAME}" \
  --expect-injections "${EXPECT_N}"

MIN_INJECTIONS="${TWIRL_PEAK_MIN_INJECTIONS:-${EXPECT_N}}"
if [[ "${MIN_INJECTIONS}" == "0" ]]; then
  MIN_INJECTIONS=1
fi
MIN_CANDIDATE_ROWS="${TWIRL_PEAK_MIN_CANDIDATE_ROWS:-$(( MIN_INJECTIONS * 5 ))}"
VERIFY_ARGS=(
  scripts/stage5_validation/verify_injection_peak_training_table.py
  --peak-table "${FINAL_TABLE}"
  --out-json "${VERIFY_JSON}"
  --min-injections "${MIN_INJECTIONS}"
  --min-candidate-rows "${MIN_CANDIDATE_ROWS}"
  --min-apertures "${MIN_APERTURES}"
  --min-positive-peak-ranks "${MIN_POSITIVE_PEAK_RANKS}"
  --require-cadence-diagnostics
)
if [[ "${REQUIRE_SIGNAL_PEAKS}" != "1" ]]; then
  VERIFY_ARGS+=(--no-require-signal-peaks)
fi
log "verifying merged table"
"${PYTHON_BIN}" "${VERIFY_ARGS[@]}"

log "complete"
log "final_table=${FINAL_TABLE}"
log "verification_json=${VERIFY_JSON}"
