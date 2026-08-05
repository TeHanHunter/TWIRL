#!/usr/bin/env bash
# Run teacher-v1 real-candidate enrichment from one A2v1 sector on PDO.
set -euo pipefail

REPO="${TWIRL_PDO_TEACHER_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
A2V1_ROOT="${TWIRL_A2V1_ROOT:-/pdo/users/tehan/tglc-gpu-production-A2v1}"
SECTOR="${TWIRL_TEACHER_SECTOR:-56}"
SECTOR_TAG="s$(printf '%04d' "${SECTOR}")"
SECTOR_SHORT="s${SECTOR}"
if [[ "${SECTOR}" != "56" ]]; then
  echo "Blocked: this teacher-search runner implements the S56-only candidate-quality contract; sector ${SECTOR} is not supported." >&2
  exit 64
fi
case "${SECTOR}" in
  56) DEFAULT_ORBITS="119,120" ;;
  *) DEFAULT_ORBITS="" ;;
esac
ORBIT_CSV="${TWIRL_TEACHER_ORBITS:-${DEFAULT_ORBITS}}"
[[ -n "${ORBIT_CSV}" ]] || {
  echo "Set TWIRL_TEACHER_ORBITS for sector ${SECTOR}" >&2
  exit 2
}
IFS=',' read -r -a ORBITS <<< "${ORBIT_CSV}"
HLSP_ROOT="${TWIRL_A2V1_HLSP_ROOT:-${A2V1_ROOT}/hlsp_${SECTOR_TAG}_A2v1}"
OUT_ROOT="${TWIRL_TEACHER_SEARCH_ROOT:-/pdo/users/tehan/twirl_stage5/${SECTOR_SHORT}_A2v1_teacher_search_v2}"
BASE_PYTHON="${TWIRL_PDO_BASE_PYTHON:-/sw/qlp-environment/.venv/bin/python}"
TORCH_PYTHON="${TWIRL_PDO_TORCH_PYTHON:-/pdo/users/tehan/envs/twirl-teacher-pdo-v2/bin/python}"
COMPACT_LC="${OUT_ROOT}/inputs/${SECTOR_SHORT}_A2v1_adp_pair.h5"
CADENCE_REFERENCE="${TWIRL_A2V1_CADENCE_REFERENCE:-${OUT_ROOT}/inputs/${SECTOR_SHORT}_A2v1_cadence_reference.csv}"
CADENCE_REFERENCE_MANIFEST="${TWIRL_A2V1_CADENCE_REFERENCE_MANIFEST:-${OUT_ROOT}/inputs/${SECTOR_SHORT}_A2v1_cadence_reference_manifest.json}"
TIER1_GATE_JSON="${TWIRL_A2V1_TIER1_GATE_JSON:-${OUT_ROOT}/inputs/${SECTOR_SHORT}_A2v1_tier1_gate.json}"
TIER1_TARGET_ELIGIBILITY="${TWIRL_A2V1_TIER1_TARGET_ELIGIBILITY:-${OUT_ROOT}/inputs/${SECTOR_SHORT}_A2v1_tier1_target_eligibility.csv}"
BLS_ROOT="${OUT_ROOT}/adp_bls"
CANDIDATE_ROOT="${OUT_ROOT}/candidates"
RAW_SOURCE="${OUT_ROOT}/inputs/${SECTOR_SHORT}_A2v1_tglc_raw_sources.h5"
NATIVE_H5="${OUT_ROOT}/inputs/${SECTOR_SHORT}_A2v1_adp_raw_pair_v2.h5"
CHECKPOINT_ROOT="${TWIRL_PDO_TEACHER_CHECKPOINT_ROOT:-/pdo/users/tehan/twirl_stage5/s56_A2v1_teacher_search_v2/checkpoints/shape_plus_periodogram_bls}"
CHECKPOINT_MANIFEST="${TWIRL_PDO_TEACHER_CHECKPOINT_MANIFEST:-$(dirname "${CHECKPOINT_ROOT}")/selected_checkpoint_manifest.json}"
SCORE_ROOT="${OUT_ROOT}/teacher_scores"
TRANSFER_LABELS="${TWIRL_PDO_TRANSFER_LABELS:-${OUT_ROOT}/inputs/human_vetting_training_table_adjudicated.csv}"
TRANSFER_ROOT="${OUT_ROOT}/transfer_gate"
PRODUCT_QA_SUMMARY="${TWIRL_A2V1_PHOTOMETRIC_QA_SUMMARY:-${OUT_ROOT}/inputs/${SECTOR_SHORT}_A2v1_photometric_qa_gate.json}"
PRODUCT_SCHEMA_SUMMARY="${TWIRL_A2V1_SCHEMA_SUMMARY:-${A2V1_ROOT}/twirl_logs/${SECTOR_SHORT}_A2v1_validation_full_schema.json}"
QA_ROOT="${OUT_ROOT}/photometric_qa"
REFERENCE_RAW_ROOT="${TWIRL_REFERENCE_RAW_ROOT:-/pdo/users/tehan/tglc-gpu-production}"
S56_REFERENCE_QA="${TWIRL_S56_REFERENCE_QA:-/pdo/users/tehan/twirl_stage5/s56_A2v1_teacher_search_v2/inputs/s56_A2v1_photometric_qa_gate.json}"
BATCH_ROOT="${OUT_ROOT}/active_learning/batch01"
EXCLUSION_ROOT="${OUT_ROOT}/inputs/exclusions"
SHEET_ROOT="${BATCH_ROOT}/twirl_vet_sheets"
BLS_WORKERS="${TWIRL_PDO_BLS_WORKERS:-32}"
METADATA_WORKERS="${TWIRL_PDO_METADATA_WORKERS:-32}"
N_PERIODS="${TWIRL_PDO_BLS_N_PERIODS:-50000}"
NATIVE_SHARDS="${TWIRL_PDO_NATIVE_SHARDS:-16}"
NATIVE_PARALLEL="${TWIRL_PDO_NATIVE_PARALLEL:-4}"
QA_WORKERS="${TWIRL_PDO_QA_WORKERS:-8}"
QA_SAMPLE_SIZE="${TWIRL_PDO_QA_SAMPLE_SIZE:-256}"
COMMAND="${1:-status}"
(( NATIVE_PARALLEL >= 1 )) || {
  echo "TWIRL_PDO_NATIVE_PARALLEL must be at least 1" >&2
  exit 2
}

export LD_LIBRARY_PATH="/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

mkdir -p "${OUT_ROOT}/inputs" "${BLS_ROOT}" "${CANDIDATE_ROOT}" "${SCORE_ROOT}" "${TRANSFER_ROOT}" "${BATCH_ROOT}" "${EXCLUSION_ROOT}" "${OUT_ROOT}/logs"
if [[ "${COMMAND}" != "status" && "${COMMAND}" != "probe" ]]; then
  exec 9>"${OUT_ROOT}/.teacher-search.lock"
  if ! flock -n 9; then
    echo "Another ${SECTOR_SHORT} A2v1 teacher-search command already holds the shared output lock" >&2
    exit 9
  fi
fi
cd "${REPO}"

peak_table() {
  if [[ -s "${BLS_ROOT}/real_adp_bls_peaks.parquet" ]]; then
    printf '%s\n' "${BLS_ROOT}/real_adp_bls_peaks.parquet"
  elif [[ -s "${BLS_ROOT}/real_adp_bls_peaks.csv" ]]; then
    printf '%s\n' "${BLS_ROOT}/real_adp_bls_peaks.csv"
  else
    return 1
  fi
}

candidate_table() {
  if [[ -s "${CANDIDATE_ROOT}/teacher_v1_scoring_candidates.parquet" ]]; then
    printf '%s\n' "${CANDIDATE_ROOT}/teacher_v1_scoring_candidates.parquet"
  elif [[ -s "${CANDIDATE_ROOT}/teacher_v1_scoring_candidates.csv" ]]; then
    printf '%s\n' "${CANDIDATE_ROOT}/teacher_v1_scoring_candidates.csv"
  else
    return 1
  fi
}

score_table() {
  if [[ -s "${SCORE_ROOT}/teacher_v1_real_candidate_scores.parquet" ]]; then
    printf '%s\n' "${SCORE_ROOT}/teacher_v1_real_candidate_scores.parquet"
  elif [[ -s "${SCORE_ROOT}/teacher_v1_real_candidate_scores.csv" ]]; then
    printf '%s\n' "${SCORE_ROOT}/teacher_v1_real_candidate_scores.csv"
  else
    return 1
  fi
}

run_export() {
  if [[ -s "${COMPACT_LC}" ]]; then
    echo "[teacher-search] reuse compact A2v1 export: ${COMPACT_LC}"
    return
  fi
  "${BASE_PYTHON}" scripts/stage3_injections/export_s56_lc_training_set.py \
    --hlsp-root "${HLSP_ROOT}" \
    --out-h5 "${COMPACT_LC}" \
    --sector "${SECTOR}" \
    --columns DET_FLUX_ADP_SML DET_FLUX_ADP \
    --progress-every 500
}

run_bls() {
  [[ -s "${CADENCE_REFERENCE}" ]] || {
    echo "Missing authoritative cadence reference: ${CADENCE_REFERENCE}" >&2
    exit 13
  }
  [[ -s "${CADENCE_REFERENCE_MANIFEST}" ]] || {
    echo "Missing cadence-reference manifest: ${CADENCE_REFERENCE_MANIFEST}" >&2
    exit 14
  }
  "${BASE_PYTHON}" scripts/stage5_validation/build_s56_adp_real_bls_peaks.py \
    --compact-lc "${COMPACT_LC}" \
    --cadence-reference "${CADENCE_REFERENCE}" \
    --cadence-reference-manifest "${CADENCE_REFERENCE_MANIFEST}" \
    --out-dir "${BLS_ROOT}" \
    --workers "${BLS_WORKERS}" \
    --n-periods "${N_PERIODS}" \
    --n-peaks 10 \
    --source-product-tag A2v1 \
    --progress-every 100 \
    --resume
}

require_tier1_evidence() {
  [[ -s "${TIER1_GATE_JSON}" ]] || {
    echo "Missing Tier-1 enrichment gate: ${TIER1_GATE_JSON}" >&2
    exit 15
  }
  [[ -s "${TIER1_TARGET_ELIGIBILITY}" ]] || {
    echo "Missing Tier-1 target eligibility: ${TIER1_TARGET_ELIGIBILITY}" >&2
    exit 16
  }
  "${BASE_PYTHON}" - "${TIER1_GATE_JSON}" <<'PY'
import json, sys
summary = json.load(open(sys.argv[1]))
if not (
    summary.get("status") == "pass"
    and summary.get("passed") is True
    and summary.get("enrichment_ready") is True
):
    raise SystemExit(f"Tier-1 enrichment gate is not passed/ready: {summary}")
PY
}

run_candidates() {
  require_tier1_evidence
  local existing=""
  if existing="$(candidate_table 2>/dev/null)" && [[ -s "${CANDIDATE_ROOT}/summary.json" ]]; then
    if "${BASE_PYTHON}" - \
      "${existing}" "${CANDIDATE_ROOT}/summary.json" \
      "$(peak_table)" "${BLS_ROOT}/summary.json" "${COMPACT_LC}" \
      "${TIER1_TARGET_ELIGIBILITY}" "${TIER1_GATE_JSON}" \
      "${CADENCE_REFERENCE}" "${CADENCE_REFERENCE_MANIFEST}" <<'PY'
import hashlib, json, sys

candidate, summary_path, peaks, peaks_summary, compact, eligibility, gate, cadence, cadence_manifest = sys.argv[1:]
summary = json.load(open(summary_path))

def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

valid = (
    summary.get("provenance_contract_version")
        == "s56_a2v1_teacher_scoring_candidates_v1"
    and summary.get("candidate_table_sha256") == sha256(candidate)
    and summary.get("adp_peaks_sha256") == sha256(peaks)
    and summary.get("adp_peaks_summary_sha256") == sha256(peaks_summary)
    and summary.get("compact_lc_sha256") == sha256(compact)
    and summary.get("tier1_target_eligibility_sha256") == sha256(eligibility)
    and summary.get("tier1_gate_json_sha256") == sha256(gate)
    and summary.get("cadence_reference_sha256") == sha256(cadence)
    and summary.get("cadence_reference_manifest_sha256") == sha256(cadence_manifest)
    and summary.get("bls_evidence", {}).get("status") == "pass"
    and summary.get("tier1_gate", {}).get("enrichment_ready") is True
    and summary.get("bls_search_contract_version")
        == summary.get("bls_evidence", {}).get("bls_search_contract_version")
    and summary.get("bls_config_sha256")
        == summary.get("bls_evidence", {}).get("bls_config_sha256")
)
raise SystemExit(0 if valid else 1)
PY
    then
      echo "[teacher-search] reuse verified scoring candidates: ${existing}"
      return
    fi
    echo "[teacher-search] cached scoring candidates failed provenance checks; rebuilding" >&2
  fi
  "${BASE_PYTHON}" scripts/stage5_validation/build_s56_teacher_v1_scoring_table.py \
    --adp-peaks "$(peak_table)" \
    --adp-peaks-summary "${BLS_ROOT}/summary.json" \
    --compact-lc "${COMPACT_LC}" \
    --cadence-reference-table "${CADENCE_REFERENCE}" \
    --cadence-reference-manifest "${CADENCE_REFERENCE_MANIFEST}" \
    --tier1-target-eligibility "${TIER1_TARGET_ELIGIBILITY}" \
    --tier1-gate-json "${TIER1_GATE_JSON}" \
    --out-dir "${CANDIDATE_ROOT}" \
    --sector "${SECTOR}" \
    --small-peaks-per-tic 3 \
    --workers "${METADATA_WORKERS}" \
    --progress-every 250
}

verify_native_cache() {
  local path="$1"
  local shard_index="${2:--1}"
  [[ -s "${path}" ]] || return 1
  "${BASE_PYTHON}" - \
    "${path}" "$(candidate_table)" "${CANDIDATE_ROOT}/summary.json" \
    "${shard_index}" "${NATIVE_SHARDS}" <<'PY'
import json, sys
from pathlib import Path
from twirl.vetting.harmonic_inputs import verify_native_candidate_binding

native, candidates, summary, shard_index, n_shards = sys.argv[1:]
shard_index = int(shard_index)
kwargs = {}
if shard_index >= 0:
    kwargs.update(
        expected_shard_index=shard_index,
        expected_n_shards=int(n_shards),
    )
verification = verify_native_candidate_binding(
    Path(native),
    candidate_table=Path(candidates),
    candidate_summary=Path(summary),
    require_periodograms=True,
    expected_periodogram_n=4096,
    **kwargs,
)
if not verification["passed"]:
    print(json.dumps(verification, indent=2, sort_keys=True), file=sys.stderr)
raise SystemExit(0 if verification["passed"] else 1)
PY
}

run_qa() {
  if [[ -s "${PRODUCT_QA_SUMMARY}" ]]; then
    if "${BASE_PYTHON}" - "${PRODUCT_QA_SUMMARY}" <<'PY'
import json, sys
raise SystemExit(0 if json.load(open(sys.argv[1])).get("passed", False) else 1)
PY
    then
      echo "[teacher-search] reuse passed photometric QA: ${PRODUCT_QA_SUMMARY}"
      return
    fi
  fi
  [[ -s "${PRODUCT_SCHEMA_SUMMARY}" ]] || {
    echo "Missing ${SECTOR_SHORT} A2v1 schema/coverage summary: ${PRODUCT_SCHEMA_SUMMARY}" >&2
    exit 11
  }
  local reference_args=()
  if [[ "${SECTOR}" != 56 ]]; then
    [[ -s "${S56_REFERENCE_QA}" ]] || {
      echo "Missing passed S56 reference QA for ${SECTOR_SHORT}: ${S56_REFERENCE_QA}" >&2
      exit 12
    }
    reference_args+=(--reference-qa-summary "${S56_REFERENCE_QA}")
  fi
  "${BASE_PYTHON}" scripts/stage1_lightcurves/audit_a2v1_photometric_qa.py \
    --sector "${SECTOR}" \
    --orbits "${ORBITS[@]}" \
    --compact-lc "${COMPACT_LC}" \
    --schema-summary "${PRODUCT_SCHEMA_SUMMARY}" \
    --bls-peaks "$(peak_table)" \
    --a2v1-root "${A2V1_ROOT}" \
    --reference-raw-root "${REFERENCE_RAW_ROOT}" \
    --out-dir "${QA_ROOT}" \
    --gate-json "${PRODUCT_QA_SUMMARY}" \
    --sample-size "${QA_SAMPLE_SIZE}" \
    --workers "${QA_WORKERS}" \
    "${reference_args[@]}"
}

run_raw_export() {
  if [[ -s "${RAW_SOURCE}" ]]; then
    echo "[teacher-search] reuse compact raw/error export: ${RAW_SOURCE}"
    return
  fi
  "${BASE_PYTHON}" scripts/stage5_validation/export_s56_harmonic_raw_sources.py \
    --training-table "$(candidate_table)" \
    --raw-root "${A2V1_ROOT}" \
    --out-h5 "${RAW_SOURCE}" \
    --compact-adp-h5 "${COMPACT_LC}" \
    --orbits "${ORBIT_CSV}"
}

run_native() {
  if verify_native_cache "${NATIVE_H5}"; then
    echo "[teacher-search] reuse verified native model input: ${NATIVE_H5}"
    return
  fi
  if [[ -s "${NATIVE_H5}" ]]; then
    echo "[teacher-search] cached native model input failed verification; rebuilding" >&2
  fi
  local shard_root="${OUT_ROOT}/native_shards"
  mkdir -p "${shard_root}"
  local pids=()
  local failed=0
  for ((shard=0; shard<NATIVE_SHARDS; shard++)); do
    local path="${shard_root}/native_${shard}.h5"
    if verify_native_cache "${path}" "${shard}"; then
      echo "[teacher-search] reuse verified native shard ${shard}: ${path}"
      continue
    fi
    if [[ -s "${path}" ]]; then
      echo "[teacher-search] cached native shard ${shard} failed verification; rebuilding" >&2
    fi
    "${BASE_PYTHON}" scripts/stage5_validation/build_s56_harmonic_native_inputs.py \
      --training-table "$(candidate_table)" \
      --training-summary "${CANDIDATE_ROOT}/summary.json" \
      --raw-source-h5 "${RAW_SOURCE}" \
      --compact-adp-h5 "${COMPACT_LC}" \
      --cadence-reference-table "${CADENCE_REFERENCE}" \
      --cadence-reference-manifest "${CADENCE_REFERENCE_MANIFEST}" \
      --sector "${SECTOR}" \
      --out-h5 "${path}" \
      --repo-root "${REPO}" \
      --n-periods 4096 \
      --shard-index "${shard}" \
      --n-shards "${NATIVE_SHARDS}" \
      >"${OUT_ROOT}/logs/native_${shard}.out" \
      2>"${OUT_ROOT}/logs/native_${shard}.err" &
    pids+=("$!")
    if (( ${#pids[@]} >= NATIVE_PARALLEL )); then
      for pid in "${pids[@]}"; do
        if ! wait "${pid}"; then
          failed=1
        fi
      done
      pids=()
      if [[ "${failed}" != 0 ]]; then
        break
      fi
    fi
  done
  # PDO still runs a Bash version where expanding an empty array under
  # `set -u` raises "unbound variable". The final batch can be empty when all
  # native shards were reused, so guard the expansion explicitly.
  if (( ${#pids[@]} > 0 )); then
    for pid in "${pids[@]}"; do
      if ! wait "${pid}"; then
        failed=1
      fi
    done
  fi
  if [[ "${failed}" != 0 ]]; then
    echo "[teacher-search] at least one native-input shard failed" >&2
    exit 5
  fi
  local shards=("${shard_root}"/native_*.h5)
  "${BASE_PYTHON}" scripts/stage5_validation/merge_s56_harmonic_native_shards.py \
    --training-table "$(candidate_table)" \
    --training-summary "${CANDIDATE_ROOT}/summary.json" \
    --shards "${shards[@]}" \
    --out-h5 "${NATIVE_H5}"
}

checkpoint_args() {
  local paths=()
  for fold in 0 1 2 3 4; do
    local path="${CHECKPOINT_ROOT}/fold_${fold}/teacher.pt"
    [[ -s "${path}" ]] || { echo "Missing checkpoint: ${path}" >&2; return 1; }
    paths+=("${path}")
  done
  printf '%s\n' "${paths[@]}"
}

verify_checkpoint_ensemble() {
  [[ -s "${CHECKPOINT_MANIFEST}" ]] || {
    echo "Missing native-v2 checkpoint manifest: ${CHECKPOINT_MANIFEST}" >&2
    return 1
  }
  "${TORCH_PYTHON}" \
    scripts/stage5_validation/verify_s56_harmonic_teacher_checkpoint_manifest.py \
    --manifest "${CHECKPOINT_MANIFEST}" \
    --checkpoint-root "${CHECKPOINT_ROOT}" \
    --expected-profile "shape_plus_periodogram_bls"
}

require_product_qa() {
  [[ -s "${PRODUCT_QA_SUMMARY}" ]] || {
    echo "Missing ${SECTOR_SHORT} A2v1 photometric-QA gate: ${PRODUCT_QA_SUMMARY}" >&2
    exit 10
  }
  "${BASE_PYTHON}" - "${PRODUCT_QA_SUMMARY}" <<'PY'
import json, sys
summary = json.load(open(sys.argv[1]))
if not summary.get("passed", False):
    raise SystemExit(f"A2v1 photometric-QA gate failed: {summary}")
PY
}

run_score() {
  run_qa
  [[ -x "${TORCH_PYTHON}" ]] || { echo "Missing PDO Torch Python: ${TORCH_PYTHON}" >&2; exit 6; }
  verify_checkpoint_ensemble || exit 7
  mapfile -t checkpoints < <(checkpoint_args)
  [[ "${#checkpoints[@]}" == 5 ]] || exit 7
  if score_table >/dev/null 2>&1 && [[ -s "${SCORE_ROOT}/summary.json" ]]; then
    if "${BASE_PYTHON}" - \
      "${SCORE_ROOT}/summary.json" "$(candidate_table)" \
      "${CANDIDATE_ROOT}/summary.json" "${NATIVE_H5}" \
      "${checkpoints[@]}" <<'PY'
import json, sys
from pathlib import Path
from twirl.vetting.harmonic_inference import verify_teacher_score_cache

summary, candidates, candidate_summary, native, *checkpoints = sys.argv[1:]
verification = verify_teacher_score_cache(
    summary_path=Path(summary),
    candidates_path=Path(candidates),
    candidate_summary_path=Path(candidate_summary),
    native_h5=Path(native),
    checkpoint_paths=[Path(path) for path in checkpoints],
)
if not verification["passed"]:
    print(json.dumps(verification, indent=2, sort_keys=True), file=sys.stderr)
raise SystemExit(0 if verification["passed"] else 1)
PY
    then
      echo "[teacher-search] reuse hash-verified teacher scores: $(score_table)"
      return
    fi
    echo "[teacher-search] cached teacher scores failed provenance checks; rebuilding" >&2
  fi
  CUDA_VISIBLE_DEVICES="${TWIRL_PDO_TEACHER_GPU:-4}" "${TORCH_PYTHON}" \
    scripts/stage5_validation/score_s56_teacher_v1_candidates.py \
    --candidates "$(candidate_table)" \
    --candidate-summary "${CANDIDATE_ROOT}/summary.json" \
    --native-h5 "${NATIVE_H5}" \
    --checkpoints "${checkpoints[@]}" \
    --out-dir "${SCORE_ROOT}" \
    --batch-size "${TWIRL_PDO_TEACHER_BATCH_SIZE:-32}" \
    --workers "${TWIRL_PDO_TEACHER_WORKERS:-4}"
}

run_transfer() {
  [[ -s "${TRANSFER_LABELS}" ]] || {
    echo "Missing adjudicated transfer labels: ${TRANSFER_LABELS}" >&2
    exit 8
  }
  "${BASE_PYTHON}" scripts/stage5_validation/audit_s56_a2v1_teacher_transfer.py \
    --labels "${TRANSFER_LABELS}" \
    --adp-peaks "$(peak_table)" \
    --teacher-scores "$(score_table)" \
    --out-dir "${TRANSFER_ROOT}"
}

run_batch() {
  require_product_qa
  [[ -s "${TRANSFER_ROOT}/summary.json" ]] || run_transfer
  "${BASE_PYTHON}" - "${TRANSFER_ROOT}/summary.json" <<'PY'
import json, sys
summary = json.load(open(sys.argv[1]))
gate = summary.get("scored_transfer", {}).get("transfer_gate", {})
if not gate.get("passed", False):
    raise SystemExit(f"A2v1 transfer gate failed: {gate}")
PY
  local exclusions=()
  while IFS= read -r path; do
    exclusions+=("${path}")
  done < <(find "${EXCLUSION_ROOT}" -maxdepth 1 -type f -name '*.csv' -print | sort)
  "${BASE_PYTHON}" scripts/stage5_validation/build_s56_teacher_active_learning_batch.py \
    --scores "$(score_table)" \
    --out-dir "${BATCH_ROOT}" \
    --batch-index 1 \
    --exclude "${TRANSFER_LABELS}" "${exclusions[@]}"
}

run_render() {
  local queue="${BATCH_ROOT}/review_queue_batch01_1k.csv"
  [[ -s "${queue}" ]] || run_batch
  "${BASE_PYTHON}" scripts/stage5_validation/render_two_aperture_vet_sheets.py \
    --queue-csv "${queue}" \
    --lc-export-h5 "${COMPACT_LC}" \
    --out-dir "${SHEET_ROOT}" \
    --metrics-csv "${BATCH_ROOT}/twirl_vet_metrics.csv" \
    --summary-json "${BATCH_ROOT}/twirl_two_aperture_vet_summary.json" \
    --branch-name current_adp \
    --apertures DET_FLUX_ADP_SML DET_FLUX_ADP \
    --workers "${TWIRL_PDO_RENDER_WORKERS:-32}" \
    --n-periods "${TWIRL_PDO_RENDER_N_PERIODS:-20000}" \
    --n-peaks 10 \
    --use-row-ephemeris \
    --no-pdf \
    --harmonic-factors 0.25 0.5 1 2 4
}

show_status() {
  echo "repo=${REPO}"
  echo "sector=${SECTOR}"
  echo "orbits=${ORBIT_CSV}"
  echo "a2v1=${A2V1_ROOT}"
  echo "out=${OUT_ROOT}"
  for path in "${COMPACT_LC}" "${CADENCE_REFERENCE}" "${CADENCE_REFERENCE_MANIFEST}" "${TIER1_GATE_JSON}" "${TIER1_TARGET_ELIGIBILITY}" "${RAW_SOURCE}" "${NATIVE_H5}"; do
    if [[ -s "${path}" ]]; then ls -lh "${path}"; else echo "MISSING ${path}"; fi
  done
  peak_table 2>/dev/null || echo "MISSING peak table"
  candidate_table 2>/dev/null || echo "MISSING candidate table"
  score_table 2>/dev/null || echo "MISSING score table"
  if [[ -s "${PRODUCT_QA_SUMMARY}" ]]; then
    echo "photometric_qa=${PRODUCT_QA_SUMMARY}"
  else
    echo "MISSING photometric QA gate"
  fi
  checkpoint_args 2>/dev/null || echo "MISSING one or more checkpoints"
}

case "${COMMAND}" in
  probe) "${BASE_PYTHON}" -c 'import astropy,h5py,numpy,pandas; print("PDO base environment OK")'; show_status ;;
  export) run_export ;;
  bls) run_export; run_bls ;;
  candidates) run_export; run_bls; run_qa; run_candidates ;;
  qa) run_export; run_bls; run_qa ;;
  prep) run_export; run_bls; run_qa; run_candidates; run_raw_export; run_native ;;
  score) run_export; run_bls; run_qa; run_candidates; run_raw_export; run_native; run_score ;;
  transfer) run_export; run_bls; run_qa; run_candidates; run_raw_export; run_native; run_score; run_transfer ;;
  batch) run_export; run_bls; run_qa; run_candidates; run_raw_export; run_native; run_score; run_transfer; run_batch ;;
  render) run_export; run_bls; run_qa; run_candidates; run_raw_export; run_native; run_score; run_transfer; run_batch; run_render ;;
  all) run_export; run_bls; run_qa; run_candidates; run_raw_export; run_native; run_score; run_transfer; run_batch; run_render ;;
  status) show_status ;;
  *) echo "Usage: $0 {probe|export|bls|candidates|qa|prep|score|transfer|batch|render|all|status}" >&2; exit 2 ;;
esac
