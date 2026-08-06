#!/usr/bin/env bash
# Fail-closed PDO stages for sealed prospective S63 Teacher-v3 enrichment.
set -euo pipefail

AUTHORIZATION_CONTRACT="s63_teacher_v3_prospective_v1"
PDO_USER_ROOT="/pdo/users/tehan"
SECTOR=63
ORBIT_ONE=133
ORBIT_TWO=134
REPO="${TWIRL_ROOT:-/pdo/users/tehan/TWIRL}"
QLP_PYTHON="${TWIRL_QLP_PYTHON:-${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}}"
STAGE2_PYTHON="${TWIRL_STAGE2_PYTHON:-/pdo/users/tehan/twirl-stage2-venv/bin/python}"
A2V1_ROOT="${TWIRL_A2V1_ROOT:-/pdo/users/tehan/tglc-gpu-production-A2v1}"
HLSP_ROOT="${TWIRL_A2V1_HLSP_ROOT:-${A2V1_ROOT}/hlsp_s0063_A2v1}"
OUT_ROOT="${TWIRL_S63_TEACHER_V3_ROOT:-/pdo/users/tehan/twirl_stage5/s63_teacher_v3_prospective_v1}"
SOURCE_VALIDATION="${TWIRL_S63_ACCEPTED_VALIDATION:-${A2V1_ROOT}/twirl_logs/s63_A2v1_validation_full_schema.json}"
VALIDATION="${TWIRL_S63_ATTESTED_VALIDATION:-${OUT_ROOT}/stage1_receipt/s63_A2v1_validation_full_schema.attested.json}"
CADENCE_DIR="${TWIRL_S63_CADENCE_DIR:-${OUT_ROOT}/cadence}"
CADENCE_TABLE="${TWIRL_S63_CADENCE_TABLE:-${CADENCE_DIR}/cadence_reference.csv}"
CADENCE_MANIFEST="${TWIRL_S63_CADENCE_MANIFEST:-${CADENCE_DIR}/cadence_reference.manifest.json}"
COMPACT_H5="${TWIRL_S63_COMPACT_H5:-${OUT_ROOT}/compact/s63_A2v1_adp_pair.h5}"
COMPACT_MANIFEST="${COMPACT_H5%.h5}.manifest.json"
PLAN="${TWIRL_S63_PROSPECTIVE_PLAN:-${REPO}/reports/stage5_validation/teacher_v3_s63_prospective_v1/preregistered/prospective_plan_v1.json}"
SELECTION_POLICY="${TWIRL_S63_SELECTION_POLICY:-${REPO}/reports/stage5_validation/teacher_v3_s63_prospective_v1/preregistered/selection_policy_v1.json}"
RESERVED_TICS="${TWIRL_S63_RESERVED_TICS:-${REPO}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp/preregistered/s63_reserved_tics.txt}"
TEACHER_V3_CORPUS="${TWIRL_TEACHER_V3_CORPUS:-${REPO}/reports/stage5_validation/teacher_v3_s56_s62_a2v1_current_adp/frozen/observation_morphology_corpus.csv}"
MODEL_READY_ALLOWLIST="${TWIRL_S63_MODEL_READY_ALLOWLIST:-${OUT_ROOT}/prospective/model_ready_tics.txt}"
MODEL_READY_SUMMARY="${TWIRL_S63_MODEL_READY_SUMMARY:-${OUT_ROOT}/prospective/model_ready_tics.summary.json}"
COHORT_DIR="${TWIRL_S63_COHORT_DIR:-${OUT_ROOT}/prospective/cohorts}"
PRIMARY_COHORT="${TWIRL_S63_PRIMARY_COHORT:-${COHORT_DIR}/primary_teacher_v3_disjoint_tics.txt}"
REPEATED_COHORT="${TWIRL_S63_REPEATED_COHORT:-${COHORT_DIR}/repeated_host_tics.txt}"
RESERVED_NOT_READY="${TWIRL_S63_RESERVED_NOT_READY:-${COHORT_DIR}/reserved_not_model_ready_tics.txt}"
COHORT_SUMMARY="${TWIRL_S63_COHORT_SUMMARY:-${COHORT_DIR}/cohort_summary.json}"
BLS_SHARD_DIR="${TWIRL_S63_BLS_SHARD_DIR:-${OUT_ROOT}/bls/shards}"
BLS_MERGED="${TWIRL_S63_BLS_MERGED:-${OUT_ROOT}/bls/real_adp_bls_peaks.parquet}"
CANDIDATES="${TWIRL_S63_CANDIDATES:-${OUT_ROOT}/candidates/s63_teacher_v3_rank1_candidates.parquet}"
CANDIDATE_SUMMARY="${CANDIDATES%.parquet}.summary.json"
RAW_SOURCE_H5="${TWIRL_S63_RAW_SOURCE_H5:-${OUT_ROOT}/raw/s63_teacher_v3_tglc_raw_sources.h5}"
NATIVE_H5="${TWIRL_S63_NATIVE_H5:-${OUT_ROOT}/native/s63_teacher_v3_native_v2.h5}"
NATIVE_SUMMARY="${TWIRL_S63_NATIVE_SUMMARY:-${NATIVE_H5%.h5}.summary.json}"
LAUNCH_MANIFEST="${TWIRL_S63_LAUNCH_MANIFEST:-${OUT_ROOT}/launch/preprocessing_launch_manifest.json}"
COMMAND="${1:-plan}"

require_authorization() {
  if [[ "${TWIRL_TEACHER_V3_PROSPECTIVE_CONTRACT:-}" != "${AUTHORIZATION_CONTRACT}" ]]; then
    echo "S63 preprocessing requires TWIRL_TEACHER_V3_PROSPECTIVE_CONTRACT=${AUTHORIZATION_CONTRACT}." >&2
    exit 64
  fi
}

require_pdo_output_root() {
  local canonical_user_root canonical_path output_var
  canonical_user_root="$(realpath -e -- "${PDO_USER_ROOT}")"
  REPO="$(realpath -e -- "${REPO}")"
  if [[ "${REPO}" != "${canonical_user_root}"/* ]]; then
    echo "Pinned PDO checkout must remain under ${canonical_user_root}: ${REPO}" >&2
    exit 2
  fi
  local output_vars=(
    OUT_ROOT VALIDATION CADENCE_DIR CADENCE_TABLE CADENCE_MANIFEST
    COMPACT_H5 COMPACT_MANIFEST MODEL_READY_ALLOWLIST MODEL_READY_SUMMARY
    COHORT_DIR PRIMARY_COHORT REPEATED_COHORT RESERVED_NOT_READY COHORT_SUMMARY
    BLS_SHARD_DIR BLS_MERGED CANDIDATES CANDIDATE_SUMMARY RAW_SOURCE_H5
    NATIVE_H5 NATIVE_SUMMARY LAUNCH_MANIFEST
  )
  for output_var in "${output_vars[@]}"; do
    canonical_path="$(realpath -m -- "${!output_var}")"
    if [[ "${canonical_path}" != "${canonical_user_root}"/* \
       || "${canonical_path}" == "${REPO}" \
       || "${canonical_path}" == "${REPO}"/* ]]; then
      echo "PDO output ${output_var} must be canonical, under ${canonical_user_root}, and outside the checkout: ${canonical_path}" >&2
      exit 2
    fi
    printf -v "${output_var}" '%s' "${canonical_path}"
  done
}

require_clean_pinned_checkout() {
  if [[ -z "${TWIRL_EXPECTED_GIT_SHA:-}" ]]; then
    echo "Set TWIRL_EXPECTED_GIT_SHA to the deployed full Git SHA." >&2
    exit 2
  fi
  "${REPO}/scripts/assert_clean_checkout.sh" "${REPO}"
}

preflight_stage1() {
  (
    export PYTHONNOUSERSITE=1
    unset LD_LIBRARY_PATH
    "${STAGE2_PYTHON}" \
      "${REPO}/scripts/stage5_validation/build_s63_teacher_v3_launch_manifest.py" \
      --accepted-validation "${SOURCE_VALIDATION}" \
      --preflight-only
  )
}

run_stage1_receipt() {
  prepare_mutating_stage stage1_receipt stage2
  require_file "${SOURCE_VALIDATION}" "immutable accepted Stage-1 validation"
  require_absent "${VALIDATION}" "attested Stage-1 receipt copy"
  local pending_validation="${VALIDATION}.pending"
  require_absent "${pending_validation}" "pending Stage-1 receipt attestation"
  mkdir -p "$(dirname "${VALIDATION}")"
  cp -p -- "${SOURCE_VALIDATION}" "${pending_validation}"
  local source_sha
  source_sha="$(sha256_file "${SOURCE_VALIDATION}")"
  "${STAGE2_PYTHON}" scripts/stage5_validation/attest_s63_teacher_v3_json.py \
    --json "${pending_validation}" \
    --expected-sha256 "${source_sha}" \
    --producer-git-sha "${TWIRL_EXPECTED_GIT_SHA}" \
    --source-json "${SOURCE_VALIDATION}" \
    --expected-source-sha256 "${source_sha}"
  mv -- "${pending_validation}" "${VALIDATION}"
}

require_file() {
  local path=$1
  local label=$2
  if [[ ! -s "${path}" ]]; then
    echo "Missing ${label}: ${path}" >&2
    exit 3
  fi
}

require_absent() {
  local path=$1
  local label=$2
  if [[ -e "${path}" ]]; then
    echo "Refusing to overwrite existing ${label}: ${path}" >&2
    exit 4
  fi
}

sha256_file() {
  local path=$1
  local digest remainder
  read -r digest remainder < <(sha256sum "${path}")
  if [[ ! "${digest}" =~ ^[0-9a-f]{64}$ ]]; then
    echo "Could not hash ${path}" >&2
    exit 3
  fi
  printf '%s\n' "${digest}"
}

claim_stage() {
  local name=$1
  local claim_root="${OUT_ROOT}/claims"
  ACTIVE_CLAIM="${claim_root}/${name}.claim"
  mkdir -p "${claim_root}"
  if ! mkdir "${ACTIVE_CLAIM}" 2>/dev/null; then
    echo "Stage claim already exists; audit before retry: ${ACTIVE_CLAIM}" >&2
    exit 9
  fi
  trap 'if [[ -n "${ACTIVE_CLAIM:-}" ]]; then rmdir "${ACTIVE_CLAIM}" 2>/dev/null || true; fi' EXIT
}

prepare_mutating_stage() {
  local name=$1
  local runtime=$2
  require_authorization
  require_pdo_output_root
  require_clean_pinned_checkout
  cd "${REPO}"
  export PYTHONPATH="${REPO}/src"
  export OMP_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1
  export MKL_NUM_THREADS=1
  export VECLIB_MAXIMUM_THREADS=1
  export NUMEXPR_NUM_THREADS=1
  export HDF5_USE_FILE_LOCKING=FALSE
  export PYTHONNOUSERSITE=1
  if [[ "${runtime}" == "qlp" ]]; then
    export LD_LIBRARY_PATH="/sw/python-versions/python-3.11.9/lib:/pdo/app/anaconda/anaconda2-4.4.0/lib:${LD_LIBRARY_PATH:-}"
    "${QLP_PYTHON}" -c 'import astropy, h5py, numpy, pandas; print("S63 QLP Python import preflight OK")'
  elif [[ "${runtime}" == "stage2" ]]; then
    unset LD_LIBRARY_PATH
    "${STAGE2_PYTHON}" -c 'import astropy, h5py, numpy, pandas, pyarrow; print("S63 Stage2/Parquet import preflight OK")'
  else
    echo "Unknown S63 runtime class: ${runtime}" >&2
    exit 2
  fi
  preflight_stage1 >/dev/null
  claim_stage "${name}"
}

run_cadence() {
  prepare_mutating_stage cadence qlp
  require_file "${VALIDATION}" "attested Stage-1 receipt copy"
  require_absent "${CADENCE_TABLE}" "cadence table"
  require_absent "${CADENCE_MANIFEST}" "cadence manifest"
  TWIRL_TEACHER_V3_PROSPECTIVE_CONTRACT="${AUTHORIZATION_CONTRACT}" \
    TWIRL_ROOT="${REPO}" \
    TWIRL_EXPECTED_GIT_SHA="${TWIRL_EXPECTED_GIT_SHA}" \
    PYTHON_BIN="${QLP_PYTHON}" \
    bash scripts/stage1_lightcurves/run_a2v1_cadence_reference_pdo.sh \
      "${SECTOR}" "${ORBIT_ONE}" "${ORBIT_TWO}" "${CADENCE_DIR}"
}

run_compact() {
  prepare_mutating_stage compact qlp
  require_file "${VALIDATION}" "attested Stage-1 receipt copy"
  require_absent "${COMPACT_H5}" "compact HDF5"
  require_absent "${COMPACT_MANIFEST}" "compact manifest"
  local pending_h5="${COMPACT_H5%.h5}.pending.h5"
  local pending_manifest="${COMPACT_H5%.h5}.pending.manifest.json"
  require_absent "${pending_h5}" "pending compact HDF5"
  require_absent "${pending_manifest}" "pending compact manifest"
  mkdir -p "$(dirname "${COMPACT_H5}")"
  "${QLP_PYTHON}" scripts/stage3_injections/export_s56_lc_training_set.py \
    --hlsp-root "${HLSP_ROOT}" \
    --out-h5 "${pending_h5}" \
    --sector "${SECTOR}" \
    --columns DET_FLUX_ADP_SML DET_FLUX_ADP \
    --progress-every 500
  local compact_h5_pre_sha compact_manifest_pre_sha validation_sha
  compact_h5_pre_sha="$(sha256_file "${pending_h5}")"
  compact_manifest_pre_sha="$(sha256_file "${pending_manifest}")"
  validation_sha="$(sha256_file "${VALIDATION}")"
  "${QLP_PYTHON}" scripts/stage5_validation/attest_teacher_v3_s63_compact.py \
    --accepted-validation "${VALIDATION}" \
    --expected-accepted-validation-sha256 "${validation_sha}" \
    --compact-h5 "${pending_h5}" \
    --expected-compact-h5-sha256 "${compact_h5_pre_sha}" \
    --compact-manifest "${pending_manifest}" \
    --expected-compact-manifest-sha256 "${compact_manifest_pre_sha}" \
    --producer-git-sha "${TWIRL_EXPECTED_GIT_SHA}" \
    --published-h5-path "${COMPACT_H5}"
  mv -- "${pending_h5}" "${COMPACT_H5}"
  mv -- "${pending_manifest}" "${COMPACT_MANIFEST}"
}

run_model_ready() {
  prepare_mutating_stage model_ready stage2
  require_file "${PLAN}" "prospective plan"
  require_file "${VALIDATION}" "accepted Stage-1 validation"
  require_file "${COMPACT_H5}" "compact HDF5"
  require_file "${COMPACT_MANIFEST}" "compact manifest"
  require_file "${RESERVED_TICS}" "sealed S63 TIC inventory"
  require_absent "${MODEL_READY_ALLOWLIST}" "model-ready allowlist"
  require_absent "${MODEL_READY_SUMMARY}" "model-ready summary"
  mkdir -p "$(dirname "${MODEL_READY_ALLOWLIST}")"
  "${STAGE2_PYTHON}" scripts/stage5_validation/build_teacher_v3_s63_model_ready_allowlist.py \
    --plan "${PLAN}" \
    --expected-plan-sha256 "$(sha256_file "${PLAN}")" \
    --accepted-validation "${VALIDATION}" \
    --expected-accepted-validation-sha256 "$(sha256_file "${VALIDATION}")" \
    --compact-h5 "${COMPACT_H5}" \
    --expected-compact-h5-sha256 "$(sha256_file "${COMPACT_H5}")" \
    --compact-manifest "${COMPACT_MANIFEST}" \
    --expected-compact-manifest-sha256 "$(sha256_file "${COMPACT_MANIFEST}")" \
    --reserved-tics "${RESERVED_TICS}" \
    --expected-reserved-tics-sha256 "$(sha256_file "${RESERVED_TICS}")" \
    --producer-git-sha "${TWIRL_EXPECTED_GIT_SHA}" \
    --out "${MODEL_READY_ALLOWLIST}" \
    --summary "${MODEL_READY_SUMMARY}"
}

run_cohorts() {
  prepare_mutating_stage cohorts stage2
  require_file "${PLAN}" "prospective plan"
  require_file "${MODEL_READY_ALLOWLIST}" "model-ready allowlist"
  require_file "${TEACHER_V3_CORPUS}" "frozen Teacher-v3 morphology corpus"
  require_file "${RESERVED_TICS}" "sealed S63 TIC inventory"
  require_absent "${PRIMARY_COHORT}" "primary cohort"
  require_absent "${REPEATED_COHORT}" "repeated-host cohort"
  require_absent "${RESERVED_NOT_READY}" "reserved-not-model-ready inventory"
  require_absent "${COHORT_SUMMARY}" "cohort summary"
  "${STAGE2_PYTHON}" scripts/stage5_validation/build_teacher_v3_s63_cohorts.py \
    --plan "${PLAN}" \
    --expected-plan-sha256 "$(sha256_file "${PLAN}")" \
    --model-ready-allowlist "${MODEL_READY_ALLOWLIST}" \
    --expected-model-ready-allowlist-sha256 "$(sha256_file "${MODEL_READY_ALLOWLIST}")" \
    --teacher-v3-corpus "${TEACHER_V3_CORPUS}" \
    --reserved-tics "${RESERVED_TICS}" \
    --producer-git-sha "${TWIRL_EXPECTED_GIT_SHA}" \
    --out-dir "${COHORT_DIR}"
}

run_bls_shard() {
  local shard_index=${2:?usage: $0 bls-shard <shard-index> <n-shards>}
  local n_shards=${3:?usage: $0 bls-shard <shard-index> <n-shards>}
  if (( n_shards < 1 || shard_index < 0 || shard_index >= n_shards )); then
    echo "BLS shard index must satisfy 0 <= index < n-shards." >&2
    exit 2
  fi
  prepare_mutating_stage "bls_${shard_index}_of_${n_shards}" stage2
  require_file "${COMPACT_H5}" "compact HDF5"
  require_file "${CADENCE_TABLE}" "cadence table"
  require_file "${CADENCE_MANIFEST}" "cadence manifest"
  require_file "${MODEL_READY_ALLOWLIST}" "model-ready allowlist"
  local suffix=""
  if (( n_shards > 1 )); then
    printf -v suffix '_%03d' "${shard_index}"
  fi
  require_absent "${BLS_SHARD_DIR}/real_adp_bls_peaks${suffix}.parquet" "BLS shard"
  require_absent "${BLS_SHARD_DIR}/summary${suffix}.json" "BLS shard summary"
  mkdir -p "${BLS_SHARD_DIR}"
  "${STAGE2_PYTHON}" scripts/stage5_validation/build_s56_adp_real_bls_peaks.py \
    --sector "${SECTOR}" \
    --compact-lc "${COMPACT_H5}" \
    --expected-compact-lc-sha256 "$(sha256_file "${COMPACT_H5}")" \
    --cadence-reference "${CADENCE_TABLE}" \
    --cadence-reference-manifest "${CADENCE_MANIFEST}" \
    --target-allowlist "${MODEL_READY_ALLOWLIST}" \
    --expected-target-allowlist-sha256 "$(sha256_file "${MODEL_READY_ALLOWLIST}")" \
    --orbitid-policy strict \
    --out-dir "${BLS_SHARD_DIR}" \
    --workers "${TWIRL_S63_BLS_WORKERS:-8}" \
    --n-periods 50000 \
    --n-peaks 10 \
    --shard-index "${shard_index}" \
    --n-shards "${n_shards}" \
    --progress-every 25 \
    --source-product-tag A2v1
  local bls_summary="${BLS_SHARD_DIR}/summary${suffix}.json"
  "${STAGE2_PYTHON}" scripts/stage5_validation/attest_s63_teacher_v3_json.py \
    --json "${bls_summary}" \
    --expected-sha256 "$(sha256_file "${bls_summary}")" \
    --producer-git-sha "${TWIRL_EXPECTED_GIT_SHA}"
}

run_bls_merge() {
  local n_shards=${2:?usage: $0 bls-merge <n-shards>}
  if (( n_shards < 1 )); then
    echo "n-shards must be positive." >&2
    exit 2
  fi
  prepare_mutating_stage "bls_merge_${n_shards}" stage2
  require_absent "${BLS_MERGED}" "merged BLS table"
  require_absent "${BLS_MERGED%.parquet}.summary.json" "merged BLS summary"
  mkdir -p "$(dirname "${BLS_MERGED}")"
  "${STAGE2_PYTHON}" scripts/stage5_validation/merge_adp_real_bls_shards.py \
    --shard-dir "${BLS_SHARD_DIR}" \
    --out-path "${BLS_MERGED}" \
    --n-shards "${n_shards}"
  local merged_summary="${BLS_MERGED%.parquet}.summary.json"
  "${STAGE2_PYTHON}" scripts/stage5_validation/attest_s63_teacher_v3_json.py \
    --json "${merged_summary}" \
    --expected-sha256 "$(sha256_file "${merged_summary}")" \
    --producer-git-sha "${TWIRL_EXPECTED_GIT_SHA}"
}

run_candidates() {
  prepare_mutating_stage candidates stage2
  require_file "${PLAN}" "prospective plan"
  require_file "${MODEL_READY_ALLOWLIST}" "model-ready allowlist"
  require_file "${BLS_MERGED}" "merged BLS table"
  require_file "${BLS_MERGED%.parquet}.summary.json" "merged BLS summary"
  require_absent "${CANDIDATES}" "rank-one candidate table"
  require_absent "${CANDIDATE_SUMMARY}" "rank-one candidate summary"
  mkdir -p "$(dirname "${CANDIDATES}")"
  "${STAGE2_PYTHON}" scripts/stage5_validation/build_teacher_v3_s63_candidates.py \
    --plan "${PLAN}" \
    --expected-plan-sha256 "$(sha256_file "${PLAN}")" \
    --model-ready-allowlist "${MODEL_READY_ALLOWLIST}" \
    --expected-model-ready-allowlist-sha256 "$(sha256_file "${MODEL_READY_ALLOWLIST}")" \
    --bls-peaks "${BLS_MERGED}" \
    --expected-bls-peaks-sha256 "$(sha256_file "${BLS_MERGED}")" \
    --bls-summary "${BLS_MERGED%.parquet}.summary.json" \
    --expected-bls-summary-sha256 "$(sha256_file "${BLS_MERGED%.parquet}.summary.json")" \
    --producer-git-sha "${TWIRL_EXPECTED_GIT_SHA}" \
    --out "${CANDIDATES}"
}

run_raw() {
  prepare_mutating_stage raw_source qlp
  require_file "${CANDIDATES}" "rank-one S63 candidate table"
  require_file "${COMPACT_H5}" "compact HDF5"
  require_absent "${RAW_SOURCE_H5}" "raw-source HDF5"
  require_absent "${RAW_SOURCE_H5%.h5}.summary.json" "raw-source summary"
  mkdir -p "$(dirname "${RAW_SOURCE_H5}")"
  TWIRL_TEACHER_V3_PROSPECTIVE_CONTRACT="${AUTHORIZATION_CONTRACT}" \
    TWIRL_TEACHER_V3_SECTOR_TABLE="${CANDIDATES}" \
    TWIRL_ROOT="${REPO}" \
    TWIRL_EXPECTED_GIT_SHA="${TWIRL_EXPECTED_GIT_SHA}" \
    PYTHON_BIN="${QLP_PYTHON}" \
    bash scripts/stage5_validation/run_teacher_v3_raw_export_pdo.sh \
      "${SECTOR}" "${ORBIT_ONE}" "${ORBIT_TWO}" "${COMPACT_H5}" \
      "${A2V1_ROOT}" "${RAW_SOURCE_H5}"
}

run_launch() {
  prepare_mutating_stage launch_manifest stage2
  for path in \
    "${PLAN}" "${SELECTION_POLICY}" "${VALIDATION}" "${CADENCE_TABLE}" "${CADENCE_MANIFEST}" \
    "${COMPACT_H5}" "${COMPACT_MANIFEST}" "${RESERVED_TICS}" \
    "${TEACHER_V3_CORPUS}" \
    "${MODEL_READY_ALLOWLIST}" "${MODEL_READY_SUMMARY}" \
    "${PRIMARY_COHORT}" "${REPEATED_COHORT}" \
    "${RESERVED_NOT_READY}" "${COHORT_SUMMARY}" "${BLS_MERGED}" \
    "${BLS_MERGED%.parquet}.summary.json" "${CANDIDATES}" \
    "${CANDIDATE_SUMMARY}" "${RAW_SOURCE_H5}" \
    "${RAW_SOURCE_H5%.h5}.summary.json" "${NATIVE_H5}" "${NATIVE_SUMMARY}"; do
    require_file "${path}" "launch input"
  done
  require_absent "${LAUNCH_MANIFEST}" "launch manifest"
  mkdir -p "$(dirname "${LAUNCH_MANIFEST}")"
  "${STAGE2_PYTHON}" scripts/stage5_validation/build_s63_teacher_v3_launch_manifest.py \
    --preregistered-contract "${PLAN}" \
    --selection-policy "${SELECTION_POLICY}" \
    --accepted-validation "${VALIDATION}" \
    --cadence-table "${CADENCE_TABLE}" \
    --cadence-manifest "${CADENCE_MANIFEST}" \
    --compact-h5 "${COMPACT_H5}" \
    --compact-manifest "${COMPACT_MANIFEST}" \
    --reserved-tics "${RESERVED_TICS}" \
    --teacher-v3-corpus "${TEACHER_V3_CORPUS}" \
    --model-ready-allowlist "${MODEL_READY_ALLOWLIST}" \
    --model-ready-summary "${MODEL_READY_SUMMARY}" \
    --primary-cohort "${PRIMARY_COHORT}" \
    --repeated-host-cohort "${REPEATED_COHORT}" \
    --reserved-not-model-ready "${RESERVED_NOT_READY}" \
    --cohort-summary "${COHORT_SUMMARY}" \
    --bls-peaks "${BLS_MERGED}" \
    --bls-summary "${BLS_MERGED%.parquet}.summary.json" \
    --candidates "${CANDIDATES}" \
    --candidate-summary "${CANDIDATE_SUMMARY}" \
    --raw-source-h5 "${RAW_SOURCE_H5}" \
    --raw-source-summary "${RAW_SOURCE_H5%.h5}.summary.json" \
    --native-h5 "${NATIVE_H5}" \
    --native-summary "${NATIVE_SUMMARY}" \
    --repo "${REPO}" \
    --expected-git-sha "${TWIRL_EXPECTED_GIT_SHA}" \
    --out "${LAUNCH_MANIFEST}"
}

show_plan() {
  cat <<EOF
Prospective contract: ${AUTHORIZATION_CONTRACT}
Pinned checkout:      TWIRL_EXPECTED_GIT_SHA=<full SHA> (required to execute)
Accepted Stage 1:     ${VALIDATION}
Immutable source:     ${SOURCE_VALIDATION} (never modified; copied + attested)
Cadence authority:    S63, orbits 133/134 -> ${CADENCE_DIR}
Compact ADP pair:     DET_FLUX_ADP_SML + DET_FLUX_ADP -> ${COMPACT_H5}
Model-ready allowlist:validated compact /targets keys -> ${MODEL_READY_ALLOWLIST}
Prospective cohorts:  ${COHORT_DIR}
Locked full BLS:      ORCD 16-task array, 50,000 periods, 10 peaks,
                      two ADP apertures, strict orbit IDs
Merged BLS:           ${BLS_MERGED}
Rank-one candidates:  ${CANDIDATES} (built by the prospective candidate CLI)
Raw source:           ${RAW_SOURCE_H5}
Native:               ORCD 32-task array + merge, strict orbit IDs, 4096 periods
Launch freeze:         ${LAUNCH_MANIFEST} (requires staged native output)

No stage is submitted by this controller. Run one explicit stage at a time:
  $0 preflight
  $0 stage1-receipt
  $0 cadence
  $0 compact
  $0 model-ready
  $0 cohorts
  $0 orcd-plan
  $0 bls-shard <index> <n-shards>   # PDO compatibility/smoke only
  $0 bls-merge <n-shards>           # PDO compatibility/smoke only
  $0 candidates
  $0 raw
  $0 launch
EOF
}

show_orcd_plan() {
  cat <<'EOF'
Full S63 downstream execution is intentionally staged on ORCD and never
submitted by this PDO controller. Export the exact TWIRL_S63_* input paths and
SHA-256 variables required by each wrapper, then submit in this order:

  sbatch scripts/orcd/slurm_teacher_v3_s63_bls_cpu.sbatch
  sbatch --dependency=afterok:<bls-array-jobid> scripts/orcd/slurm_teacher_v3_s63_bls_merge_cpu.sbatch
  # Build rank-one candidates from the merged BLS product.
  # Export/stage raw sources from PDO with this controller's `raw` stage.
  sbatch scripts/orcd/slurm_teacher_v3_s63_native_cpu.sbatch
  sbatch --dependency=afterok:<native-array-jobid> scripts/orcd/slurm_teacher_v3_s63_native_merge_cpu.sbatch
  # Stage the merged native result where `launch` can see it, then freeze launch.
  # Stage the frozen launch manifest back to ORCD and only then submit:
  sbatch scripts/orcd/slurm_teacher_v3_s63_score_h200.sbatch

The BLS array is locked to 16 shards (8 concurrent, 8 CPUs each). The native
array is locked to 32 shards (16 concurrent) with a separate 128-GiB merge.
Scoring is a distinct one-H200 job with CUDA required by default. Every wrapper
refuses existing outputs and requires the exact fully clean Git SHA.
EOF
}

show_status() {
  for path in \
    "${VALIDATION}" "${CADENCE_TABLE}" "${CADENCE_MANIFEST}" \
    "${COMPACT_H5}" "${COMPACT_MANIFEST}" "${MODEL_READY_ALLOWLIST}" \
    "${MODEL_READY_SUMMARY}" "${PRIMARY_COHORT}" "${REPEATED_COHORT}" \
    "${RESERVED_NOT_READY}" "${COHORT_SUMMARY}" \
    "${BLS_MERGED}" "${BLS_MERGED%.parquet}.summary.json" \
    "${CANDIDATES}" "${CANDIDATE_SUMMARY}" "${RAW_SOURCE_H5}" \
    "${RAW_SOURCE_H5%.h5}.summary.json" "${NATIVE_H5}" "${NATIVE_SUMMARY}" \
    "${LAUNCH_MANIFEST}"; do
    if [[ -s "${path}" ]]; then
      ls -lh "${path}"
    else
      echo "MISSING ${path}"
    fi
  done
}

case "${COMMAND}" in
  plan) show_plan ;;
  status) show_status ;;
  orcd-plan) show_orcd_plan ;;
  preflight)
    require_authorization
    export PYTHONNOUSERSITE=1
    unset LD_LIBRARY_PATH
    "${STAGE2_PYTHON}" -c 'import astropy, h5py, numpy, pandas, pyarrow; print("S63 Stage2/Parquet import preflight OK")'
    preflight_stage1
    ;;
  stage1-receipt) run_stage1_receipt ;;
  cadence) run_cadence ;;
  compact) run_compact ;;
  model-ready) run_model_ready ;;
  cohorts) run_cohorts ;;
  bls-shard) run_bls_shard "$@" ;;
  bls-merge) run_bls_merge "$@" ;;
  candidates) run_candidates ;;
  raw) run_raw ;;
  launch) run_launch ;;
  *)
    echo "Usage: $0 {plan|status|orcd-plan|preflight|stage1-receipt|cadence|compact|model-ready|cohorts|bls-shard I N|bls-merge N|candidates|raw|launch}" >&2
    exit 2
    ;;
esac
