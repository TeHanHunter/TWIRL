#!/usr/bin/env bash
# Guarded local controller for the approved Teacher v4-SSL full-pool v3 flow.
#
# This is an explicit, non-chained state machine.  It never starts password,
# Duo, or keyboard-interactive authentication, never resubmits a failed stage,
# and never admits native-v1/native-v2 products into the fresh v3 release.
set -euo pipefail

DRY_RUN=1
ACTION=""

usage() {
  cat <<'EOF'
Usage: TWIRL_EXPECTED_GIT_SHA=<40-hex> run_teacher_ssl_fullpool_v3_orcd.sh [--run] ACTION

Actions:
  probe                    Verify the user-opened ORCD control socket.
  deploy                   Create/reuse a clean detached worktree at the exact SHA.
  stage-eligibility        Copy the locked eligibility authority byte-identically
                           from v2 into the wholly fresh v3 run root.
  submit-native-canaries   Submit S56 shard 3 exclusion and S60 shard 4 ADP canaries.
  submit-native-canary-audit
                           Audit every eligible S60 shard-4 row through the exact
                           production dataset/collator/numerical-envelope path.
  submit-native            Gate both canaries, then submit seven 0-15%4 CPU arrays.
  submit-native-registry   Require 112 COMPLETED tasks, products, hashes, and empty
                           stderr files before submitting the native registry.
  submit-ssl-registry      Submit the 175,366-row audit / 175,347-row include registry.
  submit-numeric-audit     Submit the exact 112-task model-input numeric audit.
  submit-numeric-release   Merge all 112 passed numeric reports.
  submit-smoke             Submit fold-2/4,096-row/one-epoch/one-H200 smoke.
  submit-folds             Gate the smoke, then submit five one-H200 folds 0-4%4.
  validate-folds           Require five clean completed fold tasks, then submit
                           independent CPU validation/completion-release freeze.
  verify-completion        Read-only final proof: gate the validation job and
                           revalidate the exact immutable completion release.
  status                   Show v3 jobs and gate artifacts read-only.

Default mode is dry-run.  Every mutating action requires --run.  A failed job
is diagnosed separately; this controller does not cancel, retry, resubmit,
resume, weaken, or bypass a gate.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run)
      DRY_RUN=0
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      if [[ -n "${ACTION}" ]]; then
        echo "Only one action may be supplied." >&2
        exit 2
      fi
      ACTION="$1"
      shift
      ;;
  esac
done
if [[ -z "${ACTION}" ]]; then
  usage >&2
  exit 2
fi

: "${TWIRL_EXPECTED_GIT_SHA:?set TWIRL_EXPECTED_GIT_SHA to the exact reviewed 40-hex commit}"
EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}"
if [[ ! "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]]; then
  echo "TWIRL_EXPECTED_GIT_SHA must be a full lowercase Git SHA." >&2
  exit 2
fi

LOCAL_REPO="${LOCAL_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
readonly LOCKED_ORCD_HOST="tehan@orcd-login.mit.edu"
readonly LOCKED_CONTROL_PATH="/Users/tehan/.ssh/cm/%r@%h:%p"
if [[ -n "${ORCD_HOST+x}" && "${ORCD_HOST}" != "${LOCKED_ORCD_HOST}" ]]; then
  echo "ORCD_HOST differs from the locked noninteractive endpoint." >&2
  exit 2
fi
if [[ -n "${ORCD_CONTROL_PATH+x}" && "${ORCD_CONTROL_PATH}" != "${LOCKED_CONTROL_PATH}" ]]; then
  echo "ORCD_CONTROL_PATH differs from the locked user-opened socket." >&2
  exit 2
fi
readonly ORCD_HOST="${LOCKED_ORCD_HOST}"
readonly CONTROL_PATH="${LOCKED_CONTROL_PATH}"
readonly TWIRL_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly REMOTE_SOURCE="${TWIRL_ROOT}/code/TWIRL"
readonly REMOTE_REPO="${TWIRL_ROOT}/code/TWIRL_teacher_ssl_fullpool_v3_${EXPECTED_SHA:0:12}"
readonly GIT_BRANCH="codex/teacher-ssl-effective-quality-adp-v3"
readonly SOURCE_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp"
readonly HISTORICAL_V2_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v2_s56_s62_a2v1_current_adp_bls_eligible"
readonly RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v3_s56_s62_a2v1_effective_quality_adp_bls_eligible"
readonly CLAIM_PREFIX="${TWIRL_ROOT}/reports/stage5_validation/.teacher_ssl_fullpool_v3_s56_s62_a2v1_effective_quality_adp_bls_eligible"
readonly V2_ELIGIBILITY_DIR="${HISTORICAL_V2_ROOT}/frozen/native_model_eligibility"
readonly ELIGIBILITY_DIR="${RUN_ROOT}/frozen/native_model_eligibility"
readonly ELIGIBILITY_EXCLUSIONS="${ELIGIBILITY_DIR}/native_model_exclusions.csv"
readonly ELIGIBILITY_SUMMARY="${ELIGIBILITY_DIR}/native_model_eligibility.summary.json"
readonly CANARY_ROOT="${RUN_ROOT}/native_canary"
readonly S56_CANARY_H5="${CANARY_ROOT}/s56/native_3.h5"
readonly S56_CANARY_SUMMARY="${CANARY_ROOT}/s56/native_3.summary.json"
readonly S60_CANARY_H5="${CANARY_ROOT}/s60/native_4.h5"
readonly S60_CANARY_SUMMARY="${CANARY_ROOT}/s60/native_4.summary.json"
readonly S60_CANARY_NUMERIC_REPORT="${CANARY_ROOT}/s60/native_4.model_facing_numeric.json"
readonly LAUNCH_DIR="${RUN_ROOT}/launch"
readonly CANARY_JOBS="${LAUNCH_DIR}/native_canaries.tsv"
readonly CANARY_NUMERIC_JOB="${LAUNCH_DIR}/native_canary_numeric.job"
readonly NATIVE_JOBS="${LAUNCH_DIR}/native_arrays.tsv"
readonly NATIVE_REGISTRY_JOB="${LAUNCH_DIR}/native_registry.job"
readonly SSL_REGISTRY_JOB="${LAUNCH_DIR}/ssl_registry.job"
readonly NUMERIC_AUDIT_JOB="${LAUNCH_DIR}/numeric_audit.job"
readonly NUMERIC_RELEASE_JOB="${LAUNCH_DIR}/numeric_release.job"
readonly SMOKE_JOB="${LAUNCH_DIR}/smoke.job"
readonly FOLDS_JOB="${LAUNCH_DIR}/folds.job"
readonly FOLD_VALIDATION_JOB="${LAUNCH_DIR}/fold_validation.job"
readonly NATIVE_REGISTRY_DIR="${RUN_ROOT}/frozen/native_registry"
readonly NATIVE_REGISTRY="${NATIVE_REGISTRY_DIR}/native_input_registry.parquet"
readonly NATIVE_REGISTRY_SUMMARY="${NATIVE_REGISTRY_DIR}/native_input_registry.summary.json"
readonly NATIVE_RELEASE_SUMMARY="${NATIVE_REGISTRY_DIR}/native_fullpool_release.summary.json"
readonly SSL_REGISTRY_DIR="${RUN_ROOT}/frozen/ssl_registry"
readonly SSL_REGISTRY="${SSL_REGISTRY_DIR}/teacher_ssl_fullpool_registry.parquet"
readonly SSL_REGISTRY_SUMMARY="${SSL_REGISTRY_DIR}/teacher_ssl_fullpool_registry.summary.json"
readonly NUMERIC_GATE_DIR="${RUN_ROOT}/frozen/model_input_numeric_gate_v2"
readonly NUMERIC_SHARD_DIR="${NUMERIC_GATE_DIR}/shard_reports"
readonly NUMERIC_AUDIT="${NUMERIC_GATE_DIR}/model_input_numeric_audit.parquet"
readonly NUMERIC_RELEASE="${NUMERIC_GATE_DIR}/model_input_numeric_gate.release.json"
readonly MODEL_RUN_ROOT="${RUN_ROOT}/model_runs/effective_quality_adp_v1"
readonly SMOKE_SUMMARY="${MODEL_RUN_ROOT}/smoke/one_epoch/encoder_pretraining/fold_2/summary.json"
readonly FOLD_COMPLETION_RELEASE="${RUN_ROOT}/frozen/model/teacher_ssl_fullpool_v3_five_fold.complete.json"
readonly NATIVE_CONTRACT="twirl_teacher_ssl_fullpool_real_native_v3_effective_quality_adp"
readonly RELEASE_BINDING="teacher_ssl_fullpool_v3_s56_s62_a2v1_effective_quality_adp_bls_eligible"
readonly DETREND_CONFIG_SHA256="2fcded6ae0dc1ada429c37ea1f39f68d71feba143d47925af3e04e38e1ecdec6"
readonly FULL_IDENTITY_SHA256="8e9e9c12a24d5ebc7be94b03a4e35411cd10066d62a87d921a8443b06cc188d1"
readonly ELIGIBLE_IDENTITY_SHA256="6ddc8e57bb5fb938ce05389c1629221c73e0e73ac3bf40da47a2019e1a5660e6"
readonly EXCLUDED_IDENTITY_SHA256="b9f536144265e54a70bff17c782e0668ddbd96efcdf6c223ebc58f46edb7d976"
readonly S56_CANARY_EXCLUDED_SHA256="ddda9b053eddc744e5032ba350598d58d74cea4f4cc5cd705932ccf6e41022ab"
readonly S60_CANARY_TIC="704538814"
readonly SMOKE_OBSERVATION_ID="s0060-tic0000000722078603"
readonly NATIVE_SHARDS=16

SSH=(
  ssh
  -o BatchMode=yes
  -o PasswordAuthentication=no
  -o KbdInteractiveAuthentication=no
  -o PubkeyAuthentication=no
  -o HostbasedAuthentication=no
  -o GSSAPIAuthentication=no
  -o NumberOfPasswordPrompts=0
  -o ConnectTimeout=15
  -o ConnectionAttempts=1
  -o ControlMaster=no
  -o ControlPath="${CONTROL_PATH}"
)

print_cmd() {
  printf '+'
  printf ' %q' "$@"
  printf '\n'
}

remote() {
  if [[ "${DRY_RUN}" == "1" ]]; then
    print_cmd "${SSH[@]}" "${ORCD_HOST}" "$1"
  else
    "${SSH[@]}" "${ORCD_HOST}" "$1"
  fi
}

require_socket() {
  if [[ "${DRY_RUN}" == "1" ]]; then
    return
  fi
  if ! "${SSH[@]}" -O check "${ORCD_HOST}" >/dev/null 2>&1; then
    cat >&2 <<EOF
The documented ORCD control socket is unavailable at:
  ${CONTROL_PATH}

Open it from a user terminal as documented in doc/orcd_h200_usage.md.
This controller will not fall back to password or Duo authentication.
EOF
    exit 4
  fi
}

base_export() {
  printf '%s' \
    "TWIRL_ORCD_REPO=${REMOTE_REPO}" \
    ",TWIRL_SSL_FULLPOOL_V1_RUN_ROOT=${SOURCE_RUN_ROOT}" \
    ",TWIRL_SSL_FULLPOOL_V3_RUN_ROOT=${RUN_ROOT}" \
    ",TWIRL_EXPECTED_GIT_SHA=${EXPECTED_SHA}"
}

case "${ACTION}" in
  probe)
    require_socket
    remote "hostname; whoami; sinfo -h -p pg_mki_aryeh -o '%P %D %t %G' | head -10"
    ;;

  deploy)
    require_socket
    remote "
      set -euo pipefail
      test -d '${REMOTE_SOURCE}/.git'
      mkdir -- '${CLAIM_PREFIX}.deploy.claim'
      GIT_TERMINAL_PROMPT=0 git -C '${REMOTE_SOURCE}' fetch origin '${GIT_BRANCH}'
      git -C '${REMOTE_SOURCE}' cat-file -e '${EXPECTED_SHA}^{commit}'
      if [[ -e '${REMOTE_REPO}' ]]; then
        test -e '${REMOTE_REPO}/.git'
        [[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]
      else
        git -C '${REMOTE_SOURCE}' worktree add --detach '${REMOTE_REPO}' '${EXPECTED_SHA}'
      fi
      [[ -z \$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all) ]]
      git -C '${REMOTE_REPO}' rev-parse HEAD
    "
    ;;

  stage-eligibility)
    require_socket
    remote "
      set -euo pipefail
      [[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]
      [[ -z \$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all) ]]
      source_exclusions='${V2_ELIGIBILITY_DIR}/native_model_exclusions.csv'
      source_summary='${V2_ELIGIBILITY_DIR}/native_model_eligibility.summary.json'
      test -s \"\${source_exclusions}\"
      test -s \"\${source_summary}\"
      [[ ! -e '${RUN_ROOT}' ]]
      staging='${RUN_ROOT}.staging-eligibility'
      [[ ! -e \"\${staging}\" ]]
      mkdir -- '${CLAIM_PREFIX}.stage-eligibility.claim'
      mkdir -p \"\${staging}/frozen/native_model_eligibility\"
      cp -- \"\${source_exclusions}\" \"\${staging}/frozen/native_model_eligibility/native_model_exclusions.csv\"
      cp -- \"\${source_summary}\" \"\${staging}/frozen/native_model_eligibility/native_model_eligibility.summary.json\"
      cmp -s -- \"\${source_exclusions}\" \"\${staging}/frozen/native_model_eligibility/native_model_exclusions.csv\"
      cmp -s -- \"\${source_summary}\" \"\${staging}/frozen/native_model_eligibility/native_model_eligibility.summary.json\"
      '${TWIRL_ROOT}/envs/twirl-s56/bin/python' -c 'import json,sys; d=json.load(open(sys.argv[1])); c=d[\"counts\"]; h=d[\"identity_hashes\"]; a=d[\"partition_audit\"]; assert d.get(\"passed\") is True; assert c[\"full\"][\"n_observations\"] == 175366; assert c[\"eligible\"][\"n_observations\"] == 175347; assert c[\"excluded\"][\"n_observations\"] == 19; assert h[\"full_observations_sha256\"] == sys.argv[2]; assert h[\"eligible_observations_sha256\"] == sys.argv[3]; assert h[\"excluded_observations_sha256\"] == sys.argv[4]; assert a[\"eligible_excluded_disjoint\"] is True; assert a[\"eligible_excluded_union_equals_full\"] is True; assert a[\"production_lock_passed\"] is True' \"\${staging}/frozen/native_model_eligibility/native_model_eligibility.summary.json\" '${FULL_IDENTITY_SHA256}' '${ELIGIBLE_IDENTITY_SHA256}' '${EXCLUDED_IDENTITY_SHA256}'
      chmod a-w \"\${staging}/frozen/native_model_eligibility/native_model_exclusions.csv\" \"\${staging}/frozen/native_model_eligibility/native_model_eligibility.summary.json\"
      mv -- \"\${staging}\" '${RUN_ROOT}'
      cmp -s -- \"\${source_exclusions}\" '${ELIGIBILITY_EXCLUSIONS}'
      cmp -s -- \"\${source_summary}\" '${ELIGIBILITY_SUMMARY}'
      sha256sum \"\${source_exclusions}\" '${ELIGIBILITY_EXCLUSIONS}' \"\${source_summary}\" '${ELIGIBILITY_SUMMARY}'
    "
    ;;

  submit-native-canaries)
    require_socket
    remote "
      set -euo pipefail
      [[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]
      [[ -z \$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all) ]]
      cmp -s '${V2_ELIGIBILITY_DIR}/native_model_exclusions.csv' '${ELIGIBILITY_EXCLUSIONS}'
      cmp -s '${V2_ELIGIBILITY_DIR}/native_model_eligibility.summary.json' '${ELIGIBILITY_SUMMARY}'
      for output in \
        '${S56_CANARY_H5}' '${S56_CANARY_H5}.sha256' '${S56_CANARY_SUMMARY}' \
        '${S60_CANARY_H5}' '${S60_CANARY_H5}.sha256' '${S60_CANARY_SUMMARY}' \
        '${S60_CANARY_NUMERIC_REPORT}' '${S60_CANARY_NUMERIC_REPORT}.sha256' \
        '${CANARY_JOBS}' '${CANARY_JOBS}.submitting' \
        '${CANARY_NUMERIC_JOB}' '${CANARY_NUMERIC_JOB}.tmp'; do
        [[ ! -e \"\${output}\" ]]
      done
      mkdir -- '${CLAIM_PREFIX}.submit-native-canaries.claim'
      mkdir -p '${LAUNCH_DIR}'
      : > '${CANARY_JOBS}.submitting'
      cd '${REMOTE_REPO}'
      job=\$(sbatch --parsable --job-name=twirl-ssl-v3-canary-s56 \
        --output='${TWIRL_ROOT}/logs/twirl-ssl-v3-canary-s56-%j.out' \
        --error='${TWIRL_ROOT}/logs/twirl-ssl-v3-canary-s56-%j.err' \
        --export='$(base_export),TWIRL_FULLPOOL_SECTOR=56,TWIRL_FULLPOOL_NATIVE_SHARDS=${NATIVE_SHARDS},TWIRL_FULLPOOL_SHARD_INDEX=3,TWIRL_FULLPOOL_NATIVE_OUT=${S56_CANARY_H5}' \
        scripts/orcd/slurm_teacher_ssl_fullpool_v3_native_cpu.sbatch)
      job=\${job%%;*}
      printf '56\\t3\\t%s\\ttwirl-ssl-v3-canary-s56\\n' \"\${job}\" >> '${CANARY_JOBS}.submitting'
      job=\$(sbatch --parsable --job-name=twirl-ssl-v3-canary-s60 \
        --output='${TWIRL_ROOT}/logs/twirl-ssl-v3-canary-s60-%j.out' \
        --error='${TWIRL_ROOT}/logs/twirl-ssl-v3-canary-s60-%j.err' \
        --export='$(base_export),TWIRL_FULLPOOL_SECTOR=60,TWIRL_FULLPOOL_NATIVE_SHARDS=${NATIVE_SHARDS},TWIRL_FULLPOOL_SHARD_INDEX=4,TWIRL_FULLPOOL_NATIVE_OUT=${S60_CANARY_H5}' \
        scripts/orcd/slurm_teacher_ssl_fullpool_v3_native_cpu.sbatch)
      job=\${job%%;*}
      printf '60\\t4\\t%s\\ttwirl-ssl-v3-canary-s60\\n' \"\${job}\" >> '${CANARY_JOBS}.submitting'
      [[ \$(wc -l < '${CANARY_JOBS}.submitting') -eq 2 ]]
      mv '${CANARY_JOBS}.submitting' '${CANARY_JOBS}'
      cat '${CANARY_JOBS}'
    "
    ;;

  submit-native-canary-audit)
    require_socket
    remote "
      set -euo pipefail
      [[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]
      [[ -z \$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all) ]]
      test -s '${CANARY_JOBS}'
      [[ \$(wc -l < '${CANARY_JOBS}') -eq 2 ]]
      record=\$(awk -F'\\t' '\$1 == 60 && \$2 == 4 {print; found += 1} END {if (found != 1) exit 1}' '${CANARY_JOBS}')
      IFS=\$'\\t' read -r sector shard job name <<< \"\${record}\"
      [[ \"\${sector}\" == 60 && \"\${shard}\" == 4 ]]
      state=\$(sacct -X -n -P -j \"\${job}\" -o JobIDRaw,State,ExitCode | awk -F'|' -v j=\"\${job}\" '\$1==j {print \$2 \"|\" \$3}')
      [[ \"\${state}\" == 'COMPLETED|0:0' ]]
      err='${TWIRL_ROOT}/logs/'\"\${name}\"'-'\"\${job}\"'.err'
      test -e \"\${err}\"
      test ! -s \"\${err}\"
      for input in '${S60_CANARY_H5}' '${S60_CANARY_H5}.sha256' '${S60_CANARY_SUMMARY}'; do
        test -s \"\${input}\"
      done
      sha256sum -c '${S60_CANARY_H5}.sha256'
      for output in \
        '${S60_CANARY_NUMERIC_REPORT}' '${S60_CANARY_NUMERIC_REPORT}.sha256' \
        '${CANARY_NUMERIC_JOB}' '${CANARY_NUMERIC_JOB}.tmp'; do
        [[ ! -e \"\${output}\" ]]
      done
      mkdir -- '${CLAIM_PREFIX}.submit-native-canary-audit.claim'
      cd '${REMOTE_REPO}'
      job=\$(sbatch --parsable --export='$(base_export)' \
        scripts/orcd/slurm_teacher_ssl_fullpool_v3_canary_numeric_cpu.sbatch)
      job=\${job%%;*}
      printf '%s\\n' \"\${job}\" > '${CANARY_NUMERIC_JOB}.tmp'
      mv '${CANARY_NUMERIC_JOB}.tmp' '${CANARY_NUMERIC_JOB}'
      printf 'native_canary_numeric=%s\\n' \"\${job}\"
    "
    ;;

  submit-native)
    require_socket
    remote "
      set -euo pipefail
      [[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]
      [[ -z \$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all) ]]
      test -s '${CANARY_NUMERIC_JOB}'
      numeric_job=\$(< '${CANARY_NUMERIC_JOB}')
      state=\$(sacct -X -n -P -j \"\${numeric_job}\" -o JobIDRaw,State,ExitCode | awk -F'|' -v j=\"\${numeric_job}\" '\$1==j {print \$2 \"|\" \$3}')
      [[ \"\${state}\" == 'COMPLETED|0:0' ]]
      test -e '${TWIRL_ROOT}/logs/twirl-ssl-v3-canary-numeric-'\"\${numeric_job}\"'.err'
      test ! -s '${TWIRL_ROOT}/logs/twirl-ssl-v3-canary-numeric-'\"\${numeric_job}\"'.err'
      test -s '${S60_CANARY_NUMERIC_REPORT}'
      test -s '${S60_CANARY_NUMERIC_REPORT}.sha256'
      (
        cd '${CANARY_ROOT}/s60'
        sha256sum -c 'native_4.model_facing_numeric.json.sha256'
      )
      native_sha=\$(sha256sum '${S60_CANARY_H5}' | awk '{print \$1}')
      '${TWIRL_ROOT}/envs/twirl-s56/bin/python' -c 'import json,sys; d=json.load(open(sys.argv[1])); c=d[\"counts\"]; i=d[\"identities\"]; assert d[\"schema_version\"]==\"twirl_teacher_ssl_fullpool_v3_s60_shard4_model_facing_canary_v1\"; assert d[\"code_revision\"]==sys.argv[2]; assert d[\"sector\"]==60 and d[\"shard_index\"]==4 and d[\"n_shards\"]==16; assert d[\"required_tic\"]==704538814 and d[\"required_tic_passed\"] is True; assert d[\"passed\"] is True and d[\"n_failures\"]==0; assert c[\"eligible_observations\"]==c[\"audited_observations\"]==c[\"passed_observations\"]; assert c[\"failed_observations\"]==0; assert i[\"full_observations_sha256\"]==sys.argv[3]; assert i[\"eligible_observations_sha256\"]==sys.argv[4]; assert d[\"native_h5\"][\"sha256\"]==sys.argv[5]; assert d[\"native_contract_version\"]==sys.argv[6]; assert d[\"native_release_binding\"]==sys.argv[7]; assert d[\"action\"]==\"audit_only_no_clip_no_exclusion\"; assert all(r[\"passed\"] is True and r[\"n_failures\"]==0 and r[\"action\"]==\"audit_only_no_clip_no_exclusion\" for r in d[\"rows\"])' '${S60_CANARY_NUMERIC_REPORT}' '${EXPECTED_SHA}' '${FULL_IDENTITY_SHA256}' '${ELIGIBLE_IDENTITY_SHA256}' \"\${native_sha}\" '${NATIVE_CONTRACT}' '${RELEASE_BINDING}'
      test -s '${CANARY_JOBS}'
      [[ \$(wc -l < '${CANARY_JOBS}') -eq 2 ]]
      while IFS=\$'\\t' read -r sector shard job name; do
        state=\$(sacct -X -n -P -j \"\${job}\" -o JobIDRaw,State,ExitCode | awk -F'|' -v j=\"\${job}\" '\$1==j {print \$2 \"|\" \$3}')
        [[ \"\${state}\" == 'COMPLETED|0:0' ]]
        err='${TWIRL_ROOT}/logs/'\"\${name}\"'-'\"\${job}\"'.err'
        test -e \"\${err}\"
        test ! -s \"\${err}\"
      done < '${CANARY_JOBS}'
      for product in \
        '${S56_CANARY_H5}' '${S56_CANARY_H5}.sha256' '${S56_CANARY_SUMMARY}' \
        '${S60_CANARY_H5}' '${S60_CANARY_H5}.sha256' '${S60_CANARY_SUMMARY}'; do
        test -s \"\${product}\"
      done
      sha256sum -c '${S56_CANARY_H5}.sha256'
      sha256sum -c '${S60_CANARY_H5}.sha256'
      PYTHONPATH='${REMOTE_REPO}/src' '${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python' -c 'import h5py,json,sys; s56=json.load(open(sys.argv[1])); s60=json.load(open(sys.argv[2])); assert s56[\"sector\"]==56 and s56[\"shard_index\"]==3 and s56[\"n_shards\"]==16; assert s56[\"native_contract_version\"]==sys.argv[5]; assert s56[\"n_shard_excluded_observations\"]==1; assert s56[\"shard_excluded_observation_identity_sha256\"]==sys.argv[6]; assert s56[\"verification\"][\"passed\"] is True; assert s60[\"sector\"]==60 and s60[\"shard_index\"]==4 and s60[\"n_shards\"]==16; assert s60[\"native_contract_version\"]==sys.argv[5]; assert s60[\"verification\"][\"passed\"] is True; h=h5py.File(sys.argv[3],\"r\"); assert h.attrs[\"contract_version\"]==sys.argv[5]; assert h.attrs[\"release_binding\"]==sys.argv[7]; assert h.attrs[\"builder_contract_version\"]==\"twirl_teacher_ssl_fullpool_effective_quality_adp_builder_v1\"; assert h.attrs[\"detrend_contract_version\"]==\"twirl_fs_adp03q_effective_quality_v1\"; assert h.attrs[\"detrend_config_sha256\"]==sys.argv[8]; assert h.attrs[\"detrend_quality_source\"]==\"final_effective_quality\"; assert int(h.attrs[\"raw_photometry_only\"])==1; assert int(h.attrs[\"compact_adp_photometry_reused\"])==0; assert int(h.attrs[\"compact_adp_flux_reused\"])==0; h.close(); h=h5py.File(sys.argv[4],\"r\"); assert h.attrs[\"contract_version\"]==sys.argv[5]; assert h.attrs[\"release_binding\"]==sys.argv[7]; g=h[\"targets/0000000704538814\"]; assert int(g.attrs[\"raw_compact_internal_quality_agreement\"])==1; assert int(h.attrs[\"raw_photometry_only\"])==1; assert int(h.attrs[\"compact_adp_photometry_reused\"])==0; assert int(h.attrs[\"compact_adp_flux_reused\"])==0; h.close()' '${S56_CANARY_SUMMARY}' '${S60_CANARY_SUMMARY}' '${S56_CANARY_H5}' '${S60_CANARY_H5}' '${NATIVE_CONTRACT}' '${S56_CANARY_EXCLUDED_SHA256}' '${RELEASE_BINDING}' '${DETREND_CONFIG_SHA256}'
      for sector in {56..62}; do
        test ! -e '${RUN_ROOT}/native/s'\${sector}
      done
      [[ ! -e '${NATIVE_JOBS}' && ! -e '${NATIVE_JOBS}.submitting' ]]
      mkdir -- '${CLAIM_PREFIX}.submit-native.claim'
      : > '${NATIVE_JOBS}.submitting'
      cd '${REMOTE_REPO}'
      for sector in {56..62}; do
        name='twirl-ssl-v3-native-s'\${sector}
        job=\$(sbatch --parsable --job-name=\"\${name}\" --array=0-15%4 \
          --export='$(base_export),TWIRL_FULLPOOL_NATIVE_SHARDS=${NATIVE_SHARDS},TWIRL_FULLPOOL_SECTOR='\${sector} \
          scripts/orcd/slurm_teacher_ssl_fullpool_v3_native_cpu.sbatch)
        job=\${job%%;*}
        printf '%s\\t%s\\t%s\\n' \"\${sector}\" \"\${job}\" \"\${name}\" >> '${NATIVE_JOBS}.submitting'
      done
      [[ \$(wc -l < '${NATIVE_JOBS}.submitting') -eq 7 ]]
      mv '${NATIVE_JOBS}.submitting' '${NATIVE_JOBS}'
      cat '${NATIVE_JOBS}'
    "
    ;;

  submit-native-registry)
    require_socket
    remote "
      set -euo pipefail
      test -s '${NATIVE_JOBS}'
      [[ \$(wc -l < '${NATIVE_JOBS}') -eq 7 ]]
      total_tasks=0
      while IFS=\$'\\t' read -r sector job name; do
        mapfile -t states < <(sacct -n -P -j \"\${job}\" -o JobIDRaw,State,ExitCode | awk -F'|' -v p=\"\${job}_\" '\$1 ~ (\"^\" p \"[0-9]+$\") {print \$2 \"|\" \$3}')
        [[ \${#states[@]} -eq 16 ]]
        for state in \"\${states[@]}\"; do [[ \"\${state}\" == 'COMPLETED|0:0' ]]; done
        shopt -s nullglob
        errs=('${TWIRL_ROOT}/logs/'\"\${name}\"'-'\"\${job}\"_*.err)
        [[ \${#errs[@]} -eq 16 ]]
        for err in \"\${errs[@]}\"; do test ! -s \"\${err}\"; done
        total_tasks=\$((total_tasks + 16))
      done < '${NATIVE_JOBS}'
      [[ \${total_tasks} -eq 112 ]]
      products=0
      for sector in {56..62}; do
        for shard in {0..15}; do
          h5='${RUN_ROOT}/native/s'\${sector}'/native_'\${shard}'.h5'
          summary='${RUN_ROOT}/native/s'\${sector}'/native_'\${shard}'.summary.json'
          test -s \"\${h5}\"
          test -s \"\${summary}\"
          test -s \"\${h5}.sha256\"
          sha256sum -c \"\${h5}.sha256\"
          products=\$((products + 1))
        done
      done
      [[ \${products} -eq 112 ]]
      PYTHONPATH='${REMOTE_REPO}/src' '${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python' -c 'import json,pathlib,sys; root=pathlib.Path(sys.argv[1]); summaries=sorted(root.glob(\"s*/native_*.summary.json\")); assert len(summaries)==112; ds=[json.load(open(p)) for p in summaries]; assert {(int(d[\"sector\"]),int(d[\"shard_index\"])) for d in ds}=={(s,i) for s in range(56,63) for i in range(16)}; assert all(d[\"native_contract_version\"]==sys.argv[2] and d[\"release_binding\"]==sys.argv[3] and d[\"builder_contract_version\"]==\"twirl_teacher_ssl_fullpool_effective_quality_adp_builder_v1\" and d[\"detrend_contract_version\"]==\"twirl_fs_adp03q_effective_quality_v1\" and d[\"detrend_config_sha256\"]==sys.argv[4] and d[\"expected_git_sha\"]==sys.argv[5] and d[\"verification\"][\"passed\"] is True for d in ds); assert sum(int(d[\"n_source_shard_observations\"]) for d in ds)==175366; assert sum(int(d[\"n_shard_observations\"]) for d in ds)==175347; assert sum(int(d[\"n_shard_excluded_observations\"]) for d in ds)==19' '${RUN_ROOT}/native' '${NATIVE_CONTRACT}' '${RELEASE_BINDING}' '${DETREND_CONFIG_SHA256}' '${EXPECTED_SHA}'
      for output in '${NATIVE_REGISTRY}' '${NATIVE_REGISTRY_SUMMARY}' '${NATIVE_RELEASE_SUMMARY}'; do [[ ! -e \"\${output}\" ]]; done
      [[ ! -e '${NATIVE_REGISTRY_JOB}' && ! -e '${NATIVE_REGISTRY_JOB}.tmp' ]]
      mkdir -- '${CLAIM_PREFIX}.submit-native-registry.claim'
      cd '${REMOTE_REPO}'
      job=\$(sbatch --parsable --export='$(base_export),TWIRL_FULLPOOL_NATIVE_SHARDS=${NATIVE_SHARDS}' \
        scripts/orcd/slurm_teacher_ssl_fullpool_v3_native_registry_cpu.sbatch)
      job=\${job%%;*}
      printf '%s\\n' \"\${job}\" > '${NATIVE_REGISTRY_JOB}.tmp'
      mv '${NATIVE_REGISTRY_JOB}.tmp' '${NATIVE_REGISTRY_JOB}'
      printf 'native_registry=%s\\n' \"\${job}\"
    "
    ;;

  submit-ssl-registry)
    require_socket
    remote "
      set -euo pipefail
      test -s '${NATIVE_REGISTRY_JOB}'
      job=\$(< '${NATIVE_REGISTRY_JOB}')
      state=\$(sacct -X -n -P -j \"\${job}\" -o JobIDRaw,State,ExitCode | awk -F'|' -v j=\"\${job}\" '\$1==j {print \$2 \"|\" \$3}')
      [[ \"\${state}\" == 'COMPLETED|0:0' ]]
      test -e '${TWIRL_ROOT}/logs/twirl-ssl-v3-native-registry-'\"\${job}\"'.err'
      test ! -s '${TWIRL_ROOT}/logs/twirl-ssl-v3-native-registry-'\"\${job}\"'.err'
      for input in '${NATIVE_REGISTRY}' '${NATIVE_REGISTRY_SUMMARY}' '${NATIVE_RELEASE_SUMMARY}'; do test -s \"\${input}\"; done
      for output in '${SSL_REGISTRY}' '${SSL_REGISTRY_SUMMARY}'; do [[ ! -e \"\${output}\" ]]; done
      [[ ! -e '${SSL_REGISTRY_JOB}' && ! -e '${SSL_REGISTRY_JOB}.tmp' ]]
      [[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]
      [[ -z \$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all) ]]
      mkdir -- '${CLAIM_PREFIX}.submit-ssl-registry.claim'
      cd '${REMOTE_REPO}'
      job=\$(sbatch --parsable --export='$(base_export)' \
        scripts/orcd/slurm_teacher_ssl_fullpool_v3_registry_cpu.sbatch)
      job=\${job%%;*}
      printf '%s\\n' \"\${job}\" > '${SSL_REGISTRY_JOB}.tmp'
      mv '${SSL_REGISTRY_JOB}.tmp' '${SSL_REGISTRY_JOB}'
      printf 'ssl_registry=%s\\n' \"\${job}\"
    "
    ;;

  submit-numeric-audit)
    require_socket
    remote "
      set -euo pipefail
      test -s '${SSL_REGISTRY_JOB}'
      job=\$(< '${SSL_REGISTRY_JOB}')
      state=\$(sacct -X -n -P -j \"\${job}\" -o JobIDRaw,State,ExitCode | awk -F'|' -v j=\"\${job}\" '\$1==j {print \$2 \"|\" \$3}')
      [[ \"\${state}\" == 'COMPLETED|0:0' ]]
      test -e '${TWIRL_ROOT}/logs/twirl-ssl-v3-registry-'\"\${job}\"'.err'
      test ! -s '${TWIRL_ROOT}/logs/twirl-ssl-v3-registry-'\"\${job}\"'.err'
      for input in \
        '${NATIVE_REGISTRY}' '${NATIVE_REGISTRY_SUMMARY}' '${NATIVE_RELEASE_SUMMARY}' \
        '${SSL_REGISTRY}' '${SSL_REGISTRY_SUMMARY}'; do test -s \"\${input}\"; done
      [[ ! -e '${NUMERIC_GATE_DIR}' ]]
      [[ ! -e '${NUMERIC_AUDIT_JOB}' && ! -e '${NUMERIC_AUDIT_JOB}.tmp' ]]
      [[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]
      [[ -z \$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all) ]]
      mkdir -- '${CLAIM_PREFIX}.submit-numeric-audit.claim'
      cd '${REMOTE_REPO}'
      job=\$(sbatch --parsable --export='$(base_export)' \
        scripts/orcd/slurm_teacher_ssl_fullpool_v3_numeric_audit_cpu.sbatch)
      job=\${job%%;*}
      printf '%s\\n' \"\${job}\" > '${NUMERIC_AUDIT_JOB}.tmp'
      mv '${NUMERIC_AUDIT_JOB}.tmp' '${NUMERIC_AUDIT_JOB}'
      printf 'numeric_audit=%s\\n' \"\${job}\"
    "
    ;;

  submit-numeric-release)
    require_socket
    remote "
      set -euo pipefail
      test -s '${NUMERIC_AUDIT_JOB}'
      job=\$(< '${NUMERIC_AUDIT_JOB}')
      mapfile -t states < <(sacct -n -P -j \"\${job}\" -o JobIDRaw,State,ExitCode | awk -F'|' -v p=\"\${job}_\" '\$1 ~ (\"^\" p \"[0-9]+$\") {print \$2 \"|\" \$3}')
      [[ \${#states[@]} -eq 112 ]]
      for state in \"\${states[@]}\"; do [[ \"\${state}\" == 'COMPLETED|0:0' ]]; done
      shopt -s nullglob
      errs=('${TWIRL_ROOT}/logs/twirl-ssl-v3-numeric-audit-'\"\${job}\"_*.err)
      [[ \${#errs[@]} -eq 112 ]]
      for err in \"\${errs[@]}\"; do test ! -s \"\${err}\"; done
      reports=0
      for sector in {56..62}; do
        for shard in {0..15}; do
          report='${NUMERIC_SHARD_DIR}/s'\${sector}'/native_'\${shard}'.numeric.json'
          test -s \"\${report}\"
          test -s \"\${report}.sha256\"
          reports=\$((reports + 1))
        done
      done
      [[ \${reports} -eq 112 ]]
      for output in '${NUMERIC_AUDIT}' '${NUMERIC_AUDIT}.sha256' '${NUMERIC_RELEASE}' '${NUMERIC_RELEASE}.sha256'; do [[ ! -e \"\${output}\" ]]; done
      [[ ! -e '${NUMERIC_RELEASE_JOB}' && ! -e '${NUMERIC_RELEASE_JOB}.tmp' ]]
      mkdir -- '${CLAIM_PREFIX}.submit-numeric-release.claim'
      cd '${REMOTE_REPO}'
      job=\$(sbatch --parsable --export='$(base_export)' \
        scripts/orcd/slurm_teacher_ssl_fullpool_v3_numeric_release_cpu.sbatch)
      job=\${job%%;*}
      printf '%s\\n' \"\${job}\" > '${NUMERIC_RELEASE_JOB}.tmp'
      mv '${NUMERIC_RELEASE_JOB}.tmp' '${NUMERIC_RELEASE_JOB}'
      printf 'numeric_release=%s\\n' \"\${job}\"
    "
    ;;

  submit-smoke)
    require_socket
    remote "
      set -euo pipefail
      test -s '${NUMERIC_RELEASE_JOB}'
      job=\$(< '${NUMERIC_RELEASE_JOB}')
      state=\$(sacct -X -n -P -j \"\${job}\" -o JobIDRaw,State,ExitCode | awk -F'|' -v j=\"\${job}\" '\$1==j {print \$2 \"|\" \$3}')
      [[ \"\${state}\" == 'COMPLETED|0:0' ]]
      test -e '${TWIRL_ROOT}/logs/twirl-ssl-v3-numeric-release-'\"\${job}\"'.err'
      test ! -s '${TWIRL_ROOT}/logs/twirl-ssl-v3-numeric-release-'\"\${job}\"'.err'
      test -s '${NUMERIC_RELEASE}'
      test -s '${NUMERIC_RELEASE}.sha256'
      test ! -e '${MODEL_RUN_ROOT}'
      [[ ! -e '${SMOKE_JOB}' && ! -e '${SMOKE_JOB}.tmp' ]]
      PYTHONPATH='${REMOTE_REPO}/src' '${TWIRL_ROOT}/envs/twirl-s56/bin/python' -c 'from pathlib import Path; import sys; from twirl.vetting.ssl_full_pool_numeric import validate_numeric_gate_release; validate_numeric_gate_release(Path(sys.argv[1]), expected_code_revision=sys.argv[2])' '${NUMERIC_RELEASE}' '${EXPECTED_SHA}'
      mkdir -- '${CLAIM_PREFIX}.submit-smoke.claim'
      cd '${REMOTE_REPO}'
      job=\$(sbatch --parsable --export='$(base_export)' \
        scripts/orcd/slurm_teacher_ssl_fullpool_v3_smoke_h200.sbatch)
      job=\${job%%;*}
      printf '%s\\n' \"\${job}\" > '${SMOKE_JOB}.tmp'
      mv '${SMOKE_JOB}.tmp' '${SMOKE_JOB}'
      printf 'smoke=%s\\n' \"\${job}\"
    "
    ;;

  submit-folds)
    require_socket
    remote "
      set -euo pipefail
      test -s '${SMOKE_JOB}'
      job=\$(< '${SMOKE_JOB}')
      state=\$(sacct -X -n -P -j \"\${job}\" -o JobIDRaw,State,ExitCode | awk -F'|' -v j=\"\${job}\" '\$1==j {print \$2 \"|\" \$3}')
      [[ \"\${state}\" == 'COMPLETED|0:0' ]]
      test -e '${TWIRL_ROOT}/logs/twirl-ssl-v3-smoke-'\"\${job}\"'.err'
      test ! -s '${TWIRL_ROOT}/logs/twirl-ssl-v3-smoke-'\"\${job}\"'.err'
      test -s '${SMOKE_SUMMARY}'
      test ! -e '${MODEL_RUN_ROOT}/training/five_fold'
      [[ ! -e '${FOLDS_JOB}' && ! -e '${FOLDS_JOB}.tmp' ]]
      cd '${REMOTE_REPO}'
      '${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python' \
        scripts/stage5_validation/validate_teacher_ssl_fullpool_v3_smoke.py \
        --summary '${SMOKE_SUMMARY}' \
        --expected-code-revision '${EXPECTED_SHA}' \
        --eligibility-exclusions '${ELIGIBILITY_EXCLUSIONS}' \
        --eligibility-summary '${ELIGIBILITY_SUMMARY}' \
        --native-registry '${NATIVE_REGISTRY}' \
        --native-registry-summary '${NATIVE_REGISTRY_SUMMARY}' \
        --native-release-summary '${NATIVE_RELEASE_SUMMARY}' \
        --registry '${SSL_REGISTRY}' \
        --registry-summary '${SSL_REGISTRY_SUMMARY}' \
        --numeric-gate-release '${NUMERIC_RELEASE}' \
        --expected-fold 2
      mkdir -- '${CLAIM_PREFIX}.submit-folds.claim'
      job=\$(sbatch --parsable --array=0-4%4 --export='$(base_export)' \
        scripts/orcd/slurm_teacher_ssl_fullpool_v3_fold_h200.sbatch)
      job=\${job%%;*}
      printf '%s\\n' \"\${job}\" > '${FOLDS_JOB}.tmp'
      mv '${FOLDS_JOB}.tmp' '${FOLDS_JOB}'
      printf 'folds=%s\\n' \"\${job}\"
    "
    ;;

  validate-folds)
    require_socket
    remote "
      set -euo pipefail
      test -s '${FOLDS_JOB}'
      job=\$(< '${FOLDS_JOB}')
      declare -A fold_states=()
      while IFS='|' read -r job_id state exit_code; do
        [[ -n \"\${job_id}\" ]]
        fold=\${job_id#\${job}_}
        [[ \"\${fold}\" =~ ^[0-4]$ ]]
        [[ -z \"\${fold_states[\${fold}]+x}\" ]]
        fold_states[\${fold}]=\"\${state}|\${exit_code}\"
      done < <(sacct -n -P -j \"\${job}\" -o JobIDRaw,State,ExitCode | awk -F'|' -v p=\"\${job}_\" '\$1 ~ (\"^\" p \"[0-4]$\") {print \$1 \"|\" \$2 \"|\" \$3}')
      [[ \${#fold_states[@]} -eq 5 ]]
      for fold in {0..4}; do
        [[ \"\${fold_states[\${fold}]}\" == 'COMPLETED|0:0' ]]
      done
      shopt -s nullglob
      errs=('${TWIRL_ROOT}/logs/twirl-ssl-v3-fold-'\"\${job}\"_*.err)
      [[ \${#errs[@]} -eq 5 ]]
      for err in \"\${errs[@]}\"; do test ! -s \"\${err}\"; done
      [[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]
      [[ -z \$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all) ]]
      for output in '${FOLD_COMPLETION_RELEASE}' '${FOLD_COMPLETION_RELEASE}.sha256'; do
        [[ ! -e \"\${output}\" ]]
      done
      [[ ! -e '${FOLD_VALIDATION_JOB}' && ! -e '${FOLD_VALIDATION_JOB}.tmp' ]]
      mkdir -- '${CLAIM_PREFIX}.validate-folds.claim'
      cd '${REMOTE_REPO}'
      validation_job=\$(sbatch --parsable --export='$(base_export)' \
        scripts/orcd/slurm_teacher_ssl_fullpool_v3_validate_folds_cpu.sbatch)
      validation_job=\${validation_job%%;*}
      printf '%s\\n' \"\${validation_job}\" > '${FOLD_VALIDATION_JOB}.tmp'
      mv '${FOLD_VALIDATION_JOB}.tmp' '${FOLD_VALIDATION_JOB}'
      printf 'fold_validation=%s\\n' \"\${validation_job}\"
    "
    ;;

  verify-completion)
    require_socket
    remote "
      set -euo pipefail
      test -s '${FOLD_VALIDATION_JOB}'
      job=\$(< '${FOLD_VALIDATION_JOB}')
      state=\$(sacct -X -n -P -j \"\${job}\" -o JobIDRaw,State,ExitCode | awk -F'|' -v j=\"\${job}\" '\$1==j {print \$2 \"|\" \$3}')
      [[ \"\${state}\" == 'COMPLETED|0:0' ]]
      test -e '${TWIRL_ROOT}/logs/twirl-ssl-v3-validate-folds-'\"\${job}\"'.err'
      test ! -s '${TWIRL_ROOT}/logs/twirl-ssl-v3-validate-folds-'\"\${job}\"'.err'
      test -s '${FOLD_COMPLETION_RELEASE}'
      test -s '${FOLD_COMPLETION_RELEASE}.sha256'
      (
        cd '${RUN_ROOT}/frozen/model'
        sha256sum -c 'teacher_ssl_fullpool_v3_five_fold.complete.json.sha256'
      )
      [[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]
      [[ -z \$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all) ]]
      before=\$(stat -c '%s:%Y:%i' '${FOLD_COMPLETION_RELEASE}' '${FOLD_COMPLETION_RELEASE}.sha256' | sha256sum | awk '{print \$1}')
      cd '${REMOTE_REPO}'
      '${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python' \
        scripts/stage5_validation/validate_teacher_ssl_fullpool_v3_training.py \
        --model-root '${MODEL_RUN_ROOT}/training/five_fold' \
        --expected-code-revision '${EXPECTED_SHA}' \
        --output-release '${FOLD_COMPLETION_RELEASE}' >/dev/null
      after=\$(stat -c '%s:%Y:%i' '${FOLD_COMPLETION_RELEASE}' '${FOLD_COMPLETION_RELEASE}.sha256' | sha256sum | awk '{print \$1}')
      [[ \"\${after}\" == \"\${before}\" ]]
      (
        cd '${RUN_ROOT}/frozen/model'
        sha256sum -c 'teacher_ssl_fullpool_v3_five_fold.complete.json.sha256'
      )
      '${TWIRL_ROOT}/envs/twirl-s56/bin/python' -c 'import json,sys; d=json.load(open(sys.argv[1])); assert d[\"passed\"] is True; assert d[\"schema_version\"]==\"twirl_teacher_ssl_fullpool_v3_five_fold_completion_release_v1\"; assert d[\"release_binding\"]==\"teacher_ssl_fullpool_v3_effective_quality_adp_five_fold_complete_v1\"; assert d[\"expected_code_revision\"]==sys.argv[2]; assert d[\"model_namespace\"]==\"effective_quality_adp_v1\"; assert d[\"training_hyperparameters\"][\"folds\"]==[0,1,2,3,4]; assert d[\"counts\"]=={\"folds\":5,\"completed_epochs\":100}; assert [f[\"fold\"] for f in d[\"folds\"]]==[0,1,2,3,4]' '${FOLD_COMPLETION_RELEASE}' '${EXPECTED_SHA}'
      printf 'verified_completion=%s\\n' '${FOLD_COMPLETION_RELEASE}'
    "
    ;;

  status)
    require_socket
    remote "
      squeue -u tehan -o '%.18i %.36j %.2t %.10M %.4D %R'
      find '${RUN_ROOT}' -maxdepth 8 -type f \\( -name '*summary.json' -o -name '*.sha256' -o -name '*.release.json' -o -name '*.complete.json' -o -name '*.tsv' \\) 2>/dev/null | sort
    "
    ;;

  *)
    usage >&2
    exit 2
    ;;
esac
