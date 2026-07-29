#!/usr/bin/env bash
# Guarded local controller for the approved Teacher v4-SSL full-pool v2 flow.
#
# Every submit action is an explicit gate.  This controller intentionally does
# not construct a dependency chain: a later action is available only after the
# previous stage has published and validated its required release artifacts.
set -euo pipefail

DRY_RUN=1
ACTION=""

usage() {
  cat <<'EOF'
Usage: run_teacher_ssl_fullpool_v2_orcd.sh [--run] ACTION

Actions:
  probe                   Verify the user-opened ORCD control socket.
  deploy                  Create/reuse a clean detached worktree at the exact Git SHA.
  submit-eligibility      Freeze the 175,347/19 native-model eligibility authority.
  submit-native-canary    Build isolated S56 shard 3 (eligible plus excluded inputs).
  submit-native           Submit seven independent 16-shard native-v2 CPU arrays.
  submit-native-registry  Freeze exact 112-shard native coverage and release summaries.
  submit-ssl-registry     Build the 175,366-row audit / 175,347-row include registry.
  submit-numeric-audit    Audit all 112 native shards under the model-facing transform.
  submit-numeric-release  Merge 112 passed audits into one immutable numeric release.
  submit-smoke            Run one bounded one-epoch smoke on one H200.
  submit-folds            Launch five one-H200 folds as array 0-4%4.
  status                  Show v2 jobs and published gate artifacts.

Default mode is dry-run.  No action starts password, Duo, or
keyboard-interactive authentication.  The native wrapper resolves each compact
ADP path from the checksum-bound frozen-pool summary and uses the already
staged seven-sector cadence-authority root, validating every byte before use.
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

LOCAL_REPO="${LOCAL_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
ORCD_HOST="${ORCD_HOST:-tehan@orcd-login.mit.edu}"
CONTROL_PATH="${ORCD_CONTROL_PATH:-${HOME}/.ssh/cm/%r@%h:%p}"
EXPECTED_TWIRL_ROOT="/orcd/data/mki_aryeh/001/twirl"
TWIRL_ROOT="${TWIRL_ORCD_ROOT:-${EXPECTED_TWIRL_ROOT}}"
if [[ "${TWIRL_ROOT}" != "${EXPECTED_TWIRL_ROOT}" ]]; then
  echo "TWIRL_ORCD_ROOT differs from the fixed Teacher v4-SSL root." >&2
  exit 2
fi
REMOTE_SOURCE="${TWIRL_ORCD_SOURCE_REPO:-${TWIRL_ROOT}/code/TWIRL}"
EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA:-$(git -C "${LOCAL_REPO}" rev-parse HEAD)}"
REMOTE_REPO="${TWIRL_ORCD_REPO:-${TWIRL_ROOT}/code/TWIRL_teacher_ssl_fullpool_v2_${EXPECTED_SHA:0:12}}"
GIT_BRANCH="${TWIRL_GIT_BRANCH:-codex/teacher-ssl-quality-mask}"
SOURCE_RUN_ROOT="${TWIRL_SSL_FULLPOOL_V1_RUN_ROOT:-${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp}"
RUN_ROOT="${TWIRL_SSL_FULLPOOL_V2_RUN_ROOT:-${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v2_s56_s62_a2v1_current_adp_bls_eligible}"
ELIGIBILITY_DIR="${RUN_ROOT}/frozen/native_model_eligibility"
ELIGIBILITY_EXCLUSIONS="${ELIGIBILITY_DIR}/native_model_exclusions.csv"
ELIGIBILITY_SUMMARY="${ELIGIBILITY_DIR}/native_model_eligibility.summary.json"
CANARY_ROOT="${RUN_ROOT}/native_canary"
CANARY_SUMMARY="${CANARY_ROOT}/s56/native_3.summary.json"
NATIVE_REGISTRY_DIR="${RUN_ROOT}/frozen/native_registry"
NATIVE_REGISTRY="${NATIVE_REGISTRY_DIR}/native_input_registry.parquet"
NATIVE_REGISTRY_SUMMARY="${NATIVE_REGISTRY_DIR}/native_input_registry.summary.json"
NATIVE_RELEASE_SUMMARY="${NATIVE_REGISTRY_DIR}/native_fullpool_release.summary.json"
SSL_REGISTRY_DIR="${RUN_ROOT}/frozen/ssl_registry"
SSL_REGISTRY="${SSL_REGISTRY_DIR}/teacher_ssl_fullpool_registry.parquet"
SSL_REGISTRY_SUMMARY="${SSL_REGISTRY_DIR}/teacher_ssl_fullpool_registry.summary.json"
NUMERIC_GATE_DIR="${RUN_ROOT}/frozen/model_input_numeric_gate_v1"
NUMERIC_SHARD_DIR="${NUMERIC_GATE_DIR}/shard_reports"
NUMERIC_AUDIT="${NUMERIC_GATE_DIR}/model_input_numeric_audit.parquet"
NUMERIC_RELEASE="${NUMERIC_GATE_DIR}/model_input_numeric_gate.release.json"
MODEL_RUN_ROOT="${RUN_ROOT}/model_runs/effective_quality_mask_v1"
SMOKE_ROOT="${MODEL_RUN_ROOT}/smoke/one_epoch"
readonly SMOKE_FOLD=2
SMOKE_SUMMARY="${SMOKE_ROOT}/encoder_pretraining/fold_${SMOKE_FOLD}/summary.json"
NATIVE_SHARDS=16
FULL_POOL_IDENTITY_SHA256="8e9e9c12a24d5ebc7be94b03a4e35411cd10066d62a87d921a8443b06cc188d1"
ELIGIBLE_IDENTITY_SHA256="6ddc8e57bb5fb938ce05389c1629221c73e0e73ac3bf40da47a2019e1a5660e6"
EXCLUDED_IDENTITY_SHA256="b9f536144265e54a70bff17c782e0668ddbd96efcdf6c223ebc58f46edb7d976"
CANARY_EXCLUDED_IDENTITY_SHA256="ddda9b053eddc744e5032ba350598d58d74cea4f4cc5cd705932ccf6e41022ab"

if [[ -n "${TWIRL_SSL_FULLPOOL_SMOKE_FOLD+x}" ]]; then
  echo "Controller forbids a smoke-fold environment override." >&2
  exit 2
fi

EXPECTED_SOURCE_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp"
EXPECTED_RUN_ROOT="${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v2_s56_s62_a2v1_current_adp_bls_eligible"
if [[ "${SOURCE_RUN_ROOT}" != "${EXPECTED_SOURCE_RUN_ROOT}" || "${RUN_ROOT}" != "${EXPECTED_RUN_ROOT}" ]]; then
  echo "Controller roots differ from the fixed Teacher v4-SSL v2 contract." >&2
  exit 2
fi

SSH=(
  ssh
  -o BatchMode=yes
  -o PasswordAuthentication=no
  -o KbdInteractiveAuthentication=no
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
The documented ORCD control socket is not available at:
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
    ",TWIRL_SSL_FULLPOOL_V2_RUN_ROOT=${RUN_ROOT}" \
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
      source_repo='${REMOTE_SOURCE}'
      target_repo='${REMOTE_REPO}'
      expected_sha='${EXPECTED_SHA}'
      test -d \"\${source_repo}/.git\"
      GIT_TERMINAL_PROMPT=0 git -C \"\${source_repo}\" fetch origin '${GIT_BRANCH}'
      git -C \"\${source_repo}\" cat-file -e \"\${expected_sha}^{commit}\"
      if [[ -e \"\${target_repo}\" ]]; then
        test -e \"\${target_repo}/.git\"
        [[ \$(git -C \"\${target_repo}\" rev-parse HEAD) == \"\${expected_sha}\" ]]
      else
        git -C \"\${source_repo}\" worktree add --detach \"\${target_repo}\" \"\${expected_sha}\"
      fi
      [[ -z \$(git -C \"\${target_repo}\" status --porcelain=v1 --untracked-files=all) ]]
      git -C \"\${target_repo}\" rev-parse HEAD
    "
    ;;
  submit-eligibility)
    require_socket
    remote "
      set -euo pipefail
      repo='${REMOTE_REPO}'
      [[ \$(git -C \"\${repo}\" rev-parse HEAD) == '${EXPECTED_SHA}' ]]
      [[ -z \$(git -C \"\${repo}\" status --porcelain=v1 --untracked-files=all) ]]
      [[ ! -e '${ELIGIBILITY_EXCLUSIONS}' && ! -e '${ELIGIBILITY_SUMMARY}' ]]
      cd \"\${repo}\"
      sbatch --parsable --export='$(base_export)' \
        scripts/orcd/slurm_teacher_ssl_fullpool_native_eligibility_cpu.sbatch
    "
    ;;
  submit-native-canary)
    require_socket
    # With the locked 16-way shard function, S56 shard 3 contains excluded
    # TIC 2019898202 plus ordinary eligible rows.  The isolated canary thus
    # exercises subtraction without publishing into the production tree.
    remote "
      set -euo pipefail
      repo='${REMOTE_REPO}'
      [[ \$(git -C \"\${repo}\" rev-parse HEAD) == '${EXPECTED_SHA}' ]]
      [[ -z \$(git -C \"\${repo}\" status --porcelain=v1 --untracked-files=all) ]]
      test -s '${ELIGIBILITY_EXCLUSIONS}'
      test -s '${ELIGIBILITY_SUMMARY}'
      '${TWIRL_ROOT}/envs/twirl-s56/bin/python' -c 'import json,sys; d=json.load(open(sys.argv[1])); c=d[\"counts\"]; h=d[\"identity_hashes\"]; a=d[\"partition_audit\"]; assert d.get(\"passed\") is True; assert c[\"full\"][\"n_observations\"] == 175366; assert c[\"eligible\"][\"n_observations\"] == 175347; assert c[\"excluded\"][\"n_observations\"] == 19; assert h[\"full_observations_sha256\"] == sys.argv[2]; assert h[\"eligible_observations_sha256\"] == sys.argv[3]; assert h[\"excluded_observations_sha256\"] == sys.argv[4]; assert a[\"eligible_excluded_disjoint\"] is True; assert a[\"eligible_excluded_union_equals_full\"] is True; assert a[\"production_lock_passed\"] is True' '${ELIGIBILITY_SUMMARY}' '${FULL_POOL_IDENTITY_SHA256}' '${ELIGIBLE_IDENTITY_SHA256}' '${EXCLUDED_IDENTITY_SHA256}'
      test ! -e '${CANARY_ROOT}/s56/native_3.h5'
      test ! -e '${CANARY_SUMMARY}'
      cd \"\${repo}\"
      sbatch --parsable --export='$(base_export),TWIRL_FULLPOOL_SECTOR=56,TWIRL_FULLPOOL_NATIVE_SHARDS=${NATIVE_SHARDS},TWIRL_FULLPOOL_SHARD_INDEX=3,TWIRL_FULLPOOL_NATIVE_OUT=${CANARY_ROOT}/s56/native_3.h5' \
        scripts/orcd/slurm_teacher_ssl_fullpool_native_cpu.sbatch
    "
    ;;
  submit-native)
    require_socket
    if [[ "${DRY_RUN}" == "1" ]]; then
      for sector in {56..62}; do
        remote "cd '${REMOTE_REPO}' && sbatch --array=0-15%4 --export='$(base_export),TWIRL_FULLPOOL_SECTOR=${sector},TWIRL_FULLPOOL_NATIVE_SHARDS=${NATIVE_SHARDS}' scripts/orcd/slurm_teacher_ssl_fullpool_native_cpu.sbatch"
      done
      exit 0
    fi
    "${SSH[@]}" "${ORCD_HOST}" "
      set -euo pipefail
      test -s '${CANARY_SUMMARY}'
      '${TWIRL_ROOT}/envs/twirl-s56/bin/python' -c 'import json,sys; d=json.load(open(sys.argv[1])); v=d.get(\"verification\", {}); assert d.get(\"sector\") == 56; assert d.get(\"shard_index\") == 3; assert d.get(\"n_shards\") == 16; assert d.get(\"native_contract_version\") == \"twirl_teacher_ssl_fullpool_real_native_v2\"; assert d.get(\"n_shard_excluded_observations\") == 1; assert d.get(\"shard_excluded_observation_identity_sha256\") == sys.argv[2]; assert d.get(\"n_source_shard_observations\") == d.get(\"n_shard_observations\") + 1; assert v.get(\"passed\") is True; assert v.get(\"counts\", {}).get(\"targets\") == d.get(\"n_shard_observations\")' '${CANARY_SUMMARY}' '${CANARY_EXCLUDED_IDENTITY_SHA256}'
      [[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]
      [[ -z \$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all) ]]
    "
    for sector in {56..62}; do
      job="$("${SSH[@]}" "${ORCD_HOST}" "
        set -euo pipefail
        cd '${REMOTE_REPO}'
        sbatch --parsable --array=0-15%4 --export='$(base_export),TWIRL_FULLPOOL_SECTOR=${sector},TWIRL_FULLPOOL_NATIVE_SHARDS=${NATIVE_SHARDS}' \
          scripts/orcd/slurm_teacher_ssl_fullpool_native_cpu.sbatch
      ")"
      printf 'S%s_native=%s\n' "${sector}" "${job}"
    done
    ;;
  submit-native-registry)
    require_socket
    remote "
      set -euo pipefail
      [[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]
      [[ -z \$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all) ]]
      for sector in {56..62}; do
        for shard in {0..15}; do
          test -s '${RUN_ROOT}/native/s'\${sector}'/native_'\${shard}'.h5'
          test -s '${RUN_ROOT}/native/s'\${sector}'/native_'\${shard}'.summary.json'
        done
      done
      cd '${REMOTE_REPO}'
      sbatch --parsable --export='$(base_export),TWIRL_FULLPOOL_NATIVE_SHARDS=${NATIVE_SHARDS}' \
        scripts/orcd/slurm_teacher_ssl_fullpool_native_registry_cpu.sbatch
    "
    ;;
  submit-ssl-registry)
    require_socket
    remote "
      set -euo pipefail
      test -s '${NATIVE_REGISTRY}'
      test -s '${NATIVE_REGISTRY_SUMMARY}'
      test -s '${NATIVE_RELEASE_SUMMARY}'
      [[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]
      [[ -z \$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all) ]]
      cd '${REMOTE_REPO}'
      sbatch --parsable --export='$(base_export)' \
        scripts/orcd/slurm_teacher_ssl_fullpool_registry_cpu.sbatch
    "
    ;;
  submit-numeric-audit)
    require_socket
    remote "
      set -euo pipefail
      test -s '${ELIGIBILITY_EXCLUSIONS}'
      test -s '${ELIGIBILITY_SUMMARY}'
      test -s '${NATIVE_REGISTRY}'
      test -s '${NATIVE_REGISTRY_SUMMARY}'
      test -s '${NATIVE_RELEASE_SUMMARY}'
      test -s '${SSL_REGISTRY}'
      test -s '${SSL_REGISTRY_SUMMARY}'
      for sector in {56..62}; do
        for shard in {0..15}; do
          report='${NUMERIC_SHARD_DIR}/s'\${sector}'/native_'\${shard}'.numeric.json'
          if [[ -e "\${report}.sha256" ]]; then
            test -s "\${report}"
          fi
        done
      done
      test ! -e '${NUMERIC_AUDIT}'
      test ! -e '${NUMERIC_AUDIT}.sha256'
      test ! -e '${NUMERIC_RELEASE}'
      test ! -e '${NUMERIC_RELEASE}.sha256'
      [[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]
      [[ -z \$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all) ]]
      cd '${REMOTE_REPO}'
      sbatch --parsable --export='$(base_export)' \
        scripts/orcd/slurm_teacher_ssl_fullpool_numeric_audit_cpu.sbatch
    "
    ;;
  submit-numeric-release)
    require_socket
    remote "
      set -euo pipefail
      for sector in {56..62}; do
        for shard in {0..15}; do
          test -s '${NUMERIC_SHARD_DIR}/s'\${sector}'/native_'\${shard}'.numeric.json'
          test -s '${NUMERIC_SHARD_DIR}/s'\${sector}'/native_'\${shard}'.numeric.json.sha256'
        done
      done
      [[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]
      [[ -z \$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all) ]]
      cd '${REMOTE_REPO}'
      sbatch --parsable --export='$(base_export)' \
        scripts/orcd/slurm_teacher_ssl_fullpool_numeric_release_cpu.sbatch
    "
    ;;
  submit-smoke)
    require_socket
    remote "
      set -euo pipefail
      test -s '${ELIGIBILITY_EXCLUSIONS}'
      test -s '${ELIGIBILITY_SUMMARY}'
      test -s '${NATIVE_REGISTRY}'
      test -s '${NATIVE_REGISTRY_SUMMARY}'
      test -s '${NATIVE_RELEASE_SUMMARY}'
      test -s '${SSL_REGISTRY}'
      test -s '${SSL_REGISTRY_SUMMARY}'
      test -s '${NUMERIC_RELEASE}'
      PYTHONPATH='${REMOTE_REPO}/src' '${TWIRL_ROOT}/envs/twirl-s56/bin/python' -c 'from pathlib import Path; import sys; from twirl.vetting.ssl_full_pool_numeric import validate_numeric_gate_release; validate_numeric_gate_release(Path(sys.argv[1]), expected_code_revision=sys.argv[2])' '${NUMERIC_RELEASE}' '${EXPECTED_SHA}'
      [[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]
      [[ -z \$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all) ]]
      cd '${REMOTE_REPO}'
      sbatch --parsable --export='$(base_export)' \
        scripts/orcd/slurm_teacher_ssl_fullpool_smoke_h200.sbatch
    "
    ;;
  submit-folds)
    require_socket
    remote "
      set -euo pipefail
      test -s '${SMOKE_SUMMARY}'
      test -s '${ELIGIBILITY_EXCLUSIONS}'
      test -s '${ELIGIBILITY_SUMMARY}'
      test -s '${NATIVE_REGISTRY}'
      test -s '${NATIVE_REGISTRY_SUMMARY}'
      test -s '${NATIVE_RELEASE_SUMMARY}'
      test -s '${SSL_REGISTRY}'
      test -s '${SSL_REGISTRY_SUMMARY}'
      test -s '${NUMERIC_RELEASE}'
      PYTHONPATH='${REMOTE_REPO}/src' '${TWIRL_ROOT}/envs/twirl-s56/bin/python' -c 'from pathlib import Path; import sys; from twirl.vetting.ssl_full_pool_numeric import validate_numeric_gate_release; validate_numeric_gate_release(Path(sys.argv[1]), expected_code_revision=sys.argv[2])' '${NUMERIC_RELEASE}' '${EXPECTED_SHA}'
      [[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]
      [[ -z \$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all) ]]
      cd '${REMOTE_REPO}'
      '${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python' \
        scripts/stage5_validation/validate_teacher_ssl_fullpool_v2_smoke.py \
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
        --expected-fold '${SMOKE_FOLD}'
      sbatch --parsable --export='$(base_export)' \
        scripts/orcd/slurm_teacher_ssl_fullpool_fold_h200.sbatch
    "
    ;;
  status)
    require_socket
    remote "
      squeue -u tehan -o '%.18i %.32j %.2t %.10M %.4D %R'
      find '${RUN_ROOT}' -maxdepth 8 -type f \\( -name '*summary.json' -o -name '*.sha256' -o -name '*.release.json' \\) 2>/dev/null | sort
    "
    ;;
  *)
    usage >&2
    exit 2
    ;;
esac
