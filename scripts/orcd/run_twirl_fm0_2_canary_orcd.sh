#!/usr/bin/env bash
# Fail-closed controller for the frozen TWIRL-FM0.2.1 objective canary.
# Dry-run by default; never opens an interactive SSH session.
set -euo pipefail
umask 027

DRY_RUN=1
ACTION=""
readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly PARTITION="pg_mki_aryeh"
readonly RECLAIMABLE_PARTITION="mit_preemptable"
readonly GPU_NODE="node4900"
readonly LOCKED_HOST="tehan@orcd-login.mit.edu"
readonly CONTROL_PATH="/Users/tehan/.ssh/cm/%r@%h:%p"
readonly CONFIG_SHA="394a0340f9f8518423e39a1c83bc32acde4e2ebff95c4837b0edb5399bd50801"
readonly FREEZE_SHA="900264edabe13aa77dfffe3c87f5fe0c4967cbcf693b9a68020ab29ff85e691e"
readonly REUSE_SHA="46faed954e6758dfea71f7754670611fab03af63613dad3be433af71dd38e641"
readonly EVENT_SHA="eb84f27866e42bbf942256abe8fa7a532eaa8584e6800ce012aa0552c08a6630"
readonly EVALUATOR_SHA="ed8004180b28c7e2bd85911d0128f9180eddd3c35bdaae299e2f7edba6e5baac"
readonly ORCD_CONFIG_SHA="a100514e8c83a4d17bab7924a5733015fddbf2a97cadf96cafcfab22cf92cc61"
readonly INPUT_RECEIPT="${PROJECT_ROOT}/reports/stage5_validation/twirl_fm0_1_s56_s67_poc/ece8619fd72f/frozen/input_release_validation/input_release.receipt.json"
readonly INPUT_RECEIPT_SHA="a0908e43e92cae0f3e832382e986319e5942cce12b90bf3f90db17a03eb792f7"

usage() {
  cat <<'EOF'
Usage: TWIRL_EXPECTED_GIT_SHA=<40-hex> run_twirl_fm0_2_canary_orcd.sh [--execute] ACTION

Actions:
  probe                    Read the existing socket, queue, and H200 node state.
  deploy                   Create/reuse a clean detached exact-SHA worktree.
  submit-loader-restart    Submit the 4-CPU/16-GiB exact-code gate; no GPU.
  submit-fp32-smoke        Fresh admission + one-H200 eight-step FP32 smoke.
  submit-post-validation   CPU-validate one exact smoke/canary milestone.
  submit-canary            Fresh admission + one-H200 FM0.2.1 milestone.
  status                   Read this exact-SHA run root and matching jobs.

submit-fp32-smoke and submit-canary require:
  TWIRL_FM0_LOADER_RESTART_RECEIPT and _SHA256

submit-canary also requires:
  TWIRL_FM0_STOP_AFTER_STEP (64, 500, 1000, or 2000)
  TWIRL_FM0_SMOKE_POST_VALIDATION and _SHA256
For stops after 64 it additionally requires the immediately preceding
milestone checkpoint and post-validation receipt bindings used by the wrapper.

submit-post-validation requires TWIRL_FM0_RUN_KIND, TWIRL_FM0_STOP_AFTER_STEP,
TWIRL_FM0_RUN_OUTPUT, and TWIRL_FM0_ADMISSION_RECEIPT plus its SHA256.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --execute) DRY_RUN=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) [[ -z "${ACTION}" ]] || { echo "Only one action is allowed." >&2; exit 2; }; ACTION="$1"; shift ;;
  esac
done
[[ -n "${ACTION}" ]] || { usage >&2; exit 2; }
: "${TWIRL_EXPECTED_GIT_SHA:?set the exact reviewed 40-hex commit}"
readonly EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}"
[[ "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]] || { echo "Expected SHA must be full lowercase hex." >&2; exit 2; }

readonly LOCAL_REPO="${LOCAL_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
readonly REMOTE_SOURCE="${PROJECT_ROOT}/code/TWIRL"
readonly REMOTE_REPO="${PROJECT_ROOT}/code/TWIRL_fm0_2_${EXPECTED_SHA:0:12}"
readonly RUN_ROOT="${PROJECT_ROOT}/reports/stage5_validation/twirl_fm0_2_s56_s64_objective_canary/${EXPECTED_SHA:0:12}"
readonly LAUNCH_DIR="${RUN_ROOT}/launch"
readonly REUSE_RECEIPT="${REMOTE_REPO}/reports/stage5_validation/twirl_fm0_2_design_freeze_v1/input_reuse.receipt.json"
readonly ENV_PREFIX="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1"

SSH=(
  ssh -o BatchMode=yes -o PasswordAuthentication=no
  -o KbdInteractiveAuthentication=no -o PubkeyAuthentication=no
  -o HostbasedAuthentication=no -o GSSAPIAuthentication=no
  -o NumberOfPasswordPrompts=0 -o ConnectTimeout=15
  -o ConnectionAttempts=1 -o ControlMaster=no -o "ControlPath=${CONTROL_PATH}"
)
print_cmd() { printf '+'; printf ' %q' "$@"; printf '\n'; }
remote() { if [[ "${DRY_RUN}" == 1 ]]; then print_cmd "${SSH[@]}" "${LOCKED_HOST}" "$1"; else "${SSH[@]}" "${LOCKED_HOST}" "$1"; fi; }
capture() { "${SSH[@]}" "${LOCKED_HOST}" "$1"; }
require_socket() {
  [[ "${DRY_RUN}" == 1 ]] && return
  "${SSH[@]}" -O check "${LOCKED_HOST}" >/dev/null 2>&1 || {
    echo "ORCD control socket is absent/stale; open it yourself." >&2; exit 4;
  }
}
safe_remote_path() { [[ "$1" == "${PROJECT_ROOT}"/* && "$1" != *".."* && "$1" != *","* && "$1" != *$'\n'* ]]; }
require_hash() { [[ "$2" =~ ^[0-9a-f]{${1}}$ ]] || { echo "Invalid ${1}-hex hash." >&2; exit 2; }; }
base_export() { printf '%s' "ALL,TWIRL_ORCD_REPO=${REMOTE_REPO},TWIRL_EXPECTED_GIT_SHA=${EXPECTED_SHA},TWIRL_FM0_RUN_ROOT=${RUN_ROOT}"; }

verify_remote_repo="
  [[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]
  [[ -z \$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all) ]]
  [[ \$(sha256sum '${REMOTE_REPO}/configs/models/twirl_fm0_2_s56_s64_objective_canary.yaml' | awk '{print \$1}') == '${CONFIG_SHA}' ]]
  [[ \$(sha256sum '${REMOTE_REPO}/reports/stage5_validation/twirl_fm0_2_design_freeze_v1/freeze.json' | awk '{print \$1}') == '${FREEZE_SHA}' ]]
  [[ \$(sha256sum '${REUSE_RECEIPT}' | awk '{print \$1}') == '${REUSE_SHA}' ]]
  [[ \$(sha256sum '${REMOTE_REPO}/reports/stage5_validation/twirl_fm0_2_design_freeze_v1/development_event_retention.contract.json' | awk '{print \$1}') == '${EVENT_SHA}' ]]
  [[ \$(sha256sum '${REMOTE_REPO}/reports/stage5_validation/twirl_fm0_2_design_freeze_v1/development_evaluator_v2.contract.json' | awk '{print \$1}') == '${EVALUATOR_SHA}' ]]
  [[ \$(sha256sum '${REMOTE_REPO}/configs/orcd/twirl_fm0_2_s56_s64_objective_canary.json' | awk '{print \$1}') == '${ORCD_CONFIG_SHA}' ]]
  [[ \$(sha256sum '${INPUT_RECEIPT}' | awk '{print \$1}') == '${INPUT_RECEIPT_SHA}' ]]
  test -x '${ENV_PREFIX}/bin/python'
"

submit_once() {
  local record=$1 name=$2 wrapper=$3 exports=$4
  remote "
    set -euo pipefail
    ${verify_remote_repo}
    mkdir -p '${LAUNCH_DIR}' '${PROJECT_ROOT}/logs'
    test ! -e '${record}'
    claim='${record}.claim'; mkdir -- \"\${claim}\"
    cleanup() { rmdir \"\${claim}\" 2>/dev/null || true; rm -f -- '${record}.tmp'; }
    trap cleanup EXIT
    cd '${REMOTE_REPO}'
    job=\$(sbatch --parsable --job-name='${name}' --export='${exports}' '${wrapper}')
    job=\${job%%;*}; [[ \${job} =~ ^[0-9]+$ ]]
    printf '%s\n' \"\${job}\" > '${record}.tmp'; mv -- '${record}.tmp' '${record}'
    trap - EXIT; rmdir \"\${claim}\"; cat '${record}'
  "
}

loader_binding() {
  : "${TWIRL_FM0_LOADER_RESTART_RECEIPT:?set loader/restart receipt}"
  : "${TWIRL_FM0_LOADER_RESTART_RECEIPT_SHA256:?set loader/restart receipt hash}"
  safe_remote_path "${TWIRL_FM0_LOADER_RESTART_RECEIPT}" || { echo "Unsafe loader receipt path." >&2; exit 2; }
  require_hash 64 "${TWIRL_FM0_LOADER_RESTART_RECEIPT_SHA256}"
}

make_admission_and_submit() {
  local kind=$1 stop=$2 wrapper=$3 name=$4 output=$5 extra_exports=$6
  loader_binding; require_socket
  if [[ "${DRY_RUN}" == 1 ]]; then
    echo "+ verify exact loader/restart receipt"
    echo "+ inspect all ${PARTITION} jobs and ${GPU_NODE}; reserve four Stage-1 lanes"
    echo "+ publish a <=300-second receipt and immediately submit one H200 ${kind} job"
    return
  fi
  capture "set -euo pipefail; ${verify_remote_repo} [[ \$(sha256sum '${TWIRL_FM0_LOADER_RESTART_RECEIPT}' | awk '{print \$1}') == '${TWIRL_FM0_LOADER_RESTART_RECEIPT_SHA256}' ]]; '${ENV_PREFIX}/bin/python' -c 'import json,sys; p=json.load(open(sys.argv[1])); assert p[\"schema_version\"]==\"twirl_fm0_2_loader_restart_receipt_v1\" and p[\"passed\"] is True and p[\"expected_git_sha\"]==sys.argv[2] and p[\"sealed_test_access_count\"]==0 and p[\"gpu_requested\"] is False' '${TWIRL_FM0_LOADER_RESTART_RECEIPT}' '${EXPECTED_SHA}'" >/dev/null

  local tmp remote_epoch local_epoch skew admission admission_sha record
  tmp=$(mktemp -d "${TMPDIR:-/tmp}/twirl-fm0-2-admission.XXXXXX")
  trap 'rm -rf -- "${tmp:-}"' EXIT
  remote_epoch=$(capture "date -u +%s"); local_epoch=$(date -u +%s)
  skew=$(( local_epoch > remote_epoch ? local_epoch - remote_epoch : remote_epoch - local_epoch ))
  (( skew <= 30 )) || { echo "Local/ORCD clock skew exceeds 30 seconds." >&2; exit 5; }
  capture "scontrol show partition -o '${PARTITION}'" > "${tmp}/partition.txt"
  capture "scontrol show partition -o '${RECLAIMABLE_PARTITION}'" > "${tmp}/reclaimable.txt"
  capture "scontrol show config | grep -E '^Preempt(Type|Mode|Parameters|ExemptTime)'" > "${tmp}/preemption.txt"
  capture "scontrol show node -o '${GPU_NODE}'" > "${tmp}/node.txt"
  capture "squeue -p '${PARTITION}' -h -t RUNNING,PENDING,CONFIGURING,COMPLETING -o '%i|%u|%T|%r|%C|%m|%b|%j|%N|%E'" > "${tmp}/queue.txt"
  capture "squeue -w '${GPU_NODE}' -h -t RUNNING,CONFIGURING,COMPLETING -o '%i|%P|%u|%T|%r|%C|%m|%b|%j|%N|%E'" > "${tmp}/node-queue.txt"
  git -C "${LOCAL_REPO}" diff --quiet "${EXPECTED_SHA}" -- scripts/orcd/build_twirl_fm0_2_admission_receipt.py || {
    echo "Local admission builder differs from the exact reviewed commit." >&2; exit 6;
  }
  python3 "${LOCAL_REPO}/scripts/orcd/build_twirl_fm0_2_admission_receipt.py" \
    --partition "${tmp}/partition.txt" --reclaimable-partition "${tmp}/reclaimable.txt" \
    --preemption "${tmp}/preemption.txt" --node "${tmp}/node.txt" \
    --queue "${tmp}/queue.txt" --node-queue "${tmp}/node-queue.txt" \
    --output "${tmp}/receipt.json" --now-epoch "${remote_epoch}" \
    --expected-git-sha "${EXPECTED_SHA}" --run-kind "${kind}" --stop-after-step "${stop}" \
    --input-receipt "${INPUT_RECEIPT}" --input-receipt-sha256 "${INPUT_RECEIPT_SHA}" \
    --input-reuse-receipt "${REUSE_RECEIPT}" --input-reuse-receipt-sha256 "${REUSE_SHA}" \
    --loader-restart-receipt "${TWIRL_FM0_LOADER_RESTART_RECEIPT}" \
    --loader-restart-receipt-sha256 "${TWIRL_FM0_LOADER_RESTART_RECEIPT_SHA256}"
  admission="${RUN_ROOT}/admission/admission-${kind}-step${stop}-${remote_epoch}.json"
  admission_sha=$(shasum -a 256 "${tmp}/receipt.json" | awk '{print $1}')
  "${SSH[@]}" "${LOCKED_HOST}" "set -euo pipefail; mkdir -p '${RUN_ROOT}/admission'; test ! -e '${admission}'; tmp='${admission}.tmp'; cat > \"\${tmp}\"; chmod 0444 \"\${tmp}\"; mv -- \"\${tmp}\" '${admission}'" < "${tmp}/receipt.json"
  record="${LAUNCH_DIR}/${kind}-step${stop}-${remote_epoch}.job"
  submit_once "${record}" "${name}" "${wrapper}" "$(base_export),TWIRL_FM0_INPUT_RECEIPT=${INPUT_RECEIPT},TWIRL_FM0_INPUT_RECEIPT_SHA256=${INPUT_RECEIPT_SHA},TWIRL_FM0_INPUT_REUSE_RECEIPT=${REUSE_RECEIPT},TWIRL_FM0_INPUT_REUSE_RECEIPT_SHA256=${REUSE_SHA},TWIRL_FM0_LOADER_RESTART_RECEIPT=${TWIRL_FM0_LOADER_RESTART_RECEIPT},TWIRL_FM0_LOADER_RESTART_RECEIPT_SHA256=${TWIRL_FM0_LOADER_RESTART_RECEIPT_SHA256},TWIRL_FM0_ADMISSION_RECEIPT=${admission},TWIRL_FM0_ADMISSION_RECEIPT_SHA256=${admission_sha},${extra_exports}"
  printf 'admission=%s\nadmission_sha256=%s\nrun_output=%s\n' "${admission}" "${admission_sha}" "${output}"
  rm -rf -- "${tmp}"; trap - EXIT
}

case "${ACTION}" in
  probe)
    require_socket
    remote "hostname; whoami; sinfo -h -p '${PARTITION}' -o '%P|%D|%t|%G'; scontrol show node -o '${GPU_NODE}'; squeue -w '${GPU_NODE}' -h -o '%i|%P|%u|%T|%C|%m|%b|%j|%N'"
    ;;
  deploy)
    require_socket
    remote "set -euo pipefail; test -d '${REMOTE_SOURCE}/.git'; git -C '${REMOTE_SOURCE}' cat-file -e '${EXPECTED_SHA}^{commit}'; if [[ -e '${REMOTE_REPO}' ]]; then [[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]; else git -C '${REMOTE_SOURCE}' worktree add --detach '${REMOTE_REPO}' '${EXPECTED_SHA}'; fi; if git -C '${REMOTE_REPO}' sparse-checkout list >/dev/null 2>&1; then git -C '${REMOTE_REPO}' sparse-checkout add reports/stage5_validation/twirl_fm0_2_design_freeze_v1 reports/stage5_validation/twirl_fm0_1_s56_s64_development_comparison_v1; fi; ${verify_remote_repo} git -C '${REMOTE_REPO}' rev-parse HEAD"
    ;;
  submit-loader-restart)
    require_socket
    submit_once "${LAUNCH_DIR}/loader-restart.job" "twirl-fm0-2-loader" \
      "scripts/orcd/slurm_twirl_fm0_2_loader_restart_cpu.sbatch" \
      "$(base_export),TWIRL_FM0_INPUT_RECEIPT=${INPUT_RECEIPT},TWIRL_FM0_INPUT_RECEIPT_SHA256=${INPUT_RECEIPT_SHA},TWIRL_FM0_INPUT_REUSE_RECEIPT=${REUSE_RECEIPT},TWIRL_FM0_INPUT_REUSE_RECEIPT_SHA256=${REUSE_SHA}"
    ;;
  submit-fp32-smoke)
    epoch=$(date -u +%s); output="${RUN_ROOT}/model_runs/fp32_synthetic_smoke/${epoch}"
    make_admission_and_submit synthetic_fp32 8 scripts/orcd/slurm_twirl_fm0_2_fp32_smoke_h200.sbatch twirl-fm0-2-fp32 "${output}" "TWIRL_FM0_SMOKE_OUTPUT=${output}"
    ;;
  submit-canary)
    : "${TWIRL_FM0_STOP_AFTER_STEP:?set authorized stop-after step}"
    case "${TWIRL_FM0_STOP_AFTER_STEP}" in 64|500|1000|2000) ;; *) echo "Unauthorized stop." >&2; exit 2 ;; esac
    : "${TWIRL_FM0_SMOKE_POST_VALIDATION:?set smoke post-validation receipt}"
    : "${TWIRL_FM0_SMOKE_POST_VALIDATION_SHA256:?set smoke post-validation hash}"
    safe_remote_path "${TWIRL_FM0_SMOKE_POST_VALIDATION}" || { echo "Unsafe smoke post-validation path." >&2; exit 2; }
    require_hash 64 "${TWIRL_FM0_SMOKE_POST_VALIDATION_SHA256}"
    output="${RUN_ROOT}/model_runs/fm0_2_1_canary/seed560067"
    extra="TWIRL_FM0_CANARY_OUTPUT=${output},TWIRL_FM0_STOP_AFTER_STEP=${TWIRL_FM0_STOP_AFTER_STEP},TWIRL_FM0_SMOKE_POST_VALIDATION=${TWIRL_FM0_SMOKE_POST_VALIDATION},TWIRL_FM0_SMOKE_POST_VALIDATION_SHA256=${TWIRL_FM0_SMOKE_POST_VALIDATION_SHA256}"
    if [[ "${TWIRL_FM0_STOP_AFTER_STEP}" != 64 ]]; then
      : "${TWIRL_FM0_RESUME_CHECKPOINT:?set preceding milestone checkpoint}"
      : "${TWIRL_FM0_PREVIOUS_POST_VALIDATION:?set preceding post-validation receipt}"
      : "${TWIRL_FM0_PREVIOUS_POST_VALIDATION_SHA256:?set preceding post-validation hash}"
      extra="${extra},TWIRL_FM0_RESUME_CHECKPOINT=${TWIRL_FM0_RESUME_CHECKPOINT},TWIRL_FM0_PREVIOUS_POST_VALIDATION=${TWIRL_FM0_PREVIOUS_POST_VALIDATION},TWIRL_FM0_PREVIOUS_POST_VALIDATION_SHA256=${TWIRL_FM0_PREVIOUS_POST_VALIDATION_SHA256}"
    fi
    make_admission_and_submit real_canary "${TWIRL_FM0_STOP_AFTER_STEP}" scripts/orcd/slurm_twirl_fm0_2_canary_h200.sbatch twirl-fm0-2-canary "${output}" "${extra}"
    ;;
  submit-post-validation)
    : "${TWIRL_FM0_RUN_KIND:?set synthetic_fp32 or real_canary}"
    : "${TWIRL_FM0_STOP_AFTER_STEP:?set exact stop-after step}"
    : "${TWIRL_FM0_RUN_OUTPUT:?set exact run output}"
    : "${TWIRL_FM0_ADMISSION_RECEIPT:?set exact admission receipt}"
    : "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256:?set admission receipt hash}"
    require_socket
    submit_once "${LAUNCH_DIR}/post-${TWIRL_FM0_RUN_KIND}-step${TWIRL_FM0_STOP_AFTER_STEP}-$(date -u +%s).job" "twirl-fm0-2-validate" \
      "scripts/orcd/slurm_twirl_fm0_2_post_validation_cpu.sbatch" \
      "$(base_export),TWIRL_FM0_RUN_GIT_SHA=${EXPECTED_SHA},TWIRL_FM0_RUN_KIND=${TWIRL_FM0_RUN_KIND},TWIRL_FM0_STOP_AFTER_STEP=${TWIRL_FM0_STOP_AFTER_STEP},TWIRL_FM0_RUN_OUTPUT=${TWIRL_FM0_RUN_OUTPUT},TWIRL_FM0_INPUT_RECEIPT=${INPUT_RECEIPT},TWIRL_FM0_INPUT_RECEIPT_SHA256=${INPUT_RECEIPT_SHA},TWIRL_FM0_INPUT_REUSE_RECEIPT=${REUSE_RECEIPT},TWIRL_FM0_INPUT_REUSE_RECEIPT_SHA256=${REUSE_SHA},TWIRL_FM0_ADMISSION_RECEIPT=${TWIRL_FM0_ADMISSION_RECEIPT},TWIRL_FM0_ADMISSION_RECEIPT_SHA256=${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}"
    ;;
  status)
    require_socket
    remote "set -euo pipefail; echo RUN_ROOT='${RUN_ROOT}'; find '${RUN_ROOT}' -maxdepth 5 -type f \( -name '*.job' -o -name '*.receipt.json' -o -name 'summary.json' \) -print 2>/dev/null | sort || true; squeue -u \$USER -h -o '%i|%P|%T|%C|%m|%b|%j|%N' | grep -E 'twirl-fm0-2|twirl-s[0-9]+' || true"
    ;;
  *) usage >&2; exit 2 ;;
esac
