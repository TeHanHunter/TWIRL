#!/usr/bin/env bash
# Guarded local controller for the Teacher v4-SSL full-pool ORCD workflow.
set -euo pipefail

DRY_RUN=1
ACTION=""

usage() {
  cat <<'EOF'
Usage: run_teacher_ssl_fullpool_orcd.sh [--run] ACTION

Actions:
  probe              Verify the existing user-opened ORCD control socket.
  deploy             Create/reuse a clean detached worktree at the exact Git SHA.
  stage-preregister  Stage and verify the immutable S63 identity reservation.
  submit-preprocess  Submit freeze -> BLS canary -> BLS -> sector merge -> global merge.
  status             Show this workflow's queue and published summaries.

Default mode is dry-run. No action may initiate password or Duo authentication.
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
TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REMOTE_SOURCE="${TWIRL_ORCD_SOURCE_REPO:-${TWIRL_ROOT}/code/TWIRL}"
EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA:-$(git -C "${LOCAL_REPO}" rev-parse HEAD)}"
REMOTE_REPO="${TWIRL_ORCD_REPO:-${TWIRL_ROOT}/code/TWIRL_teacher_ssl_fullpool_${EXPECTED_SHA:0:12}}"
GIT_BRANCH="${TWIRL_GIT_BRANCH:-codex/s56-tier1-qa}"
RUN_ROOT="${TWIRL_SSL_FULLPOOL_RUN_ROOT:-${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp}"
LOCAL_PREREG="${LOCAL_REPO}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp/preregistered"
REMOTE_PREREG="${RUN_ROOT}/preregistered"
S63_TICS_SHA256="0560961c44e73d24f4466f6873c6308b89b9c3ac1d32df5cee7ee87141f6c98b"
S63_SUMMARY_SHA256="b94d94556b036b6819a77596eba9c2514de264a545587ae90628fe4da395eb58"

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
SCP=(
  scp
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
This controller will not fall back to a password/Duo-capable SSH attempt.
EOF
    exit 4
  fi
}

verify_local_preregistration() {
  local tics="${LOCAL_PREREG}/s63_reserved_tics.txt"
  local summary="${LOCAL_PREREG}/s63_reserved_tics.summary.json"
  for input in "${tics}" "${summary}"; do
    if [[ ! -s "${input}" ]]; then
      echo "Missing local preregistration artifact: ${input}" >&2
      exit 3
    fi
  done
  if [[ "$(shasum -a 256 "${tics}" | awk '{print $1}')" != "${S63_TICS_SHA256}" ]]; then
    echo "Local S63 TIC reservation hash changed." >&2
    exit 3
  fi
  if [[ "$(shasum -a 256 "${summary}" | awk '{print $1}')" != "${S63_SUMMARY_SHA256}" ]]; then
    echo "Local S63 reservation summary hash changed." >&2
    exit 3
  fi
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
        observed=\$(git -C \"\${target_repo}\" rev-parse HEAD)
        [[ \"\${observed}\" == \"\${expected_sha}\" ]]
      else
        git -C \"\${source_repo}\" worktree add --detach \"\${target_repo}\" \"\${expected_sha}\"
      fi
      [[ -z \$(git -C \"\${target_repo}\" status --porcelain=v1 --untracked-files=all) ]]
      git -C \"\${target_repo}\" rev-parse HEAD
    "
    ;;
  stage-preregister)
    verify_local_preregistration
    require_socket
    if [[ "${DRY_RUN}" == "1" ]]; then
      remote "mkdir -p '${REMOTE_PREREG}'"
      print_cmd "${SCP[@]}" \
        "${LOCAL_PREREG}/s63_reserved_tics.txt" \
        "${LOCAL_PREREG}/s63_reserved_tics.summary.json" \
        "${ORCD_HOST}:${REMOTE_PREREG}/"
    else
      remote "mkdir -p '${REMOTE_PREREG}'"
      "${SCP[@]}" \
        "${LOCAL_PREREG}/s63_reserved_tics.txt" \
        "${LOCAL_PREREG}/s63_reserved_tics.summary.json" \
        "${ORCD_HOST}:${REMOTE_PREREG}/"
      remote "
        set -euo pipefail
        [[ \$(sha256sum '${REMOTE_PREREG}/s63_reserved_tics.txt' | awk '{print \$1}') == '${S63_TICS_SHA256}' ]]
        [[ \$(sha256sum '${REMOTE_PREREG}/s63_reserved_tics.summary.json' | awk '{print \$1}') == '${S63_SUMMARY_SHA256}' ]]
      "
    fi
    ;;
  submit-preprocess)
    require_socket
    if [[ "${DRY_RUN}" == "1" ]]; then
      echo "Would submit the five-job dependency chain from ${REMOTE_REPO} at ${EXPECTED_SHA}."
      exit 0
    fi
    jobs="$("${SSH[@]}" "${ORCD_HOST}" "
      set -euo pipefail
      repo='${REMOTE_REPO}'
      run='${RUN_ROOT}'
      expected_sha='${EXPECTED_SHA}'
      [[ \$(git -C \"\${repo}\" rev-parse HEAD) == \"\${expected_sha}\" ]]
      [[ -z \$(git -C \"\${repo}\" status --porcelain=v1 --untracked-files=all) ]]
      [[ \$(sha256sum '${REMOTE_PREREG}/s63_reserved_tics.txt' | awk '{print \$1}') == '${S63_TICS_SHA256}' ]]
      [[ \$(sha256sum '${REMOTE_PREREG}/s63_reserved_tics.summary.json' | awk '{print \$1}') == '${S63_SUMMARY_SHA256}' ]]
      mkdir -p '${TWIRL_ROOT}/logs'
      export_arg='ALL,TWIRL_ORCD_REPO=${REMOTE_REPO},TWIRL_SSL_FULLPOOL_RUN_ROOT=${RUN_ROOT}'
      freeze=\$(sbatch --parsable --export=\"\${export_arg}\" scripts/orcd/slurm_teacher_ssl_fullpool_freeze_cpu.sbatch)
      canary=\$(sbatch --parsable --dependency=afterok:\${freeze} --export=\"\${export_arg}\" scripts/orcd/slurm_teacher_ssl_fullpool_bls_canary_cpu.sbatch)
      bls=\$(sbatch --parsable --dependency=afterok:\${canary} --export=\"\${export_arg}\" scripts/orcd/slurm_teacher_ssl_fullpool_bls_cpu.sbatch)
      merge=\$(sbatch --parsable --dependency=afterok:\${bls} --export=\"\${export_arg}\" scripts/orcd/slurm_teacher_ssl_fullpool_bls_merge_cpu.sbatch)
      global=\$(sbatch --parsable --dependency=afterok:\${merge} --export=\"\${export_arg}\" scripts/orcd/slurm_teacher_ssl_fullpool_bls_global_cpu.sbatch)
      printf 'freeze=%s\ncanary=%s\nbls=%s\nmerge=%s\nglobal=%s\n' \"\${freeze}\" \"\${canary}\" \"\${bls}\" \"\${merge}\" \"\${global}\"
    ")"
    printf '%s\n' "${jobs}"
    ;;
  status)
    require_socket
    remote "
      squeue -u tehan -o '%.18i %.32j %.2t %.10M %.4D %R'
      find '${RUN_ROOT}' -maxdepth 5 -type f \\( -name '*summary.json' -o -name '*manifest.json' \\) 2>/dev/null | sort
    "
    ;;
  *)
    usage >&2
    exit 2
    ;;
esac
