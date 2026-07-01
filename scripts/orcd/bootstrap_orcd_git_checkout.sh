#!/usr/bin/env bash
# Make the ORCD TWIRL directory a real Git checkout that tracks origin/main.
#
# The ORCD tree also holds large local data/report artifacts under gitignored
# paths. This helper initializes Git around the existing directory and fetches
# the committed source tree without deleting those local products.
set -euo pipefail

DRY_RUN=1

usage() {
  cat <<'EOF'
Usage: bootstrap_orcd_git_checkout.sh [--run]

Initializes or refreshes /orcd/data/mki_aryeh/001/twirl/code/TWIRL as a Git
checkout tracking the TWIRL GitHub repository. Default mode is dry-run.

This helper requires the user-opened ORCD SSH control socket. It never prompts
for password or Duo.
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
      echo "[orcd-git] unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

ORCD_HOST="${ORCD_HOST:-tehan@orcd-login.mit.edu}"
ORCD_REPO="${ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
ORCD_CONTROL_PATH="${ORCD_CONTROL_PATH:-$HOME/.ssh/cm/%r@%h:%p}"
TWIRL_GIT_REMOTE="${TWIRL_GIT_REMOTE:-https://github.com/TeHanHunter/TWIRL.git}"
TWIRL_GIT_BRANCH="${TWIRL_GIT_BRANCH:-main}"
TWIRL_ORCD_SPARSE="${TWIRL_ORCD_SPARSE:-1}"

ORCD_SSH=(
  ssh
  -o BatchMode=yes
  -o PasswordAuthentication=no
  -o KbdInteractiveAuthentication=no
  -o NumberOfPasswordPrompts=0
  -o ControlMaster=auto
  -o ControlPath="${ORCD_CONTROL_PATH}"
)

print_cmd() {
  printf '+'
  printf ' %q' "$@"
  printf '\n'
}

remote_script=$(cat <<'EOS'
set -euo pipefail

echo "[orcd-git] repo=${ORCD_REPO}"
echo "[orcd-git] remote=${TWIRL_GIT_REMOTE}"
echo "[orcd-git] branch=${TWIRL_GIT_BRANCH}"
echo "[orcd-git] sparse=${TWIRL_ORCD_SPARSE}"

export GIT_SSH_COMMAND="${GIT_SSH_COMMAND:-ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new}"
export GIT_TERMINAL_PROMPT=0

mkdir -p "$(dirname "${ORCD_REPO}")"
mkdir -p "${ORCD_REPO}"
cd "${ORCD_REPO}"

if [[ ! -d .git ]]; then
  git init -b "${TWIRL_GIT_BRANCH}"
fi
if git remote get-url origin >/dev/null 2>&1; then
  git remote set-url origin "${TWIRL_GIT_REMOTE}"
else
  git remote add origin "${TWIRL_GIT_REMOTE}"
fi
git fetch origin "${TWIRL_GIT_BRANCH}"
git symbolic-ref HEAD "refs/heads/${TWIRL_GIT_BRANCH}"
git branch --set-upstream-to="origin/${TWIRL_GIT_BRANCH}" "${TWIRL_GIT_BRANCH}" >/dev/null 2>&1 || true

if [[ "${TWIRL_ORCD_SPARSE}" == "1" ]]; then
  git sparse-checkout init --cone
  git sparse-checkout set \
    catalogs \
    configs \
    doc \
    scripts \
    src \
    tests
fi

git reset --hard "origin/${TWIRL_GIT_BRANCH}"

git rev-parse --show-toplevel
git remote -v
git status --short --branch
EOS
)

remote_command="$(
  printf 'ORCD_REPO=%q TWIRL_GIT_REMOTE=%q TWIRL_GIT_BRANCH=%q TWIRL_ORCD_SPARSE=%q bash -lc %q' \
    "${ORCD_REPO}" \
    "${TWIRL_GIT_REMOTE}" \
    "${TWIRL_GIT_BRANCH}" \
    "${TWIRL_ORCD_SPARSE}" \
    "${remote_script}"
)"

echo "[orcd-git] orcd=${ORCD_HOST}:${ORCD_REPO}"
if [[ "${DRY_RUN}" == "1" ]]; then
  echo "[orcd-git] dry run. Re-run with --run after the ORCD control socket is open."
  print_cmd "${ORCD_SSH[@]}" "${ORCD_HOST}" "${remote_command}"
else
  "${ORCD_SSH[@]}" "${ORCD_HOST}" "${remote_command}"
fi
