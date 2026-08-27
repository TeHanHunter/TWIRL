#!/usr/bin/env bash
# Interactively open the persistent SSH and Globus sessions used by TWIRL.
#
# This helper is intentionally user-run from a real terminal. Password, Duo,
# and browser-consent prompts may appear here; TWIRL automation only reuses
# the resulting control sockets in BatchMode.

set -euo pipefail

readonly TWIRL_SOCKET_DIR="${HOME}/.ssh/cm"
readonly TWIRL_PERSIST="${TWIRL_CONTROL_PERSIST:-604800}"
readonly TWIRL_IDENTITY="${TWIRL_SSH_IDENTITY:-${HOME}/.ssh/id_tess}"
readonly TWIRL_HOSTESS_HOST="hostess3"
readonly TWIRL_BLACKHOLE_HOST="blackhole.mit.edu"
readonly TWIRL_ORCD_HOST="orcd-login.mit.edu"
readonly TWIRL_PDO_HOST="pdogpu1"
readonly TWIRL_ORCD_GLOBUS_REL=".local/twirl-globus-cli/bin/globus"
readonly TWIRL_BLACKHOLE_COLLECTION="d8ea14bc-dca1-4cbc-92b6-b76d7289b6d2"
readonly TWIRL_ORCD_COLLECTION="ec54b570-cac5-47f7-b2a1-100c2078686f"
readonly TWIRL_ORCD_DATA_ACCESS_SCOPE="urn:globus:auth:scope:transfer.api.globus.org:all[*https://auth.globus.org/scopes/${TWIRL_BLACKHOLE_COLLECTION}/data_access *https://auth.globus.org/scopes/${TWIRL_ORCD_COLLECTION}/data_access]"
readonly TWIRL_SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

TWIRL_ORCD_ONLY=0

usage() {
  cat <<'EOF'
Usage: open_twirl_control_masters.sh [--orcd-only]

Run interactively with no arguments to authenticate every current TWIRL route:
Hostess3, Blackhole, ORCD, PDO, and the ORCD Globus CLI. This command only
authenticates; it never starts or resumes production.

  --orcd-only  Open ORCD and its Globus session only.

Environment overrides:
  TWIRL_SSH_IDENTITY      SSH private key (default: ~/.ssh/id_tess)
  TWIRL_CONTROL_PERSIST   Socket lifetime in seconds (default: 604800)
EOF
}

case ${1:-} in
  "") ;;
  --orcd-only) TWIRL_ORCD_ONLY=1 ;;
  --help|-h) usage; exit 0 ;;
  *) usage >&2; exit 2 ;;
esac
if [[ $# -gt 1 ]]; then
  usage >&2
  exit 2
fi
if [[ ! -t 0 || ! -t 1 ]]; then
  printf 'Run this helper directly from an interactive terminal so authentication can prompt.\n' >&2
  exit 2
fi
if [[ ! -r ${TWIRL_IDENTITY} ]]; then
  printf 'TWIRL SSH identity is unavailable: %s\n' "${TWIRL_IDENTITY}" >&2
  exit 2
fi

mkdir -p "${TWIRL_SOCKET_DIR}"
chmod 700 "${TWIRL_SOCKET_DIR}"

control_path() {
  local host=$1
  ssh -G "${host}" 2>/dev/null | awk \
    -v socket_dir="${TWIRL_SOCKET_DIR}" \
    '$1 == "hostname" {hostname=$2}
     $1 == "user" {user=$2}
     $1 == "port" {port=$2}
     END {print socket_dir "/" user "@" hostname ":" port}'
}

master_is_live() {
  local host=$1
  local path=$2
  ssh -S "${path}" -O check -o BatchMode=yes "${host}" >/dev/null 2>&1
}

open_local_master() {
  local label=$1
  local host=$2
  local proxy_command=${3:-}
  local path
  path=$(control_path "${host}")
  if master_is_live "${host}" "${path}"; then
    printf '[already open] %s: %s\n' "${label}" "${path}"
    return
  fi

  printf '\nOpening %s (%s); answer password/Duo if requested.\n' "${label}" "${host}"
  local -a ssh_command=(
    ssh -tt -A
    -i "${TWIRL_IDENTITY}"
    -o IdentitiesOnly=yes
    -o ControlMaster=yes
    -o "ControlPath=${path}"
    -o "ControlPersist=${TWIRL_PERSIST}"
    -o ConnectTimeout=30
    -o ConnectionAttempts=1
    -o KbdInteractiveAuthentication=yes
    -o PreferredAuthentications=publickey,keyboard-interactive,password
    -o ServerAliveInterval=30
    -o ServerAliveCountMax=6
    -o TCPKeepAlive=yes
    -o ForwardAgent=yes
  )
  if [[ -n ${proxy_command} ]]; then
    ssh_command+=(-o "ProxyCommand=${proxy_command}")
  fi
  ssh_command+=("${host}" true)
  "${ssh_command[@]}"
  master_is_live "${host}" "${path}" || {
    printf 'Failed to establish %s control master.\n' "${label}" >&2
    exit 2
  }
  printf '[opened] %s: %s\n' "${label}" "${path}"
}

orcd_globus_data_access_ready() {
  local orcd_path=$1
  ssh -S "${orcd_path}" -o BatchMode=yes "${TWIRL_ORCD_HOST}" \
    "\$HOME/${TWIRL_ORCD_GLOBUS_REL} ls ${TWIRL_BLACKHOLE_COLLECTION}:/twirl/a2v1 --filter '=s0069' >/dev/null 2>&1" \
    && ssh -S "${orcd_path}" -o BatchMode=yes "${TWIRL_ORCD_HOST}" \
      "\$HOME/${TWIRL_ORCD_GLOBUS_REL} ls ${TWIRL_ORCD_COLLECTION}:/orcd/data/mki_aryeh/001/twirl --filter '=stage1' >/dev/null 2>&1"
}

open_orcd_globus_session() {
  local orcd_path
  orcd_path=$(control_path "${TWIRL_ORCD_HOST}")
  if ssh -S "${orcd_path}" -o BatchMode=yes "${TWIRL_ORCD_HOST}" \
    "\$HOME/${TWIRL_ORCD_GLOBUS_REL} session show" >/dev/null 2>&1; then
    printf '[already authenticated] ORCD Globus CLI; checking collection consent.\n'
  else
    printf '\nOpening ORCD Globus browser login. Complete it once in your browser.\n'
    ssh -tt -S "${orcd_path}" "${TWIRL_ORCD_HOST}" \
      "\$HOME/${TWIRL_ORCD_GLOBUS_REL} login"
  fi
  if ! ssh -S "${orcd_path}" -o BatchMode=yes "${TWIRL_ORCD_HOST}" \
    "\$HOME/${TWIRL_ORCD_GLOBUS_REL} session show" >/dev/null 2>&1; then
    printf 'ORCD Globus login did not complete.\n' >&2
    exit 2
  fi
  if orcd_globus_data_access_ready "${orcd_path}"; then
    printf '[authenticated] ORCD Globus CLI and collection data access\n'
    return
  fi
  printf '\nGranting the required TESS TSO and ORCD data-access scopes.\n'
  ssh -tt -S "${orcd_path}" "${TWIRL_ORCD_HOST}" \
    "\$HOME/${TWIRL_ORCD_GLOBUS_REL} session consent '${TWIRL_ORCD_DATA_ACCESS_SCOPE}'"
  orcd_globus_data_access_ready "${orcd_path}" || {
    printf 'ORCD Globus collection data access did not complete.\n' >&2
    exit 2
  }
  printf '[authenticated] ORCD Globus CLI and collection data access\n'
}

if [[ ${TWIRL_ORCD_ONLY} -eq 0 ]]; then
  open_local_master "Hostess3 jump host" "${TWIRL_HOSTESS_HOST}"
  TWIRL_HOSTESS_PATH=$(control_path "${TWIRL_HOSTESS_HOST}")
  readonly TWIRL_HOSTESS_PATH
  readonly TWIRL_HOSTESS_PROXY="ssh -i ${TWIRL_IDENTITY} -o IdentitiesOnly=yes -S ${TWIRL_HOSTESS_PATH} -o BatchMode=yes -W %h:%p ${TWIRL_HOSTESS_HOST}"
  open_local_master "Blackhole" "${TWIRL_BLACKHOLE_HOST}" "${TWIRL_HOSTESS_PROXY}"
fi

open_local_master "ORCD" "${TWIRL_ORCD_HOST}"

if [[ ${TWIRL_ORCD_ONLY} -eq 0 ]]; then
  open_local_master "PDO" "${TWIRL_PDO_HOST}" "${TWIRL_HOSTESS_PROXY}"
fi

open_orcd_globus_session

printf '\nVerifying noninteractive routes.\n'
if [[ ${TWIRL_ORCD_ONLY} -eq 1 ]]; then
  "${TWIRL_SCRIPT_DIR}/twirl_auth_preflight.sh" --orcd-only
else
  "${TWIRL_SCRIPT_DIR}/twirl_auth_preflight.sh"
fi

printf '[ok] Authentication complete; no production process was started.\n'
