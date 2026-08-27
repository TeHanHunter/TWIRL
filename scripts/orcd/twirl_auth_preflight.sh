#!/usr/bin/env bash
# Verify TWIRL control sockets without opening a fresh authenticated session.

set -euo pipefail

readonly TWIRL_CM_DIR="${HOME}/.ssh/cm"
readonly HOSTESS_HOST="hostess3"
readonly BLACKHOLE_HOST="blackhole.mit.edu"
readonly ORCD_HOST="orcd-login.mit.edu"
readonly PDO_HOST="pdogpu1"
readonly ORCD_GLOBUS_REL=".local/twirl-globus-cli/bin/globus"
readonly BLACKHOLE_COLLECTION="d8ea14bc-dca1-4cbc-92b6-b76d7289b6d2"
readonly ORCD_COLLECTION="ec54b570-cac5-47f7-b2a1-100c2078686f"

orcd_only=0
case ${1:-} in
  "") ;;
  --orcd-only) orcd_only=1 ;;
  *) printf 'Usage: %s [--orcd-only]\n' "${0##*/}" >&2; exit 2 ;;
esac

control_path() {
  local host=$1
  ssh -G "${host}" 2>/dev/null | awk \
    -v socket_dir="${TWIRL_CM_DIR}" \
    '$1 == "hostname" {hostname=$2}
     $1 == "user" {user=$2}
     $1 == "port" {port=$2}
     END {print socket_dir "/" user "@" hostname ":" port}'
}

check_master() {
  local host=$1
  local path
  path=$(control_path "${host}")
  if ssh -S "${path}" -O check -o BatchMode=yes "${host}" >/dev/null 2>&1; then
    printf '[ok] %s control master: %s\n' "${host}" "${path}"
    return
  fi
  printf '[missing or stale] %s control master: %s\n' "${host}" "${path}" >&2
  return 1
}

check_orcd_globus_session() {
  local orcd_path
  orcd_path=$(control_path "${ORCD_HOST}")
  if ssh -S "${orcd_path}" -o BatchMode=yes "${ORCD_HOST}" \
    "\$HOME/${ORCD_GLOBUS_REL} session show" >/dev/null 2>&1 \
    && ssh -S "${orcd_path}" -o BatchMode=yes "${ORCD_HOST}" \
      "\$HOME/${ORCD_GLOBUS_REL} ls ${BLACKHOLE_COLLECTION}:/twirl/a2v1 --filter '=s0069' >/dev/null 2>&1" \
    && ssh -S "${orcd_path}" -o BatchMode=yes "${ORCD_HOST}" \
      "\$HOME/${ORCD_GLOBUS_REL} ls ${ORCD_COLLECTION}:/orcd/data/mki_aryeh/001/twirl --filter '=stage1' >/dev/null 2>&1"; then
    printf '[ok] ORCD Globus CLI and collection data access\n'
    return
  fi
  printf '[missing] ORCD Globus CLI session or collection consent\n' >&2
  return 1
}

failed=0
if [[ ${orcd_only} -eq 0 ]]; then
  check_master "${HOSTESS_HOST}" || failed=1
  check_master "${BLACKHOLE_HOST}" || failed=1
fi
check_master "${ORCD_HOST}" || failed=1
if [[ ${orcd_only} -eq 0 ]]; then
  check_master "${PDO_HOST}" || failed=1
fi
if [[ ${failed} -eq 0 ]]; then
  check_orcd_globus_session || failed=1
fi

if [[ ${failed} -ne 0 ]]; then
  printf 'TWIRL authentication preflight failed; rerun open_twirl_control_masters.sh interactively.\n' >&2
  exit 2
fi
printf 'TWIRL authentication preflight passed; automation may reuse BatchMode-only sockets.\n'
