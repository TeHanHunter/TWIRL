#!/usr/bin/env bash
# Sync ORCD-rendered ADP015 two-aperture vet sheets back to local and PDO.
#
# This copies only review-sheet images/PDFs and metrics/summary CSV/JSON files.
# It does not transfer light curves, injections, or FITS trees.
set -euo pipefail

DRY_RUN=1

usage() {
  cat <<'EOF'
Usage: sync_s56_orcd_twoap_outputs.sh [--run]

Copies ORCD-rendered ADP015 two-aperture vet sheets to both the local checkout
and the PDO checkout. Default mode is dry-run. Use --run only after the ORCD
twoap-smoke or twoap-full job has completed and the ORCD control socket is open.

Useful overrides:
  TWIRL_ORCD_TWOAP_OUT_DIR
  TWIRL_ORCD_TWOAP_METRICS_CSV
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
      echo "[sync-orcd-twoap] unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

LOCAL_REPO="${LOCAL_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
ORCD_HOST="${ORCD_HOST:-tehan@orcd-login.mit.edu}"
ORCD_REPO="${ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
ORCD_CONTROL_PATH="${ORCD_CONTROL_PATH:-$HOME/.ssh/cm/%r@%h:%p}"
ORCD_CONNECT_TIMEOUT="${ORCD_CONNECT_TIMEOUT:-15}"
PDO_HOST="${PDO_HOST:-pdogpu1}"
PDO_REPO="${PDO_REPO:-/pdo/users/tehan/TWIRL}"
SOURCE_DIR="${TWIRL_ORCD_TWOAP_OUT_DIR:-reports/stage5_validation/s56_mixed_teacher_queue_pdo/twirl_vet_sheets_adp015q_orcd}"
METRICS_CSV="${TWIRL_ORCD_TWOAP_METRICS_CSV:-reports/stage5_validation/s56_mixed_teacher_queue_pdo/twirl_vet_metrics_adp015q_orcd.csv}"
SUMMARY_JSON="${TWIRL_ORCD_TWOAP_SUMMARY_JSON:-reports/stage5_validation/s56_mixed_teacher_queue_pdo/twirl_two_aperture_vet_summary_adp015q_orcd.json}"

ORCD_SSH=(
  ssh
  -o BatchMode=yes
  -o PasswordAuthentication=no
  -o KbdInteractiveAuthentication=no
  -o NumberOfPasswordPrompts=0
  -o ConnectTimeout="${ORCD_CONNECT_TIMEOUT}"
  -o ConnectionAttempts=1
  -o ServerAliveInterval=15
  -o ServerAliveCountMax=1
  -o ControlMaster=auto
  -o ControlPath="${ORCD_CONTROL_PATH}"
)
PDO_SSH=(
  ssh
  -o BatchMode=yes
  -o ConnectTimeout=20
  -o ControlMaster=no
)

print_cmd() {
  printf '+'
  printf ' %q' "$@"
  printf '\n'
}

echo "[sync-orcd-twoap] orcd=${ORCD_HOST}:${ORCD_REPO}/${SOURCE_DIR}"
echo "[sync-orcd-twoap] metrics=${METRICS_CSV}"
echo "[sync-orcd-twoap] summary=${SUMMARY_JSON}"
echo "[sync-orcd-twoap] local=${LOCAL_REPO}/${SOURCE_DIR}"
echo "[sync-orcd-twoap] pdo=${PDO_HOST}:${PDO_REPO}/${SOURCE_DIR}"
if [[ "${DRY_RUN}" == "1" ]]; then
  echo "[sync-orcd-twoap] dry run. Re-run with --run after the ORCD twoap job completes."
fi

required=("${SOURCE_DIR}" "${METRICS_CSV}" "${SUMMARY_JSON}")
remote_check=$(
  printf '%s\n' "${required[@]}" |
    "${ORCD_SSH[@]}" "${ORCD_HOST}" "cd '${ORCD_REPO}' && while IFS= read -r p; do [[ -e \"\$p\" ]] && printf 'OK %s\\n' \"\$p\" || printf 'MISSING %s\\n' \"\$p\"; done" \
    2>/dev/null || true
)
if [[ -z "${remote_check}" ]]; then
  echo "[sync-orcd-twoap] could not inspect ORCD outputs. Is the control socket open?" >&2
  exit 4
fi
printf '%s\n' "${remote_check}" | sed 's/^/[sync-orcd-twoap] /'
if printf '%s\n' "${remote_check}" | grep -q '^MISSING '; then
  echo "[sync-orcd-twoap] required ORCD two-aperture outputs are incomplete; not syncing." >&2
  exit 5
fi

if [[ "${DRY_RUN}" == "1" ]]; then
  print_cmd "${ORCD_SSH[@]}" "${ORCD_HOST}" "cd '${ORCD_REPO}' && tar -cf - '${SOURCE_DIR}' '${METRICS_CSV}' '${SUMMARY_JSON}'"
  print_cmd tar -C "${LOCAL_REPO}" -xf -
  print_cmd "${PDO_SSH[@]}" "${PDO_HOST}" "mkdir -p '${PDO_REPO}' && tar -C '${PDO_REPO}' -xf -"
else
  "${ORCD_SSH[@]}" "${ORCD_HOST}" "cd '${ORCD_REPO}' && tar -cf - '${SOURCE_DIR}' '${METRICS_CSV}' '${SUMMARY_JSON}'" |
    tar -C "${LOCAL_REPO}" -xf -
  "${ORCD_SSH[@]}" "${ORCD_HOST}" "cd '${ORCD_REPO}' && tar -cf - '${SOURCE_DIR}' '${METRICS_CSV}' '${SUMMARY_JSON}'" |
    "${PDO_SSH[@]}" "${PDO_HOST}" "mkdir -p '${PDO_REPO}' && tar -C '${PDO_REPO}' -xf -"
fi

echo "[sync-orcd-twoap] done"
