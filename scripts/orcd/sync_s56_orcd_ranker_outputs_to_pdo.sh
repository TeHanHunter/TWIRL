#!/usr/bin/env bash
# Sync small ORCD ranker-selected real-candidate outputs back to PDO for LEO.
#
# Run this from the local TWIRL checkout after the ORCD apply job has produced
# selected_ephemerides.csv. It reuses the user-opened ORCD control socket and
# copies only small report artifacts, not light curves or injection HDF5 files.
set -euo pipefail

ORCD_HOST="${ORCD_HOST:-tehan@orcd-login.mit.edu}"
ORCD_REPO="${ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
ORCD_CONTROL_PATH="${ORCD_CONTROL_PATH:-$HOME/.ssh/cm/%r@%h:%p}"
PDO_HOST="${PDO_HOST:-pdogpu1}"
PDO_REPO="${PDO_REPO:-/pdo/users/tehan/TWIRL}"
SOURCE_DIR="${TWIRL_ORCD_SELECTED_DIR:-reports/stage5_validation/s56_ranker_selected_real_candidates_orcd}"
DEST_DIR="${TWIRL_PDO_SELECTED_DIR:-reports/stage5_validation/s56_ranker_selected_real_candidates_orcd}"
DRY_RUN=1

usage() {
  cat <<'EOF'
Usage: sync_s56_orcd_ranker_outputs_to_pdo.sh [--run]

Copies ORCD ranker-selected real-candidate outputs to the PDO checkout so
scripts/stage5_validation/run_s56_ranker_selected_real_leo_pdo.sh can render
WD-tuned LEO sheets from:

  TWIRL_SELECTED_EPHEMERIDES=reports/stage5_validation/s56_ranker_selected_real_candidates_orcd/selected_ephemerides.csv

Default mode is dry-run. Use --run only after the ORCD apply job has completed.
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
      echo "[sync-orcd-ranker] unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

ORCD_SSH=(
  ssh
  -o BatchMode=yes
  -o PasswordAuthentication=no
  -o KbdInteractiveAuthentication=no
  -o NumberOfPasswordPrompts=0
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

run_or_echo() {
  if [[ "${DRY_RUN}" == "1" ]]; then
    print_cmd "$@"
  else
    "$@"
  fi
}

echo "[sync-orcd-ranker] orcd=${ORCD_HOST}:${ORCD_REPO}/${SOURCE_DIR}"
echo "[sync-orcd-ranker] pdo=${PDO_HOST}:${PDO_REPO}/${DEST_DIR}"
if [[ "${DRY_RUN}" == "1" ]]; then
  echo "[sync-orcd-ranker] dry run. Re-run with --run after ORCD apply completes."
fi

required=(
  selected_ephemerides.csv
  selected_ephemerides_verification.json
  summary.json
  real_peak_table_verification.json
  ranker_selection_summary/summary.json
  ranker_selection_summary/summary.md
  ranker_selection_summary/by_aperture.csv
  ranker_selection_summary/by_period_bin.csv
  ranker_selection_summary/by_tmag_bin.csv
  ranker_selection_summary/by_duration_bin.csv
  ranker_selection_summary/by_selection_rank.csv
  ranker_selection_summary/top_ranker_rows.csv
)

remote_check=$(
  printf '%s\n' "${required[@]}" |
    "${ORCD_SSH[@]}" "${ORCD_HOST}" "if [[ ! -d '${ORCD_REPO}/${SOURCE_DIR}' ]]; then echo 'MISSING_DIR ${SOURCE_DIR}'; exit 0; fi; cd '${ORCD_REPO}/${SOURCE_DIR}' && while IFS= read -r p; do [[ -f \"\$p\" ]] && printf '%s\\n' \"\$p\" || printf 'MISSING %s\\n' \"\$p\"; done" \
    2>/dev/null || true
)
if [[ -z "${remote_check}" ]]; then
  echo "[sync-orcd-ranker] could not inspect ORCD source directory. Is the control socket open?" >&2
  exit 4
fi
printf '%s\n' "${remote_check}" | sed 's/^/[sync-orcd-ranker] /'
if printf '%s\n' "${remote_check}" | grep -q '^MISSING_DIR '; then
  echo "[sync-orcd-ranker] ORCD apply output directory is not present yet; not syncing." >&2
  exit 5
fi
if printf '%s\n' "${remote_check}" | grep -q '^MISSING '; then
  echo "[sync-orcd-ranker] required ORCD outputs are not complete yet; not syncing." >&2
  exit 5
fi

if [[ "${DRY_RUN}" == "1" ]]; then
  print_cmd "${ORCD_SSH[@]}" "${ORCD_HOST}" "cd '${ORCD_REPO}' && tar -cf - '${SOURCE_DIR}'"
  print_cmd "${PDO_SSH[@]}" "${PDO_HOST}" "mkdir -p '${PDO_REPO}' && tar -C '${PDO_REPO}' -xf -"
else
  "${ORCD_SSH[@]}" "${ORCD_HOST}" "cd '${ORCD_REPO}' && tar -cf - '${SOURCE_DIR}'" |
    "${PDO_SSH[@]}" "${PDO_HOST}" "mkdir -p '${PDO_REPO}' && tar -C '${PDO_REPO}' -xf -"
fi

echo "[sync-orcd-ranker] done"
echo "[sync-orcd-ranker] next PDO LEO command:"
cat <<EOF
  cd ${PDO_REPO}
  TWIRL_SELECTED_EPHEMERIDES=${DEST_DIR}/selected_ephemerides.csv \\
    bash scripts/stage5_validation/run_s56_ranker_selected_real_leo_pdo.sh
EOF
