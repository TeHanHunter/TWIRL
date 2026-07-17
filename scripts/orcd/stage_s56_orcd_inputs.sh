#!/usr/bin/env bash
# Stage compact S56 TWIRL products and the current checkout to ORCD.
#
# This is meant to be run from the local TWIRL checkout after opening the ORCD
# SSH control socket described in doc/orcd_h200_usage.md. It streams large PDO
# artifacts through the local machine without staging raw TGLC/TICA trees.
set -euo pipefail

DRY_RUN=1
STAGE_SUBSET="${TWIRL_STAGE_SUBSET:-all}"

usage() {
  cat <<'EOF'
Usage: stage_s56_orcd_inputs.sh [--run] [--subset all|balanced|allhost]

Bootstraps/refetches the ORCD Git checkout, then stages compact PDO artifacts
to ORCD. Default mode is dry-run. Use --run only after the ORCD control socket
is open.

By default code is updated through Git, not rsync. Set
TWIRL_STAGE_RSYNC_CODE=1 only for an explicit dirty-checkout deployment smoke.

Subsets:
  all       Stage the full compact S56 downstream bundle.
  balanced  Stage the compact LC export, balanced 20k injection product, real
            S56 BLS table, and recovery summaries.
  allhost   Stage the all-host sharded injection product, all-host validation
            reports/peak chunks, and real S56 BLS table.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run)
      DRY_RUN=0
      shift
      ;;
    --subset)
      STAGE_SUBSET="${2:?missing value for --subset}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "[stage-orcd] unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

LOCAL_REPO="${LOCAL_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
PDO_HOST="${PDO_HOST:-pdogpu1}"
PDO_REPO="${PDO_REPO:-/pdo/users/tehan/TWIRL}"
ORCD_HOST="${ORCD_HOST:-tehan@orcd-login.mit.edu}"
ORCD_REPO="${ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
ORCD_CONTROL_PATH="${ORCD_CONTROL_PATH:-$HOME/.ssh/cm/%r@%h:%p}"
ORCD_CONNECT_TIMEOUT="${ORCD_CONNECT_TIMEOUT:-15}"
TWIRL_STAGE_RSYNC_CODE="${TWIRL_STAGE_RSYNC_CODE:-0}"

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
ORCD_SSH_CMD="ssh -o BatchMode=yes -o PasswordAuthentication=no -o KbdInteractiveAuthentication=no -o NumberOfPasswordPrompts=0 -o ConnectTimeout=${ORCD_CONNECT_TIMEOUT} -o ConnectionAttempts=1 -o ServerAliveInterval=15 -o ServerAliveCountMax=1 -o ControlMaster=auto -o ControlPath=${ORCD_CONTROL_PATH}"
PDO_SSH=(
  ssh
  -o BatchMode=yes
  -o ConnectTimeout=20
  -o ControlMaster=no
)

BALANCED_PDO_PATHS=(
  data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_lc_export_pdo.h5
  data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_lc_export_pdo.manifest.json
  data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp_lc_export_pdo.h5
  data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp_lc_export_pdo.manifest.json
  data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp015_lc_export_pdo.h5
  data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp015_lc_export_pdo.manifest.json
  data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_20k_predetrend_batman_periodradius_bright_balanced/injected_lightcurves.h5
  data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_20k_predetrend_batman_periodradius_bright_balanced/injection_manifest.csv
  data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_20k_predetrend_batman_periodradius_bright_balanced/injection_labels.csv
  data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_20k_predetrend_batman_periodradius_bright_balanced/summary.json
  data_local/stage2/bls_first_pass_v2/sector_0056/candidates.parquet
  data_local/stage2/bls_first_pass_v2/sector_0056/candidates.csv
  reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/small_pair_200k/injection_bls_recoveries.csv
  reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/small_pair_200k/recovery_mode_summary/summary.json
  reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/duration_aware/summary.json
  reports/stage5_validation/s56_mixed_teacher_queue_pdo/mixed_teacher_pool.csv
  reports/stage5_validation/s56_mixed_teacher_queue_pdo/review_queue_1k.csv
  reports/stage5_validation/s56_mixed_teacher_queue_pdo/review_queue_1k_pre_leo.csv
  reports/stage5_validation/s56_mixed_teacher_queue_pdo/summary.json
  reports/stage5_validation/s56_mixed_teacher_queue_pdo/verification.json
)

ALLHOST_PDO_PATHS=(
  data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_allhost_predetrend_batman_periodradius_grid/injected_lightcurves.h5
  data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_allhost_predetrend_batman_periodradius_grid/injection_manifest.csv
  data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_allhost_predetrend_batman_periodradius_grid/injection_labels.csv
  data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_allhost_predetrend_batman_periodradius_grid/summary.json
  data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_allhost_predetrend_batman_periodradius_grid_sharded
  reports/stage5_validation/s56_allhost_predetrend_batman_periodradius_grid_sharded_pdo
  data_local/stage2/bls_first_pass_v2/sector_0056/candidates.parquet
  data_local/stage2/bls_first_pass_v2/sector_0056/candidates.csv
)

case "${STAGE_SUBSET}" in
  all)
    PDO_PATHS=("${BALANCED_PDO_PATHS[@]}" "${ALLHOST_PDO_PATHS[@]}")
    ;;
  balanced)
    PDO_PATHS=("${BALANCED_PDO_PATHS[@]}")
    ;;
  allhost)
    PDO_PATHS=("${ALLHOST_PDO_PATHS[@]}")
    ;;
  *)
    echo "[stage-orcd] unknown subset: ${STAGE_SUBSET}" >&2
    usage >&2
    exit 2
    ;;
esac

run_or_echo() {
  if [[ "${DRY_RUN}" == "1" ]]; then
    printf '+'
    printf ' %q' "$@"
    printf '\n'
  else
    "$@"
  fi
}

echo "[stage-orcd] local_repo=${LOCAL_REPO}"
echo "[stage-orcd] pdo=${PDO_HOST}:${PDO_REPO}"
echo "[stage-orcd] orcd=${ORCD_HOST}:${ORCD_REPO}"
echo "[stage-orcd] subset=${STAGE_SUBSET}"

if [[ "${DRY_RUN}" == "1" ]]; then
  echo "[stage-orcd] dry run. Re-run with --run after the ORCD control socket is open."
fi

echo "[stage-orcd] ensuring ORCD code checkout is Git-tracked"
if [[ "${DRY_RUN}" == "1" ]]; then
  run_or_echo "${LOCAL_REPO}/scripts/orcd/bootstrap_orcd_git_checkout.sh" --run
else
  "${LOCAL_REPO}/scripts/orcd/bootstrap_orcd_git_checkout.sh" --run
fi

RSYNC_ARGS=(
  -az
  --delete
  --exclude .git/
  --exclude .venv/
  --exclude __pycache__/
  --exclude .pytest_cache/
  --exclude .matplotlib_cache/
  --exclude data_local/
  --exclude reports/
  --exclude logs/
  --exclude outputs/current_keynote_edit/
)
if [[ "${DRY_RUN}" == "1" ]]; then
  RSYNC_ARGS=(-n "${RSYNC_ARGS[@]}")
fi

if [[ "${TWIRL_STAGE_RSYNC_CODE}" == "1" ]]; then
  echo "[stage-orcd] syncing current local code into the Git checkout"
  run_or_echo rsync "${RSYNC_ARGS[@]}" -e "${ORCD_SSH_CMD}" "${LOCAL_REPO}/" "${ORCD_HOST}:${ORCD_REPO}/"
else
  echo "[stage-orcd] code sync skipped; ORCD code is managed by Git"
fi

echo "[stage-orcd] checking compact PDO artifact list"
if [[ "${DRY_RUN}" == "1" ]]; then
  printf '%s\n' "${PDO_PATHS[@]}" | sed 's/^/  /'
else
  printf '%s\n' "${PDO_PATHS[@]}" |
    "${PDO_SSH[@]}" "${PDO_HOST}" "cd '${PDO_REPO}' && while IFS= read -r p; do [[ -e \"\$p\" ]] && printf '%s\\0' \"\$p\" || echo \"[stage-orcd] missing optional PDO artifact: \$p\" >&2; done" |
    "${PDO_SSH[@]}" "${PDO_HOST}" "cd '${PDO_REPO}' && tar --null -T - -cf -" |
    "${ORCD_SSH[@]}" "${ORCD_HOST}" "mkdir -p '${ORCD_REPO}' && tar -C '${ORCD_REPO}' -xf -"
fi

echo "[stage-orcd] done"
echo "[stage-orcd] next ORCD check:"
echo "  ssh -o BatchMode=yes -o PasswordAuthentication=no -o KbdInteractiveAuthentication=no -o NumberOfPasswordPrompts=0 -o ControlMaster=auto -o ControlPath=${ORCD_CONTROL_PATH} ${ORCD_HOST} 'cd ${ORCD_REPO} && ls -lh data_local/stage3_injections/s56_twirlfs_v2_lc_export || true'"
