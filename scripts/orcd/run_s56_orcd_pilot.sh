#!/usr/bin/env bash
# Local helper for staging/submitting the S56 downstream ORCD pilot.
#
# Run this after opening the ORCD control socket described in
# doc/orcd_h200_usage.md. By default it only prints the commands it would run;
# pass --run to execute them.
set -euo pipefail

LOCAL_REPO="${LOCAL_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
ORCD_HOST="${ORCD_HOST:-tehan@orcd-login.mit.edu}"
ORCD_REPO="${ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
ORCD_CONTROL_PATH="${ORCD_CONTROL_PATH:-$HOME/.ssh/cm/%r@%h:%p}"
ORCD_CONNECT_TIMEOUT="${ORCD_CONNECT_TIMEOUT:-15}"
PARTITION="${TWIRL_ORCD_PARTITION:-pg_mki_aryeh}"
DRY_RUN=1
COMMANDS=()

usage() {
  cat <<'EOF'
Usage: run_s56_orcd_pilot.sh [--run] COMMAND [COMMAND ...]

Commands:
  probe          Verify the ORCD control socket and partition are reachable.
  sync-code      Rsync the current local code checkout to ORCD, preserving data.
  stage          Stage the local checkout and compact S56 PDO artifacts.
  stage-balanced Stage only the compact balanced-grid S56 artifacts.
  stage-allhost  Stage only the all-host S56 artifacts and reports.
  smoke          Submit a 100-injection CPU peak-table smoke job.
  peak-full      Submit the full restartable CPU peak-table job.
  ranker         Submit peak-ranker training plus the post-BLS gate.
  short-branch   Submit short-period branch comparison after the standard table exists.
  apply          Submit ranker application to real S56 BLS peaks.
  allhost-peak   Submit all-host sharded injected peak-table build.
  allhost-ranker Submit all-host peak-ranker training plus post-BLS gate.
  allhost-apply  Submit real-candidate ranker application using all-host model.
  predetrend-bls-smoke
                 Submit a CPU raw-flux detrending-strength BLS audit smoke.
  predetrend-bls-full
                 Submit the full CPU raw-flux detrending-strength BLS audit.
  adpplus-smoke  Submit a CPU ADP+ two-aperture audit smoke.
  adpplus-full   Submit the full CPU ADP+ two-aperture audit.
  twoap-smoke    Submit a CPU ADP015 two-aperture vet-sheet smoke.
  twoap-full     Submit the full CPU ADP015 two-aperture vet-sheet render.
  sync-twoap-smoke
                 Sync ORCD twoap-smoke outputs back to local/PDO.
  sync-twoap     Sync full ORCD two-aperture outputs back to local/PDO.
  allhost-sync-apply
                 Sync verified all-host ORCD selected outputs back to PDO.
  allhost-monitor-apply
                 Wait for all-host ORCD apply outputs, then sync to PDO.
  sync-apply     Sync verified ORCD ranker-selected outputs back to PDO.
  monitor-apply  Wait for ORCD apply outputs, then sync them back to PDO.
  monitor-apply-leo-smoke
                 Wait for ORCD apply outputs, sync to PDO, then render a
                 bounded 50-row LEO queue on PDO.
  h200-smoke     Submit a one-H200 smoke job.
  h200-tensor-smoke
                 Submit a one-H200 candidate-centered tensor data smoke.
  h200-torch-tensor-smoke
                 Submit tensor smoke with the torch env and require torch/CUDA.
  stage-labels   Stage human-label teacher/readiness products to ORCD.
  h200-tensor-train-smoke
                 Submit one-H200 tensor training smoke using real teacher labels.
  h200-tensor-train-synthetic-smoke
                 Submit one-H200 tensor training smoke with explicit synthetic labels.
  status         Show active TWIRL jobs and recent ORCD log files.
  smoke-pilot    Run probe, stage, smoke, and h200-smoke.

Default mode is dry-run. Use --run after the ORCD control socket is open.

Example:
  scripts/orcd/run_s56_orcd_pilot.sh probe
  scripts/orcd/run_s56_orcd_pilot.sh --run smoke-pilot
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
      COMMANDS+=("$1")
      shift
      ;;
  esac
done

if [[ "${#COMMANDS[@]}" -eq 0 ]]; then
  usage >&2
  exit 2
fi

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

orcd() {
  run_or_echo "${ORCD_SSH[@]}" "${ORCD_HOST}" "$1"
}

require_socket() {
  if [[ "${DRY_RUN}" == "1" ]]; then
    return 0
  fi
  if ! "${ORCD_SSH[@]}" "${ORCD_HOST}" "true" >/dev/null 2>&1; then
    cat >&2 <<EOF
[orcd-pilot] ORCD control socket is not available.

Open it from a user terminal first:

  mkdir -p ~/.ssh/cm
  chmod 700 ~/.ssh ~/.ssh/cm
  ssh -MNf -o ControlMaster=yes -o ControlPath="\$HOME/.ssh/cm/%r@%h:%p" -o ControlPersist=8h -o ServerAliveInterval=60 -o ServerAliveCountMax=3 ${ORCD_HOST}

The agent must not retry in a way that can trigger password or Duo prompts.

EOF
    exit 4
  fi
}

submit() {
  local name="$1"
  local export_arg="$2"
  local script="$3"
  echo "[orcd-pilot] submit ${name}: ${script}"
  enforce_h200_cap "${script}"
  orcd "cd '${ORCD_REPO}' && mkdir -p /orcd/data/mki_aryeh/001/twirl/logs && sbatch ${export_arg} '${script}'"
}

enforce_h200_cap() {
  local script="$1"
  local gres n_h200
  if [[ ! -f "${LOCAL_REPO}/${script}" ]]; then
    return 0
  fi
  gres="$(grep -E '^[[:space:]]*#SBATCH[[:space:]]+--gres=gpu:h200:[0-9]+' "${LOCAL_REPO}/${script}" || true)"
  if [[ -z "${gres}" ]]; then
    return 0
  fi
  while IFS= read -r line; do
    n_h200="${line##*:}"
    if [[ "${n_h200}" -gt 2 && "${TWIRL_ALLOW_MANY_H200:-0}" != "1" ]]; then
      cat >&2 <<EOF
[orcd-pilot] Refusing to submit ${script}: it requests ${n_h200} H200s.

TWIRL's default ORCD policy is CPU-only where possible, 1 H200 for smoke tests,
and at most 2 H200s for routine development. Get explicit user approval before
running larger jobs, then set TWIRL_ALLOW_MANY_H200=1 for that specific submit.
EOF
      exit 5
    fi
  done <<<"${gres}"
}

cmd_probe() {
  echo "[orcd-pilot] probe ${ORCD_HOST}"
  orcd "hostname; whoami; sinfo -h -p '${PARTITION}' -o '%P %D %t %G' | head -10"
}

cmd_sync_code() {
  echo "[orcd-pilot] sync current local code to ORCD"
  local rsync_args=(
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
    rsync_args=(-n "${rsync_args[@]}")
  fi
  run_or_echo rsync "${rsync_args[@]}" -e "${ORCD_SSH_CMD}" "${LOCAL_REPO}/" "${ORCD_HOST}:${ORCD_REPO}/"
}

cmd_stage() {
  echo "[orcd-pilot] stage compact S56 inputs"
  if [[ "${DRY_RUN}" == "1" ]]; then
    print_cmd "${LOCAL_REPO}/scripts/orcd/stage_s56_orcd_inputs.sh" --run
  else
    "${LOCAL_REPO}/scripts/orcd/stage_s56_orcd_inputs.sh" --run
  fi
}

cmd_stage_subset() {
  local subset="$1"
  echo "[orcd-pilot] stage compact S56 inputs subset=${subset}"
  if [[ "${DRY_RUN}" == "1" ]]; then
    print_cmd "${LOCAL_REPO}/scripts/orcd/stage_s56_orcd_inputs.sh" --run --subset "${subset}"
  else
    "${LOCAL_REPO}/scripts/orcd/stage_s56_orcd_inputs.sh" --run --subset "${subset}"
  fi
}

cmd_smoke() {
  submit "peak-smoke" "--export=ALL,TWIRL_N_INJECTIONS=100,TWIRL_PEAK_OUT_ROOT=reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training_orcd_smoke100,TWIRL_FINAL_TABLE=reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training_orcd_smoke100/s56_20k_injection_bls_peaks_orcd_smoke100.csv" "scripts/orcd/slurm_s56_peak_training_cpu.sbatch"
}

cmd_peak_full() {
  submit "peak-full" "--export=ALL,TWIRL_PEAK_OUT_ROOT=reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training_orcd_full,TWIRL_FINAL_TABLE=reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training_orcd_full/s56_20k_injection_bls_peaks_orcd_full.csv" "scripts/orcd/slurm_s56_peak_training_cpu.sbatch"
}

cmd_ranker() {
  submit "ranker" "" "scripts/orcd/slurm_s56_peak_ranker_train_cpu.sbatch"
}

cmd_short_branch() {
  submit "short-branch" "" "scripts/orcd/slurm_s56_short_branch_compare_cpu.sbatch"
}

cmd_apply() {
  submit "ranker-apply" "" "scripts/orcd/slurm_s56_peak_ranker_apply_cpu.sbatch"
}

cmd_allhost_peak() {
  submit "allhost-peak" "" "scripts/orcd/slurm_s56_allhost_peak_training_cpu.sbatch"
}

cmd_allhost_ranker() {
  submit "allhost-ranker" "--export=ALL,TWIRL_PEAK_TABLE=reports/stage5_validation/s56_allhost_predetrend_batman_periodradius_grid_sharded_pdo/peak_training_orcd/s56_allhost_injection_bls_peaks_orcd.csv,TWIRL_RANKER_OUT_DIR=reports/stage5_validation/s56_allhost_predetrend_batman_periodradius_grid_sharded_pdo/peak_ranker_orcd,TWIRL_GATE_OUT_DIR=reports/stage5_validation/s56_allhost_predetrend_batman_periodradius_grid_sharded_pdo/peak_training_gate_orcd,TWIRL_INJECTION_MANIFEST=data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_allhost_predetrend_batman_periodradius_grid_sharded/injection_manifest.csv,TWIRL_INJECTION_SUMMARY=data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_allhost_predetrend_batman_periodradius_grid_sharded/summary.json" "scripts/orcd/slurm_s56_peak_ranker_train_cpu.sbatch"
}

cmd_allhost_apply() {
  submit "allhost-ranker-apply" "--export=ALL,TWIRL_RANKER_MODEL=reports/stage5_validation/s56_allhost_predetrend_batman_periodradius_grid_sharded_pdo/peak_ranker_orcd/peak_ranker_model.npz,TWIRL_RANKER_SELECTED_OUT_DIR=reports/stage5_validation/s56_allhost_ranker_selected_real_candidates_orcd" "scripts/orcd/slurm_s56_peak_ranker_apply_cpu.sbatch"
}

cmd_predetrend_bls_smoke() {
  submit "predetrend-bls-smoke" "--exclude=node4900 --cpus-per-task=24 --mem=96G --export=ALL,TWIRL_PREDET_BLS_N_INJECTIONS=100,TWIRL_PREDET_BLS_N_PERIODS=5000,TWIRL_PREDET_BLS_N_PEAKS=10,TWIRL_PREDET_BLS_WORKERS=24,TWIRL_PREDET_BLS_OUT_DIR=reports/stage5_validation/s56_predetrend_detrending_bls_audit_orcd_smoke" "scripts/orcd/slurm_s56_predetrend_detrending_bls_audit_cpu.sbatch"
}

cmd_predetrend_bls_full() {
  submit "predetrend-bls-full" "--exclude=node4900 --cpus-per-task=96 --mem=384G --export=ALL,TWIRL_PREDET_BLS_N_INJECTIONS=3000,TWIRL_PREDET_BLS_N_PERIODS=50000,TWIRL_PREDET_BLS_N_PEAKS=20,TWIRL_PREDET_BLS_WORKERS=96,TWIRL_PREDET_BLS_OUT_DIR=reports/stage5_validation/s56_predetrend_detrending_bls_audit_orcd_full" "scripts/orcd/slurm_s56_predetrend_detrending_bls_audit_cpu.sbatch"
}

cmd_adpplus_smoke() {
  submit "adpplus-smoke" "--exclude=node4900 --cpus-per-task=24 --mem=96G --export=ALL,TWIRL_ADPPLUS_N_INJECTED=100,TWIRL_ADPPLUS_N_REAL=50,TWIRL_ADPPLUS_N_PERIODS=5000,TWIRL_ADPPLUS_N_PEAKS=10,TWIRL_ADPPLUS_WORKERS=24,TWIRL_ADPPLUS_OUT_DIR=reports/stage5_validation/s56_adpplus_bls_audit_orcd_smoke" "scripts/orcd/slurm_s56_adpplus_bls_audit_cpu.sbatch"
}

cmd_adpplus_full() {
  submit "adpplus-full" "--exclude=node4900 --cpus-per-task=96 --mem=384G --export=ALL,TWIRL_ADPPLUS_N_INJECTED=3000,TWIRL_ADPPLUS_N_REAL=2000,TWIRL_ADPPLUS_N_PERIODS=50000,TWIRL_ADPPLUS_N_PEAKS=20,TWIRL_ADPPLUS_WORKERS=96,TWIRL_ADPPLUS_OUT_DIR=reports/stage5_validation/s56_adpplus_bls_audit_orcd_full" "scripts/orcd/slurm_s56_adpplus_bls_audit_cpu.sbatch"
}

cmd_twoap_smoke() {
  submit "twoap-smoke" "--exclude=node4900 --cpus-per-task=8 --mem=32G --export=ALL,TWIRL_TWOAP_LIMIT=10,TWIRL_TWOAP_WORKERS=2,TWIRL_TWOAP_N_PERIODS=5000,TWIRL_TWOAP_OUT_DIR=reports/stage5_validation/s56_mixed_teacher_queue_pdo/twirl_vet_sheets_adp015q_orcd_smoke,TWIRL_TWOAP_METRICS_CSV=reports/stage5_validation/s56_mixed_teacher_queue_pdo/twirl_vet_metrics_adp015q_orcd_smoke.csv,TWIRL_TWOAP_OVERWRITE=1" "scripts/orcd/slurm_s56_two_aperture_vet_sheets_cpu.sbatch"
}

cmd_twoap_full() {
  submit "twoap-full" "--exclude=node4900 --cpus-per-task=48 --mem=192G --export=ALL,TWIRL_TWOAP_LIMIT=0,TWIRL_TWOAP_WORKERS=48,TWIRL_TWOAP_N_PERIODS=20000,TWIRL_TWOAP_OUT_DIR=reports/stage5_validation/s56_mixed_teacher_queue_pdo/twirl_vet_sheets_adp015q_orcd,TWIRL_TWOAP_METRICS_CSV=reports/stage5_validation/s56_mixed_teacher_queue_pdo/twirl_vet_metrics_adp015q_orcd.csv" "scripts/orcd/slurm_s56_two_aperture_vet_sheets_cpu.sbatch"
}

cmd_sync_twoap_smoke() {
  echo "[orcd-pilot] sync ORCD twoap-smoke outputs to local/PDO"
  if [[ "${DRY_RUN}" == "1" ]]; then
    print_cmd env \
      TWIRL_ORCD_TWOAP_OUT_DIR=reports/stage5_validation/s56_mixed_teacher_queue_pdo/twirl_vet_sheets_adp015q_orcd_smoke \
      TWIRL_ORCD_TWOAP_METRICS_CSV=reports/stage5_validation/s56_mixed_teacher_queue_pdo/twirl_vet_metrics_adp015q_orcd_smoke.csv \
      "${LOCAL_REPO}/scripts/orcd/sync_s56_orcd_twoap_outputs.sh" --run
  else
    env \
      TWIRL_ORCD_TWOAP_OUT_DIR=reports/stage5_validation/s56_mixed_teacher_queue_pdo/twirl_vet_sheets_adp015q_orcd_smoke \
      TWIRL_ORCD_TWOAP_METRICS_CSV=reports/stage5_validation/s56_mixed_teacher_queue_pdo/twirl_vet_metrics_adp015q_orcd_smoke.csv \
      bash "${LOCAL_REPO}/scripts/orcd/sync_s56_orcd_twoap_outputs.sh" --run
  fi
}

cmd_sync_twoap() {
  echo "[orcd-pilot] sync full ORCD two-aperture outputs to local/PDO"
  if [[ "${DRY_RUN}" == "1" ]]; then
    print_cmd "${LOCAL_REPO}/scripts/orcd/sync_s56_orcd_twoap_outputs.sh" --run
  else
    bash "${LOCAL_REPO}/scripts/orcd/sync_s56_orcd_twoap_outputs.sh" --run
  fi
}

cmd_sync_apply() {
  echo "[orcd-pilot] sync ORCD ranker-selected outputs to PDO"
  if [[ "${DRY_RUN}" == "1" ]]; then
    print_cmd "${LOCAL_REPO}/scripts/orcd/sync_s56_orcd_ranker_outputs_to_pdo.sh" --run
  else
    bash "${LOCAL_REPO}/scripts/orcd/sync_s56_orcd_ranker_outputs_to_pdo.sh" --run
  fi
}

cmd_monitor_apply() {
  echo "[orcd-pilot] monitor ORCD apply outputs and sync to PDO"
  if [[ "${DRY_RUN}" == "1" ]]; then
    print_cmd bash "${LOCAL_REPO}/scripts/orcd/monitor_s56_orcd_apply_to_pdo.sh"
  else
    bash "${LOCAL_REPO}/scripts/orcd/monitor_s56_orcd_apply_to_pdo.sh"
  fi
}

cmd_monitor_apply_leo_smoke() {
  echo "[orcd-pilot] monitor ORCD apply outputs, sync to PDO, and render a bounded PDO LEO smoke queue"
  if [[ "${DRY_RUN}" == "1" ]]; then
    print_cmd env \
      TWIRL_RUN_PDO_LEO=1 \
      TWIRL_PDO_LEO_N_REAL=50 \
      TWIRL_PDO_LEO_MAX_REPORTS=50 \
      TWIRL_PDO_LEO_OUT_DIR=reports/stage5_validation/s56_ranker_selected_real_leo_queue_orcd50_pdo \
      bash "${LOCAL_REPO}/scripts/orcd/monitor_s56_orcd_apply_to_pdo.sh"
  else
    env \
      TWIRL_RUN_PDO_LEO=1 \
      TWIRL_PDO_LEO_N_REAL=50 \
      TWIRL_PDO_LEO_MAX_REPORTS=50 \
      TWIRL_PDO_LEO_OUT_DIR=reports/stage5_validation/s56_ranker_selected_real_leo_queue_orcd50_pdo \
      bash "${LOCAL_REPO}/scripts/orcd/monitor_s56_orcd_apply_to_pdo.sh"
  fi
}

cmd_allhost_sync_apply() {
  echo "[orcd-pilot] sync all-host ORCD ranker-selected outputs to PDO"
  if [[ "${DRY_RUN}" == "1" ]]; then
    print_cmd env \
      TWIRL_ORCD_SELECTED_DIR=reports/stage5_validation/s56_allhost_ranker_selected_real_candidates_orcd \
      TWIRL_PDO_SELECTED_DIR=reports/stage5_validation/s56_allhost_ranker_selected_real_candidates_orcd \
      "${LOCAL_REPO}/scripts/orcd/sync_s56_orcd_ranker_outputs_to_pdo.sh" --run
  else
    env \
      TWIRL_ORCD_SELECTED_DIR=reports/stage5_validation/s56_allhost_ranker_selected_real_candidates_orcd \
      TWIRL_PDO_SELECTED_DIR=reports/stage5_validation/s56_allhost_ranker_selected_real_candidates_orcd \
      bash "${LOCAL_REPO}/scripts/orcd/sync_s56_orcd_ranker_outputs_to_pdo.sh" --run
  fi
}

cmd_allhost_monitor_apply() {
  echo "[orcd-pilot] monitor all-host ORCD apply outputs and sync to PDO"
  if [[ "${DRY_RUN}" == "1" ]]; then
    print_cmd env \
      TWIRL_ORCD_SELECTED_DIR=reports/stage5_validation/s56_allhost_ranker_selected_real_candidates_orcd \
      TWIRL_PDO_SELECTED_DIR=reports/stage5_validation/s56_allhost_ranker_selected_real_candidates_orcd \
      bash "${LOCAL_REPO}/scripts/orcd/monitor_s56_orcd_apply_to_pdo.sh"
  else
    env \
      TWIRL_ORCD_SELECTED_DIR=reports/stage5_validation/s56_allhost_ranker_selected_real_candidates_orcd \
      TWIRL_PDO_SELECTED_DIR=reports/stage5_validation/s56_allhost_ranker_selected_real_candidates_orcd \
      bash "${LOCAL_REPO}/scripts/orcd/monitor_s56_orcd_apply_to_pdo.sh"
  fi
}

cmd_h200_smoke() {
  submit "h200-smoke" "" "scripts/orcd/slurm_h200_smoke.sbatch"
}

cmd_h200_tensor_smoke() {
  submit "h200-tensor-smoke" "" "scripts/orcd/slurm_s56_tensor_smoke_h200.sbatch"
}

cmd_h200_torch_tensor_smoke() {
  submit "h200-torch-tensor-smoke" "--export=ALL,TWIRL_ORCD_PYTHON=/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56-torch/bin/python,TWIRL_REQUIRE_TORCH=1,TWIRL_TENSOR_OUT_DIR=reports/stage5_validation/s56_candidate_tensor_smoke_h200_torch" "scripts/orcd/slurm_s56_tensor_smoke_h200.sbatch"
}

cmd_stage_labels() {
  echo "[orcd-pilot] stage human-label products to ORCD"
  if [[ "${DRY_RUN}" == "1" ]]; then
    print_cmd "${LOCAL_REPO}/scripts/orcd/stage_s56_human_label_products_to_orcd.sh" --run
  else
    "${LOCAL_REPO}/scripts/orcd/stage_s56_human_label_products_to_orcd.sh" --run
  fi
}

cmd_h200_tensor_train_smoke() {
  submit "h200-tensor-train-smoke" "" "scripts/orcd/slurm_s56_tensor_train_smoke_h200.sbatch"
}

cmd_h200_tensor_train_synthetic_smoke() {
  submit "h200-tensor-train-synthetic-smoke" "--export=ALL,TWIRL_TENSOR_SYNTHETIC_LABEL_SMOKE=1,TWIRL_TENSOR_TRAIN_OUT_DIR=reports/stage5_validation/s56_candidate_tensor_train_synthetic_smoke_h200,TWIRL_TENSOR_TRAIN_EPOCHS=8" "scripts/orcd/slurm_s56_tensor_train_smoke_h200.sbatch"
}

cmd_status() {
  echo "[orcd-pilot] status"
  orcd "squeue -u \"\$(whoami)\" -p '${PARTITION}' || true; echo '--- logs ---'; ls -lt /orcd/data/mki_aryeh/001/twirl/logs 2>/dev/null | head -20 || true"
}

expand_command() {
  case "$1" in
    smoke-pilot)
      printf '%s\n' probe stage smoke h200-smoke
      ;;
    *)
      printf '%s\n' "$1"
      ;;
  esac
}

expanded=()
for command in "${COMMANDS[@]}"; do
  while IFS= read -r item; do
    expanded+=("${item}")
  done < <(expand_command "${command}")
done

echo "[orcd-pilot] local_repo=${LOCAL_REPO}"
echo "[orcd-pilot] orcd=${ORCD_HOST}:${ORCD_REPO}"
if [[ "${DRY_RUN}" == "1" ]]; then
  echo "[orcd-pilot] dry run. Pass --run to execute."
fi

require_socket

for command in "${expanded[@]}"; do
  case "${command}" in
    probe) cmd_probe ;;
    sync-code) cmd_sync_code ;;
    stage) cmd_stage ;;
    stage-balanced) cmd_stage_subset balanced ;;
    stage-allhost) cmd_stage_subset allhost ;;
    smoke) cmd_smoke ;;
    peak-full) cmd_peak_full ;;
    ranker) cmd_ranker ;;
    short-branch) cmd_short_branch ;;
    apply) cmd_apply ;;
    allhost-peak) cmd_allhost_peak ;;
    allhost-ranker) cmd_allhost_ranker ;;
    allhost-apply) cmd_allhost_apply ;;
    predetrend-bls-smoke) cmd_predetrend_bls_smoke ;;
    predetrend-bls-full) cmd_predetrend_bls_full ;;
    adpplus-smoke) cmd_adpplus_smoke ;;
    adpplus-full) cmd_adpplus_full ;;
    twoap-smoke) cmd_twoap_smoke ;;
    twoap-full) cmd_twoap_full ;;
    sync-twoap-smoke) cmd_sync_twoap_smoke ;;
    sync-twoap) cmd_sync_twoap ;;
    allhost-sync-apply) cmd_allhost_sync_apply ;;
    allhost-monitor-apply) cmd_allhost_monitor_apply ;;
    sync-apply) cmd_sync_apply ;;
    monitor-apply) cmd_monitor_apply ;;
    monitor-apply-leo-smoke) cmd_monitor_apply_leo_smoke ;;
    h200-smoke) cmd_h200_smoke ;;
    h200-tensor-smoke) cmd_h200_tensor_smoke ;;
    h200-torch-tensor-smoke) cmd_h200_torch_tensor_smoke ;;
    stage-labels) cmd_stage_labels ;;
    h200-tensor-train-smoke) cmd_h200_tensor_train_smoke ;;
    h200-tensor-train-synthetic-smoke) cmd_h200_tensor_train_synthetic_smoke ;;
    status) cmd_status ;;
    *)
      echo "[orcd-pilot] unknown command: ${command}" >&2
      usage >&2
      exit 2
      ;;
  esac
done
