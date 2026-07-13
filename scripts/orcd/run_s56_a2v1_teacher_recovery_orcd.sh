#!/usr/bin/env bash
# Submit and inspect the gated S56 A2v1 injection/BLS/Teacher-v1 workflow.
set -euo pipefail

ACTION="${1:-status}"
HOST="${ORCD_HOST:-tehan@orcd-login.mit.edu}"
CONTROL_PATH="${ORCD_CONTROL_PATH:-${HOME}/.ssh/cm/%r@%h:%p}"
REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
OUT_ROOT="${TWIRL_RECOVERY_ROOT:-${REPO}/reports/stage5_validation/s56_A2v1_teacher_v1_recovery}"
INPUT_VERIFICATION="${TWIRL_RECOVERY_INPUT_VERIFICATION:-/orcd/data/mki_aryeh/001/twirl/data_local/stage1_lightcurves/s56_A2v1/input_verification.json}"
CHECKPOINT_ROOT="${TWIRL_TEACHER_CHECKPOINT_ROOT:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL/reports/stage5_validation/s56_label_adjudication_real343/retrained/s56_harmonic_cnn_v1/shape_plus_periodogram_bls}"
SBATCH_EXPORT="ALL,TWIRL_ORCD_REPO=${REPO},TWIRL_RECOVERY_ROOT=${OUT_ROOT},TWIRL_TEACHER_CHECKPOINT_ROOT=${CHECKPOINT_ROOT}"
SSH=(
  ssh
  -o BatchMode=yes
  -o PasswordAuthentication=no
  -o KbdInteractiveAuthentication=no
  -o NumberOfPasswordPrompts=0
  -o ControlMaster=no
  -o ControlPath="${CONTROL_PATH}"
  "${HOST}"
)

probe() {
  "${SSH[@]}" "hostname"
}

require_remote_file() {
  local path="$1"
  "${SSH[@]}" "test -s '${path}'"
}

case "${ACTION}" in
  probe)
    probe
    "${SSH[@]}" "sinfo -h -p pg_mki_aryeh -o '%P %D %t %G'"
    ;;
  prepare)
    probe >/dev/null
    "${SSH[@]}" "cd '${REPO}' && sbatch --parsable --export='${SBATCH_EXPORT}' scripts/orcd/slurm_s56_a2v1_recovery_prepare_cpu.sbatch"
    ;;
  smoke)
    probe >/dev/null
    require_remote_file "${OUT_ROOT}/schedule/adp_roundtrip_parity_summary.json"
    "${SSH[@]}" "grep -q '\"passed\": true' '${OUT_ROOT}/schedule/adp_roundtrip_parity_summary.json'"
    "${SSH[@]}" "cd '${REPO}' && sbatch --parsable --array=0 --export='${SBATCH_EXPORT},TWIRL_RECOVERY_MODE=smoke' scripts/orcd/slurm_s56_a2v1_recovery_injection_bls_cpu.sbatch"
    ;;
  smoke-gpu)
    probe >/dev/null
    require_remote_file "${OUT_ROOT}/smoke100/peak_shards/peaks_00_verification.json"
    "${SSH[@]}" "grep -q '\"passed\": true' '${OUT_ROOT}/smoke100/peak_shards/peaks_00_verification.json'"
    require_remote_file "${OUT_ROOT}/smoke100/teacher_inputs/native_00.teacher_summary.json"
    "${SSH[@]}" "cd '${REPO}' && sbatch --parsable --export='${SBATCH_EXPORT},TWIRL_RECOVERY_MODE=smoke' scripts/orcd/slurm_s56_a2v1_recovery_teacher_h200.sbatch"
    ;;
  full)
    probe >/dev/null
    require_remote_file "${OUT_ROOT}/smoke100/teacher_scores_128.summary.json"
    "${SSH[@]}" "cd '${REPO}' && /orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python scripts/stage5_validation/verify_s56_a2v1_recovery_smoke.py --recovery-root '${OUT_ROOT}' --input-verification '${INPUT_VERIFICATION}' --out-json '${OUT_ROOT}/smoke100/smoke_acceptance.json'"
    INJECT_JOB="$("${SSH[@]}" "cd '${REPO}' && sbatch --parsable --export='${SBATCH_EXPORT},TWIRL_RECOVERY_MODE=full' scripts/orcd/slurm_s56_a2v1_recovery_injection_bls_cpu.sbatch")"
    NATIVE_JOB="$("${SSH[@]}" "cd '${REPO}' && sbatch --parsable --dependency=afterok:${INJECT_JOB} --export='${SBATCH_EXPORT}' scripts/orcd/slurm_s56_a2v1_recovery_teacher_prepare_cpu.sbatch")"
    MERGE_JOB="$("${SSH[@]}" "cd '${REPO}' && sbatch --parsable --dependency=afterok:${NATIVE_JOB} --export='${SBATCH_EXPORT}' scripts/orcd/slurm_s56_a2v1_recovery_merge_cpu.sbatch")"
    TEACHER_JOB="$("${SSH[@]}" "cd '${REPO}' && sbatch --parsable --dependency=afterok:${MERGE_JOB} --export='${SBATCH_EXPORT},TWIRL_RECOVERY_MODE=full' scripts/orcd/slurm_s56_a2v1_recovery_teacher_h200.sbatch")"
    FINAL_JOB="$("${SSH[@]}" "cd '${REPO}' && sbatch --parsable --dependency=afterok:${TEACHER_JOB} --export='${SBATCH_EXPORT}' scripts/orcd/slurm_s56_a2v1_recovery_finalize_cpu.sbatch")"
    printf 'injection_bls=%s\nnative=%s\nmerge=%s\nteacher=%s\nfinal=%s\n' \
      "${INJECT_JOB}" "${NATIVE_JOB}" "${MERGE_JOB}" "${TEACHER_JOB}" "${FINAL_JOB}"
    ;;
  status)
    probe >/dev/null
    "${SSH[@]}" "squeue -u tehan -o '%.18i %.28j %.2t %.10M %.4D %R'"
    "${SSH[@]}" "find '${OUT_ROOT}' -maxdepth 4 -type f -name '*summary.json' -o -name '*verification.json' 2>/dev/null | sort"
    ;;
  *)
    echo "usage: $0 {probe|prepare|smoke|smoke-gpu|full|status}" >&2
    exit 2
    ;;
esac
