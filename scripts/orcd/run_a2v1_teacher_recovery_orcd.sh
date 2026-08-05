#!/usr/bin/env bash
# Submit and inspect a gated sector A2v1 injection/BLS/Teacher-v1 workflow.
set -euo pipefail

ACTION="${1:-status}"
SECTOR="${TWIRL_RECOVERY_SECTOR:-56}"
SECTOR_TAG="s${SECTOR}"
HOST="${ORCD_HOST:-tehan@orcd-login.mit.edu}"
CONTROL_PATH="${ORCD_CONTROL_PATH:-${HOME}/.ssh/cm/%r@%h:%p}"
REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
INPUT_ROOT="${TWIRL_A2V1_INPUT_ROOT:-/orcd/data/mki_aryeh/001/twirl/data_local/stage1_lightcurves/${SECTOR_TAG}_A2v1}"
OUT_ROOT="${TWIRL_RECOVERY_ROOT:-${REPO}/reports/stage5_validation/${SECTOR_TAG}_A2v1_teacher_v1_recovery}"
CONFIG="${TWIRL_RECOVERY_CONFIG:-${REPO}/configs/injections/${SECTOR_TAG}_a2v1_teacher_v1_recovery_v1.yaml}"
INPUT_VERIFICATION="${TWIRL_RECOVERY_INPUT_VERIFICATION:-${INPUT_ROOT}/input_verification.json}"
CADENCE_REFERENCE="${TWIRL_A2V1_CADENCE_REFERENCE:-${INPUT_ROOT}/${SECTOR_TAG}_A2v1_cadence_reference.csv}"
CADENCE_REFERENCE_MANIFEST="${TWIRL_A2V1_CADENCE_REFERENCE_MANIFEST:-${INPUT_ROOT}/${SECTOR_TAG}_A2v1_cadence_reference_manifest.json}"
CHECKPOINT_ROOT="${TWIRL_TEACHER_CHECKPOINT_ROOT:-${REPO}/reports/stage5_validation/s56_label_adjudication_real343/retrained/s56_harmonic_cnn_v1_native_v2/shape_plus_periodogram_bls}"
CHECKPOINT_MANIFEST="${TWIRL_TEACHER_CHECKPOINT_MANIFEST:-$(dirname "${CHECKPOINT_ROOT}")/selected_checkpoint_manifest.json}"
LOG_ROOT="${TWIRL_RECOVERY_LOG_ROOT:-/orcd/data/mki_aryeh/001/twirl/logs/${SECTOR_TAG}_a2v1_recovery}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
PRIOR_SCHEDULE="${TWIRL_RECOVERY_PRIOR_SCHEDULE:-}"
HOST_OVERLAP_AUDIT_TABLE="${TWIRL_RECOVERY_HOST_OVERLAP_AUDIT_TABLE:-}"
if (( SECTOR > 56 )) && [[ -z "${HOST_OVERLAP_AUDIT_TABLE}" ]]; then
  HOST_OVERLAP_AUDIT_TABLE="${REPO}/reports/stage5_validation/s56_A2v1_teacher_v1_recovery/schedule/injection_schedule.parquet"
fi
SBATCH_EXPORT="ALL,TWIRL_ORCD_REPO=${REPO},TWIRL_RECOVERY_SECTOR=${SECTOR},TWIRL_A2V1_INPUT_ROOT=${INPUT_ROOT},TWIRL_RECOVERY_ROOT=${OUT_ROOT},TWIRL_RECOVERY_CONFIG=${CONFIG},TWIRL_RECOVERY_LOG_ROOT=${LOG_ROOT},TWIRL_TEACHER_CHECKPOINT_ROOT=${CHECKPOINT_ROOT},TWIRL_TEACHER_CHECKPOINT_MANIFEST=${CHECKPOINT_MANIFEST},TWIRL_RECOVERY_PRIOR_SCHEDULE=${PRIOR_SCHEDULE},TWIRL_RECOVERY_HOST_OVERLAP_AUDIT_TABLE=${HOST_OVERLAP_AUDIT_TABLE},TWIRL_A2V1_CADENCE_REFERENCE=${CADENCE_REFERENCE},TWIRL_A2V1_CADENCE_REFERENCE_MANIFEST=${CADENCE_REFERENCE_MANIFEST}"
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

if (( SECTOR < 56 )); then
  echo "A2v1 recovery sectors must be >= 56" >&2
  exit 2
fi

if (( SECTOR == 57 )); then
  echo "Blocked: the legacy S57 recovery contract is retired. Build the observation-keyed, per-sector quality-aware S56-S62 teacher contract instead." >&2
  exit 64
fi

probe() {
  "${SSH[@]}" "hostname"
}

prepare_submit_root() {
  "${SSH[@]}" "mkdir -p '${LOG_ROOT}'"
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
  rebuild-input)
    probe >/dev/null
    prepare_submit_root
    "${SSH[@]}" "cd '${REPO}' && sbatch --parsable --job-name='${SECTOR_TAG}-a2rec-compact' --output='${LOG_ROOT}/%x-%j.out' --error='${LOG_ROOT}/%x-%j.err' --export='${SBATCH_EXPORT}' scripts/orcd/slurm_s56_a2v1_recovery_rebuild_compact_cpu.sbatch"
    ;;
  prepare)
    probe >/dev/null
    prepare_submit_root
    "${SSH[@]}" "cd '${REPO}' && sbatch --parsable --job-name='${SECTOR_TAG}-a2rec-prep' --output='${LOG_ROOT}/%x-%j.out' --error='${LOG_ROOT}/%x-%j.err' --export='${SBATCH_EXPORT}' scripts/orcd/slurm_s56_a2v1_recovery_prepare_cpu.sbatch"
    ;;
  smoke)
    probe >/dev/null
    prepare_submit_root
    require_remote_file "${OUT_ROOT}/schedule/adp_roundtrip_parity_summary.json"
    require_remote_file "${CADENCE_REFERENCE}"
    require_remote_file "${CADENCE_REFERENCE_MANIFEST}"
    "${SSH[@]}" "grep -q '\"passed\": true' '${OUT_ROOT}/schedule/adp_roundtrip_parity_summary.json'"
    "${SSH[@]}" "cd '${REPO}' && sbatch --parsable --array=0 --job-name='${SECTOR_TAG}-a2rec-bls' --output='${LOG_ROOT}/%x-%A_%a.out' --error='${LOG_ROOT}/%x-%A_%a.err' --export='${SBATCH_EXPORT},TWIRL_RECOVERY_MODE=smoke' scripts/orcd/slurm_s56_a2v1_recovery_injection_bls_cpu.sbatch"
    ;;
  smoke-gpu)
    probe >/dev/null
    prepare_submit_root
    require_remote_file "${OUT_ROOT}/smoke100/peak_shards/peaks_00_verification.json"
    require_remote_file "${CHECKPOINT_MANIFEST}"
    "${SSH[@]}" "grep -q '\"passed\": true' '${OUT_ROOT}/smoke100/peak_shards/peaks_00_verification.json'"
    require_remote_file "${OUT_ROOT}/smoke100/teacher_inputs/native_00.teacher_summary.json"
    "${SSH[@]}" "cd '${REPO}' && sbatch --parsable --job-name='${SECTOR_TAG}-a2rec-teacher' --output='${LOG_ROOT}/%x-%j.out' --error='${LOG_ROOT}/%x-%j.err' --export='${SBATCH_EXPORT},TWIRL_RECOVERY_MODE=smoke' scripts/orcd/slurm_s56_a2v1_recovery_teacher_h200.sbatch"
    ;;
  full)
    probe >/dev/null
    prepare_submit_root
    require_remote_file "${OUT_ROOT}/smoke100/teacher_scores_128.summary.json"
    require_remote_file "${CHECKPOINT_MANIFEST}"
    require_remote_file "${CADENCE_REFERENCE}"
    require_remote_file "${CADENCE_REFERENCE_MANIFEST}"
    "${SSH[@]}" "cd '${REPO}' && '${PYTHON}' scripts/stage5_validation/verify_s56_a2v1_recovery_smoke.py --recovery-root '${OUT_ROOT}' --input-verification '${INPUT_VERIFICATION}' --checkpoint-manifest '${CHECKPOINT_MANIFEST}' --out-json '${OUT_ROOT}/smoke100/smoke_acceptance.json'"
    "${SSH[@]}" "cd '${REPO}' && git rev-parse HEAD > '${OUT_ROOT}/schedule/producer_git_sha.txt'"
    INJECT_JOB="$("${SSH[@]}" "cd '${REPO}' && sbatch --parsable --job-name='${SECTOR_TAG}-a2rec-bls' --output='${LOG_ROOT}/%x-%A_%a.out' --error='${LOG_ROOT}/%x-%A_%a.err' --export='${SBATCH_EXPORT},TWIRL_RECOVERY_MODE=full' scripts/orcd/slurm_s56_a2v1_recovery_injection_bls_cpu.sbatch")"
    NATIVE_JOB="$("${SSH[@]}" "cd '${REPO}' && sbatch --parsable --job-name='${SECTOR_TAG}-a2rec-native' --output='${LOG_ROOT}/%x-%A_%a.out' --error='${LOG_ROOT}/%x-%A_%a.err' --dependency=afterok:${INJECT_JOB} --export='${SBATCH_EXPORT}' scripts/orcd/slurm_s56_a2v1_recovery_teacher_prepare_cpu.sbatch")"
    MERGE_JOB="$("${SSH[@]}" "cd '${REPO}' && sbatch --parsable --job-name='${SECTOR_TAG}-a2rec-merge' --output='${LOG_ROOT}/%x-%j.out' --error='${LOG_ROOT}/%x-%j.err' --dependency=afterok:${NATIVE_JOB} --export='${SBATCH_EXPORT}' scripts/orcd/slurm_s56_a2v1_recovery_merge_cpu.sbatch")"
    TEACHER_JOB="$("${SSH[@]}" "cd '${REPO}' && sbatch --parsable --job-name='${SECTOR_TAG}-a2rec-teacher' --output='${LOG_ROOT}/%x-%j.out' --error='${LOG_ROOT}/%x-%j.err' --dependency=afterok:${MERGE_JOB} --export='${SBATCH_EXPORT},TWIRL_RECOVERY_MODE=full' scripts/orcd/slurm_s56_a2v1_recovery_teacher_h200.sbatch")"
    FINAL_JOB="$("${SSH[@]}" "cd '${REPO}' && sbatch --parsable --job-name='${SECTOR_TAG}-a2rec-final' --output='${LOG_ROOT}/%x-%j.out' --error='${LOG_ROOT}/%x-%j.err' --dependency=afterok:${TEACHER_JOB} --export='${SBATCH_EXPORT}' scripts/orcd/slurm_s56_a2v1_recovery_finalize_cpu.sbatch")"
    printf 'injection_bls=%s\nnative=%s\nmerge=%s\nteacher=%s\nfinal=%s\n' \
      "${INJECT_JOB}" "${NATIVE_JOB}" "${MERGE_JOB}" "${TEACHER_JOB}" "${FINAL_JOB}"
    ;;
  status)
    probe >/dev/null
    "${SSH[@]}" "squeue -u tehan -o '%.18i %.28j %.2t %.10M %.4D %R'"
    "${SSH[@]}" "find '${OUT_ROOT}' -maxdepth 4 -type f \( -name '*summary.json' -o -name '*verification.json' \) 2>/dev/null | sort"
    ;;
  *)
    echo "usage: $0 {probe|rebuild-input|prepare|smoke|smoke-gpu|full|status}" >&2
    exit 2
    ;;
esac
