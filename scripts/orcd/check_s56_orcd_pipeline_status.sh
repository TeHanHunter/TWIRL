#!/usr/bin/env bash
# Read-only status check for the S56 ORCD peak-table/ranker chain.
set -euo pipefail

ORCD_HOST="${ORCD_HOST:-tehan@orcd-login.mit.edu}"
ORCD_REPO="${ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
ORCD_CONTROL_PATH="${ORCD_CONTROL_PATH:-$HOME/.ssh/cm/%r@%h:%p}"
PEAK_JOB="${TWIRL_ORCD_PEAK_JOB:-16869261}"
RANKER_JOB="${TWIRL_ORCD_RANKER_JOB:-16869487}"
APPLY_JOB="${TWIRL_ORCD_APPLY_JOB:-16869492}"

ORCD_SSH=(
  ssh
  -o BatchMode=yes
  -o PasswordAuthentication=no
  -o KbdInteractiveAuthentication=no
  -o NumberOfPasswordPrompts=0
  -o ConnectTimeout=10
  -o ConnectionAttempts=1
  -o ControlMaster=auto
  -o ControlPath="${ORCD_CONTROL_PATH}"
)
ORCD_SOCKET_CHECK=(
  ssh
  -O check
  -o BatchMode=yes
  -o PasswordAuthentication=no
  -o KbdInteractiveAuthentication=no
  -o NumberOfPasswordPrompts=0
  -o ConnectTimeout=5
  -o ConnectionAttempts=1
  -S "${ORCD_CONTROL_PATH}"
)

jobs_csv="${PEAK_JOB},${RANKER_JOB},${APPLY_JOB}"

if ! socket_check_err="$("${ORCD_SOCKET_CHECK[@]}" "${ORCD_HOST}" 2>&1 >/dev/null)"; then
  printf '[orcd-status] existing ORCD control socket is not usable; refusing to start a new interactive SSH attempt.\n' >&2
  printf '%s\n' "${socket_check_err}" >&2
  exit 2
fi

"${ORCD_SSH[@]}" "${ORCD_HOST}" "
set -euo pipefail
echo '[orcd-status] host='\"\$(hostname)\"' date='\"\$(date -Is)\"
echo '[orcd-status] jobs=${jobs_csv}'
echo '--- queue ---'
squeue -j '${jobs_csv}' -o '%.18i %.28j %.9P %.2t %.10M %.10l %.6D %.40R' || true
echo '--- accounting ---'
sacct -j '${jobs_csv}' --format=JobID,JobName%28,State,ExitCode,Elapsed,ReqTRES%80,AllocTRES%80 -P 2>/dev/null || true
echo '--- gpu/gres guardrail ---'
tres_report=\$(sacct -j '${jobs_csv}' --format=JobID,ReqTRES%80,AllocTRES%80 -P 2>/dev/null || true)
if printf '%s\n' \"\${tres_report}\" | grep -Eiq 'gres/gpu|gpu:h200'; then
  echo 'TRACKED_JOB_GPU_GRES_PRESENT'
  printf '%s\n' \"\${tres_report}\" | grep -Ei 'gres/gpu|gpu:h200' || true
else
  echo 'tracked_jobs_gpu_gres=none'
fi
cd '${ORCD_REPO}'
out_root='reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training_orcd_full'
ranker_dir='reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_ranker_orcd'
gate_dir='reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training_gate_orcd'
apply_dir='reports/stage5_validation/s56_ranker_selected_real_candidates_orcd'
echo '--- peak progress ---'
cat \"\${out_root}/chunk_ids/manifest.txt\" 2>/dev/null || true
n_chunks=\$(awk -F= '\$1 == \"n_chunks\" {print \$2}' \"\${out_root}/chunk_ids/manifest.txt\" 2>/dev/null || true)
done_chunks=\$(find \"\${out_root}/chunks\" -maxdepth 2 -name injection_bls_peaks_summary.json 2>/dev/null | wc -l | tr -d ' ')
echo \"done_chunks=\${done_chunks}\"
if [[ -n \"\${n_chunks:-}\" && \"\${n_chunks}\" != \"0\" ]]; then
  awk -v done=\"\${done_chunks}\" -v total=\"\${n_chunks}\" 'BEGIN {printf \"done_fraction=%.1f%%\\n\", 100.0 * done / total}'
fi
if [[ -d \"\${out_root}/chunks\" ]]; then
  active_chunks=\$(find \"\${out_root}/chunks\" -mindepth 1 -maxdepth 1 -type d \
    -name 'chunk_*' ! -exec test -e '{}/injection_bls_peaks_summary.json' \\; -print 2>/dev/null |
    sort |
    while read -r chunk_dir; do
      if [[ -f \"\${chunk_dir}/run.log\" ]]; then
        basename \"\${chunk_dir}\"
      fi
    done |
    paste -sd, -)
  echo \"active_or_incomplete_chunks=\${active_chunks:-none}\"
  latest_chunk=\$(find \"\${out_root}/chunks\" -maxdepth 2 -name run.log -printf '%T@ %p\n' 2>/dev/null | sort -n | tail -1 | cut -d' ' -f2-)
  if [[ -n \"\${latest_chunk:-}\" ]]; then
    echo \"latest_chunk_log=\${latest_chunk}\"
    tail -n 20 \"\${latest_chunk}\" || true
  fi
fi
echo '--- expected artifacts ---'
for p in \
  \"\${out_root}/s56_20k_injection_bls_peaks_orcd_full.csv\" \
  \"\${out_root}/s56_20k_injection_bls_peaks_orcd_full_verification.json\" \
  \"\${ranker_dir}/summary.json\" \
  \"\${gate_dir}/host_coverage/summary.json\" \
  \"\${gate_dir}/summary.json\" \
  \"\${gate_dir}/failure_modes/summary.json\" \
  \"\${apply_dir}/real_peak_table_verification.json\" \
  \"\${apply_dir}/selected_ephemerides.csv\" \
  \"\${apply_dir}/selected_ephemerides_verification.json\" \
  \"\${apply_dir}/ranker_selection_summary/summary.json\"; do
  if [[ -e \"\${p}\" ]]; then
    stat -c 'ARTIFACT %s %n' \"\${p}\"
  else
    echo \"MISSING \${p}\"
  fi
done
echo '--- handoff audit ---'
\"\${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}\" \
  scripts/stage5_validation/audit_s56_peak_handoff_status.py \
  --layout orcd \
  --out-dir reports/stage5_validation/s56_peak_handoff_audit_orcd \
  2>/dev/null | tail -n 18 || true
"
