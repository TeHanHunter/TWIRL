#!/usr/bin/env bash
# Copy only the selected five teacher-v1 checkpoints from ORCD to pdogpu6.
set -euo pipefail

ORCD_HOST="${ORCD_HOST:-tehan@orcd-login.mit.edu}"
ORCD_CONTROL_PATH="${ORCD_CONTROL_PATH:-$HOME/.ssh/cm/%r@%h:%p}"
ORCD_REPO="${ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PDO_HOST="${PDO_HOST:-pdogpu6}"
SOURCE_ROOT="${TWIRL_TEACHER_SOURCE_ROOT:-reports/stage5_validation/s56_label_adjudication_real343/retrained/s56_harmonic_cnn_v1/shape_plus_periodogram_bls}"
PDO_DEST="${TWIRL_PDO_TEACHER_CHECKPOINT_ROOT:-/pdo/users/tehan/twirl_stage5/s56_A2v1_teacher_search_v1/checkpoints/shape_plus_periodogram_bls}"

ORCD_SSH=(
  ssh
  -o BatchMode=yes
  -o PasswordAuthentication=no
  -o KbdInteractiveAuthentication=no
  -o NumberOfPasswordPrompts=0
  -o ConnectTimeout=10
  -o ControlMaster=auto
  -o ControlPath="${ORCD_CONTROL_PATH}"
)
PDO_SSH=(
  ssh
  -o BatchMode=yes
  -o PasswordAuthentication=no
  -o KbdInteractiveAuthentication=no
  -o NumberOfPasswordPrompts=0
  -o ConnectTimeout=10
)

"${ORCD_SSH[@]}" "${ORCD_HOST}" true >/dev/null
"${PDO_SSH[@]}" "${PDO_HOST}" true >/dev/null

files=()
for fold in 0 1 2 3 4; do
  files+=("${SOURCE_ROOT}/fold_${fold}/teacher.pt")
done
quoted_files=()
for path in "${files[@]}"; do
  printf -v quoted '%q' "${path}"
  quoted_files+=("${quoted}")
done

"${PDO_SSH[@]}" "${PDO_HOST}" "mkdir -p '${PDO_DEST}'"
"${ORCD_SSH[@]}" "${ORCD_HOST}" \
  "cd '${ORCD_REPO}' && for p in ${quoted_files[*]}; do test -s \"\$p\" || { echo \"missing \$p\" >&2; exit 4; }; done && tar -cf - ${quoted_files[*]}" \
  | "${PDO_SSH[@]}" "${PDO_HOST}" \
      "tmp='/pdo/users/tehan/twirl_stage5/s56_A2v1_teacher_search_v1/checkpoint_stage'; mkdir -p \"\$tmp\"; tar -C \"\$tmp\" -xf -; for fold in 0 1 2 3 4; do src=\"\$tmp/${SOURCE_ROOT}/fold_\$fold/teacher.pt\"; dst='${PDO_DEST}'/fold_\$fold/teacher.pt; mkdir -p \"\$(dirname \"\$dst\")\"; if test -e \"\$dst\"; then cmp -s \"\$src\" \"\$dst\" || { echo \"refusing checkpoint mismatch: \$dst\" >&2; exit 6; }; else cp \"\$src\" \"\$dst\"; fi; done"

"${PDO_SSH[@]}" "${PDO_HOST}" \
  "for fold in 0 1 2 3 4; do test -s '${PDO_DEST}'/fold_\$fold/teacher.pt || exit 5; sha256sum '${PDO_DEST}'/fold_\$fold/teacher.pt; done"
