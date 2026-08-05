#!/usr/bin/env bash
# Copy the selected native-v2 teacher-v1 ensemble and its provenance manifest.
set -euo pipefail

ORCD_HOST="${ORCD_HOST:-tehan@orcd-login.mit.edu}"
ORCD_CONTROL_PATH="${ORCD_CONTROL_PATH:-$HOME/.ssh/cm/%r@%h:%p}"
ORCD_REPO="${ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
ORCD_PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56-torch/bin/python}"
PDO_HOST="${PDO_HOST:-pdogpu6}"
PDO_REPO="${TWIRL_PDO_TEACHER_REPO:-/pdo/users/tehan/TWIRL_teacher_search}"
PDO_PYTHON="${TWIRL_PDO_TORCH_PYTHON:-/pdo/users/tehan/envs/twirl-teacher-pdo-v2/bin/python}"
SOURCE_MODEL_ROOT="${TWIRL_TEACHER_SOURCE_MODEL_ROOT:-reports/stage5_validation/s56_label_adjudication_real343/retrained/s56_harmonic_cnn_v1_native_v2}"
PROFILE="${TWIRL_TEACHER_PROFILE:-shape_plus_periodogram_bls}"
SOURCE_ROOT="${SOURCE_MODEL_ROOT}/${PROFILE}"
SOURCE_MANIFEST="${SOURCE_MODEL_ROOT}/selected_checkpoint_manifest.json"
TRAINING_TABLE="${TWIRL_TEACHER_TRAINING_TABLE:-reports/stage5_validation/s56_label_adjudication_real343/adjudicated_training_table/human_vetting_training_table_adjudicated.csv}"
NATIVE_H5="${TWIRL_HARMONIC_NATIVE_H5:-reports/stage5_validation/s56_label_adjudication_real343/native_inputs/s56_adp_raw_pair_v2.h5}"
PDO_MODEL_DEST="${TWIRL_PDO_TEACHER_MODEL_ROOT:-/pdo/users/tehan/twirl_stage5/s56_A2v1_teacher_search_v2/checkpoints}"
PDO_DEST="${TWIRL_PDO_TEACHER_CHECKPOINT_ROOT:-${PDO_MODEL_DEST}/${PROFILE}}"

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
"${PDO_SSH[@]}" "${PDO_HOST}" \
  "test -x '${PDO_PYTHON}' && test -f '${PDO_REPO}/scripts/stage5_validation/verify_s56_harmonic_teacher_checkpoint_manifest.py'"

"${ORCD_SSH[@]}" "${ORCD_HOST}" \
  "cd '${ORCD_REPO}' && PYTHONPATH=src '${ORCD_PYTHON}' scripts/stage5_validation/verify_s56_harmonic_teacher_checkpoint_manifest.py --manifest '${SOURCE_MANIFEST}' --training-table '${TRAINING_TABLE}' --native-h5 '${NATIVE_H5}' --expected-profile '${PROFILE}' >/dev/null"

files=("${SOURCE_MANIFEST}")
for fold in 0 1 2 3 4; do
  files+=("${SOURCE_ROOT}/fold_${fold}/teacher.pt")
done
quoted_files=()
for path in "${files[@]}"; do
  printf -v quoted '%q' "${path}"
  quoted_files+=("${quoted}")
done

"${PDO_SSH[@]}" "${PDO_HOST}" "mkdir -p '${PDO_MODEL_DEST}' '${PDO_DEST}'"
"${ORCD_SSH[@]}" "${ORCD_HOST}" \
  "cd '${ORCD_REPO}' && for p in ${quoted_files[*]}; do test -s \"\$p\" || { echo \"missing \$p\" >&2; exit 4; }; done && tar -cf - ${quoted_files[*]}" \
  | "${PDO_SSH[@]}" "${PDO_HOST}" \
      "tmp='/pdo/users/tehan/twirl_stage5/s56_A2v1_teacher_search_v2/checkpoint_stage'; mkdir -p \"\$tmp\"; tar -C \"\$tmp\" -xf -; manifest_src=\"\$tmp/${SOURCE_MANIFEST}\"; manifest_dst='${PDO_MODEL_DEST}'/selected_checkpoint_manifest.json; if test -e \"\$manifest_dst\"; then cmp -s \"\$manifest_src\" \"\$manifest_dst\" || { echo \"refusing manifest mismatch: \$manifest_dst\" >&2; exit 6; }; else cp \"\$manifest_src\" \"\$manifest_dst\"; fi; for fold in 0 1 2 3 4; do src=\"\$tmp/${SOURCE_ROOT}/fold_\$fold/teacher.pt\"; dst='${PDO_DEST}'/fold_\$fold/teacher.pt; mkdir -p \"\$(dirname \"\$dst\")\"; if test -e \"\$dst\"; then cmp -s \"\$src\" \"\$dst\" || { echo \"refusing checkpoint mismatch: \$dst\" >&2; exit 6; }; else cp \"\$src\" \"\$dst\"; fi; done"

"${PDO_SSH[@]}" "${PDO_HOST}" \
  "cd '${PDO_REPO}' && PYTHONPATH=src '${PDO_PYTHON}' scripts/stage5_validation/verify_s56_harmonic_teacher_checkpoint_manifest.py --manifest '${PDO_MODEL_DEST}/selected_checkpoint_manifest.json' --checkpoint-root '${PDO_DEST}' --expected-profile '${PROFILE}' >/dev/null"

"${PDO_SSH[@]}" "${PDO_HOST}" \
  "test -s '${PDO_MODEL_DEST}'/selected_checkpoint_manifest.json || exit 5; sha256sum '${PDO_MODEL_DEST}'/selected_checkpoint_manifest.json; for fold in 0 1 2 3 4; do test -s '${PDO_DEST}'/fold_\$fold/teacher.pt || exit 5; sha256sum '${PDO_DEST}'/fold_\$fold/teacher.pt; done"
