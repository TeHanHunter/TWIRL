#!/bin/bash
# Run the S56 recovery50 teacher/student smoke pipeline on ORCD.
set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"

QUEUE="${TWIRL_RECOVERY50_QUEUE:-reports/stage5_validation/s56_recovery50_teacher_queue/review_queue_1k.csv}"
LABELS="${TWIRL_RECOVERY50_LABELS:-reports/stage5_validation/s56_recovery50_teacher_queue/human_labels_vetted.csv}"
ROOT="${TWIRL_RECOVERY50_ROOT:-reports/stage5_validation/s56_recovery50_teacher_queue}"
COMPACT_LC="${TWIRL_RECOVERY50_COMPACT_LC:-data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp_lc_export_pdo.h5}"
HLSP_ROOT="${TWIRL_RECOVERY50_HLSP_ROOT:-data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare}"
INJECTION_H5="${TWIRL_RECOVERY50_INJECTION_H5:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/recovery50_adp_pair_subset/injected_lightcurves.h5}"
MIN_CLASS_COUNT="${TWIRL_RECOVERY50_MIN_CLASS_COUNT:-40}"
MAX_SHAPE_ROWS="${TWIRL_RECOVERY50_MAX_SHAPE_ROWS:-}"

echo "[recovery50] host=$(hostname)"
echo "[recovery50] date=$(date -Is)"
echo "[recovery50] repo=${REPO}"
echo "[recovery50] python=${PYTHON}"
echo "[recovery50] queue=${QUEUE}"
echo "[recovery50] labels=${LABELS}"

cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

"${PYTHON}" scripts/stage5_validation/audit_s56_recovery50_teacher_labels.py \
  --queue-csv "${QUEUE}" \
  --labels-csv "${LABELS}" \
  --out-dir "${ROOT}/human_label_audit"

"${PYTHON}" scripts/stage5_validation/build_s56_recovery50_training_table.py \
  --queue-csv "${QUEUE}" \
  --labels-csv "${LABELS}" \
  --out-dir "${ROOT}/human_training_table" \
  --min-multiclass-count "${MIN_CLASS_COUNT}"

SHAPE_ARGS=(
  scripts/stage5_validation/build_s56_recovery50_shape_features.py
  --queue-csv "${QUEUE}"
  --out-dir "${ROOT}/folded_shape_features"
  --compact-lc-h5 "${COMPACT_LC}"
  --hlsp-root "${HLSP_ROOT}"
  --injection-h5-override "${INJECTION_H5}"
)
if [[ -n "${MAX_SHAPE_ROWS}" ]]; then
  SHAPE_ARGS+=(--max-rows "${MAX_SHAPE_ROWS}")
fi
"${PYTHON}" "${SHAPE_ARGS[@]}"

"${PYTHON}" scripts/stage5_validation/train_s56_recovery50_teacher_student.py \
  --training-table "${ROOT}/human_training_table/human_vetting_training_table.csv" \
  --shape-features "${ROOT}/folded_shape_features/folded_shape_features.csv" \
  --metrics-table "${ROOT}/twirl_vet_metrics_real_fullphase_binmatch.csv" \
  --metrics-table "${ROOT}/twirl_vet_metrics_injected_fullphase_binmatch.csv" \
  --out-dir "${ROOT}/teacher_smoke" \
  --min-class-count "${MIN_CLASS_COUNT}"

mkdir -p "${ROOT}/student_smoke"
cp "${ROOT}/teacher_smoke/summary.json" "${ROOT}/student_smoke/teacher_student_summary.json"

echo "[recovery50] complete $(date -Is)"
