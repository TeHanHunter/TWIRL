#!/bin/bash
# Build recovery50 training products and CNN tensors on ORCD, then submit the H200 teacher job.
set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
ROOT="${TWIRL_RECOVERY50_ROOT:-reports/stage5_validation/s56_recovery50_teacher_queue}"
QUEUE="${TWIRL_RECOVERY50_QUEUE:-${ROOT}/review_queue_1k.csv}"
LABELS="${TWIRL_RECOVERY50_LABELS:-${ROOT}/human_labels_vetted.csv}"
COMPACT_LC="${TWIRL_RECOVERY50_COMPACT_LC:-data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp_lc_export_pdo.h5}"
HLSP_ROOT="${TWIRL_RECOVERY50_HLSP_ROOT:-data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare}"
INJECTION_H5="${TWIRL_RECOVERY50_INJECTION_H5:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/recovery50_adp_pair_subset/injected_lightcurves.h5}"
MIN_CLASS_COUNT="${TWIRL_RECOVERY50_MIN_CLASS_COUNT:-40}"
SUBMIT_H200="${TWIRL_RECOVERY50_SUBMIT_H200:-1}"
APERTURES="${TWIRL_RECOVERY50_APERTURES:-DET_FLUX_ADP_SML,DET_FLUX_ADP}"

cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

echo "[recovery50-cnn] host=$(hostname)"
echo "[recovery50-cnn] date=$(date -Is)"
echo "[recovery50-cnn] repo=${REPO}"
echo "[recovery50-cnn] python=${PYTHON}"
echo "[recovery50-cnn] queue=${QUEUE}"
echo "[recovery50-cnn] labels=${LABELS}"
echo "[recovery50-cnn] apertures=${APERTURES}"

"${PYTHON}" scripts/stage5_validation/audit_s56_recovery50_teacher_labels.py \
  --queue-csv "${QUEUE}" \
  --labels-csv "${LABELS}" \
  --out-dir "${ROOT}/human_label_audit"

"${PYTHON}" scripts/stage5_validation/build_s56_recovery50_training_table.py \
  --queue-csv "${QUEUE}" \
  --labels-csv "${LABELS}" \
  --out-dir "${ROOT}/human_training_table" \
  --min-multiclass-count "${MIN_CLASS_COUNT}"

"${PYTHON}" scripts/stage5_validation/build_s56_recovery50_cnn_tensors.py \
  --training-table "${ROOT}/human_training_table/human_vetting_training_table.csv" \
  --out-dir "${ROOT}/cnn_tensors" \
  --compact-lc-h5 "${COMPACT_LC}" \
  --hlsp-root "${HLSP_ROOT}" \
  --injection-h5-override "${INJECTION_H5}" \
  --apertures "${APERTURES}"

if [[ "${SUBMIT_H200}" == "1" ]]; then
  mkdir -p /orcd/data/mki_aryeh/001/twirl/logs
  job_id="$(
    sbatch --parsable \
      --export=ALL,TWIRL_RECOVERY50_ROOT="${ROOT}",TWIRL_RECOVERY50_MIN_CLASS_COUNT="${MIN_CLASS_COUNT}" \
      scripts/orcd/slurm_s56_recovery50_cnn_teacher_h200.sbatch
  )"
  echo "[recovery50-cnn] submitted H200 teacher job ${job_id}"
else
  echo "[recovery50-cnn] tensor build complete; H200 submit skipped"
fi

echo "[recovery50-cnn] complete $(date -Is)"
