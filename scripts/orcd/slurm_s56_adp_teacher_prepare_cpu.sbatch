#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-adp-teacher-prep
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=160G
#SBATCH -t 06:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
ROOT="${TWIRL_RECOVERY50_ROOT:-reports/stage5_validation/s56_recovery50_teacher_queue_2k}"
TRAINING_TABLE="${TWIRL_RECOVERY50_TRAINING_TABLE:-${ROOT}/human_training_table_adp_only/human_vetting_training_table.csv}"
TENSOR_DIR="${TWIRL_RECOVERY50_TENSOR_DIR:-${ROOT}/cnn_tensors_adp_only}"
COMPACT_LC="${TWIRL_RECOVERY50_COMPACT_LC:-data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp_lc_export_pdo.h5}"
HLSP_ROOT="${TWIRL_RECOVERY50_HLSP_ROOT:-data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[adp-teacher-prep] host=$(hostname)"
echo "[adp-teacher-prep] date=$(date -Is)"
echo "[adp-teacher-prep] training_table=${TRAINING_TABLE}"
echo "[adp-teacher-prep] tensor_dir=${TENSOR_DIR}"

"${PYTHON}" scripts/stage5_validation/apply_s56_recovery50_revisit_labels.py

"${PYTHON}" scripts/stage5_validation/audit_build_s56_2k_adp_only.py

"${PYTHON}" scripts/stage5_validation/build_s56_recovery50_cnn_tensors.py \
  --training-table "${TRAINING_TABLE}" \
  --out-dir "${TENSOR_DIR}" \
  --compact-lc-h5 "${COMPACT_LC}" \
  --hlsp-root "${HLSP_ROOT}" \
  --use-row-injection-h5

echo "[adp-teacher-prep] complete $(date -Is)"
