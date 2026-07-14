#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=s56-a2rec-prep
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/s56_a2v1_recovery/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/s56_a2v1_recovery/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
SECTOR="${TWIRL_RECOVERY_SECTOR:-56}"
SECTOR_TAG="s${SECTOR}"
LOG_ROOT="${TWIRL_RECOVERY_LOG_ROOT:-/orcd/data/mki_aryeh/001/twirl/logs/${SECTOR_TAG}_a2v1_recovery}"
INPUT_ROOT="${TWIRL_A2V1_INPUT_ROOT:-/orcd/data/mki_aryeh/001/twirl/data_local/stage1_lightcurves/${SECTOR_TAG}_A2v1}"
OUT_ROOT="${TWIRL_RECOVERY_ROOT:-${REPO}/reports/stage5_validation/${SECTOR_TAG}_A2v1_teacher_v1_recovery}"
CONFIG="${TWIRL_RECOVERY_CONFIG:-${REPO}/configs/injections/${SECTOR_TAG}_a2v1_teacher_v1_recovery_v1.yaml}"
ADP_H5="${TWIRL_A2V1_ADP_H5:-${INPUT_ROOT}/${SECTOR_TAG}_A2v1_adp_pair_rebuilt.h5}"
PRIOR_SCHEDULE="${TWIRL_RECOVERY_PRIOR_SCHEDULE:-}"
OVERLAP_AUDIT_TABLE="${TWIRL_RECOVERY_HOST_OVERLAP_AUDIT_TABLE:-}"
EXCLUSION_ARGS=()
OVERLAP_AUDIT_ARGS=()
if [[ -n "${PRIOR_SCHEDULE}" ]]; then
  test -s "${PRIOR_SCHEDULE}"
  EXCLUSION_ARGS=(--additional-exclusion-table "${PRIOR_SCHEDULE}")
fi
if [[ -n "${OVERLAP_AUDIT_TABLE}" ]]; then
  test -s "${OVERLAP_AUDIT_TABLE}"
  OVERLAP_AUDIT_ARGS=(--host-overlap-audit-table "${OVERLAP_AUDIT_TABLE}")
fi

mkdir -p "${LOG_ROOT}" "${OUT_ROOT}/schedule"
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE

test -s "${INPUT_ROOT}/input_verification.json"
test -s "${ADP_H5}"
grep -q '"passed": true' "${INPUT_ROOT}/input_verification.json"

"${PYTHON}" scripts/stage3_injections/build_s56_a2v1_recovery_schedule.py \
  --config "${CONFIG}" \
  --raw-h5 "${INPUT_ROOT}/${SECTOR_TAG}_A2v1_tglc_raw_sources.h5" \
  --adp-h5 "${ADP_H5}" \
  --teacher-table "${INPUT_ROOT}/human_vetting_training_table_adjudicated.csv" \
  --out-dir "${OUT_ROOT}/schedule" \
  "${EXCLUSION_ARGS[@]}" \
  "${OVERLAP_AUDIT_ARGS[@]}"

echo "[${SECTOR_TAG}-a2rec-prep] complete $(date -Is)"
