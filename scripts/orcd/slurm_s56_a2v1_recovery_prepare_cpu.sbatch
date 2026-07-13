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
INPUT_ROOT="${TWIRL_A2V1_INPUT_ROOT:-/orcd/data/mki_aryeh/001/twirl/data_local/stage1_lightcurves/s56_A2v1}"
OUT_ROOT="${TWIRL_RECOVERY_ROOT:-${REPO}/reports/stage5_validation/s56_A2v1_teacher_v1_recovery}"
CONFIG="${TWIRL_RECOVERY_CONFIG:-${REPO}/configs/injections/s56_a2v1_teacher_v1_recovery_v1.yaml}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs/s56_a2v1_recovery "${OUT_ROOT}/schedule"
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE

test -s "${INPUT_ROOT}/input_verification.json"
grep -q '"passed": true' "${INPUT_ROOT}/input_verification.json"

"${PYTHON}" scripts/stage3_injections/build_s56_a2v1_recovery_schedule.py \
  --config "${CONFIG}" \
  --raw-h5 "${INPUT_ROOT}/s56_A2v1_tglc_raw_sources.h5" \
  --adp-h5 "${INPUT_ROOT}/s56_A2v1_adp_pair.h5" \
  --teacher-table "${INPUT_ROOT}/human_vetting_training_table_adjudicated.csv" \
  --out-dir "${OUT_ROOT}/schedule"

echo "[s56-a2rec-prep] complete $(date -Is)"
