#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=s56-a2rec-compact
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/s56_a2v1_recovery/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/s56_a2v1_recovery/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
SECTOR="${TWIRL_RECOVERY_SECTOR:-56}"
SECTOR_TAG="s${SECTOR}"
SECTOR_PADDED="$(printf 's%04d' "${SECTOR}")"
LOG_ROOT="${TWIRL_RECOVERY_LOG_ROOT:-/orcd/data/mki_aryeh/001/twirl/logs/${SECTOR_TAG}_a2v1_recovery}"
INPUT_ROOT="${TWIRL_A2V1_INPUT_ROOT:-/orcd/data/mki_aryeh/001/twirl/data_local/stage1_lightcurves/${SECTOR_TAG}_A2v1}"
FITS_ROOT="${INPUT_ROOT}/hlsp_${SECTOR_PADDED}_A2v1"
ACTIVE_ADP="${INPUT_ROOT}/${SECTOR_TAG}_A2v1_adp_pair_rebuilt.h5"

mkdir -p "${LOG_ROOT}"
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE

"${PYTHON}" scripts/stage3_injections/export_s56_lc_training_set.py \
  --hlsp-root "${FITS_ROOT}" \
  --out-h5 "${ACTIVE_ADP}" \
  --sector "${SECTOR}" \
  --columns DET_FLUX_ADP_SML DET_FLUX_ADP \
  --compression lzf \
  --progress-every 500 \
  --overwrite

"${PYTHON}" scripts/orcd/verify_s56_a2v1_recovery_inputs.py \
  --input-root "${INPUT_ROOT}" \
  --sector "${SECTOR}" \
  --out-json "${INPUT_ROOT}/input_verification.json" \
  --hash-large-files

echo "[${SECTOR_TAG}-a2rec-compact] complete $(date -Is)"
