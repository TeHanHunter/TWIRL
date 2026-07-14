#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=s56-a2rec-final
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -t 08:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/s56_a2v1_recovery/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/s56_a2v1_recovery/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
SECTOR="${TWIRL_RECOVERY_SECTOR:-56}"
SECTOR_TAG="s${SECTOR}"
OUT_ROOT="${TWIRL_RECOVERY_ROOT:-${REPO}/reports/stage5_validation/${SECTOR_TAG}_A2v1_teacher_v1_recovery}"
FULL="${OUT_ROOT}/full"
INPUT_ROOT="${TWIRL_A2V1_INPUT_ROOT:-/orcd/data/mki_aryeh/001/twirl/data_local/stage1_lightcurves/${SECTOR_TAG}_A2v1}"
CONFIG="${TWIRL_RECOVERY_CONFIG:-${REPO}/configs/injections/${SECTOR_TAG}_a2v1_teacher_v1_recovery_v1.yaml}"
CHECKPOINT_ROOT="${TWIRL_TEACHER_CHECKPOINT_ROOT:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL/reports/stage5_validation/s56_label_adjudication_real343/retrained/s56_harmonic_cnn_v1/shape_plus_periodogram_bls}"
PRODUCER_SHA_FILE="${OUT_ROOT}/schedule/producer_git_sha.txt"

mkdir -p "${FULL}/final"
cd "${REPO}"
export PYTHONPATH=src
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
test -s "${PRODUCER_SHA_FILE}"

"${PYTHON}" scripts/stage5_validation/audit_a2v1_recovery_products.py \
  --config "${CONFIG}" \
  --schedule "${OUT_ROOT}/schedule/injection_schedule.parquet" \
  --shards "${FULL}"/injection_shards/injections_*.h5 \
  --raw-h5 "${INPUT_ROOT}/${SECTOR_TAG}_A2v1_tglc_raw_sources.h5" \
  --adp-h5 "${INPUT_ROOT}/${SECTOR_TAG}_A2v1_adp_pair_rebuilt.h5" \
  --out-json "${FULL}/final/injection_product_audit.json"

"${PYTHON}" scripts/stage5_validation/finalize_s56_a2v1_teacher_recovery.py \
  --manifests "${FULL}"/injection_shards/injections_*_manifest.parquet \
  --peak-tables "${FULL}"/peak_shards/peaks_*.parquet \
  --teacher-scores "${FULL}/teacher_scores.parquet" \
  --out-dir "${FULL}/final" \
  --expect-injections 20000 \
  --preserve-threshold 0.5

"${PYTHON}" scripts/stage5_validation/write_s56_a2v1_recovery_run_manifest.py \
  --repo "${REPO}" \
  --config "${CONFIG}" \
  --input-verification "${INPUT_ROOT}/input_verification.json" \
  --schedule-summary "${OUT_ROOT}/schedule/schedule_summary.json" \
  --teacher-summary "${FULL}/teacher_scores.summary.json" \
  --final-summary "${FULL}/final/summary.json" \
  --product-audit "${FULL}/final/injection_product_audit.json" \
  --producer-git-sha "$(tr -d '[:space:]' < "${PRODUCER_SHA_FILE}")" \
  --checkpoints \
    "${CHECKPOINT_ROOT}/fold_0/teacher.pt" \
    "${CHECKPOINT_ROOT}/fold_1/teacher.pt" \
    "${CHECKPOINT_ROOT}/fold_2/teacher.pt" \
    "${CHECKPOINT_ROOT}/fold_3/teacher.pt" \
    "${CHECKPOINT_ROOT}/fold_4/teacher.pt" \
  --out-json "${FULL}/final/run_manifest.json"

"${PYTHON}" scripts/stage5_validation/verify_a2v1_teacher_recovery_final.py \
  --config "${CONFIG}" \
  --recovery-root "${OUT_ROOT}" \
  --input-verification "${INPUT_ROOT}/input_verification.json" \
  --teacher-table "${INPUT_ROOT}/human_vetting_training_table_adjudicated.csv" \
  --out-json "${FULL}/final/acceptance.json"

echo "[${SECTOR_TAG}-a2rec-final] complete $(date -Is)"
