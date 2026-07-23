#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=s56-a2rec-bls
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=64G
#SBATCH -t 2-00:00:00
#SBATCH --array=0-39%8
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/s56_a2v1_recovery/%x-%A_%a.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/s56_a2v1_recovery/%x-%A_%a.err

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
CADENCE_REFERENCE="${TWIRL_A2V1_CADENCE_REFERENCE:?set TWIRL_A2V1_CADENCE_REFERENCE}"
CADENCE_REFERENCE_MANIFEST="${TWIRL_A2V1_CADENCE_REFERENCE_MANIFEST:?set TWIRL_A2V1_CADENCE_REFERENCE_MANIFEST}"
MODE="${TWIRL_RECOVERY_MODE:-full}"
SHARD="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
SHARD_TAG="$(printf '%02d' "${SHARD}")"

if [[ "${MODE}" == "smoke" ]]; then
  [[ "${SHARD}" == "0" ]] || { echo "smoke mode requires array task 0" >&2; exit 2; }
  STAGE_ROOT="${OUT_ROOT}/smoke100"
  LIMIT_ARGS=(--limit 100 --selection-mode parameter_spanning --overwrite)
  MIN_INJECTIONS=100
  MIN_CANDIDATES=1000
else
  STAGE_ROOT="${OUT_ROOT}/full"
  LIMIT_ARGS=()
  MIN_INJECTIONS=500
  MIN_CANDIDATES=5000
fi

PAIR_DIR="${STAGE_ROOT}/injection_shards"
PEAK_DIR="${STAGE_ROOT}/peak_shards"
PAIR_H5="${PAIR_DIR}/injections_${SHARD_TAG}.h5"
PEAK_TABLE="${PEAK_DIR}/peaks_${SHARD_TAG}.parquet"
mkdir -p "${LOG_ROOT}" "${PAIR_DIR}" "${PEAK_DIR}"
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE

"${PYTHON}" scripts/stage3_injections/run_s56_a2v1_injection_shard.py \
  --config "${CONFIG}" \
  --raw-h5 "${INPUT_ROOT}/${SECTOR_TAG}_A2v1_tglc_raw_sources.h5" \
  --adp-h5 "${ADP_H5}" \
  --cadence-reference-table "${CADENCE_REFERENCE}" \
  --cadence-reference-manifest "${CADENCE_REFERENCE_MANIFEST}" \
  --schedule "${OUT_ROOT}/schedule/injection_schedule.parquet" \
  --shard-index "${SHARD}" \
  --out-h5 "${PAIR_H5}" \
  "${LIMIT_ARGS[@]}"

"${PYTHON}" scripts/stage5_validation/build_injection_peak_training_table.py \
  --injection-h5 "${PAIR_H5}" \
  --out-table "${PEAK_TABLE}" \
  --summary-json "${PEAK_DIR}/peaks_${SHARD_TAG}_summary.json" \
  --apertures DET_FLUX_ADP_SML DET_FLUX_ADP \
  --n-injections 0 \
  --workers "${SLURM_CPUS_PER_TASK}" \
  --n-periods 200000 \
  --n-peaks 10 \
  --period-tol 0.02 \
  --min-window-overlap-fraction 0.50 \
  --harmonic-factors 0.25,0.3333333333333333,0.5,1,2,3,4 \
  --search-branch "S${SECTOR}_A2v1_ADP_small_anchor_v1" \
  --overwrite

"${PYTHON}" scripts/stage5_validation/verify_injection_peak_training_table.py \
  --peak-table "${PEAK_TABLE}" \
  --out-json "${PEAK_DIR}/peaks_${SHARD_TAG}_verification.json" \
  --min-injections "${MIN_INJECTIONS}" \
  --min-candidate-rows "${MIN_CANDIDATES}" \
  --min-apertures 2 \
  --min-positive-peak-ranks 10 \
  --require-cadence-diagnostics \
  --require-signal-peaks

if [[ "${MODE}" == "smoke" ]]; then
  mkdir -p "${STAGE_ROOT}/teacher_inputs"
  "${PYTHON}" scripts/stage5_validation/prepare_s56_a2v1_recovery_teacher_shard.py \
    --peak-table "${PEAK_TABLE}" \
    --pair-h5 "${PAIR_H5}" \
    --raw-h5 "${INPUT_ROOT}/${SECTOR_TAG}_A2v1_tglc_raw_sources.h5" \
    --adp-h5 "${ADP_H5}" \
    --cadence-reference-table "${CADENCE_REFERENCE}" \
    --cadence-reference-manifest "${CADENCE_REFERENCE_MANIFEST}" \
    --sector "${SECTOR}" \
    --out-candidates "${STAGE_ROOT}/teacher_inputs/candidates_00.csv" \
    --out-native-h5 "${STAGE_ROOT}/teacher_inputs/native_00.h5" \
    --overwrite
fi

echo "[${SECTOR_TAG}-a2rec-bls] mode=${MODE} shard=${SHARD} complete $(date -Is)"
