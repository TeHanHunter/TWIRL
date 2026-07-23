#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=s56-t1-canary-v2
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH --array=0,13,26,39%4
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/s56_A2v1_tier1_canary_v2/%x-%A_%a.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/s56_A2v1_tier1_canary_v2/%x-%A_%a.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:?set TWIRL_ORCD_REPO to an exact versioned checkout}"
PRODUCER_GIT_SHA="${TWIRL_PRODUCER_GIT_SHA:?set TWIRL_PRODUCER_GIT_SHA}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
INPUT_ROOT="${TWIRL_A2V1_INPUT_ROOT:-/orcd/data/mki_aryeh/001/twirl/data_local/stage1_lightcurves/s56_A2v1}"
EVIDENCE_ROOT="${TWIRL_TIER1_EVIDENCE_ROOT:-${INPUT_ROOT}/tier1_active_search_pair_v1/evidence}"
RECOVERY_ROOT="${TWIRL_RECOVERY_ROOT:-/orcd/data/mki_aryeh/001/twirl/retired_worktrees/20260717_cleanup/TWIRL_s56_a2v1_recovery_2b4f57bb/reports/stage5_validation/s56_A2v1_teacher_v1_recovery}"
OUT_ROOT="${TWIRL_TIER1_CANARY_ROOT:-/orcd/data/mki_aryeh/001/twirl/data_local/stage3_injections/s56_A2v1_tier1_canary_v2_20260723}"
CONFIG="${TWIRL_RECOVERY_CONFIG:-${REPO}/configs/injections/s56_a2v1_teacher_v1_recovery_v1.yaml}"
SHARD="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
SHARD_TAG="$(printf '%02d' "${SHARD}")"
CADENCE_TABLE="${EVIDENCE_ROOT}/s56_A2v1_cadence_reference.csv"
CADENCE_MANIFEST="${EVIDENCE_ROOT}/s56_A2v1_cadence_reference_manifest.json"
SCHEDULE="${RECOVERY_ROOT}/schedule/injection_schedule.parquet"

verify_sha256() {
  local expected="$1"
  local path="$2"
  local observed
  observed="$(sha256sum "${path}" | awk '{print $1}')"
  if [[ "${observed}" != "${expected}" ]]; then
    echo "SHA-256 mismatch for ${path}: ${observed} != ${expected}" >&2
    exit 2
  fi
}

mkdir -p "${OUT_ROOT}/injection_shards"
cd "${REPO}"
OBSERVED_GIT_SHA="$(git rev-parse HEAD)"
if [[ "${OBSERVED_GIT_SHA}" != "${PRODUCER_GIT_SHA}" ]]; then
  echo "producer revision mismatch: ${OBSERVED_GIT_SHA} != ${PRODUCER_GIT_SHA}" >&2
  exit 2
fi
verify_sha256 \
  "6da30db4feff76b52575f83b6acedd39173e29b8f7328c364f2927329a167f4a" \
  "${CADENCE_TABLE}"
verify_sha256 \
  "a55bd499b46305eaa55ab6be03979d7c15253cd1a1b6d1e928782a2d23db5b62" \
  "${CADENCE_MANIFEST}"
verify_sha256 \
  "e22cc0dc7d6ac4e59d65ab0bc55ae11b1c39e2b85302d7f6813b2b03952262dd" \
  "${SCHEDULE}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

"${PYTHON}" scripts/stage3_injections/run_s56_a2v1_injection_shard.py \
  --config "${CONFIG}" \
  --raw-h5 "${INPUT_ROOT}/s56_A2v1_tglc_raw_sources.h5" \
  --adp-h5 "${INPUT_ROOT}/s56_A2v1_adp_pair_rebuilt.h5" \
  --cadence-reference-table "${CADENCE_TABLE}" \
  --cadence-reference-manifest "${CADENCE_MANIFEST}" \
  --schedule "${SCHEDULE}" \
  --shard-index "${SHARD}" \
  --out-h5 "${OUT_ROOT}/injection_shards/injections_${SHARD_TAG}.h5"

echo "[s56-tier1-canary-v2] shard=${SHARD} complete $(date -Is)"
