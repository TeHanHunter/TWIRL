#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-eb-h200
#SBATCH -N 1
#SBATCH -c 48
#SBATCH --mem=320G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 12:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56-torch/bin/python}"
ROOT="${TWIRL_EB_MINER_ROOT:-reports/stage5_validation/s56_eb_miner_adp_only}"
TRAINING_TABLE="${TWIRL_EB_TRAINING_TABLE:-${ROOT}/training/eb_miner_training_table.csv}"
CANDIDATE_POOL="${TWIRL_EB_CANDIDATE_POOL:-${ROOT}/candidate_pool/eb_miner_candidate_pool.csv}"
OUT_DIR="${TWIRL_EB_CNN_OUT_DIR:-${ROOT}/eb_cnn_miner}"
COMPACT_LC="${TWIRL_EB_COMPACT_LC:-data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp_lc_export_pdo.h5}"
HLSP_ROOT="${TWIRL_EB_HLSP_ROOT:-data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare}"
CHUNK_SIZE="${TWIRL_EB_CHUNK_SIZE:-512}"
SCORE_CHUNK_WORKERS="${TWIRL_EB_SCORE_CHUNK_WORKERS:-8}"
ENSEMBLE_SIZE="${TWIRL_EB_ENSEMBLE_SIZE:-5}"
EPOCHS="${TWIRL_EB_EPOCHS:-80}"
BATCH_SIZE="${TWIRL_EB_BATCH_SIZE:-32}"
REUSE_MODELS="${TWIRL_EB_REUSE_MODELS:-0}"
RESUME_SCORING="${TWIRL_EB_RESUME_SCORING:-1}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs
cd "${REPO}"

export PYTHONPATH=src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[eb-h200] host=$(hostname)"
echo "[eb-h200] date=$(date -Is)"
echo "[eb-h200] repo=${REPO}"
echo "[eb-h200] python=${PYTHON}"
echo "[eb-h200] root=${ROOT}"
echo "[eb-h200] training_table=${TRAINING_TABLE}"
echo "[eb-h200] candidate_pool=${CANDIDATE_POOL}"
echo "[eb-h200] out_dir=${OUT_DIR}"
echo "[eb-h200] ensemble=${ENSEMBLE_SIZE} epochs=${EPOCHS} chunk_size=${CHUNK_SIZE} score_chunk_workers=${SCORE_CHUNK_WORKERS}"

module load cuda/13.1.0 || module load cuda/12.9.1 || true
nvidia-smi

extra_args=()
if [[ "${REUSE_MODELS}" == "1" ]]; then
  extra_args+=(--reuse-existing-models)
fi
if [[ "${RESUME_SCORING}" != "1" ]]; then
  extra_args+=(--no-resume-scoring)
fi

"${PYTHON}" scripts/stage5_validation/train_score_s56_eb_miner.py \
  --training-table "${TRAINING_TABLE}" \
  --candidate-pool "${CANDIDATE_POOL}" \
  --out-dir "${OUT_DIR}" \
  --compact-lc-h5 "${COMPACT_LC}" \
  --hlsp-root "${HLSP_ROOT}" \
  --ensemble-size "${ENSEMBLE_SIZE}" \
  --epochs "${EPOCHS}" \
  --batch-size "${BATCH_SIZE}" \
  --candidate-chunk-size "${CHUNK_SIZE}" \
  --score-chunk-workers "${SCORE_CHUNK_WORKERS}" \
  "${extra_args[@]}"

echo "[eb-h200] complete $(date -Is)"
