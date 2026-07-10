#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-adjud-retrain
#SBATCH -N 1
#SBATCH -c 48
#SBATCH --mem=384G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 16:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56-torch/bin/python}"
ROOT="${TWIRL_ADJUDICATION_ROOT:-reports/stage5_validation/s56_label_adjudication_real343}"
TRAINING="${ROOT}/adjudicated_training_table/human_vetting_training_table_adjudicated.csv"
TENSORS="${ROOT}/retrained/cnn_tensors"
CNN_OUT="${ROOT}/retrained/cnn_teacher"
EB_ROOT="${ROOT}/retrained/eb_miner"
COMPACT_LC="${TWIRL_ADJUDICATION_COMPACT_LC:-data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp_lc_export_pdo.h5}"
HLSP_ROOT="${TWIRL_ADJUDICATION_HLSP_ROOT:-data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

module load cuda/13.1.0 || module load cuda/12.9.1 || true
nvidia-smi

"${PYTHON}" scripts/stage5_validation/audit_s56_real_label_adjudication.py \
  --root "${ROOT}" \
  --require-complete

"${PYTHON}" scripts/stage5_validation/build_s56_recovery50_cnn_tensors.py \
  --training-table "${TRAINING}" \
  --out-dir "${TENSORS}" \
  --compact-lc-h5 "${COMPACT_LC}" \
  --hlsp-root "${HLSP_ROOT}" \
  --use-row-injection-h5 \
  --apertures DET_FLUX_ADP_SML,DET_FLUX_ADP

"${PYTHON}" scripts/stage5_validation/train_s56_recovery50_cnn_teacher.py \
  --tensor-npz "${TENSORS}/recovery50_cnn_tensors.npz" \
  --tensor-rows "${TENSORS}/recovery50_cnn_tensor_rows.csv" \
  --training-table "${TRAINING}" \
  --out-dir "${CNN_OUT}" \
  --epochs 80 \
  --batch-size 32 \
  --min-class-count 40

OLD_ROOT="reports/stage5_validation/s56_recovery50_teacher_queue_2k/cnn_teacher_harmonic_raw_multiview"
for profile in cnn_shape_only cnn_shape_plus_bls; do
  "${PYTHON}" scripts/stage5_validation/score_s56_recovery50_cnn_teacher.py \
    --tensor-npz "${TENSORS}/recovery50_cnn_tensors.npz" \
    --tensor-rows "${TENSORS}/recovery50_cnn_tensor_rows.csv" \
    --model-path "${OLD_ROOT}/${profile}/cnn_teacher.pt" \
    --feature-table "${TRAINING}" \
    --out-dir "${ROOT}/retrained/old_checkpoint_scores/${profile}" \
    --require-cuda
done

"${PYTHON}" scripts/stage5_validation/compare_s56_adjudicated_cnn_checkpoints.py \
  --new-root "${CNN_OUT}" \
  --old-score-root "${ROOT}/retrained/old_checkpoint_scores" \
  --out-dir "${ROOT}/retrained/checkpoint_comparison"

"${PYTHON}" scripts/stage5_validation/build_s56_adjudicated_eb_miner_tables.py \
  --training-table "${TRAINING}" \
  --out-root "${EB_ROOT}"

"${PYTHON}" scripts/stage5_validation/train_score_s56_eb_miner.py \
  --training-table "${EB_ROOT}/training/eb_miner_training_table.csv" \
  --candidate-pool "${EB_ROOT}/candidate_pool/eb_miner_candidate_pool.csv" \
  --out-dir "${EB_ROOT}/eb_cnn_miner" \
  --compact-lc-h5 "${COMPACT_LC}" \
  --hlsp-root "${HLSP_ROOT}" \
  --ensemble-size 5 \
  --epochs 80 \
  --batch-size 32 \
  --candidate-chunk-size 512 \
  --score-chunk-workers 8 \
  --no-resume-scoring

echo "[adjudicated-retrain] complete $(date -Is)"
