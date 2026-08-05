#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-teacher-v1-native-v2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 48
#SBATCH --mem=384G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 24:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56-torch/bin/python}"
ROOT="${TWIRL_ADJUDICATION_ROOT:-reports/stage5_validation/s56_label_adjudication_real343}"
TRAINING="${ROOT}/adjudicated_training_table/human_vetting_training_table_adjudicated.csv"
NATIVE_H5="${TWIRL_HARMONIC_NATIVE_H5:-${ROOT}/native_inputs/s56_adp_raw_pair_v2.h5}"
OUT="${ROOT}/retrained/s56_harmonic_cnn_v1_native_v2"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs "${OUT}"
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

"${PYTHON}" scripts/stage5_validation/train_s56_harmonic_cnn_v1.py \
  --training-table "${TRAINING}" \
  --native-h5 "${NATIVE_H5}" \
  --out-dir "${OUT}" \
  --profiles metadata_only,single_period_native_fold,seven_harmonic_shape,shape_plus_raw_chronology,shape_plus_periodogram_bls,full_combined \
  --epochs 100 \
  --pretrain-epochs 20 \
  --batch-size 32 \
  --workers 16 \
  --learning-rate 3e-4 \
  --weight-decay 1e-4 \
  --patience 12

echo "[harmonic-cnn-v1-native-v2] complete $(date -Is): ${OUT}"
