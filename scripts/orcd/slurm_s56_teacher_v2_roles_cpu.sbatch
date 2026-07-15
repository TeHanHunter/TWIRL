#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=s56-t2-roles
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -t 08:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_teacher_v2}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
OUT_ROOT="${TWIRL_TEACHER_V2_ROOT:-${REPO}/reports/stage5_validation/s56_A2v1_teacher_v2}"
S56_RECOVERY="${TWIRL_S56_RECOVERY_ROOT:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_s56_a2v1_recovery_2b4f57bb/reports/stage5_validation/s56_A2v1_teacher_v1_recovery}"
S57_RECOVERY="${TWIRL_S57_RECOVERY_ROOT:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_a2v1_recovery_1090f011/reports/stage5_validation/s57_A2v1_teacher_v1_recovery}"
N_SHARDS="${TWIRL_TEACHER_V2_NATIVE_SHARDS:-32}"
FINAL_NATIVE="${OUT_ROOT}/real_native/s56_a2v1_real_native.h5"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2 "${OUT_ROOT}/roles"
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

shards=()
for ((index=0; index<N_SHARDS; index++)); do
  shards+=("${OUT_ROOT}/real_native/shards/native_${index}.h5")
done
"${PYTHON}" scripts/stage5_validation/merge_teacher_v2_native_shards.py \
  --shards "${shards[@]}" \
  --out-h5 "${FINAL_NATIVE}" \
  --candidate-table "${OUT_ROOT}/human_transfer/s56_a2v1_real_candidates_top5.parquet"

"${PYTHON}" scripts/stage5_validation/prepare_s56_teacher_v2_roles.py \
  --s56-schedule "${S56_RECOVERY}/schedule/injection_schedule.parquet" \
  --s56-candidates "${S56_RECOVERY}/full/teacher_inputs/candidates_top5.csv" \
  --s56-injection-native-h5 "${S56_RECOVERY}/full/teacher_inputs/native_top5.h5" \
  --human-table "${OUT_ROOT}/human_transfer/human_a2v1_transferred_training_rows.csv" \
  --human-native-h5 "${FINAL_NATIVE}" \
  --s56-real-candidates "${OUT_ROOT}/human_transfer/s56_a2v1_real_candidates_top5.parquet" \
  --human-prior-splits "${OUT_ROOT}/source_snapshots/training_rows_with_fixed_splits.csv" \
  --franklin-labels "${OUT_ROOT}/source_snapshots/franklin_labels_vetted.csv" \
  --s57-schedule "${S57_RECOVERY}/schedule/injection_schedule.parquet" \
  --s57-candidates "${S57_RECOVERY}/full/teacher_inputs/candidates_top5.csv" \
  --s57-native-h5 "${S57_RECOVERY}/full/teacher_inputs/native_top5.h5" \
  --out-dir "${OUT_ROOT}/roles" \
  --allow-incomplete-s57

"${PYTHON}" scripts/stage5_validation/seed_teacher_v2_native_verification_cache.py \
  --native-summary \
    "${S56_RECOVERY}/full/teacher_inputs/native_top5.h5" \
    "${S56_RECOVERY}/full/teacher_inputs/native_top5.merge_summary.json" \
  --native-summary \
    "${FINAL_NATIVE}" \
    "${FINAL_NATIVE%.h5}.summary.json" \
  --out-cache "${OUT_ROOT}/roles/native_contract_verification_cache.json"

"${PYTHON}" -c \
  'import json,sys; p=json.load(open(sys.argv[1])); assert p["s57_external"]["passed"], p["s57_external"]' \
  "${OUT_ROOT}/roles/role_summary.json"
echo "[s56-t2-roles] complete=$(date -Is)"
