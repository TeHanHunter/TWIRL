#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=s56-t2-finalize
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
FROZEN="${OUT_ROOT}/frozen/frozen_selection.json"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2 \
  "${OUT_ROOT}/s56_locked_holdout_recovery" "${OUT_ROOT}/s57_external_recovery" \
  "${OUT_ROOT}/human_holdout"
cd "${REPO}"
export PYTHONPATH=src
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

V1_THRESHOLD=$("${PYTHON}" -c \
  'import json,sys; print(json.load(open(sys.argv[1]))["workload_threshold"]["threshold"])' \
  "${OUT_ROOT}/frozen/teacher_v1_workload_threshold.json")
METADATA_THRESHOLD=$("${PYTHON}" -c \
  'import json,sys; print(json.load(open(sys.argv[1]))["profile_details"]["metadata_only"]["workload_threshold"]["threshold"])' \
  "${FROZEN}")

"${PYTHON}" scripts/stage5_validation/finalize_teacher_v2_external_recovery.py \
  --manifest "${OUT_ROOT}/roles/s56_injection_locked_holdout_manifest.parquet" \
  --candidates "${S56_RECOVERY}/full/teacher_inputs/candidates_top5.csv" \
  --teacher-v2-scores "${OUT_ROOT}/frozen_scores/s56_holdout_teacher_v2.parquet" \
  --teacher-v1-scores "${OUT_ROOT}/frozen_scores/s56_holdout_teacher_v1.parquet" \
  --teacher-v1-threshold "${V1_THRESHOLD}" \
  --metadata-scores "${OUT_ROOT}/frozen_scores/s56_holdout_metadata_only.parquet" \
  --metadata-threshold "${METADATA_THRESHOLD}" \
  --frozen-selection "${FROZEN}" \
  --expect-min-injections 5000 \
  --out-dir "${OUT_ROOT}/s56_locked_holdout_recovery"

"${PYTHON}" scripts/stage5_validation/finalize_teacher_v2_external_recovery.py \
  --manifest "${OUT_ROOT}/roles/s57_external_role_manifest.parquet" \
  --candidates "${S57_RECOVERY}/full/teacher_inputs/candidates_top5.csv" \
  --teacher-v2-scores "${OUT_ROOT}/frozen_scores/s57_external_teacher_v2.parquet" \
  --teacher-v1-scores "${OUT_ROOT}/frozen_scores/s57_external_teacher_v1.parquet" \
  --teacher-v1-threshold "${V1_THRESHOLD}" \
  --metadata-scores "${OUT_ROOT}/frozen_scores/s57_external_metadata_only.parquet" \
  --metadata-threshold "${METADATA_THRESHOLD}" \
  --frozen-selection "${FROZEN}" \
  --expect-min-injections 10000 \
  --out-dir "${OUT_ROOT}/s57_external_recovery"

human_status=0
"${PYTHON}" scripts/stage5_validation/evaluate_s56_teacher_v2_human_holdout.py \
  --human-role-table "${OUT_ROOT}/roles/s56_human_role_manifest.parquet" \
  --teacher-v2-scores "${OUT_ROOT}/frozen_scores/s56_human_teacher_v2.parquet" \
  --teacher-v1-scores "${OUT_ROOT}/frozen_scores/s56_human_teacher_v1.parquet" \
  --out-dir "${OUT_ROOT}/human_holdout" || human_status=$?

"${PYTHON}" scripts/stage5_validation/summarize_franklin_teacher_v2_audit.py \
  --scores "${OUT_ROOT}/frozen_scores/franklin_current_a2v1_teacher_v2.parquet" \
  --frozen-selection "${FROZEN}" \
  --out-dir "${OUT_ROOT}/franklin_score_audit"

"${PYTHON}" scripts/stage5_validation/verify_teacher_v2_recovery_final.py \
  --root "${OUT_ROOT}/s56_locked_holdout_recovery" \
  --role-summary "${OUT_ROOT}/roles/role_summary.json" \
  --expected-sector 56 \
  --minimum-injections 5000
"${PYTHON}" scripts/stage5_validation/verify_teacher_v2_recovery_final.py \
  --root "${OUT_ROOT}/s57_external_recovery" \
  --role-summary "${OUT_ROOT}/roles/role_summary.json" \
  --expected-sector 57 \
  --minimum-injections 10000
echo "[s56-t2-finalize] complete=$(date -Is) human_status=${human_status}"
exit "${human_status}"
