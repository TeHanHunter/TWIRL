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
OUT_ROOT="${TWIRL_RECOVERY_ROOT:-${REPO}/reports/stage5_validation/s56_A2v1_teacher_v1_recovery}"
FULL="${OUT_ROOT}/full"

mkdir -p "${FULL}/final"
cd "${REPO}"
export PYTHONPATH=src
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

"${PYTHON}" scripts/stage5_validation/finalize_s56_a2v1_teacher_recovery.py \
  --manifests "${FULL}"/injection_shards/injections_*_manifest.parquet \
  --peak-tables "${FULL}"/peak_shards/peaks_*.parquet \
  --teacher-scores "${FULL}/teacher_scores.parquet" \
  --out-dir "${FULL}/final" \
  --expect-injections 20000 \
  --preserve-threshold 0.5

echo "[s56-a2rec-final] complete $(date -Is)"
