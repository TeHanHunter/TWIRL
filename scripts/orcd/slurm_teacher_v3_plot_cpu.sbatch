#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-tv3-plot
#SBATCH --exclude=node4900
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:-${TWIRL_ROOT}/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56/bin/python}"
RUN_DIR="${TWIRL_TEACHER_V3_FULL_RUN_DIR:-${TWIRL_ROOT}/results/models/teacher_v3_s56_s62_a2v1_current_adp/full}"
OUT_DIR="${TWIRL_TEACHER_V3_PERFORMANCE_DIR:-${RUN_DIR}/performance}"
SUMMARY="${RUN_DIR}/summary.json"

if [[ ! -s "${SUMMARY}" ]]; then
  echo "Missing completed Teacher-v3 full summary: ${SUMMARY}" >&2
  exit 3
fi

mkdir -p "${TWIRL_ROOT}/logs"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg

plot_args=(
  --run-dir "${RUN_DIR}"
  --expected-bootstrap-draws 2000
  --expected-bootstrap-seed 560063
)

echo "[teacher-v3-plot] validate=$(date -Is) run=${RUN_DIR}"
"${PYTHON}" scripts/stage5_validation/plot_teacher_v3_performance.py \
  "${plot_args[@]}" \
  --check-only

echo "[teacher-v3-plot] render=$(date -Is) output=${OUT_DIR}"
"${PYTHON}" scripts/stage5_validation/plot_teacher_v3_performance.py \
  "${plot_args[@]}" \
  --out-dir "${OUT_DIR}"
echo "[teacher-v3-plot] complete=$(date -Is) output=${OUT_DIR}"
