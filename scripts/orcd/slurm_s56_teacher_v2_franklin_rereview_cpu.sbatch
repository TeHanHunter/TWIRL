#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=s56-t2-franklin-qa
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH -t 08:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_teacher_v2}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
DATA_ROOT="${TWIRL_A2V1_S56_ROOT:-/orcd/data/mki_aryeh/001/twirl/data_local/stage1_lightcurves/s56_A2v1}"
OUT_ROOT="${TWIRL_TEACHER_V2_ROOT:-${REPO}/reports/stage5_validation/s56_A2v1_teacher_v2}"
QA_ROOT="${OUT_ROOT}/franklin_rereview"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs/s56_teacher_v2 \
  /orcd/data/mki_aryeh/001/twirl/mpl-cache "${QA_ROOT}/twirl_vet_sheets"
cd "${REPO}"
export PYTHONPATH=src
export MPLCONFIGDIR=/orcd/data/mki_aryeh/001/twirl/mpl-cache
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

"${PYTHON}" scripts/stage5_validation/render_two_aperture_vet_sheets.py \
  --queue-csv "${OUT_ROOT}/roles/franklin_a2v1_rereview_queue.csv" \
  --lc-export-h5 "${DATA_ROOT}/s56_A2v1_adp_pair_rebuilt.h5" \
  --out-dir "${QA_ROOT}/twirl_vet_sheets" \
  --metrics-csv "${QA_ROOT}/twirl_vet_metrics.csv" \
  --summary-json "${QA_ROOT}/render_summary.json" \
  --branch-name current_a2v1_adp \
  --apertures DET_FLUX_ADP_SML DET_FLUX_ADP \
  --workers 16 \
  --n-periods 200000 \
  --n-peaks 10 \
  --use-row-ephemeris \
  --no-pdf \
  --harmonic-factors 0.25 0.5 1 2 4 \
  --overwrite

"${PYTHON}" scripts/stage5_validation/verify_franklin_a2v1_rereview.py \
  --queue "${OUT_ROOT}/roles/franklin_a2v1_rereview_queue.csv" \
  --private-manifest "${OUT_ROOT}/roles/franklin_a2v1_rereview_private_manifest.parquet" \
  --render-summary "${QA_ROOT}/render_summary.json" \
  --sheet-dir "${QA_ROOT}/twirl_vet_sheets" \
  --out-json "${QA_ROOT}/verification.json"
echo "[s56-t2-franklin-qa] complete=$(date -Is)"
