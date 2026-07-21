#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-franklin-render
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 48
#SBATCH --mem=192G
#SBATCH -t 12:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/teacher_enrichment/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/teacher_enrichment/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_franklin_s57_s59_20260717}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
QUEUE="${TWIRL_FRANKLIN_QUEUE:?set TWIRL_FRANKLIN_QUEUE}"
LC_EXPORT="${TWIRL_FRANKLIN_LC_EXPORT:?set TWIRL_FRANKLIN_LC_EXPORT}"
OUT_ROOT="${TWIRL_FRANKLIN_OUT_ROOT:?set TWIRL_FRANKLIN_OUT_ROOT}"
WORKERS="${TWIRL_FRANKLIN_WORKERS:-${SLURM_CPUS_PER_TASK:-48}}"
N_PERIODS="${TWIRL_FRANKLIN_N_PERIODS:-50000}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs/teacher_enrichment \
  /orcd/data/mki_aryeh/001/twirl/mpl-cache "${OUT_ROOT}/vet_sheets"
cd "${REPO}"
export PYTHONPATH=src
export MPLCONFIGDIR=/orcd/data/mki_aryeh/001/twirl/mpl-cache
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

"${PYTHON}" scripts/stage5_validation/render_two_aperture_vet_sheets.py \
  --queue-csv "${QUEUE}" \
  --lc-export-h5 "${LC_EXPORT}" \
  --out-dir "${OUT_ROOT}/vet_sheets" \
  --metrics-csv "${OUT_ROOT}/twirl_vet_metrics.csv" \
  --summary-json "${OUT_ROOT}/render_summary.json" \
  --branch-name current_adp \
  --apertures DET_FLUX_ADP_SML DET_FLUX_ADP \
  --workers "${WORKERS}" \
  --n-periods "${N_PERIODS}" \
  --n-peaks 10 \
  --use-row-ephemeris \
  --no-pdf \
  --harmonic-factors 0.25 0.5 1 2 4

echo "[franklin-sector-render] queue=${QUEUE} complete=$(date -Is)"
