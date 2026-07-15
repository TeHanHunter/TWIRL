#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-mixed-render
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 48
#SBATCH --mem=192G
#SBATCH -t 18:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/teacher_enrichment/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/teacher_enrichment/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL_teacher_v2}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
S56_SCORE_ROOT="${TWIRL_S56_SCORE_ROOT:-${REPO}/reports/stage5_validation/s56_A2v1_teacher_v2/real_scores}"
S57_SCORE_ROOT="${TWIRL_S57_SCORE_ROOT:?set TWIRL_S57_SCORE_ROOT}"
S56_DATA_ROOT="${TWIRL_S56_DATA_ROOT:-/orcd/data/mki_aryeh/001/twirl/data_local/stage1_lightcurves/s56_A2v1}"
S57_DATA_ROOT="${TWIRL_S57_DATA_ROOT:-/orcd/data/mki_aryeh/001/twirl/data_local/stage1_lightcurves/s57_A2v1}"
OUT_ROOT="${TWIRL_MIXED_ENRICH_OUT_ROOT:?set TWIRL_MIXED_ENRICH_OUT_ROOT}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/mpl-cache "${OUT_ROOT}/twirl_vet_sheets"
cd "${REPO}"
export PYTHONPATH=src
export MPLCONFIGDIR=/orcd/data/mki_aryeh/001/twirl/mpl-cache
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

build_args=(
  --sector-score-set 56 "${S56_SCORE_ROOT}/shape_plus_raw_chronology.parquet" "${S56_SCORE_ROOT}/shape_plus_periodogram_bls.parquet"
  --sector-score-set 57 "${S57_SCORE_ROOT}/shape_plus_raw_chronology.parquet" "${S57_SCORE_ROOT}/shape_plus_periodogram_bls.parquet"
  --out-dir "${OUT_ROOT}"
  --batch-index 0
)
if [[ -n "${TWIRL_MIXED_ENRICH_EXCLUDE:-}" ]]; then
  read -r -a exclude_paths <<< "${TWIRL_MIXED_ENRICH_EXCLUDE}"
  build_args+=(--exclude "${exclude_paths[@]}")
fi
"${PYTHON}" scripts/stage5_validation/build_teacher_v2_mixed_enrichment_batch.py "${build_args[@]}"

for sector in 56 57; do
  data_root="${S56_DATA_ROOT}"
  if [[ "${sector}" == "57" ]]; then
    data_root="${S57_DATA_ROOT}"
  fi
  "${PYTHON}" scripts/stage5_validation/render_two_aperture_vet_sheets.py \
    --queue-csv "${OUT_ROOT}/sector_$(printf '%04d' "${sector}")_review_queue_500.csv" \
    --lc-export-h5 "${data_root}/s${sector}_A2v1_adp_pair_rebuilt.h5" \
    --out-dir "${OUT_ROOT}/twirl_vet_sheets" \
    --metrics-csv "${OUT_ROOT}/twirl_vet_metrics_s${sector}.csv" \
    --summary-json "${OUT_ROOT}/render_summary_s${sector}.json" \
    --branch-name current_adp \
    --apertures DET_FLUX_ADP_SML DET_FLUX_ADP \
    --workers "${SLURM_CPUS_PER_TASK:-48}" \
    --n-periods 200000 \
    --n-peaks 10 \
    --use-row-ephemeris \
    --no-pdf \
    --harmonic-factors 0.25 0.5 1 2 4 \
    --overwrite
done

"${PYTHON}" scripts/stage5_validation/verify_teacher_v2_enrichment_batch.py \
  --queue "${OUT_ROOT}/review_queue_1k.csv" \
  --overlap "${OUT_ROOT}/double_review_queue_100.csv" \
  --hidden "${OUT_ROOT}/hidden_selection_provenance.parquet" \
  --sheet-dir "${OUT_ROOT}/twirl_vet_sheets" \
  --out-json "${OUT_ROOT}/verification.json"
