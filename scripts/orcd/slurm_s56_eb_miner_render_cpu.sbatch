#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-eb-render
#SBATCH -N 1
#SBATCH -c 48
#SBATCH --mem=160G
#SBATCH -t 06:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
ROOT="${TWIRL_EB_MINER_ROOT:-reports/stage5_validation/s56_eb_miner_adp_only}"
QUEUE_DIR="${TWIRL_EB_QUEUE_DIR:-${ROOT}/review_queue_eb_priority_pilot100}"
SCORED="${TWIRL_EB_SCORED:-${ROOT}/eb_cnn_miner/scored_candidates_ensemble.csv}"
MINER_SUMMARY="${TWIRL_EB_MINER_SUMMARY:-${ROOT}/eb_cnn_miner/summary.json}"
COMPACT_LC="${TWIRL_EB_COMPACT_LC:-data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp_lc_export_pdo.h5}"
HLSP_ROOT="${TWIRL_EB_HLSP_ROOT:-data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare}"
WORKERS="${TWIRL_EB_RENDER_WORKERS:-${SLURM_CPUS_PER_TASK:-8}}"
N_PERIODS="${TWIRL_EB_N_PERIODS:-20000}"
N_PEAKS="${TWIRL_EB_N_PEAKS:-10}"
N_REVIEW="${TWIRL_EB_N_REVIEW:-100}"
OVERWRITE="${TWIRL_EB_OVERWRITE:-0}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs
mkdir -p /orcd/data/mki_aryeh/001/twirl/mpl-cache
cd "${REPO}"

export PYTHONPATH=src
export MPLCONFIGDIR=/orcd/data/mki_aryeh/001/twirl/mpl-cache
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[eb-render] host=$(hostname)"
echo "[eb-render] date=$(date -Is)"
echo "[eb-render] repo=${REPO}"
echo "[eb-render] python=${PYTHON}"
echo "[eb-render] queue_dir=${QUEUE_DIR}"
echo "[eb-render] scored=${SCORED}"

"${PYTHON}" scripts/stage5_validation/build_s56_eb_priority_queue.py \
  --scored-candidates "${SCORED}" \
  --miner-summary "${MINER_SUMMARY}" \
  --out-dir "${QUEUE_DIR}" \
  --n-review "${N_REVIEW}"

QUEUE_CSV="${QUEUE_DIR}/review_queue_eb_priority_${N_REVIEW}.csv"

mkdir -p "${QUEUE_DIR}/twirl_vet_sheets"
find "${QUEUE_DIR}/twirl_vet_sheets" -maxdepth 1 -type f -name "*.pdf" -delete

render_args=(
  --queue-csv "${QUEUE_CSV}"
  --hlsp-root "${HLSP_ROOT}"
  --lc-export-h5 "${COMPACT_LC}"
  --out-dir "${QUEUE_DIR}/twirl_vet_sheets"
  --metrics-csv "${QUEUE_DIR}/twirl_vet_metrics_current_adp.csv"
  --summary-json "${QUEUE_DIR}/twirl_two_aperture_vet_summary_current_adp.json"
  --branch-name current_adp
  --apertures DET_FLUX_ADP_SML DET_FLUX_ADP
  --workers "${WORKERS}"
  --n-periods "${N_PERIODS}"
  --n-peaks "${N_PEAKS}"
  --use-row-ephemeris
  --no-pdf
)
if [[ "${OVERWRITE}" == "1" ]]; then
  render_args+=(--overwrite)
fi

"${PYTHON}" scripts/stage5_validation/render_two_aperture_vet_sheets.py "${render_args[@]}"

"${PYTHON}" scripts/stage5_validation/franklin_vetting_app.py \
  --queue "${QUEUE_CSV}" \
  --labels-out "${QUEUE_DIR}/human_labels_vetted.csv" \
  --sheet-root "${QUEUE_DIR}/twirl_vet_sheets" \
  --check-only

echo "[eb-render] complete $(date -Is)"
