#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-label-adjud
#SBATCH -N 1
#SBATCH -c 48
#SBATCH --mem=192G
#SBATCH -t 08:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
OUT_DIR="${TWIRL_ADJUDICATION_ROOT:-reports/stage5_validation/s56_label_adjudication_real343}"
COMPACT_LC="${TWIRL_ADJUDICATION_COMPACT_LC:-data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp_lc_export_pdo.h5}"
HLSP_ROOT="${TWIRL_ADJUDICATION_HLSP_ROOT:-data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare}"
WORKERS="${TWIRL_ADJUDICATION_WORKERS:-${SLURM_CPUS_PER_TASK:-48}}"
N_PERIODS="${TWIRL_ADJUDICATION_N_PERIODS:-20000}"
OVERWRITE="${TWIRL_ADJUDICATION_OVERWRITE:-1}"

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

echo "[label-adjudication] host=$(hostname) date=$(date -Is)"
echo "[label-adjudication] repo=${REPO} out=${OUT_DIR} workers=${WORKERS}"

"${PYTHON}" scripts/stage5_validation/build_s56_real_label_adjudication_queue.py --out-dir "${OUT_DIR}"

mkdir -p "${OUT_DIR}/twirl_vet_sheets"
find "${OUT_DIR}/twirl_vet_sheets" -maxdepth 1 -type f -name "*.pdf" -delete

render_args=(
  --queue-csv "${OUT_DIR}/review_queue_real343.csv"
  --hlsp-root "${HLSP_ROOT}"
  --lc-export-h5 "${COMPACT_LC}"
  --out-dir "${OUT_DIR}/twirl_vet_sheets"
  --metrics-csv "${OUT_DIR}/twirl_vet_metrics_current_adp.csv"
  --summary-json "${OUT_DIR}/twirl_two_aperture_vet_summary_current_adp.json"
  --branch-name current_adp
  --apertures DET_FLUX_ADP_SML DET_FLUX_ADP
  --workers "${WORKERS}"
  --n-periods "${N_PERIODS}"
  --n-peaks 10
  --use-row-ephemeris
  --no-pdf
  --harmonic-factors 0.25 0.5 1 2 4
)
if [[ "${OVERWRITE}" == "1" ]]; then
  render_args+=(--overwrite)
fi

"${PYTHON}" scripts/stage5_validation/render_two_aperture_vet_sheets.py "${render_args[@]}"

"${PYTHON}" scripts/stage5_validation/verify_s56_real_label_adjudication.py \
  --root "${OUT_DIR}" \
  --require-sheets

echo "[label-adjudication] complete $(date -Is)"
