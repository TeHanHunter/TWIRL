#!/bin/bash
#SBATCH -p mit_normal
#SBATCH --job-name=twirl-rec50-inj-vet
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
WORKERS="${TWIRL_REC50_VET_WORKERS:-${SLURM_CPUS_PER_TASK:-8}}"
QUEUE_DIR="${TWIRL_REC50_QUEUE_DIR:-reports/stage5_validation/s56_recovery50_teacher_queue}"
INJ_QUEUE_CSV="${TWIRL_REC50_INJ_QUEUE_CSV:-${QUEUE_DIR}/review_queue_1k_injected.csv}"
INJECTION_H5="${TWIRL_REC50_INJECTION_H5:-data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_20k_predetrend_batman_periodradius_bright_balanced/injected_lightcurves.h5}"
HLSP_ROOT="${TWIRL_REC50_HLSP_ROOT:-data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare}"
OUT_DIR="${TWIRL_REC50_OUT_DIR:-${QUEUE_DIR}/twirl_vet_sheets_rerendered}"
METRICS_CSV="${TWIRL_REC50_INJ_METRICS_CSV:-${QUEUE_DIR}/twirl_vet_metrics_injected_rerendered.csv}"
SUMMARY_JSON="${TWIRL_REC50_INJ_SUMMARY_JSON:-${QUEUE_DIR}/twirl_two_aperture_vet_summary_injected_rerendered.json}"
INJ_APERTURES="${TWIRL_REC50_INJ_APERTURES:-DET_FLUX_ADP_SML DET_FLUX_ADP}"
N_PERIODS="${TWIRL_REC50_N_PERIODS:-20000}"
N_PEAKS="${TWIRL_REC50_N_PEAKS:-10}"
LIMIT="${TWIRL_REC50_LIMIT:-0}"

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

echo "[orcd-rec50-inj-vet] host=$(hostname)"
echo "[orcd-rec50-inj-vet] repo=${REPO}"
echo "[orcd-rec50-inj-vet] python=${PYTHON}"
echo "[orcd-rec50-inj-vet] inj_queue=${INJ_QUEUE_CSV}"
echo "[orcd-rec50-inj-vet] injection_h5=${INJECTION_H5}"
echo "[orcd-rec50-inj-vet] out_dir=${OUT_DIR}"
echo "[orcd-rec50-inj-vet] apertures=${INJ_APERTURES}"
echo "[orcd-rec50-inj-vet] workers=${WORKERS} n_periods=${N_PERIODS} limit=${LIMIT}"

read -r -a inj_aperture_args <<< "${INJ_APERTURES}"

args=(
  --queue-csv "${INJ_QUEUE_CSV}"
  --hlsp-root "${HLSP_ROOT}"
  --injection-h5 "${INJECTION_H5}"
  --out-dir "${OUT_DIR}"
  --metrics-csv "${METRICS_CSV}"
  --summary-json "${SUMMARY_JSON}"
  --branch-name current_adp
  --apertures "${inj_aperture_args[@]}"
  --workers "${WORKERS}"
  --n-periods "${N_PERIODS}"
  --n-peaks "${N_PEAKS}"
  --overwrite
)

if [[ "${LIMIT}" != "0" ]]; then
  args+=(--limit "${LIMIT}")
fi

"${PYTHON}" scripts/stage5_validation/render_two_aperture_vet_sheets.py "${args[@]}"

echo "[orcd-rec50-inj-vet] complete $(date -Is)"
