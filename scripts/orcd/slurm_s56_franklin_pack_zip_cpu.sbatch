#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-franklin-zip
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
QUEUE_DIR="${TWIRL_FRANKLIN_QUEUE_DIR:-reports/stage5_validation/s56_franklin_real5k_handoff}"
EXPECTED_ROWS="${TWIRL_FRANKLIN_EXPECTED_ROWS:-5000}"
OUT_ZIP="${TWIRL_FRANKLIN_OUT_ZIP:-reports/stage5_validation/s56_franklin_real5k_handoff_png.zip}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs
cd "${REPO}"

export PYTHONPATH=src

echo "[franklin-zip] host=$(hostname)"
echo "[franklin-zip] queue_dir=${QUEUE_DIR}"
echo "[franklin-zip] out_zip=${OUT_ZIP}"

"${PYTHON}" scripts/stage5_validation/franklin_vetting_app.py \
  --queue "${QUEUE_DIR}/franklin_review_queue_5k_real.csv" \
  --labels-out "${QUEUE_DIR}/franklin_labels_vetted.csv" \
  --sheet-root "${QUEUE_DIR}/vet_sheets" \
  --check-only

"${PYTHON}" scripts/stage5_validation/package_s56_franklin_real_handoff.py \
  --queue-dir "${QUEUE_DIR}" \
  --expected-rows "${EXPECTED_ROWS}" \
  --out-zip "${OUT_ZIP}"

echo "[franklin-zip] complete $(date -Is)"
