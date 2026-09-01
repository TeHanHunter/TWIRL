#!/bin/bash
# Freeze the label-blind S66--S77 development panel; no shard payload is opened.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-temporal-panel
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 027

readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
readonly ADMISSION="${PROJECT_ROOT}/reports/stage5_validation/twirl_fm0_s66_s77_full_visit_preparation_v1/admission_v2/receipt.json"
readonly BASELINE_MANIFEST="${PROJECT_ROOT}/exports/fm0/input_releases/twirl_fm0_1_s56_s64_v2/manifest.csv"
readonly BASELINE_MANIFEST_SHA256="f9c62ba02ea40ad7b2f457cd4dab55526ee7a453ead8a1a553c4de9f2969598a"
: "${TWIRL_ORCD_REPO:?missing exact evaluator repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?missing exact evaluator Git SHA}"
: "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256:?missing admission receipt SHA-256}"
: "${TWIRL_FM0_TEMPORAL_PANEL_DIR:?missing temporal-panel output directory}"

readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly OUTPUT="$(realpath -m -- "${TWIRL_FM0_TEMPORAL_PANEL_DIR}")"
readonly CONFIG="${REPO}/configs/models/twirl_fm0_2_s66_s77_temporal_panel_v1.yaml"
[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || { echo "invalid Git SHA" >&2; exit 2; }
[[ "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}" =~ ^[0-9a-f]{64}$ ]] || { echo "invalid admission SHA-256" >&2; exit 2; }
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" ]] || { echo "repository revision differs" >&2; exit 2; }
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || { echo "repository is dirty" >&2; exit 2; }
[[ "$(sha256sum "${ADMISSION}" | awk '{print $1}')" == "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}" ]] || { echo "admission receipt differs" >&2; exit 2; }
[[ "$(sha256sum "${BASELINE_MANIFEST}" | awk '{print $1}')" == "${BASELINE_MANIFEST_SHA256}" ]] || { echo "baseline manifest differs" >&2; exit 2; }
[[ ! -e "${OUTPUT}" && ! -e "${OUTPUT}.partial" ]] || { echo "refusing to overwrite temporal panel" >&2; exit 2; }

cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
readonly CONFIG_SHA256="$(sha256sum "${CONFIG}" | awk '{print $1}')"
"${PYTHON}" scripts/stage5_validation/freeze_twirl_fm0_s66_s77_temporal_panel.py \
  --config "${CONFIG}" \
  --config-sha256 "${CONFIG_SHA256}" \
  --admission-receipt "${ADMISSION}" \
  --admission-receipt-sha256 "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}" \
  --baseline-manifest "${BASELINE_MANIFEST}" \
  --baseline-manifest-sha256 "${BASELINE_MANIFEST_SHA256}" \
  --producer-git-sha "${TWIRL_EXPECTED_GIT_SHA}" \
  --output-dir "${OUTPUT}"

echo "[fm0-temporal-panel] complete=${OUTPUT}"
