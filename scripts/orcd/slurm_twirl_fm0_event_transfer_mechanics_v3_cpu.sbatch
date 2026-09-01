#!/bin/bash
# Fast raw/quality cadence-probe mechanics control; no FM checkpoint or GPU.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-mechanics-v3
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 027

readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
readonly PANEL_SHA256="78c370e10c556472c5997c20cfe95207a0b334bafe7f024bf7ba4fc7ec4de624"
readonly CONFIG_SHA256="c820743391645bc626effae097c3e7e8e8810c8376148105910b60b14e7c78b2"

: "${TWIRL_ORCD_REPO:?missing exact evaluator repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?missing exact evaluator Git SHA}"
: "${TWIRL_FM0_TEMPORAL_PANEL_DIR:?missing frozen temporal-panel directory}"
: "${TWIRL_FM0_TEMPORAL_PANEL_RECEIPT_SHA256:?missing panel receipt SHA-256}"
: "${TWIRL_FM0_EVENT_TRANSFER_MECHANICS_OUTPUT:?missing mechanics output path}"

readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly PANEL="$(realpath -m -- "${TWIRL_FM0_TEMPORAL_PANEL_DIR}")"
readonly CONFIG="${REPO}/configs/models/twirl_fm0_2_s66_s77_event_transfer_mechanics_v3.yaml"
readonly OUTPUT="$(realpath -m -- "${TWIRL_FM0_EVENT_TRANSFER_MECHANICS_OUTPUT}")"

[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || { echo "invalid Git SHA" >&2; exit 2; }
[[ "${TWIRL_FM0_TEMPORAL_PANEL_RECEIPT_SHA256}" == "${PANEL_SHA256}" ]] || { echo "unexpected panel SHA-256" >&2; exit 2; }
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" ]] || { echo "repository revision differs" >&2; exit 2; }
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || { echo "repository is dirty" >&2; exit 2; }
[[ "$(sha256sum "${CONFIG}" | awk '{print $1}')" == "${CONFIG_SHA256}" ]] || { echo "mechanics config differs" >&2; exit 2; }
[[ "$(sha256sum "${PANEL}/receipt.json" | awk '{print $1}')" == "${PANEL_SHA256}" ]] || { echo "panel receipt differs" >&2; exit 2; }
[[ "${OUTPUT}" == "${PROJECT_ROOT}/reports/stage5_validation/"* ]] || { echo "output escaped report root" >&2; exit 2; }
[[ ! -e "${OUTPUT}" && ! -e "${OUTPUT}.sha256" ]] || { echo "refusing to overwrite mechanics output" >&2; exit 2; }

mkdir -p "$(dirname "${OUTPUT}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
"${PYTHON}" scripts/stage5_validation/evaluate_twirl_fm0_event_transfer_mechanics.py \
  --config "${CONFIG}" \
  --temporal-panel-dir "${PANEL}" \
  --temporal-panel-receipt-sha256 "${PANEL_SHA256}" \
  --output "${OUTPUT}" \
  --expected-git-sha "${TWIRL_EXPECTED_GIT_SHA}"

echo "[fm0-mechanics-v3] complete=${OUTPUT}; raw/quality center-cadence mechanics only"
