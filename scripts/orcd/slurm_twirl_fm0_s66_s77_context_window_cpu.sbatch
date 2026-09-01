#!/bin/bash
# Fixed-exposure 256/512/1024/2048 context diagnostic on S66--S77 development data.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-context-window
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 027

readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
readonly RUN_DIR="${PROJECT_ROOT}/reports/stage5_validation/twirl_fm0_2_s56_s64_objective_canary/ddf442aafb8f/model_runs/fm0_2_1_canary/seed560067"
readonly STEP0="${RUN_DIR}/checkpoint_step_00000000.pt"
readonly STEP2000="${RUN_DIR}/checkpoint_step_00002000.pt"
readonly STEP0_SHA256="92463070381486f0c6053c190d62da8d0c5c0d31be8072e2bdcd677329ac792c"
readonly STEP2000_SHA256="976b5053c857c38b9fbf7a35c9d0605f0023318b2c9bd37a88995992e5aa7bd2"
: "${TWIRL_ORCD_REPO:?missing exact evaluator repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?missing exact evaluator Git SHA}"
: "${TWIRL_FM0_TEMPORAL_PANEL_DIR:?missing frozen temporal-panel directory}"
: "${TWIRL_FM0_TEMPORAL_PANEL_RECEIPT_SHA256:?missing temporal-panel receipt SHA-256}"
: "${TWIRL_FM0_CONTEXT_WINDOW_OUTPUT:?missing context diagnostic output path}"

readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly PANEL="$(realpath -m -- "${TWIRL_FM0_TEMPORAL_PANEL_DIR}")"
readonly OUTPUT="$(realpath -m -- "${TWIRL_FM0_CONTEXT_WINDOW_OUTPUT}")"
[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || { echo "invalid Git SHA" >&2; exit 2; }
[[ "${TWIRL_FM0_TEMPORAL_PANEL_RECEIPT_SHA256}" =~ ^[0-9a-f]{64}$ ]] || { echo "invalid panel SHA-256" >&2; exit 2; }
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" ]] || { echo "repository revision differs" >&2; exit 2; }
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || { echo "repository is dirty" >&2; exit 2; }
[[ "$(sha256sum "${PANEL}/receipt.json" | awk '{print $1}')" == "${TWIRL_FM0_TEMPORAL_PANEL_RECEIPT_SHA256}" ]] || { echo "panel receipt differs" >&2; exit 2; }
[[ "$(sha256sum "${STEP0}" | awk '{print $1}')" == "${STEP0_SHA256}" ]] || { echo "step-0 checkpoint differs" >&2; exit 2; }
[[ "$(sha256sum "${STEP2000}" | awk '{print $1}')" == "${STEP2000_SHA256}" ]] || { echo "step-2000 checkpoint differs" >&2; exit 2; }
[[ ! -e "${OUTPUT}" && ! -e "${OUTPUT}.sha256" ]] || { echo "refusing to overwrite context diagnostic" >&2; exit 2; }

mkdir -p "$(dirname "${OUTPUT}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
"${PYTHON}" scripts/stage5_validation/evaluate_twirl_fm0_context_window_diagnostic.py \
  --run-dir "${RUN_DIR}" \
  --step0-checkpoint "${STEP0}" \
  --step2000-checkpoint "${STEP2000}" \
  --temporal-panel-dir "${PANEL}" \
  --temporal-panel-receipt-sha256 "${TWIRL_FM0_TEMPORAL_PANEL_RECEIPT_SHA256}" \
  --output "${OUTPUT}" \
  --max-repeated-components 64 \
  --max-new-components 64 \
  --batch-size 32 \
  --expected-git-sha "${TWIRL_EXPECTED_GIT_SHA}"

echo "[fm0-context-window] complete=${OUTPUT}; no training, event, label, or sealed access"
