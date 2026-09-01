#!/bin/bash
# Paired step-0/step-2000 evaluation on the frozen S66--S77 development panel.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-temporal-eval
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 02:00:00
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
readonly BASELINE_MANIFEST="${PROJECT_ROOT}/exports/fm0/input_releases/twirl_fm0_1_s56_s64_v2/manifest.csv"
readonly BASELINE_MANIFEST_SHA256="f9c62ba02ea40ad7b2f457cd4dab55526ee7a453ead8a1a553c4de9f2969598a"
: "${TWIRL_ORCD_REPO:?missing exact evaluator repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?missing exact evaluator Git SHA}"
: "${TWIRL_FM0_TEMPORAL_PANEL_DIR:?missing frozen temporal-panel directory}"
: "${TWIRL_FM0_TEMPORAL_PANEL_RECEIPT_SHA256:?missing temporal-panel receipt SHA-256}"
: "${TWIRL_FM0_TEMPORAL_EVALUATION_OUTPUT:?missing zero-shot output path}"

readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly PANEL="$(realpath -m -- "${TWIRL_FM0_TEMPORAL_PANEL_DIR}")"
readonly OUTPUT="$(realpath -m -- "${TWIRL_FM0_TEMPORAL_EVALUATION_OUTPUT}")"
[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || { echo "invalid Git SHA" >&2; exit 2; }
[[ "${TWIRL_FM0_TEMPORAL_PANEL_RECEIPT_SHA256}" =~ ^[0-9a-f]{64}$ ]] || { echo "invalid panel SHA-256" >&2; exit 2; }
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" ]] || { echo "repository revision differs" >&2; exit 2; }
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || { echo "repository is dirty" >&2; exit 2; }
[[ "$(sha256sum "${PANEL}/receipt.json" | awk '{print $1}')" == "${TWIRL_FM0_TEMPORAL_PANEL_RECEIPT_SHA256}" ]] || { echo "panel receipt differs" >&2; exit 2; }
[[ "$(sha256sum "${STEP0}" | awk '{print $1}')" == "${STEP0_SHA256}" ]] || { echo "step-0 checkpoint differs" >&2; exit 2; }
[[ "$(sha256sum "${STEP2000}" | awk '{print $1}')" == "${STEP2000_SHA256}" ]] || { echo "step-2000 checkpoint differs" >&2; exit 2; }
[[ "$(sha256sum "${BASELINE_MANIFEST}" | awk '{print $1}')" == "${BASELINE_MANIFEST_SHA256}" ]] || { echo "baseline manifest differs" >&2; exit 2; }
[[ ! -e "${OUTPUT}" && ! -e "${OUTPUT}.sha256" ]] || { echo "refusing to overwrite temporal evaluation" >&2; exit 2; }

mkdir -p "$(dirname "${OUTPUT}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
"${PYTHON}" scripts/stage5_validation/evaluate_twirl_fm0_temporal_zero_shot.py \
  --run-dir "${RUN_DIR}" \
  --step0-checkpoint "${STEP0}" \
  --step2000-checkpoint "${STEP2000}" \
  --temporal-panel-dir "${PANEL}" \
  --temporal-panel-receipt-sha256 "${TWIRL_FM0_TEMPORAL_PANEL_RECEIPT_SHA256}" \
  --baseline-manifest "${BASELINE_MANIFEST}" \
  --baseline-manifest-sha256 "${BASELINE_MANIFEST_SHA256}" \
  --output "${OUTPUT}" \
  --device cpu \
  --max-repeated-components 256 \
  --max-new-components 256 \
  --max-baseline-train-components 512 \
  --pca-components 32 \
  --batch-size 8 \
  --expected-git-sha "${TWIRL_EXPECTED_GIT_SHA}"

echo "[fm0-temporal-eval] complete=${OUTPUT}; no training, event, or sealed access"
