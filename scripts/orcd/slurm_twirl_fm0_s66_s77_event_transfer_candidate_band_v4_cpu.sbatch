#!/bin/bash
# Fixed-band, BLS-free cadence-transfer screen on S66--S77 development.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-band-v4
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=192G
#SBATCH -t 08:00:00
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
readonly PANEL_SHA256="78c370e10c556472c5997c20cfe95207a0b334bafe7f024bf7ba4fc7ec4de624"
readonly CONFIG_SHA256="732f525e1afb4aa9bf55966f8d55bfa807393784fd13bb5ffd5127c4fc6e4b94"
readonly MECHANICS_RESULT="${PROJECT_ROOT}/reports/stage5_validation/twirl_fm0_2_s66_s77_event_transfer_mechanics_v3/7ae0059479fd/result.json"
readonly MECHANICS_SHA256="10b78bf083595ad8f4ecacb0428d99b939abc9dae5e3f361d02559659d1e2833"

: "${TWIRL_ORCD_REPO:?missing exact evaluator repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?missing exact evaluator Git SHA}"
: "${TWIRL_FM0_TEMPORAL_PANEL_DIR:?missing frozen temporal-panel directory}"
: "${TWIRL_FM0_TEMPORAL_PANEL_RECEIPT_SHA256:?missing panel receipt SHA-256}"
: "${TWIRL_FM0_EVENT_TRANSFER_BAND_OUTPUT:?missing candidate-band output path}"

readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly PANEL="$(realpath -m -- "${TWIRL_FM0_TEMPORAL_PANEL_DIR}")"
readonly CONFIG="${REPO}/configs/models/twirl_fm0_2_s66_s77_event_transfer_candidate_band_v4.yaml"
readonly OUTPUT="$(realpath -m -- "${TWIRL_FM0_EVENT_TRANSFER_BAND_OUTPUT}")"

[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || { echo "invalid Git SHA" >&2; exit 2; }
[[ "${TWIRL_FM0_TEMPORAL_PANEL_RECEIPT_SHA256}" == "${PANEL_SHA256}" ]] || { echo "unexpected panel SHA-256" >&2; exit 2; }
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" ]] || { echo "repository revision differs" >&2; exit 2; }
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || { echo "repository is dirty" >&2; exit 2; }
[[ "$(sha256sum "${CONFIG}" | awk '{print $1}')" == "${CONFIG_SHA256}" ]] || { echo "candidate-band config differs" >&2; exit 2; }
[[ "$(sha256sum "${PANEL}/receipt.json" | awk '{print $1}')" == "${PANEL_SHA256}" ]] || { echo "panel receipt differs" >&2; exit 2; }
[[ "$(sha256sum "${MECHANICS_RESULT}" | awk '{print $1}')" == "${MECHANICS_SHA256}" ]] || { echo "v3 mechanics result differs" >&2; exit 2; }
[[ "$(sha256sum "${STEP0}" | awk '{print $1}')" == "${STEP0_SHA256}" ]] || { echo "step-0 checkpoint differs" >&2; exit 2; }
[[ "$(sha256sum "${STEP2000}" | awk '{print $1}')" == "${STEP2000_SHA256}" ]] || { echo "step-2000 checkpoint differs" >&2; exit 2; }
[[ "${OUTPUT}" == "${PROJECT_ROOT}/reports/stage5_validation/"* ]] || { echo "output escaped report root" >&2; exit 2; }
[[ ! -e "${OUTPUT}" && ! -e "${OUTPUT}.sha256" ]] || { echo "refusing to overwrite candidate-band output" >&2; exit 2; }

mkdir -p "$(dirname "${OUTPUT}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export LD_LIBRARY_PATH="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=32 MKL_NUM_THREADS=32
"${PYTHON}" scripts/stage5_validation/evaluate_twirl_fm0_event_transfer_candidate_band.py \
  --config "${CONFIG}" \
  --mechanics-result "${MECHANICS_RESULT}" \
  --run-dir "${RUN_DIR}" \
  --step0-checkpoint "${STEP0}" \
  --step2000-checkpoint "${STEP2000}" \
  --temporal-panel-dir "${PANEL}" \
  --temporal-panel-receipt-sha256 "${PANEL_SHA256}" \
  --output "${OUTPUT}" \
  --batch-size 32 \
  --expected-git-sha "${TWIRL_EXPECTED_GIT_SHA}"

echo "[fm0-band-v4] complete=${OUTPUT}; fixed indices 36--92; no cadence merging, FM training, BLS, GPU, sealed access, or architecture selection"
