#!/bin/bash
# Training-free 128-cadence TCN/stride-1-Conformer mechanics screen.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm03-mechanics
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 027

readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
readonly CONFIG_SHA256="6e16b2d1ab65fcfd3995d38903ac94ad9281fa768b233b9ab7a64510d7798612"

: "${TWIRL_ORCD_REPO:?missing exact evaluator repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?missing exact evaluator Git SHA}"
: "${TWIRL_FM0_MECHANICS_OUTPUT:?missing mechanics output path}"

readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly CONFIG="${REPO}/configs/models/twirl_fm0_3_stride1_mechanics_v1.yaml"
readonly OUTPUT="$(realpath -m -- "${TWIRL_FM0_MECHANICS_OUTPUT}")"

[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || { echo "invalid Git SHA" >&2; exit 2; }
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" ]] || { echo "repository revision differs" >&2; exit 2; }
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || { echo "repository is dirty" >&2; exit 2; }
[[ "$(sha256sum "${CONFIG}" | awk '{print $1}')" == "${CONFIG_SHA256}" ]] || { echo "mechanics config differs" >&2; exit 2; }
[[ "${OUTPUT}" == "${PROJECT_ROOT}/reports/stage5_validation/"* ]] || { echo "output escaped report root" >&2; exit 2; }
[[ ! -e "${OUTPUT}" && ! -e "${OUTPUT}.sha256" ]] || { echo "refusing to overwrite mechanics output" >&2; exit 2; }

mkdir -p "$(dirname "${OUTPUT}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
"${PYTHON}" scripts/stage5_validation/evaluate_twirl_fm0_3_stride1_mechanics.py \
  --config "${CONFIG}" \
  --output "${OUTPUT}" \
  --expected-git-sha "${TWIRL_EXPECTED_GIT_SHA}"

echo "[fm03-mechanics] complete=${OUTPUT}; zero training steps and one token per cadence"
