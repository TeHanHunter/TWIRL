#!/bin/bash
# Development-only FM0.1 representation health; never opens the sealed test.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-rep-health
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 027
readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
: "${TWIRL_ORCD_REPO:?missing detached repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?missing exact evaluator SHA}"
: "${TWIRL_FM0_RUN_ROOT:?missing run root}"
: "${TWIRL_FM0_TRAIN_OUTPUT:?missing exact real training output}"
: "${TWIRL_FM0_INPUT_RECEIPT:?missing exact input receipt}"
: "${TWIRL_FM0_INPUT_RECEIPT_SHA256:?missing input receipt hash}"
: "${TWIRL_FM0_REAL_POST_VALIDATION:?missing independent post-validation receipt}"
: "${TWIRL_FM0_REAL_POST_VALIDATION_SHA256:?missing post-validation receipt hash}"
readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly RUN_ROOT="$(realpath -m -- "${TWIRL_FM0_RUN_ROOT}")"
readonly TRAIN="$(realpath -m -- "${TWIRL_FM0_TRAIN_OUTPUT}")"
readonly INPUT_RECEIPT="$(realpath -m -- "${TWIRL_FM0_INPUT_RECEIPT}")"
readonly POST_VALIDATION="$(realpath -m -- "${TWIRL_FM0_REAL_POST_VALIDATION}")"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
readonly OUT="${RUN_ROOT}/evaluations/representation_health-$(basename "${TRAIN}").json"
die() { echo "[fm0-representation-health] ERROR: $*" >&2; exit 2; }
for path in "${TRAIN}" "${INPUT_RECEIPT}" "${POST_VALIDATION}"; do
  [[ "${path}" == "${RUN_ROOT}"/* ]] || die "authority path escapes run root: ${path}"
done
[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ && "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" =~ ^[0-9a-f]{64}$ && "${TWIRL_FM0_REAL_POST_VALIDATION_SHA256}" =~ ^[0-9a-f]{64}$ ]] || die "invalid hash"
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" && -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is not exact and clean"
[[ "$(sha256sum "${INPUT_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" ]] || die "input receipt drift"
[[ "$(sha256sum "${POST_VALIDATION}" | awk '{print $1}')" == "${TWIRL_FM0_REAL_POST_VALIDATION_SHA256}" ]] || die "post-validation receipt drift"
[[ ! -e "${OUT}" ]] || die "refusing to overwrite representation-health output"
"${PYTHON}" - "${POST_VALIDATION}" "${TRAIN}" "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" <<'PY'
import json, sys
from pathlib import Path
post_path = Path(sys.argv[1])
train_path = Path(sys.argv[2])
input_sha = sys.argv[3]
post = json.loads(post_path.read_text())
if post.get("schema_version") != "twirl_fm0_1_orcd_real_post_validation_v1" or post.get("passed") is not True:
    raise SystemExit("post-validation receipt did not pass")
if post.get("input_receipt_sha256") != input_sha:
    raise SystemExit("post-validation input binding drift")
contract = json.loads((train_path / "run_contract.json").read_text())
if contract.get("input_release", {}).get("receipt_sha256") != input_sha:
    raise SystemExit("training input binding drift")
if contract.get("synthetic_smoke") or contract.get("real_data_consumed") is not True:
    raise SystemExit("representation health requires a real-data run")
PY
mkdir -p "$(dirname "${OUT}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
"${PYTHON}" scripts/stage5_validation/evaluate_twirl_fm0_representation_health.py \
  --run-dir "${TRAIN}" \
  --output "${OUT}" \
  --device cpu \
  --max-components 256 \
  --batch-size 8 \
  --random-control-seed 0 \
  --expected-git-sha "${TWIRL_EXPECTED_GIT_SHA}"
[[ -s "${OUT}" && -s "${OUT}.sha256" ]] || die "missing representation-health output/sidecar"
echo "[fm0-representation-health] complete=${OUT}"
