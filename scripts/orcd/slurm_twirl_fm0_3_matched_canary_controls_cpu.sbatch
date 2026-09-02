#!/bin/bash
# Run only the frozen matched-FM0.3 raw and quality controls on CPU.
# This reopens the payload schedule's 1,440 shards and performs no FM inference.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-3-controls
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 027

readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
readonly FM0_ENV_LIB="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/lib"
readonly PANEL_RECEIPT_SHA="78c370e10c556472c5997c20cfe95207a0b334bafe7f024bf7ba4fc7ec4de624"

: "${TWIRL_ORCD_REPO:?missing exact detached repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?missing exact Git SHA}"
: "${TWIRL_FM0_MATCHED_CANARY_PAYLOAD_PLAN_DIR:?missing payload-plan root}"
: "${TWIRL_FM0_MATCHED_CANARY_PAYLOAD_PLAN_RECEIPT_SHA256:?missing payload-plan receipt SHA-256}"
: "${TWIRL_FM0_MATCHED_CANARY_CONTROLS_DIR:?missing controls output root}"

readonly REPO="$(realpath -e -- "${TWIRL_ORCD_REPO}")"
readonly PAYLOAD_PLAN="$(realpath -e -- "${TWIRL_FM0_MATCHED_CANARY_PAYLOAD_PLAN_DIR}")"
readonly OUTPUT="$(realpath -m -- "${TWIRL_FM0_MATCHED_CANARY_CONTROLS_DIR}")"

die() { echo "[fm0.3-controls] ERROR: $*" >&2; exit 2; }
digest() { sha256sum "$1" | awk '{print $1}'; }

[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || die "invalid Git SHA"
[[ "${TWIRL_FM0_MATCHED_CANARY_PAYLOAD_PLAN_RECEIPT_SHA256}" =~ ^[0-9a-f]{64}$ ]] || die "invalid payload-plan receipt SHA-256"
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" ]] || die "Git revision drift"
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is dirty"
[[ "${PAYLOAD_PLAN}" == "${PROJECT_ROOT}/reports/stage5_validation/"* ]] || die "payload-plan root escaped report tree"
[[ "${OUTPUT}" == "${PROJECT_ROOT}/reports/stage5_validation/"* ]] || die "output root escaped report tree"
[[ ! -e "${OUTPUT}" ]] || die "refusing to overwrite controls output"
[[ "$(digest "${PAYLOAD_PLAN}/receipt.json")" == "${TWIRL_FM0_MATCHED_CANARY_PAYLOAD_PLAN_RECEIPT_SHA256}" ]] || die "payload-plan receipt drift"
[[ "$(tr -d '\r\n' < "${PAYLOAD_PLAN}/READY")" == "${TWIRL_FM0_MATCHED_CANARY_PAYLOAD_PLAN_RECEIPT_SHA256}" ]] || die "payload-plan READY drift"

mkdir -p "$(dirname "${OUTPUT}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export LD_LIBRARY_PATH="${FM0_ENV_LIB}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"${PYTHON}" scripts/stage5_validation/evaluate_twirl_fm0_3_matched_canary.py \
  --payload-plan-root "${PAYLOAD_PLAN}" \
  --payload-plan-receipt-sha256 "${TWIRL_FM0_MATCHED_CANARY_PAYLOAD_PLAN_RECEIPT_SHA256}" \
  --output "${OUTPUT}" \
  --expected-git-sha "${TWIRL_EXPECTED_GIT_SHA}" \
  --controls-only \
  --device cpu

"${PYTHON}" - "${OUTPUT}" "${PANEL_RECEIPT_SHA}" <<'PY'
import json
import sys
from pathlib import Path

from twirl.models.fm0.matched_canary_evaluation import (
    CONTROL_ARMS,
    validate_matched_canary_evaluation,
)

root = Path(sys.argv[1])
panel_receipt_sha = sys.argv[2]
receipt_sha = (root / "READY").read_text(encoding="utf-8").strip()
result = validate_matched_canary_evaluation(
    root, expected_receipt_sha256=receipt_sha, require_read_only=True
)
receipt = result.receipt
if (
    tuple(receipt["feature_arms"]) != CONTROL_ARMS
    or receipt["gate_summary"]["evaluation_scope"] != "controls_only"
    or receipt["gate_summary"]["controls_preflight"]["passed"] is not True
    or receipt["payload_plan"]["temporal_panel_authority"]["receipt_sha256"]
    != panel_receipt_sha
    or receipt["model_checkpoint_bindings"]
    or receipt["limits"]["fm_encoder_training_executed"] is not False
    or receipt["limits"]["bls_or_search_executed"] is not False
    or receipt["limits"]["sealed_test_opened"] is not False
):
    raise SystemExit("controls-only preflight did not pass its frozen contract")
print(json.dumps({
    "receipt_sha256": result.receipt_sha256,
    "results_sha256": result.results_sha256,
    "scores_sha256": result.scores_sha256,
    "controls_preflight_passed": True,
}, sort_keys=True))
PY

echo "[fm0.3-controls] complete=${OUTPUT}; controls only; no FM inference or training"
