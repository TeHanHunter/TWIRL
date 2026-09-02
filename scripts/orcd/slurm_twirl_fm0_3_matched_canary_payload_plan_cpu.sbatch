#!/bin/bash
# Freeze exact direct-128 crops for the matched FM0.3 canary on CPU only.
# Opens only the 1,440 bound source shards; performs no inference or training.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-3-payload-plan
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 027

readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
readonly FM0_ENV_LIB="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/lib"

: "${TWIRL_ORCD_REPO:?missing exact detached repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?missing exact Git SHA}"
: "${TWIRL_FM0_MATCHED_CANARY_EVALUATION_PLAN_DIR:?missing identity-plan root}"
: "${TWIRL_FM0_MATCHED_CANARY_EVALUATION_PLAN_RECEIPT_SHA256:?missing identity-plan receipt SHA-256}"
: "${TWIRL_FM0_MATCHED_CANARY_PAYLOAD_PLAN_DIR:?missing payload-plan output root}"

readonly REPO="$(realpath -e -- "${TWIRL_ORCD_REPO}")"
readonly IDENTITY_PLAN="$(realpath -e -- "${TWIRL_FM0_MATCHED_CANARY_EVALUATION_PLAN_DIR}")"
readonly OUTPUT="$(realpath -m -- "${TWIRL_FM0_MATCHED_CANARY_PAYLOAD_PLAN_DIR}")"

die() { echo "[fm0.3-payload-plan] ERROR: $*" >&2; exit 2; }
digest() { sha256sum "$1" | awk '{print $1}'; }

[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || die "invalid Git SHA"
[[ "${TWIRL_FM0_MATCHED_CANARY_EVALUATION_PLAN_RECEIPT_SHA256}" =~ ^[0-9a-f]{64}$ ]] || die "invalid identity-plan receipt SHA-256"
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" ]] || die "Git revision drift"
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is dirty"
[[ "${IDENTITY_PLAN}" == "${PROJECT_ROOT}/reports/stage5_validation/"* ]] || die "identity-plan root escaped report tree"
[[ "${OUTPUT}" == "${PROJECT_ROOT}/reports/stage5_validation/"* ]] || die "output root escaped report tree"
[[ ! -e "${OUTPUT}" ]] || die "refusing to overwrite payload plan"
[[ "$(digest "${IDENTITY_PLAN}/receipt.json")" == "${TWIRL_FM0_MATCHED_CANARY_EVALUATION_PLAN_RECEIPT_SHA256}" ]] || die "identity-plan receipt drift"
[[ "$(tr -d '\r\n' < "${IDENTITY_PLAN}/READY")" == "${TWIRL_FM0_MATCHED_CANARY_EVALUATION_PLAN_RECEIPT_SHA256}" ]] || die "identity-plan READY drift"

mkdir -p "$(dirname "${OUTPUT}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export LD_LIBRARY_PATH="${FM0_ENV_LIB}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"${PYTHON}" scripts/stage5_validation/freeze_twirl_fm0_3_matched_canary_payload_plan.py \
  --identity-plan-root "${IDENTITY_PLAN}" \
  --identity-plan-receipt-sha256 "${TWIRL_FM0_MATCHED_CANARY_EVALUATION_PLAN_RECEIPT_SHA256}" \
  --producer-git-sha "${TWIRL_EXPECTED_GIT_SHA}" \
  --output "${OUTPUT}"

"${PYTHON}" - "${OUTPUT}" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

from twirl.models.fm0.matched_canary_payload_plan import (
    validate_matched_canary_payload_plan,
)

root = Path(sys.argv[1])
receipt_sha = hashlib.sha256((root / "receipt.json").read_bytes()).hexdigest()
result = validate_matched_canary_payload_plan(
    root, expected_receipt_sha256=receipt_sha
)
assert result.receipt["schedule"]["n_rows"] == 1440
assert result.receipt["limits"]["bound_source_payloads_opened"] is True
assert result.receipt["limits"]["events_injected"] is False
assert result.receipt["limits"]["checkpoints_loaded"] is False
assert result.receipt["limits"]["metrics_computed"] is False
assert result.receipt["limits"]["model_training_executed"] is False
print(json.dumps({
    "receipt_sha256": receipt_sha,
    "schedule_sha256": result.schedule_sha256,
    "n_schedule_rows": 1440,
}, sort_keys=True))
PY

echo "[fm0.3-payload-plan] complete=${OUTPUT}; payload-screened only; no checkpoint, evaluation, or training"
