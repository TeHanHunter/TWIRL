#!/bin/bash
# Freeze only the identity schedule for the matched FM0.3 evaluation.
# No shard payload, checkpoint, evaluator, trainer, or GPU is used.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-3-eval-plan
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 027

readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
readonly FM0_ENV_LIB="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/lib"
readonly COMPOSITE_ROOT="${PROJECT_ROOT}/exports/fm0/input_releases/twirl_fm0_3_s56_s64_s66_s77_composite_v1"
readonly PANEL_RECEIPT_SHA="78c370e10c556472c5997c20cfe95207a0b334bafe7f024bf7ba4fc7ec4de624"
readonly COMPOSITE_RECEIPT_SHA="cc5cc3bce4c24e74bef1fbf084f407855233de7893183e6bc3486e284a2f44d9"
readonly SOURCE_BINDINGS_SHA="8cbcac99409ab89fe2dee0c36687f50122e78376206495073698489ca5424f2b"
readonly ROLE_INDEX_SHA="abe9c616523f2486bf1b7be69dfcfda6193d534b40b3f4afc49d8ebf3e40ce5a"

: "${TWIRL_ORCD_REPO:?missing exact detached repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?missing exact Git SHA}"
: "${TWIRL_FM0_TEMPORAL_PANEL_DIR:?missing frozen temporal-panel root}"
: "${TWIRL_FM0_MATCHED_CANARY_EVALUATION_PLAN_DIR:?missing evaluation-plan output root}"

readonly REPO="$(realpath -e -- "${TWIRL_ORCD_REPO}")"
readonly PANEL="$(realpath -e -- "${TWIRL_FM0_TEMPORAL_PANEL_DIR}")"
readonly OUTPUT="$(realpath -m -- "${TWIRL_FM0_MATCHED_CANARY_EVALUATION_PLAN_DIR}")"

die() { echo "[fm0.3-eval-plan] ERROR: $*" >&2; exit 2; }
digest() { sha256sum "$1" | awk '{print $1}'; }

[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || die "invalid Git SHA"
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" ]] || die "Git revision drift"
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is dirty"
[[ "${PANEL}" == "${PROJECT_ROOT}/reports/stage5_validation/"* ]] || die "panel root escaped the report tree"
[[ "${OUTPUT}" == "${PROJECT_ROOT}/reports/stage5_validation/"* ]] || die "output root escaped the report tree"
[[ ! -e "${OUTPUT}" ]] || die "refusing to overwrite evaluation plan"

[[ "$(digest "${PANEL}/receipt.json")" == "${PANEL_RECEIPT_SHA}" ]] || die "temporal-panel receipt drift"
[[ "$(tr -d '\r\n' < "${PANEL}/READY")" == "${PANEL_RECEIPT_SHA}" ]] || die "temporal-panel READY drift"
for binding in \
  "receipt.json:${COMPOSITE_RECEIPT_SHA}" \
  "source_bindings.csv:${SOURCE_BINDINGS_SHA}" \
  "role_index.csv:${ROLE_INDEX_SHA}"; do
  name="${binding%%:*}"
  expected="${binding#*:}"
  [[ "$(digest "${COMPOSITE_ROOT}/${name}")" == "${expected}" ]] || die "composite ${name} drift"
done
[[ "$(tr -d '\r\n' < "${COMPOSITE_ROOT}/READY")" == "${COMPOSITE_RECEIPT_SHA}" ]] || die "composite READY drift"

mkdir -p "$(dirname "${OUTPUT}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export LD_LIBRARY_PATH="${FM0_ENV_LIB}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"${PYTHON}" scripts/stage5_validation/freeze_twirl_fm0_3_matched_canary_evaluation_plan.py \
  --temporal-panel-root "${PANEL}" \
  --temporal-panel-receipt-sha256 "${PANEL_RECEIPT_SHA}" \
  --composite-root "${COMPOSITE_ROOT}" \
  --composite-receipt-sha256 "${COMPOSITE_RECEIPT_SHA}" \
  --composite-source-bindings-sha256 "${SOURCE_BINDINGS_SHA}" \
  --composite-role-index-sha256 "${ROLE_INDEX_SHA}" \
  --producer-git-sha "${TWIRL_EXPECTED_GIT_SHA}" \
  --output "${OUTPUT}"

"${PYTHON}" - "${OUTPUT}" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

from twirl.models.fm0.matched_canary_plan import (
    SPLIT_ORDER,
    TEMPORAL_COHORTS,
    validate_matched_canary_evaluation_plan,
)

root = Path(sys.argv[1])
receipt_sha = hashlib.sha256((root / "receipt.json").read_bytes()).hexdigest()
result = validate_matched_canary_evaluation_plan(
    root, expected_receipt_sha256=receipt_sha
)
counts = result.receipt["split_contract"]["counts"]
assert result.receipt["schedule"]["n_rows"] == 1440
assert all(
    counts[split][cohort] == 240
    for split in SPLIT_ORDER
    for cohort in TEMPORAL_COHORTS
)
assert (root / "READY").read_text().strip() == receipt_sha
print(json.dumps({
    "receipt_sha256": receipt_sha,
    "schedule_sha256": result.schedule_sha256,
    "n_schedule_rows": 1440,
}, sort_keys=True))
PY

echo "[fm0.3-eval-plan] complete=${OUTPUT}; identity-only; no shard payload, evaluation, training, or GPU"
