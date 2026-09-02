#!/bin/bash
# Inspect or explicitly submit one frozen FM0.3 matched-canary segment.
# Submission is a dry run unless --execute is present.

set -euo pipefail

readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly CONTROL_PATH="${ORCD_CONTROL_PATH:-/Users/tehan/.ssh/cm/tehan@orcd-login.mit.edu:22}"
readonly ORCD_HOST="${ORCD_HOST:-tehan@orcd-login.mit.edu}"
readonly COMPOSITE_ROOT="${PROJECT_ROOT}/exports/fm0/input_releases/twirl_fm0_3_s56_s64_s66_s77_composite_v1"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
readonly FM0_ENV_LIB="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/lib"
readonly DESIGN_SHA="6ca42278e83ccd0c28eb37fe6e398e23a08e8f7f19949742c86cfe02bb17cd10"
readonly CONFIG_SHA="3a307fcb2c430a5ec7f56c381e98017d05b957c368863fde560f0938611cd337"
readonly FREEZE_SHA="340425b22fb345310cedabace02ced0e2a6f8a356428200f7764a9e22de609c8"
readonly COMPOSITE_RECEIPT_SHA="cc5cc3bce4c24e74bef1fbf084f407855233de7893183e6bc3486e284a2f44d9"
readonly SOURCE_BINDINGS_SHA="8cbcac99409ab89fe2dee0c36687f50122e78376206495073698489ca5424f2b"
readonly ROLE_INDEX_SHA="abe9c616523f2486bf1b7be69dfcfda6193d534b40b3f4afc49d8ebf3e40ce5a"
readonly PANEL_RECEIPT_SHA="78c370e10c556472c5997c20cfe95207a0b334bafe7f024bf7ba4fc7ec4de624"

usage() {
  cat <<'EOF'
Usage:
  run_twirl_fm0_3_matched_canary_orcd.sh [--run status|preflight|submit]
      [--variant TWIRL-FM0.3.1|TWIRL-FM0.3.2]
      [--stop-after-step 8|64|2000] [--execute]

The default action is read-only status. A submit action only prints the
planned sbatch command unless --execute is supplied.
EOF
}

RUN="status"
VARIANT=""
STOP_AFTER=""
EXECUTE=0
while (( $# )); do
  case "$1" in
    --run) RUN="${2:?missing --run value}"; shift 2 ;;
    --variant) VARIANT="${2:?missing --variant value}"; shift 2 ;;
    --stop-after-step) STOP_AFTER="${2:?missing stop value}"; shift 2 ;;
    --execute) EXECUTE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

readonly LOCAL_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA:-$(git -C "${LOCAL_ROOT}" rev-parse HEAD)}"
[[ "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]] || { echo "Expected Git SHA must be 40 lowercase hex characters." >&2; exit 2; }
readonly SHA12="${EXPECTED_SHA:0:12}"
readonly REMOTE_REPO="${TWIRL_ORCD_REPO:-${PROJECT_ROOT}/code/TWIRL_fm0_3_${SHA12}}"
readonly RUN_ROOT="${TWIRL_FM0_RUN_ROOT:-${PROJECT_ROOT}/reports/stage5_validation/twirl_fm0_3_matched_canary_v1/${SHA12}}"
readonly EVALUATION_PLAN="${RUN_ROOT}/evaluation/payload_plan"
[[ "${REMOTE_REPO}" == "${PROJECT_ROOT}"/* && "${RUN_ROOT}" == "${PROJECT_ROOT}"/* ]] || {
  echo "Remote paths must stay under ${PROJECT_ROOT}." >&2
  exit 2
}

SSH=(
  ssh
  -S "${CONTROL_PATH}"
  -o ControlMaster=no
  -o BatchMode=yes
  -o PasswordAuthentication=no
  -o KbdInteractiveAuthentication=no
  -o PubkeyAuthentication=no
  -o HostbasedAuthentication=no
  -o GSSAPIAuthentication=no
  -o NumberOfPasswordPrompts=0
  -o ConnectTimeout=15
  -o ConnectionAttempts=1
  "${ORCD_HOST}"
)

remote() { "${SSH[@]}" "$@"; }

preflight() {
  remote "set -euo pipefail
test -e '${REMOTE_REPO}/.git'
test \"\$(git -C '${REMOTE_REPO}' rev-parse HEAD)\" = '${EXPECTED_SHA}'
test -z \"\$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all)\"
test \"\$(sha256sum '${REMOTE_REPO}/doc/fm0_3_matched_canary_design.md' | awk '{print \$1}')\" = '${DESIGN_SHA}'
test \"\$(sha256sum '${REMOTE_REPO}/configs/models/twirl_fm0_3_s56_s64_s66_s77_matched_canary_v1.yaml' | awk '{print \$1}')\" = '${CONFIG_SHA}'
test \"\$(sha256sum '${REMOTE_REPO}/reports/stage5_validation/twirl_fm0_3_matched_canary_design_freeze_v1/freeze.json' | awk '{print \$1}')\" = '${FREEZE_SHA}'
test \"\$(sha256sum '${COMPOSITE_ROOT}/receipt.json' | awk '{print \$1}')\" = '${COMPOSITE_RECEIPT_SHA}'
test \"\$(sha256sum '${COMPOSITE_ROOT}/source_bindings.csv' | awk '{print \$1}')\" = '${SOURCE_BINDINGS_SHA}'
test \"\$(sha256sum '${COMPOSITE_ROOT}/role_index.csv' | awk '{print \$1}')\" = '${ROLE_INDEX_SHA}'
test \"\$(tr -d '\\r\\n' < '${COMPOSITE_ROOT}/READY')\" = '${COMPOSITE_RECEIPT_SHA}'
echo 'FM0.3 matched-canary preflight passed: ${EXPECTED_SHA}'"
}

validate_prior_release() {
  local kind="$1"
  local root="$2"
  local expected_step="$3"
  remote "set -euo pipefail
cd '${REMOTE_REPO}'
export PYTHONPATH='${REMOTE_REPO}/src'
export LD_LIBRARY_PATH='${FM0_ENV_LIB}'
'${PYTHON}' - '${kind}' '${root}' '${VARIANT}' '${EXPECTED_SHA}' '${expected_step}' '${DESIGN_SHA}' '${CONFIG_SHA}' '${FREEZE_SHA}'" <<'PY'
import sys
from pathlib import Path

from twirl.models.fm0.validation import (
    read_json,
    validate_real_run_release,
    validate_run_release,
)

(
    kind,
    root_arg,
    variant,
    git_sha,
    expected_step_arg,
    design_sha,
    config_sha,
    freeze_sha,
) = sys.argv[1:]
root = Path(root_arg)
expected_step = int(expected_step_arg)
expected_authorities = {
    "design_sha256": design_sha,
    "config_sha256": config_sha,
    "freeze_receipt_sha256": freeze_sha,
}

if kind == "synthetic":
    validation = validate_run_release(root, inspect_checkpoint=True)
    expected_precision = "fp32"
    expected_synthetic = True
elif kind == "real":
    validation = validate_real_run_release(root, inspect_checkpoint=True)
    expected_precision = "bf16"
    expected_synthetic = False
else:
    raise SystemExit(f"unsupported prior-release kind: {kind}")
expected_real = not expected_synthetic

summary = read_json(root / "summary.json")
contract = read_json(root / "run_contract.json")
if (
    validation.get("passed") is not True
    or validation.get("checkpoint_inspected") is not True
    or validation.get("variant") != variant
    or validation.get("global_step") != expected_step
    or summary.get("requested_stop_after_step") != expected_step
    or summary.get("precision") != expected_precision
    or contract.get("variant") != variant
    or contract.get("expected_git_sha") != git_sha
    or contract.get("authorities") != expected_authorities
    or contract.get("synthetic_smoke") is not expected_synthetic
    or contract.get("real_data_consumed") is not expected_real
    or contract.get("training_horizon_step") != 20000
):
    raise SystemExit(f"{kind} prerequisite differs from the matched canary")
PY
}

validate_evaluation_plan() {
  local receipt_sha="$1"
  [[ "${receipt_sha}" =~ ^[0-9a-f]{64}$ ]] || {
    echo "Evaluation receipt SHA-256 is malformed." >&2
    return 2
  }
  remote "set -euo pipefail
cd '${REMOTE_REPO}'
export PYTHONPATH='${REMOTE_REPO}/src'
export LD_LIBRARY_PATH='${FM0_ENV_LIB}'
test -s '${EVALUATION_PLAN}/receipt.json'
test -s '${EVALUATION_PLAN}/schedule.csv'
test -s '${EVALUATION_PLAN}/READY'
test \"\$(sha256sum '${EVALUATION_PLAN}/receipt.json' | awk '{print \$1}')\" = '${receipt_sha}'
test \"\$(tr -d '\\r\\n' < '${EVALUATION_PLAN}/READY')\" = '${receipt_sha}'
'${PYTHON}' - '${EVALUATION_PLAN}' '${receipt_sha}' '${EXPECTED_SHA}' '${PANEL_RECEIPT_SHA}' '${COMPOSITE_RECEIPT_SHA}' '${SOURCE_BINDINGS_SHA}' '${ROLE_INDEX_SHA}'" <<'PY'
import sys

from twirl.models.fm0.matched_canary_payload_plan import (
    validate_matched_canary_payload_plan,
)

(
    root,
    receipt_sha,
    git_sha,
    panel_receipt_sha,
    composite_sha,
    source_sha,
    role_sha,
) = sys.argv[1:]
result = validate_matched_canary_payload_plan(
    root, expected_receipt_sha256=receipt_sha, require_read_only=True
)
authorities = result.receipt["identity_plan"]["input_authorities"]
temporal = authorities["temporal_panel"]
composite = authorities["composite_release"]
if (
    result.receipt["producer_git_sha"] != git_sha
    or temporal["receipt_sha256"] != panel_receipt_sha
    or composite["receipt_sha256"] != composite_sha
    or composite["source_bindings_sha256"] != source_sha
    or composite["role_index_sha256"] != role_sha
):
    raise SystemExit("evaluation plan differs from this exact training campaign")
PY
}

pin_evaluation_plan() {
  local receipt_sha="$1"
  [[ "${receipt_sha}" =~ ^[0-9a-f]{64}$ ]] || {
    echo "Evaluation receipt SHA-256 is malformed." >&2
    return 2
  }
  remote "set -euo pipefail
mkdir -p '${RUN_ROOT}/launch'
pin='${RUN_ROOT}/launch/evaluation_plan.sha256'
if (set -o noclobber; printf '%s\n' '${receipt_sha}' > \"\${pin}\") 2>/dev/null; then
  chmod a-w \"\${pin}\"
fi
test -f \"\${pin}\"
test ! -L \"\${pin}\"
test ! -w \"\${pin}\"
test \"\$(tr -d '\\r\\n' < \"\${pin}\")\" = '${receipt_sha}'"
}

case "${RUN}" in
  status)
    remote "set -euo pipefail
echo 'run_root=${RUN_ROOT}'
find '${RUN_ROOT}/launch' -maxdepth 1 -type f -name '*.job' -print -exec sed -n '1p' {} \\; 2>/dev/null | sort || true
squeue -u \$USER -h -o '%i|%T|%M|%l|%j|%N' | grep 'twirl-fm0-3' || true"
    ;;
  preflight)
    preflight
    ;;
  submit)
    case "${VARIANT}" in
      TWIRL-FM0.3.1) VARIANT_SLUG="fm0_3_1_tcn"; JOB_VARIANT="031" ;;
      TWIRL-FM0.3.2) VARIANT_SLUG="fm0_3_2_conformer"; JOB_VARIANT="032" ;;
      *) echo "--variant must select FM0.3.1 or FM0.3.2." >&2; exit 2 ;;
    esac
    case "${STOP_AFTER}" in
      8) OUTPUT="${RUN_ROOT}/smokes/${VARIANT_SLUG}/seed560067" ;;
      64|2000) OUTPUT="${RUN_ROOT}/model_runs/${VARIANT_SLUG}/seed560067" ;;
      *) echo "--stop-after-step must be 8, 64, or 2000." >&2; exit 2 ;;
    esac
    WRAPPER="${REMOTE_REPO}/scripts/orcd/slurm_twirl_fm0_3_matched_canary_h200.sbatch"
    JOB_NAME="twirl-fm0-${JOB_VARIANT}-s${STOP_AFTER}"
    EXPORTS="ALL,TWIRL_ORCD_REPO=${REMOTE_REPO},TWIRL_EXPECTED_GIT_SHA=${EXPECTED_SHA},TWIRL_FM0_RUN_ROOT=${RUN_ROOT},TWIRL_FM0_VARIANT=${VARIANT},TWIRL_FM0_STOP_AFTER_STEP=${STOP_AFTER},TWIRL_FM0_OUTPUT=${OUTPUT}"
    if [[ "${STOP_AFTER}" == 2000 ]]; then
      EXPORTS+=",TWIRL_FM0_RESUME_CHECKPOINT=${OUTPUT}/checkpoint_step_00000064.pt"
    fi
    if (( ! EXECUTE )); then
      if [[ "${STOP_AFTER}" == 64 || "${STOP_AFTER}" == 2000 ]]; then
        EXPORTS+=",TWIRL_FM0_EVALUATION_PLAN_RECEIPT_SHA256=bound-at-execution"
      fi
      echo "DRY RUN: sbatch --parsable --job-name=${JOB_NAME} --export=${EXPORTS} ${WRAPPER}"
      echo "No job submitted; add --execute to submit."
      exit 0
    fi
    preflight
    EVALUATION_RECEIPT_SHA="not_applicable"
    if [[ "${STOP_AFTER}" == 64 || "${STOP_AFTER}" == 2000 ]]; then
      EVALUATION_RECEIPT_SHA="$(remote "set -euo pipefail
test -s '${EVALUATION_PLAN}/receipt.json'
sha256sum '${EVALUATION_PLAN}/receipt.json' | awk '{print \$1}'")"
      validate_evaluation_plan "${EVALUATION_RECEIPT_SHA}"
      EXPORTS+=",TWIRL_FM0_EVALUATION_PLAN_RECEIPT_SHA256=${EVALUATION_RECEIPT_SHA}"
      validate_prior_release synthetic \
        "${RUN_ROOT}/smokes/${VARIANT_SLUG}/seed560067" 8
    fi
    if [[ "${STOP_AFTER}" == 2000 ]]; then
      validate_prior_release real "${OUTPUT}" 64
    fi
    if [[ "${STOP_AFTER}" == 64 || "${STOP_AFTER}" == 2000 ]]; then
      pin_evaluation_plan "${EVALUATION_RECEIPT_SHA}"
    fi
    remote "set -euo pipefail
mkdir -p '${RUN_ROOT}/launch'
record='${RUN_ROOT}/launch/${VARIANT_SLUG}-step${STOP_AFTER}.job'
test ! -e \"\${record}\"
job=\$(sbatch --parsable --job-name='${JOB_NAME}' --export='${EXPORTS}' '${WRAPPER}')
printf '%s\\t%s\\t%s\\t%s\\t%s\\n' \"\${job}\" '${VARIANT}' '${STOP_AFTER}' '${EXPECTED_SHA}' '${EVALUATION_RECEIPT_SHA}' > \"\${record}\"
chmod a-w \"\${record}\"
printf '%s\\n' \"\${job}\""
    ;;
  *) echo "--run must be status, preflight, or submit." >&2; exit 2 ;;
esac
