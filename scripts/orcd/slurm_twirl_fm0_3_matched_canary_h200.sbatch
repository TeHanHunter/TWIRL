#!/bin/bash
# One-H200 matched FM0.3 segment. Stop 8 is synthetic FP32; stops 64 and
# 2000 are real BF16, with 2000 resuming the same variant's step-64 milestone.
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-fm0-3-canary
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 47:30:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 027

readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly COMPOSITE_ROOT="${PROJECT_ROOT}/exports/fm0/input_releases/twirl_fm0_3_s56_s64_s66_s77_composite_v1"
readonly DESIGN_SHA="6ca42278e83ccd0c28eb37fe6e398e23a08e8f7f19949742c86cfe02bb17cd10"
readonly CONFIG_SHA="3a307fcb2c430a5ec7f56c381e98017d05b957c368863fde560f0938611cd337"
readonly FREEZE_SHA="340425b22fb345310cedabace02ced0e2a6f8a356428200f7764a9e22de609c8"
readonly COMPOSITE_RECEIPT_SHA="cc5cc3bce4c24e74bef1fbf084f407855233de7893183e6bc3486e284a2f44d9"
readonly SOURCE_BINDINGS_SHA="8cbcac99409ab89fe2dee0c36687f50122e78376206495073698489ca5424f2b"
readonly ROLE_INDEX_SHA="abe9c616523f2486bf1b7be69dfcfda6193d534b40b3f4afc49d8ebf3e40ce5a"
readonly PANEL_RECEIPT_SHA="78c370e10c556472c5997c20cfe95207a0b334bafe7f024bf7ba4fc7ec4de624"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
readonly FM0_ENV_LIB="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/lib"

: "${TWIRL_ORCD_REPO:?missing exact detached repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?missing exact Git SHA}"
: "${TWIRL_FM0_RUN_ROOT:?missing run root}"
: "${TWIRL_FM0_VARIANT:?missing FM0.3 variant}"
: "${TWIRL_FM0_STOP_AFTER_STEP:?missing stop-after step}"
: "${TWIRL_FM0_OUTPUT:?missing output path}"

readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly RUN_ROOT="$(realpath -m -- "${TWIRL_FM0_RUN_ROOT}")"
readonly OUT="$(realpath -m -- "${TWIRL_FM0_OUTPUT}")"
readonly VARIANT="${TWIRL_FM0_VARIANT}"
readonly STOP_AFTER="${TWIRL_FM0_STOP_AFTER_STEP}"
readonly DESIGN="${REPO}/doc/fm0_3_matched_canary_design.md"
readonly CONFIG="${REPO}/configs/models/twirl_fm0_3_s56_s64_s66_s77_matched_canary_v1.yaml"
readonly FREEZE="${REPO}/reports/stage5_validation/twirl_fm0_3_matched_canary_design_freeze_v1/freeze.json"

die() { echo "[fm0.3-canary] ERROR: $*" >&2; exit 2; }
digest() { sha256sum "$1" | awk '{print $1}'; }

[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || die "invalid Git SHA"
[[ "${RUN_ROOT}" == "${PROJECT_ROOT}"/* && "${OUT}" == "${RUN_ROOT}"/* ]] || die "unsafe run path"
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" ]] || die "Git revision drift"
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is dirty"
[[ "$(digest "${DESIGN}")" == "${DESIGN_SHA}" ]] || die "design drift"
[[ "$(digest "${CONFIG}")" == "${CONFIG_SHA}" ]] || die "config drift"
[[ "$(digest "${FREEZE}")" == "${FREEZE_SHA}" ]] || die "freeze receipt drift"
export PYTHONPATH="${REPO}/src"
export LD_LIBRARY_PATH="${FM0_ENV_LIB}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

case "${VARIANT}" in
  TWIRL-FM0.3.1) readonly VARIANT_SLUG="fm0_3_1_tcn" ;;
  TWIRL-FM0.3.2) readonly VARIANT_SLUG="fm0_3_2_conformer" ;;
  *) die "variant must be TWIRL-FM0.3.1 or TWIRL-FM0.3.2" ;;
esac

readonly SMOKE_OUT="${RUN_ROOT}/smokes/${VARIANT_SLUG}/seed560067"
readonly REAL_OUT="${RUN_ROOT}/model_runs/${VARIANT_SLUG}/seed560067"
readonly EVALUATION_PLAN="${RUN_ROOT}/evaluation/payload_plan"
readonly EVALUATION_PLAN_PIN="${RUN_ROOT}/launch/evaluation_plan.sha256"

validate_prior_release() {
  local kind="$1"
  local root="$2"
  local expected_step="$3"
  "${PYTHON}" - "${kind}" "${root}" "${VARIANT}" "${TWIRL_EXPECTED_GIT_SHA}" "${expected_step}" "${DESIGN_SHA}" "${CONFIG_SHA}" "${FREEZE_SHA}" <<'PY'
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
  [[ -s "${EVALUATION_PLAN}/receipt.json" && -s "${EVALUATION_PLAN}/schedule.csv" && -s "${EVALUATION_PLAN}/READY" ]] || die "payload-screened evaluation plan is incomplete"
  [[ "${EVALUATION_RECEIPT_SHA}" =~ ^[0-9a-f]{64}$ ]] || die "missing externally pinned evaluation-plan receipt SHA-256"
  [[ -f "${EVALUATION_PLAN_PIN}" && ! -L "${EVALUATION_PLAN_PIN}" && ! -w "${EVALUATION_PLAN_PIN}" ]] || die "campaign evaluation-plan pin is absent, linked, or writable"
  [[ "$(tr -d '\r\n' < "${EVALUATION_PLAN_PIN}")" == "${EVALUATION_RECEIPT_SHA}" ]] || die "evaluation-plan receipt differs from campaign pin"
  [[ "$(digest "${EVALUATION_PLAN}/receipt.json")" == "${EVALUATION_RECEIPT_SHA}" ]] || die "evaluation-plan receipt differs from launch authority"
  [[ "$(tr -d '\r\n' < "${EVALUATION_PLAN}/READY")" == "${EVALUATION_RECEIPT_SHA}" ]] || die "evaluation-plan READY binding differs"
  "${PYTHON}" - "${EVALUATION_PLAN}" "${EVALUATION_RECEIPT_SHA}" "${TWIRL_EXPECTED_GIT_SHA}" "${PANEL_RECEIPT_SHA}" "${COMPOSITE_RECEIPT_SHA}" "${SOURCE_BINDINGS_SHA}" "${ROLE_INDEX_SHA}" <<'PY'
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

TRAIN_ARGS=()
EVALUATION_RECEIPT_SHA=""
case "${STOP_AFTER}" in
  8)
    [[ "${OUT}" == "${SMOKE_OUT}" ]] || die "step-8 output path differs"
    [[ ! -e "${OUT}" ]] || die "refusing to overwrite smoke output"
    readonly PRECISION="fp32"
    TRAIN_ARGS=(--synthetic-smoke)
    ;;
  64|2000)
    [[ "${OUT}" == "${REAL_OUT}" ]] || die "real output path differs"
    readonly PRECISION="bf16"
    : "${TWIRL_FM0_EVALUATION_PLAN_RECEIPT_SHA256:?real training requires the launch-pinned evaluation receipt}"
    EVALUATION_RECEIPT_SHA="${TWIRL_FM0_EVALUATION_PLAN_RECEIPT_SHA256}"
    validate_prior_release synthetic "${SMOKE_OUT}" 8
    validate_evaluation_plan
    TRAIN_ARGS=(
      --input-release "${COMPOSITE_ROOT}"
      --input-release-receipt "${COMPOSITE_ROOT}/receipt.json"
      --evaluation-plan "${EVALUATION_PLAN}"
      --evaluation-plan-receipt-sha256 "${EVALUATION_RECEIPT_SHA}"
      --prestart-smoke-run "${SMOKE_OUT}"
      --prestart-smoke-summary-sha256 "$(digest "${SMOKE_OUT}/summary.json")"
    )
    ;;
  *) die "stop-after step must be 8, 64, or 2000" ;;
esac

if [[ "${STOP_AFTER}" == 64 ]]; then
  [[ ! -e "${OUT}" ]] || die "step-64 invocation requires a fresh output"
elif [[ "${STOP_AFTER}" == 2000 ]]; then
  validate_prior_release real "${OUT}" 64
  : "${TWIRL_FM0_RESUME_CHECKPOINT:?step 2000 requires the step-64 checkpoint}"
  readonly RESUME="$(realpath -e -- "${TWIRL_FM0_RESUME_CHECKPOINT}")"
  [[ "${RESUME}" == "${OUT}/checkpoint_step_00000064.pt" ]] || die "resume is not this variant's step-64 milestone"
  [[ -s "${RESUME}.sha256" ]] || die "step-64 milestone sidecar is absent"
  (cd "${OUT}" && sha256sum -c "$(basename "${RESUME}").sha256") >/dev/null || die "step-64 checkpoint hash differs"
  TRAIN_ARGS+=(--resume-checkpoint "${RESUME}")
fi

for binding in \
  "receipt.json:${COMPOSITE_RECEIPT_SHA}" \
  "source_bindings.csv:${SOURCE_BINDINGS_SHA}" \
  "role_index.csv:${ROLE_INDEX_SHA}"; do
  name="${binding%%:*}"
  expected="${binding#*:}"
  path="${COMPOSITE_ROOT}/${name}"
  [[ -f "${path}" && ! -L "${path}" && ! -w "${path}" ]] || die "composite authority is absent, linked, or writable: ${name}"
  [[ "$(digest "${path}")" == "${expected}" ]] || die "composite authority hash differs: ${name}"
done
[[ "$(tr -d '\r\n' < "${COMPOSITE_ROOT}/READY")" == "${COMPOSITE_RECEIPT_SHA}" ]] || die "composite READY binding differs"

mkdir -p "$(dirname "${OUT}")"
cd "${REPO}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
module load cuda/13.1.0 || module load cuda/12.9.1 || true
nvidia-smi
"${PYTHON}" -c 'import torch; assert torch.cuda.is_available(); assert torch.cuda.device_count() == 1; print(torch.__version__, torch.cuda.get_device_name(0))'

"${PYTHON}" scripts/stage5_validation/train_twirl_fm0.py \
  --config "${CONFIG}" \
  --design "${DESIGN}" \
  --freeze-receipt "${FREEZE}" \
  --variant "${VARIANT}" \
  --output-dir "${OUT}" \
  --device cuda \
  --precision "${PRECISION}" \
  --stop-after-step "${STOP_AFTER}" \
  --micro-batch-windows 8 \
  --seed 560067 \
  --expected-git-sha "${TWIRL_EXPECTED_GIT_SHA}" \
  "${TRAIN_ARGS[@]}"

for path in run_contract.json checkpoint.pt summary.json; do
  [[ -s "${OUT}/${path}" && -s "${OUT}/${path}.sha256" ]] || die "missing artifact/sidecar: ${path}"
done
for step in 0 64 2000; do
  (( step <= STOP_AFTER )) || continue
  milestone="${OUT}/checkpoint_step_$(printf '%08d' "${step}").pt"
  [[ -s "${milestone}" && -s "${milestone}.sha256" ]] || die "missing milestone/sidecar: ${step}"
  chmod a-w "${milestone}" "${milestone}.sha256"
done
if [[ "${STOP_AFTER}" == 8 ]]; then
  chmod -R a-w "${OUT}"
fi
echo "[fm0.3-canary] complete=${OUT} variant=${VARIANT} stop=${STOP_AFTER} precision=${PRECISION}"
