#!/bin/bash
# One-H200 BF16 execution segment for the frozen FM0.2.1 real-data canary.
# Each invocation stops at one authorized milestone and resumes only from the
# immediately preceding immutable checkpoint after independent validation.
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-fm0-2-1-canary
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
readonly CONFIG_SHA="394a0340f9f8518423e39a1c83bc32acde4e2ebff95c4837b0edb5399bd50801"
readonly FREEZE_SHA="900264edabe13aa77dfffe3c87f5fe0c4967cbcf693b9a68020ab29ff85e691e"
readonly REUSE_SHA="46faed954e6758dfea71f7754670611fab03af63613dad3be433af71dd38e641"
readonly EVENT_SHA="eb84f27866e42bbf942256abe8fa7a532eaa8584e6800ce012aa0552c08a6630"
readonly EVALUATOR_SHA="ed8004180b28c7e2bd85911d0128f9180eddd3c35bdaae299e2f7edba6e5baac"
readonly RELEASE_ROOT="${PROJECT_ROOT}/exports/fm0/input_releases/twirl_fm0_1_s56_s64_v2"
: "${TWIRL_ORCD_REPO:?missing detached repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?missing exact clean SHA}"
: "${TWIRL_FM0_RUN_ROOT:?missing run root}"
: "${TWIRL_FM0_INPUT_RECEIPT:?missing exact upstream FM0.1 input receipt}"
: "${TWIRL_FM0_INPUT_RECEIPT_SHA256:?missing upstream input receipt hash}"
: "${TWIRL_FM0_INPUT_REUSE_RECEIPT:?missing frozen FM0.2 input-reuse receipt}"
: "${TWIRL_FM0_INPUT_REUSE_RECEIPT_SHA256:?missing input-reuse receipt hash}"
: "${TWIRL_FM0_LOADER_RESTART_RECEIPT:?missing exact-code loader/restart receipt}"
: "${TWIRL_FM0_LOADER_RESTART_RECEIPT_SHA256:?missing loader/restart receipt hash}"
: "${TWIRL_FM0_SMOKE_POST_VALIDATION:?missing independently validated eight-step smoke receipt}"
: "${TWIRL_FM0_SMOKE_POST_VALIDATION_SHA256:?missing smoke post-validation hash}"
: "${TWIRL_FM0_ADMISSION_RECEIPT:?missing fresh Stage-1-safe admission receipt}"
: "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256:?missing admission receipt hash}"
: "${TWIRL_FM0_CANARY_OUTPUT:?missing canary output path}"
: "${TWIRL_FM0_STOP_AFTER_STEP:?missing authorized stop-after step}"
readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly RUN_ROOT="$(realpath -m -- "${TWIRL_FM0_RUN_ROOT}")"
readonly INPUT_RECEIPT="$(realpath -m -- "${TWIRL_FM0_INPUT_RECEIPT}")"
readonly REUSE_RECEIPT="$(realpath -m -- "${TWIRL_FM0_INPUT_REUSE_RECEIPT}")"
readonly LOADER_RECEIPT="$(realpath -m -- "${TWIRL_FM0_LOADER_RESTART_RECEIPT}")"
readonly SMOKE_POST="$(realpath -m -- "${TWIRL_FM0_SMOKE_POST_VALIDATION}")"
readonly ADMISSION_RECEIPT="$(realpath -m -- "${TWIRL_FM0_ADMISSION_RECEIPT}")"
readonly OUT="$(realpath -m -- "${TWIRL_FM0_CANARY_OUTPUT}")"
readonly STOP_AFTER="${TWIRL_FM0_STOP_AFTER_STEP}"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
readonly CONFIG="${REPO}/configs/models/twirl_fm0_2_s56_s64_objective_canary.yaml"
readonly FREEZE="${REPO}/reports/stage5_validation/twirl_fm0_2_design_freeze_v1/freeze.json"
readonly EVENT_CONTRACT="${REPO}/reports/stage5_validation/twirl_fm0_2_design_freeze_v1/development_event_retention.contract.json"
readonly EVALUATOR_CONTRACT="${REPO}/reports/stage5_validation/twirl_fm0_2_design_freeze_v1/development_evaluator_v2.contract.json"
die() { echo "[fm0.2-canary] ERROR: $*" >&2; exit 2; }

case "${STOP_AFTER}" in
  64) readonly PREVIOUS_STEP=0 ;;
  500) readonly PREVIOUS_STEP=64 ;;
  1000) readonly PREVIOUS_STEP=500 ;;
  2000) readonly PREVIOUS_STEP=1000 ;;
  *) die "stop-after step must be one of 64, 500, 1000, or 2000" ;;
esac
[[ "${RUN_ROOT}" == "${PROJECT_ROOT}"/* && "${OUT}" == "${RUN_ROOT}"/model_runs/fm0_2_1_canary/* ]] || die "unsafe canary output"
[[ "${INPUT_RECEIPT}" == "${PROJECT_ROOT}"/* && "${REUSE_RECEIPT}" == "${REPO}"/* && "${LOADER_RECEIPT}" == "${RUN_ROOT}"/* ]] || die "unsafe input authority"
[[ "${SMOKE_POST}" == "${RUN_ROOT}"/* && "${ADMISSION_RECEIPT}" == "${RUN_ROOT}"/* ]] || die "unsafe gate authority"
for value in "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" "${TWIRL_FM0_INPUT_REUSE_RECEIPT_SHA256}" "${TWIRL_FM0_LOADER_RESTART_RECEIPT_SHA256}" "${TWIRL_FM0_SMOKE_POST_VALIDATION_SHA256}" "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}"; do
  [[ "${value}" =~ ^[0-9a-f]{64}$ ]] || die "invalid receipt hash"
done
[[ "${TWIRL_FM0_INPUT_REUSE_RECEIPT_SHA256}" == "${REUSE_SHA}" ]] || die "unexpected reuse receipt hash"
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" && -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is not exact and clean"
[[ "$(sha256sum "${CONFIG}" | awk '{print $1}')" == "${CONFIG_SHA}" ]] || die "FM0.2 config drift"
[[ "$(sha256sum "${FREEZE}" | awk '{print $1}')" == "${FREEZE_SHA}" ]] || die "FM0.2 freeze drift"
[[ "$(sha256sum "${REUSE_RECEIPT}" | awk '{print $1}')" == "${REUSE_SHA}" ]] || die "reuse receipt drift"
[[ "$(sha256sum "${EVENT_CONTRACT}" | awk '{print $1}')" == "${EVENT_SHA}" ]] || die "event contract drift"
[[ "$(sha256sum "${EVALUATOR_CONTRACT}" | awk '{print $1}')" == "${EVALUATOR_SHA}" ]] || die "evaluator contract drift"
[[ "$(sha256sum "${INPUT_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" ]] || die "upstream input receipt drift"
[[ "$(sha256sum "${LOADER_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_LOADER_RESTART_RECEIPT_SHA256}" ]] || die "loader/restart receipt drift"
[[ "$(sha256sum "${SMOKE_POST}" | awk '{print $1}')" == "${TWIRL_FM0_SMOKE_POST_VALIDATION_SHA256}" ]] || die "smoke post-validation drift"
[[ "$(sha256sum "${ADMISSION_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}" ]] || die "admission receipt drift"

readonly PREVIOUS_PADDED="$(printf '%08d' "${PREVIOUS_STEP}")"
readonly CURRENT_PADDED="$(printf '%08d' "${STOP_AFTER}")"
RESUME_ARGS=()
if [[ "${STOP_AFTER}" == 64 ]]; then
  [[ ! -e "${OUT}" ]] || die "first canary invocation requires a fresh output directory"
else
  : "${TWIRL_FM0_RESUME_CHECKPOINT:?missing immediately preceding milestone checkpoint}"
  : "${TWIRL_FM0_PREVIOUS_POST_VALIDATION:?missing previous milestone post-validation receipt}"
  : "${TWIRL_FM0_PREVIOUS_POST_VALIDATION_SHA256:?missing previous post-validation hash}"
  readonly RESUME="$(realpath -e -- "${TWIRL_FM0_RESUME_CHECKPOINT}")"
  readonly PREVIOUS_POST="$(realpath -e -- "${TWIRL_FM0_PREVIOUS_POST_VALIDATION}")"
  [[ "${RESUME}" == "${OUT}/checkpoint_step_${PREVIOUS_PADDED}.pt" ]] || die "resume checkpoint is not the immediately preceding milestone"
  [[ -s "${RESUME}.sha256" ]] || die "resume milestone sidecar is absent"
  [[ "${PREVIOUS_POST}" == "${RUN_ROOT}"/* && "${TWIRL_FM0_PREVIOUS_POST_VALIDATION_SHA256}" =~ ^[0-9a-f]{64}$ ]] || die "unsafe previous post-validation authority"
  [[ "$(sha256sum "${PREVIOUS_POST}" | awk '{print $1}')" == "${TWIRL_FM0_PREVIOUS_POST_VALIDATION_SHA256}" ]] || die "previous post-validation drift"
  RESUME_ARGS=(--resume-checkpoint "${RESUME}")
fi

"${PYTHON}" - "${REUSE_RECEIPT}" "${INPUT_RECEIPT}" "${LOADER_RECEIPT}" "${SMOKE_POST}" "${ADMISSION_RECEIPT}" "${TWIRL_EXPECTED_GIT_SHA}" "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" "${TWIRL_FM0_INPUT_REUSE_RECEIPT_SHA256}" "${TWIRL_FM0_LOADER_RESTART_RECEIPT_SHA256}" "${STOP_AFTER}" "${PREVIOUS_STEP}" "${TWIRL_FM0_PREVIOUS_POST_VALIDATION:-}" <<'PY'
import json, stat, sys, time
from pathlib import Path

reuse_path, input_path, loader_path, smoke_path, admission_path = map(Path, sys.argv[1:6])
expected_sha, upstream_sha, reuse_sha, loader_sha = sys.argv[6:10]
stop_after, previous_step = map(int, sys.argv[10:12])
previous_post_arg = sys.argv[12]
for path in (input_path, loader_path, smoke_path, admission_path):
    if path.stat().st_mode & (stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH):
        raise SystemExit(f"authority is writable: {path}")
reuse = json.loads(reuse_path.read_text())
upstream = json.loads(input_path.read_text())
loader = json.loads(loader_path.read_text())
smoke = json.loads(smoke_path.read_text())
admission = json.loads(admission_path.read_text())
assert reuse["schema_version"] == "twirl_fm0_2_input_reuse_receipt_v1"
assert reuse["passed"] is True and reuse["scientific_training_eligible"] is True
assert reuse["upstream_fm0_1_input_receipt"]["sha256"] == upstream_sha
assert reuse["sealed_test"]["access_count"] == 0
assert upstream["schema_version"] == "twirl_fm0_1_input_release_receipt_v1"
assert upstream["passed"] is True and upstream["scientific_training_eligible"] is True
assert loader["schema_version"] == "twirl_fm0_2_loader_restart_receipt_v1"
assert loader["passed"] is True and loader["expected_git_sha"] == expected_sha
assert loader["input_receipt"]["sha256"] == upstream_sha
assert loader["input_reuse_receipt"]["sha256"] == reuse_sha
assert smoke["schema_version"] == "twirl_fm0_2_milestone_post_validation_v1"
assert smoke["passed"] is True and smoke["run_kind"] == "synthetic_fp32"
assert smoke["run_git_sha"] == expected_sha and smoke["stop_after_step"] == 8
assert smoke["sealed_test_access_count"] == 0
assert admission["schema_version"] == "twirl_fm0_2_orcd_admission_receipt_v1"
assert admission["passed"] is True and admission["run_kind"] == "real_canary"
assert admission["expected_git_sha"] == expected_sha
assert admission["upstream_input_receipt"]["sha256"] == upstream_sha
assert admission["input_reuse_receipt"]["sha256"] == reuse_sha
assert admission["loader_restart_receipt"]["sha256"] == loader_sha
assert admission["stop_after_step"] == stop_after
assert admission["requested_fm_resources"] == {"h200": 1, "cpus": 4, "memory_mib": 32768}
now = int(time.time())
if now < admission["created_at_epoch"] - 30 or now > admission["expires_at_epoch"]:
    raise SystemExit("admission receipt is not fresh at job start")
if previous_step == 0:
    assert previous_post_arg == ""
else:
    previous = json.loads(Path(previous_post_arg).read_text())
    assert previous["schema_version"] == "twirl_fm0_2_milestone_post_validation_v1"
    assert previous["passed"] is True and previous["run_kind"] == "real_canary"
    assert previous["stop_after_step"] == previous_step
    assert previous["run_git_sha"] == expected_sha
    assert previous["sealed_test_access_count"] == 0
PY

cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
module load cuda/13.1.0 || module load cuda/12.9.1 || true
nvidia-smi
"${PYTHON}" -c 'import torch; assert torch.cuda.is_available(); assert torch.cuda.device_count() == 1; print(torch.__version__, torch.cuda.get_device_name(0))'
"${PYTHON}" scripts/stage5_validation/train_twirl_fm0.py \
  --config "${CONFIG}" \
  --freeze-receipt "${FREEZE}" \
  --variant TWIRL-FM0.2.1 \
  --input-release "${RELEASE_ROOT}" \
  --input-release-receipt "${INPUT_RECEIPT}" \
  --input-reuse-receipt "${REUSE_RECEIPT}" \
  --output-dir "${OUT}" \
  --device cuda \
  --precision bf16 \
  --stop-after-step "${STOP_AFTER}" \
  --micro-batch-windows 8 \
  --seed 560067 \
  --expected-git-sha "${TWIRL_EXPECTED_GIT_SHA}" \
  "${RESUME_ARGS[@]}"

for step in 0 64 500 1000 2000; do
  (( step <= STOP_AFTER )) || continue
  padded="$(printf '%08d' "${step}")"
  milestone="${OUT}/checkpoint_step_${padded}.pt"
  [[ -s "${milestone}" && -s "${milestone}.sha256" ]] || die "missing immutable milestone/sidecar: ${step}"
  chmod a-w "${milestone}" "${milestone}.sha256"
done
for path in run_contract.json checkpoint.pt summary.json; do
  [[ -s "${OUT}/${path}" && -s "${OUT}/${path}.sha256" ]] || die "missing run artifact/sidecar: ${path}"
done
readonly INVOCATION_DIR="${OUT}/invocations"
readonly INVOCATION="${INVOCATION_DIR}/step_${CURRENT_PADDED}.receipt.json"
[[ ! -e "${INVOCATION}" ]] || die "refusing to overwrite invocation receipt"
mkdir -p "${INVOCATION_DIR}"
TWIRL_FM0_INVOCATION_OUT="${INVOCATION}" TWIRL_FM0_CANARY_OUT="${OUT}" \
TWIRL_FM0_STOP_AFTER="${STOP_AFTER}" TWIRL_FM0_RUN_SHA="${TWIRL_EXPECTED_GIT_SHA}" \
TWIRL_FM0_INPUT_SHA="${TWIRL_FM0_INPUT_RECEIPT_SHA256}" \
TWIRL_FM0_REUSE_SHA="${TWIRL_FM0_INPUT_REUSE_RECEIPT_SHA256}" \
TWIRL_FM0_ADMISSION_SHA="${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}" \
"${PYTHON}" - <<'PY'
from datetime import datetime, timezone
import hashlib, json, os
from pathlib import Path

root = Path(os.environ["TWIRL_FM0_CANARY_OUT"])
step = int(os.environ["TWIRL_FM0_STOP_AFTER"])
milestone = root / f"checkpoint_step_{step:08d}.pt"
payload = {
    "schema_version": "twirl_fm0_2_orcd_canary_invocation_v1",
    "created_at_utc": datetime.now(timezone.utc).isoformat(),
    "passed": True,
    "run_kind": "real_canary",
    "variant": "TWIRL-FM0.2.1",
    "seed": 560067,
    "run_git_sha": os.environ["TWIRL_FM0_RUN_SHA"],
    "run_contract_target_step": 20000,
    "stop_after_step": step,
    "input_receipt_sha256": os.environ["TWIRL_FM0_INPUT_SHA"],
    "input_reuse_receipt_sha256": os.environ["TWIRL_FM0_REUSE_SHA"],
    "admission_receipt_sha256": os.environ["TWIRL_FM0_ADMISSION_SHA"],
    "milestone_checkpoint_sha256": hashlib.sha256(milestone.read_bytes()).hexdigest(),
    "run_contract_sha256": hashlib.sha256((root / "run_contract.json").read_bytes()).hexdigest(),
    "summary_sha256": hashlib.sha256((root / "summary.json").read_bytes()).hexdigest(),
    "sealed_test_access_count": 0,
    "production_authorized": False,
    "foundation_model_claim_authorized": False,
}
Path(os.environ["TWIRL_FM0_INVOCATION_OUT"]).write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n"
)
PY
sha256sum "${INVOCATION}" > "${INVOCATION}.sha256"
chmod a-w "${INVOCATION}" "${INVOCATION}.sha256"
echo "[fm0.2-canary] complete=${INVOCATION}; stopped_at=${STOP_AFTER}; no sealed access or model claim"
