#!/bin/bash
# CPU-only independent validation for one FM0.2 synthetic or real milestone.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-2-validate
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 027
readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
: "${TWIRL_ORCD_REPO:?missing detached validator repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?missing exact validator SHA}"
: "${TWIRL_FM0_RUN_GIT_SHA:?missing exact training SHA}"
: "${TWIRL_FM0_RUN_ROOT:?missing run root}"
: "${TWIRL_FM0_RUN_KIND:?missing synthetic_fp32 or real_canary run kind}"
: "${TWIRL_FM0_STOP_AFTER_STEP:?missing stop-after step}"
: "${TWIRL_FM0_RUN_OUTPUT:?missing run output}"
: "${TWIRL_FM0_INPUT_RECEIPT:?missing upstream input receipt}"
: "${TWIRL_FM0_INPUT_RECEIPT_SHA256:?missing upstream input receipt hash}"
: "${TWIRL_FM0_INPUT_REUSE_RECEIPT:?missing input-reuse receipt}"
: "${TWIRL_FM0_INPUT_REUSE_RECEIPT_SHA256:?missing input-reuse receipt hash}"
: "${TWIRL_FM0_ADMISSION_RECEIPT:?missing admission receipt}"
: "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256:?missing admission receipt hash}"
readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly RUN_ROOT="$(realpath -m -- "${TWIRL_FM0_RUN_ROOT}")"
readonly RUN="$(realpath -e -- "${TWIRL_FM0_RUN_OUTPUT}")"
readonly INPUT_RECEIPT="$(realpath -e -- "${TWIRL_FM0_INPUT_RECEIPT}")"
readonly REUSE_RECEIPT="$(realpath -e -- "${TWIRL_FM0_INPUT_REUSE_RECEIPT}")"
readonly ADMISSION_RECEIPT="$(realpath -e -- "${TWIRL_FM0_ADMISSION_RECEIPT}")"
readonly KIND="${TWIRL_FM0_RUN_KIND}"
readonly STOP_AFTER="${TWIRL_FM0_STOP_AFTER_STEP}"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
die() { echo "[fm0.2-validate] ERROR: $*" >&2; exit 2; }

case "${KIND}:${STOP_AFTER}" in
  synthetic_fp32:8|real_canary:64|real_canary:500|real_canary:1000|real_canary:2000) ;;
  *) die "invalid run-kind/stop-after pair" ;;
esac
[[ "${RUN_ROOT}" == "${PROJECT_ROOT}"/* && "${RUN}" == "${RUN_ROOT}"/* ]] || die "unsafe run path"
[[ "${INPUT_RECEIPT}" == "${PROJECT_ROOT}"/* && "${REUSE_RECEIPT}" == "${REPO}"/* && "${ADMISSION_RECEIPT}" == "${RUN_ROOT}"/* ]] || die "unsafe authority path"
[[ "${TWIRL_FM0_RUN_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || die "invalid training SHA"
for value in "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" "${TWIRL_FM0_INPUT_REUSE_RECEIPT_SHA256}" "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}"; do
  [[ "${value}" =~ ^[0-9a-f]{64}$ ]] || die "invalid receipt hash"
done
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" && -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "validator repository is not exact and clean"
[[ "$(sha256sum "${INPUT_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" ]] || die "upstream input receipt drift"
[[ "$(sha256sum "${REUSE_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_INPUT_REUSE_RECEIPT_SHA256}" ]] || die "input-reuse receipt drift"
[[ "$(sha256sum "${ADMISSION_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}" ]] || die "admission receipt drift"
readonly PADDED="$(printf '%08d' "${STOP_AFTER}")"
readonly VALIDATION_ROOT="${RUN_ROOT}/validations/$(basename "${RUN}")-${KIND}-step${PADDED}"
readonly RUN_VALIDATION="${VALIDATION_ROOT}/independent_run_validation.json"
readonly OUT="${VALIDATION_ROOT}/post_validation.receipt.json"
[[ ! -e "${VALIDATION_ROOT}" ]] || die "refusing to overwrite validation evidence"
mkdir -p "${VALIDATION_ROOT}"

cd "${REPO}"
export PYTHONPATH="${REPO}/src"
"${PYTHON}" scripts/stage5_validation/validate_twirl_fm0_release.py \
  --run-dir "${RUN}" \
  --output "${RUN_VALIDATION}"

TWIRL_FM0_POST_OUT="${OUT}" TWIRL_FM0_RUN_VALIDATION="${RUN_VALIDATION}" \
TWIRL_FM0_RUN="${RUN}" TWIRL_FM0_KIND="${KIND}" TWIRL_FM0_STOP_AFTER="${STOP_AFTER}" \
TWIRL_FM0_VALIDATOR_SHA="${TWIRL_EXPECTED_GIT_SHA}" TWIRL_FM0_RUN_SHA="${TWIRL_FM0_RUN_GIT_SHA}" \
TWIRL_FM0_INPUT_SHA="${TWIRL_FM0_INPUT_RECEIPT_SHA256}" \
TWIRL_FM0_REUSE_SHA="${TWIRL_FM0_INPUT_REUSE_RECEIPT_SHA256}" \
TWIRL_FM0_ADMISSION_SHA="${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}" \
"${PYTHON}" - <<'PY'
from datetime import datetime, timezone
import hashlib, json, os
from pathlib import Path

run = Path(os.environ["TWIRL_FM0_RUN"])
kind = os.environ["TWIRL_FM0_KIND"]
stop = int(os.environ["TWIRL_FM0_STOP_AFTER"])
validation_path = Path(os.environ["TWIRL_FM0_RUN_VALIDATION"])
validation = json.loads(validation_path.read_text())
contract = json.loads((run / "run_contract.json").read_text())
summary = json.loads((run / "summary.json").read_text())
if validation.get("passed") is not True or validation.get("checkpoint_inspected") is not True:
    raise SystemExit("independent checkpoint inspection did not pass")
if contract.get("campaign_id") != "twirl_fm0_2_s56_s64_objective_canary":
    raise SystemExit("wrong campaign in run contract")
if contract.get("variant") != "TWIRL-FM0.2.1" or contract.get("architecture") != "tcn":
    raise SystemExit("wrong canary model identity")
if contract.get("objective") != "masked_reconstruction_plus_same_window_vicreg":
    raise SystemExit("wrong canary objective")
if contract.get("use_vicreg") is not True or contract.get("reconstruct_second_view") is not False:
    raise SystemExit("canary objective switches drifted")
if contract.get("seed") != 560067 or contract.get("training_horizon_step") != 20000:
    raise SystemExit("seed or run-contract horizon drifted")
if contract.get("expected_git_sha") != os.environ["TWIRL_FM0_RUN_SHA"]:
    raise SystemExit("run code binding drift")
if summary.get("global_step") != stop or summary.get("requested_stop_after_step") != stop:
    raise SystemExit("summary stop-after step drifted")
if summary.get("scientific_result") is not False:
    raise SystemExit("canary claim boundary drifted")
if kind == "synthetic_fp32":
    if summary.get("synthetic_only") is not True or summary.get("real_data_consumed") is not False:
        raise SystemExit("synthetic smoke data boundary drifted")
    milestone_names = []
else:
    if summary.get("synthetic_only") is not False or summary.get("real_data_consumed") is not True:
        raise SystemExit("real canary data boundary drifted")
    if contract.get("input_release", {}).get("receipt_sha256") != os.environ["TWIRL_FM0_INPUT_SHA"]:
        raise SystemExit("run upstream input binding drifted")
    if contract.get("input_release_reuse", {}).get("sha256") != os.environ["TWIRL_FM0_REUSE_SHA"]:
        raise SystemExit("run input-reuse binding drifted")
    milestone_names = [
        f"checkpoint_step_{step:08d}.pt"
        for step in (0, 64, 500, 1000, 2000)
        if step <= stop
    ]

artifacts = ["run_contract.json", "checkpoint.pt", "summary.json", *milestone_names]
artifact_hashes = {}
for name in artifacts:
    path = run / name
    sidecar = path.with_name(path.name + ".sha256")
    if not path.is_file() or not sidecar.is_file():
        raise SystemExit(f"missing artifact or sidecar: {name}")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    fields = sidecar.read_text().strip().split()
    if len(fields) != 2 or fields[0] != digest or fields[1] != path.name:
        raise SystemExit(f"artifact sidecar mismatch: {name}")
    artifact_hashes[name] = digest

payload = {
    "schema_version": "twirl_fm0_2_milestone_post_validation_v1",
    "created_at_utc": datetime.now(timezone.utc).isoformat(),
    "passed": True,
    "run_kind": kind,
    "variant": "TWIRL-FM0.2.1",
    "seed": 560067,
    "run_git_sha": os.environ["TWIRL_FM0_RUN_SHA"],
    "validator_git_sha": os.environ["TWIRL_FM0_VALIDATOR_SHA"],
    "run_contract_target_step": 20000,
    "stop_after_step": stop,
    "input_receipt_sha256": os.environ["TWIRL_FM0_INPUT_SHA"],
    "input_reuse_receipt_sha256": os.environ["TWIRL_FM0_REUSE_SHA"],
    "admission_receipt_sha256": os.environ["TWIRL_FM0_ADMISSION_SHA"],
    "independent_run_validation_sha256": hashlib.sha256(validation_path.read_bytes()).hexdigest(),
    "artifacts_sha256": artifact_hashes,
    "sealed_test_access_count": 0,
    "scientific_result": False,
    "production_authorized": False,
    "foundation_model_claim_authorized": False,
}
Path(os.environ["TWIRL_FM0_POST_OUT"]).write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n"
)
PY
for path in "${RUN_VALIDATION}" "${OUT}"; do
  sha256sum "${path}" > "${path}.sha256"
done
chmod -R a-w "${VALIDATION_ROOT}"
echo "[fm0.2-validate] complete=${OUT}; no sealed access or model claim"
