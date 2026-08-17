#!/bin/bash
# Independently validate one immutable real-data FM0.1.1 run.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-real-validate
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 027
: "${TWIRL_ORCD_REPO:?missing detached repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?missing exact validator SHA}"
: "${TWIRL_FM0_RUN_ROOT:?missing run root}"
: "${TWIRL_FM0_TRAIN_OUTPUT:?missing real training output}"
: "${TWIRL_FM0_INPUT_RECEIPT:?missing exact input receipt}"
: "${TWIRL_FM0_INPUT_RECEIPT_SHA256:?missing exact input receipt hash}"
: "${TWIRL_FM0_ADMISSION_RECEIPT:?missing exact admission receipt}"
: "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256:?missing admission receipt hash}"
readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly RUN_ROOT="$(realpath -m -- "${TWIRL_FM0_RUN_ROOT}")"
readonly TRAIN="$(realpath -m -- "${TWIRL_FM0_TRAIN_OUTPUT}")"
readonly INPUT_RECEIPT="$(realpath -m -- "${TWIRL_FM0_INPUT_RECEIPT}")"
readonly ADMISSION_RECEIPT="$(realpath -m -- "${TWIRL_FM0_ADMISSION_RECEIPT}")"
readonly PYTHON="/orcd/data/mki_aryeh/001/twirl/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
readonly VALIDATION_ROOT="${RUN_ROOT}/validations/$(basename "${TRAIN}")"
readonly RUN_VALIDATION="${VALIDATION_ROOT}/independent_run_validation.json"
readonly OUT="${VALIDATION_ROOT}/post_validation.receipt.json"
die() { echo "[fm0-real-validate] ERROR: $*" >&2; exit 2; }
for path in "${TRAIN}" "${INPUT_RECEIPT}" "${ADMISSION_RECEIPT}"; do
  [[ "${path}" == "${RUN_ROOT}"/* ]] || die "authority path escapes run root: ${path}"
done
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" && -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "validator repository is not exact and clean"
[[ "$(sha256sum "${INPUT_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" ]] || die "input receipt drift"
[[ "$(sha256sum "${ADMISSION_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}" ]] || die "admission receipt drift"
[[ ! -e "${VALIDATION_ROOT}" ]] || die "refusing to overwrite validation evidence"
mkdir -p "${VALIDATION_ROOT}"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
"${PYTHON}" scripts/stage5_validation/validate_twirl_fm0_release.py \
  --run-dir "${TRAIN}" \
  --output "${RUN_VALIDATION}"
TWIRL_FM0_TRAIN="${TRAIN}" TWIRL_FM0_VALIDATOR_SHA="${TWIRL_EXPECTED_GIT_SHA}" \
TWIRL_FM0_INPUT_SHA="${TWIRL_FM0_INPUT_RECEIPT_SHA256}" \
TWIRL_FM0_ADMISSION_SHA="${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}" \
TWIRL_FM0_INPUT_RECEIPT="${INPUT_RECEIPT}" \
TWIRL_FM0_ADMISSION_RECEIPT="${ADMISSION_RECEIPT}" \
TWIRL_FM0_RUN_VALIDATION="${RUN_VALIDATION}" \
TWIRL_FM0_POST_OUT="${OUT}" "${PYTHON}" - <<'PY'
from datetime import datetime, timezone
import hashlib, json, os
from pathlib import Path

train = Path(os.environ["TWIRL_FM0_TRAIN"])
input_receipt = json.loads(Path(os.environ["TWIRL_FM0_INPUT_RECEIPT"]).read_text())
admission = json.loads(Path(os.environ["TWIRL_FM0_ADMISSION_RECEIPT"]).read_text())
validation_path = Path(os.environ["TWIRL_FM0_RUN_VALIDATION"])
validation = json.loads(validation_path.read_text())
contract = json.loads((train / "run_contract.json").read_text())
summary = json.loads((train / "summary.json").read_text())
if input_receipt.get("schema_version") != "twirl_fm0_1_input_release_receipt_v1" or input_receipt.get("passed") is not True:
    raise SystemExit("input receipt did not pass")
if admission.get("schema_version") != "twirl_fm0_1_orcd_admission_receipt_v1" or admission.get("passed") is not True:
    raise SystemExit("admission receipt did not pass")
if admission.get("expected_git_sha") != os.environ["TWIRL_FM0_VALIDATOR_SHA"]:
    raise SystemExit("admission code binding drift")
if admission.get("input_receipt", {}).get("sha256") != os.environ["TWIRL_FM0_INPUT_SHA"]:
    raise SystemExit("admission input binding drift")
if admission.get("requested_fm_resources") != {"h200": 1, "cpus": 4, "memory_mib": 32768}:
    raise SystemExit("admission resource drift")
if validation.get("passed") is not True or validation.get("checkpoint_inspected") is not True:
    raise SystemExit("independent checkpoint inspection did not pass")
if validation.get("real_data_consumed") is not True or validation.get("synthetic_only") is not False:
    raise SystemExit("real-data validation boundary drift")
if contract.get("expected_git_sha") != os.environ["TWIRL_FM0_VALIDATOR_SHA"]:
    raise SystemExit("run code binding drift")
if contract.get("input_release", {}).get("receipt_sha256") != os.environ["TWIRL_FM0_INPUT_SHA"]:
    raise SystemExit("run input binding drift")
if summary.get("scientific_result") is not False:
    raise SystemExit("proof-of-concept claim boundary drift")
artifact_names = ("run_contract.json", "checkpoint.pt", "summary.json")
payload = {
    "schema_version": "twirl_fm0_1_orcd_real_post_validation_v1",
    "created_at_utc": datetime.now(timezone.utc).isoformat(),
    "passed": True,
    "validator_git_sha": os.environ["TWIRL_FM0_VALIDATOR_SHA"],
    "real_data_consumed": True,
    "scientific_result": False,
    "production_authorized": False,
    "foundation_model_claim_authorized": False,
    "input_receipt_sha256": os.environ["TWIRL_FM0_INPUT_SHA"],
    "admission_receipt_sha256": os.environ["TWIRL_FM0_ADMISSION_SHA"],
    "independent_run_validation_sha256": hashlib.sha256(validation_path.read_bytes()).hexdigest(),
    "artifacts_sha256": {
        name: hashlib.sha256((train / name).read_bytes()).hexdigest()
        for name in artifact_names
    },
}
Path(os.environ["TWIRL_FM0_POST_OUT"]).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
for path in "${RUN_VALIDATION}" "${OUT}"; do
  sha256sum "${path}" > "${path}.sha256"
done
chmod -R a-w "${VALIDATION_ROOT}"
echo "[fm0-real-validate] complete=${OUT}"
