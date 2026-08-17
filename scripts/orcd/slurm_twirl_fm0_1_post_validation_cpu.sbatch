#!/bin/bash
# Independently validate the bounded synthetic H200 mechanics smoke.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-postvalidate
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
: "${TWIRL_FM0_SMOKE_OUTPUT:?missing smoke output}"
: "${TWIRL_FM0_INPUT_RECEIPT:?missing exact input receipt}"
: "${TWIRL_FM0_INPUT_RECEIPT_SHA256:?missing exact input receipt hash}"
: "${TWIRL_FM0_ADMISSION_RECEIPT:?missing exact admission receipt}"
: "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256:?missing admission receipt hash}"
readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly RUN_ROOT="$(realpath -m -- "${TWIRL_FM0_RUN_ROOT}")"
readonly SMOKE="$(realpath -m -- "${TWIRL_FM0_SMOKE_OUTPUT}")"
readonly INPUT_RECEIPT="$(realpath -m -- "${TWIRL_FM0_INPUT_RECEIPT}")"
readonly ADMISSION_RECEIPT="$(realpath -m -- "${TWIRL_FM0_ADMISSION_RECEIPT}")"
readonly PYTHON="/orcd/data/mki_aryeh/001/twirl/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
readonly OUT="${SMOKE}/post_validation.receipt.json"
readonly RUN_VALIDATION="${SMOKE}/independent_run_validation.json"
die() { echo "[fm0-postvalidate] ERROR: $*" >&2; exit 2; }
for path in "${SMOKE}" "${INPUT_RECEIPT}" "${ADMISSION_RECEIPT}"; do
  [[ "${path}" == "${RUN_ROOT}"/* ]] || die "authority path escapes run root: ${path}"
done
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" && -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "validator repository is not exact and clean"
[[ "$(sha256sum "${INPUT_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" ]] || die "input receipt drift"
[[ "$(sha256sum "${ADMISSION_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}" ]] || die "admission receipt drift"
[[ ! -e "${OUT}" && ! -e "${RUN_VALIDATION}" ]] || die "refusing to overwrite post-validation evidence"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
"${PYTHON}" scripts/stage5_validation/validate_twirl_fm0_release.py \
  --run-dir "${SMOKE}" \
  --output "${RUN_VALIDATION}"
chmod a-w "${RUN_VALIDATION}" "${RUN_VALIDATION}.sha256"
TWIRL_FM0_SMOKE="${SMOKE}" TWIRL_FM0_VALIDATOR_SHA="${TWIRL_EXPECTED_GIT_SHA}" \
TWIRL_FM0_INPUT_SHA="${TWIRL_FM0_INPUT_RECEIPT_SHA256}" \
TWIRL_FM0_ADMISSION_SHA="${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}" \
TWIRL_FM0_INPUT_RECEIPT="${INPUT_RECEIPT}" \
TWIRL_FM0_ADMISSION_RECEIPT="${ADMISSION_RECEIPT}" \
TWIRL_FM0_POST_OUT="${OUT}" "${PYTHON}" - <<'PY'
from datetime import datetime, timezone
import hashlib, json, os
from pathlib import Path
root = Path(os.environ["TWIRL_FM0_SMOKE"])
input_receipt = json.loads(Path(os.environ["TWIRL_FM0_INPUT_RECEIPT"]).read_text())
admission_receipt = json.loads(Path(os.environ["TWIRL_FM0_ADMISSION_RECEIPT"]).read_text())
if input_receipt.get("schema_version") != "twirl_fm0_1_input_release_receipt_v1" or input_receipt.get("passed") is not True:
    raise SystemExit("input receipt did not pass")
if admission_receipt.get("schema_version") != "twirl_fm0_1_orcd_admission_receipt_v1" or admission_receipt.get("passed") is not True:
    raise SystemExit("admission receipt did not pass")
execution = json.loads((root / "orcd_execution_receipt.json").read_text())
run_validation = json.loads((root / "independent_run_validation.json").read_text())
if run_validation.get("passed") is not True or run_validation.get("checkpoint_inspected") is not True:
    raise SystemExit("independent checkpoint inspection did not pass")
if execution.get("schema_version") != "twirl_fm0_1_orcd_synthetic_fp32_smoke_v1" or execution.get("passed") is not True:
    raise SystemExit("execution receipt did not pass")
if execution.get("synthetic_only") is not True or execution.get("real_input_consumed") is not False:
    raise SystemExit("smoke claim boundary drift")
if execution["input_receipt_sha256"] != os.environ["TWIRL_FM0_INPUT_SHA"] or execution["admission_receipt_sha256"] != os.environ["TWIRL_FM0_ADMISSION_SHA"]:
    raise SystemExit("execution binding drift")
for name, expected in execution["artifacts"].items():
    if hashlib.sha256((root / name).read_bytes()).hexdigest() != expected:
        raise SystemExit(f"artifact drift: {name}")
summary = json.loads((root / "summary.json").read_text())
contract = json.loads((root / "run_contract.json").read_text())
out = Path(os.environ["TWIRL_FM0_POST_OUT"])
payload = {
  "schema_version": "twirl_fm0_1_orcd_post_validation_v1",
  "created_at_utc": datetime.now(timezone.utc).isoformat(), "passed": True,
  "validator_git_sha": os.environ["TWIRL_FM0_VALIDATOR_SHA"],
  "synthetic_only": True, "real_data_poc_certified": False,
  "production_authorized": False, "foundation_model_claim_authorized": False,
  "summary_sha256": hashlib.sha256((root / "summary.json").read_bytes()).hexdigest(),
  "run_contract_sha256": hashlib.sha256((root / "run_contract.json").read_bytes()).hexdigest(),
  "independent_run_validation_sha256": hashlib.sha256((root / "independent_run_validation.json").read_bytes()).hexdigest(),
}
out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
chmod a-w "${OUT}"
sha256sum "${OUT}" > "${OUT}.sha256"
chmod a-w "${OUT}.sha256"
echo "[fm0-postvalidate] complete=${OUT} real_data_poc_certified=false"
