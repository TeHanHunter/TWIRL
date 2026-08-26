#!/bin/bash
# Dedicated development-only evaluation of the exact FM0.2 seed-560067 initialization.
# This wrapper does not train, open the sealed test, run event retention, or apply a gate.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-2-eval-init
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 06:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 027
readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly EXPECTED_TRAINING_GIT_SHA="ddf442aafb8f62966e549e2287abad3474dd556a"
readonly EXPECTED_EVALUATOR_GIT_SHA="83816d07975eebe3825d76dfe7096d22b70376f5"
readonly EXPECTED_INITIALIZATION_SHA="92463070381486f0c6053c190d62da8d0c5c0d31be8072e2bdcd677329ac792c"
readonly EXPECTED_STEP2000_POST_SHA="80d678c4621a22258c5d0e8f0be6f7e08ffa5bf93841c315abf42e9bb006b110"
readonly EXPECTED_STEP2000_REPRESENTATION_RECEIPT_SHA="616ef9a100b7b7a3a1923f81ac19a272cab8b6f4657a4a58c2815688ee6d1191"
readonly EXPECTED_STEP2000_EVALUATION_SHA="3c802741744999fe553e71bc2dccfbef6c309a29a28c30dfaac693a798bb3718"
readonly EXPECTED_AUTHORITY_SHA="bdd3e8039b312aeb21557662226a8a11b761dd1da031742be42ee9b0d1c6edc5"
readonly CONFIG_SHA="394a0340f9f8518423e39a1c83bc32acde4e2ebff95c4837b0edb5399bd50801"
readonly EVALUATOR_CONTRACT_SHA="ed8004180b28c7e2bd85911d0128f9180eddd3c35bdaae299e2f7edba6e5baac"
readonly EVALUATOR_ENTRY_SHA="2ad38165d0d89acb08970e9bdf2a07df54022bf020b6a72e310e4cb4eb3f014e"
readonly EVALUATOR_MODULE_SHA="573335ee9e2e9f3a7cc4aedf42dec76a584f0a72cdebacbf9f0e805a81d489f8"
: "${TWIRL_ORCD_REPO:?missing exact orchestration repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?missing exact orchestration commit}"
: "${TWIRL_FM0_ARTIFACT_GIT_SHA:?missing exact training-artifact commit}"
: "${TWIRL_FM0_ARTIFACT_RUN_ROOT:?missing training-artifact run root}"
: "${TWIRL_FM0_RUN_OUTPUT:?missing exact run output}"
: "${TWIRL_FM0_INITIALIZATION_CHECKPOINT:?missing exact initialization checkpoint}"
: "${TWIRL_FM0_INPUT_RECEIPT:?missing exact input receipt}"
: "${TWIRL_FM0_INPUT_RECEIPT_SHA256:?missing input receipt hash}"
: "${TWIRL_FM0_STEP2000_POST_VALIDATION:?missing step-2000 post-validation receipt}"
: "${TWIRL_FM0_STEP2000_POST_VALIDATION_SHA256:?missing step-2000 post-validation hash}"
: "${TWIRL_FM0_STEP2000_REPRESENTATION_RECEIPT:?missing step-2000 representation receipt}"
: "${TWIRL_FM0_STEP2000_REPRESENTATION_RECEIPT_SHA256:?missing step-2000 representation receipt hash}"
readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly EVALUATOR_REPO="${PROJECT_ROOT}/code/TWIRL_fm0_2_${EXPECTED_EVALUATOR_GIT_SHA:0:12}"
readonly ARTIFACT_ROOT="$(realpath -m -- "${TWIRL_FM0_ARTIFACT_RUN_ROOT}")"
readonly TRAIN="$(realpath -m -- "${TWIRL_FM0_RUN_OUTPUT}")"
readonly CHECKPOINT="$(realpath -m -- "${TWIRL_FM0_INITIALIZATION_CHECKPOINT}")"
readonly INPUT_RECEIPT="$(realpath -m -- "${TWIRL_FM0_INPUT_RECEIPT}")"
readonly POST="$(realpath -m -- "${TWIRL_FM0_STEP2000_POST_VALIDATION}")"
readonly REFERENCE_RECEIPT="$(realpath -m -- "${TWIRL_FM0_STEP2000_REPRESENTATION_RECEIPT}")"
readonly REFERENCE_OUTPUT="${ARTIFACT_ROOT}/evaluations/step_00002000/representation_health.json"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
readonly AUTHORITY_DIR="${ARTIFACT_ROOT}/evaluations/observation_sector_authority_v1"
readonly AUTHORITY_CSV="${AUTHORITY_DIR}/observation_sector_authority.csv"
readonly AUTHORITY_RECEIPT="${AUTHORITY_DIR}/observation_sector_authority.receipt.json"
readonly OUTPUT_DIR="${ARTIFACT_ROOT}/evaluations/step_00000000"
readonly OUTPUT="${OUTPUT_DIR}/representation_health.json"
readonly OUTPUT_RECEIPT="${OUTPUT_DIR}/representation_health.receipt.json"
readonly EVALUATOR_CONTRACT="${EVALUATOR_REPO}/reports/stage5_validation/twirl_fm0_2_design_freeze_v1/development_evaluator_v2.contract.json"
die() { echo "[fm0.2-initialization] ERROR: $*" >&2; exit 2; }
require_sidecar() {
  local path=$1 expected=$2 observed name extra
  read -r observed name extra < "${path}.sha256" || die "cannot read SHA-256 sidecar for ${path}"
  [[ "${observed}" == "${expected}" && "${name}" == "$(basename "${path}")" && -z "${extra}" ]] || die "SHA-256 sidecar differs for ${path}"
}

[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || die "invalid orchestration revision"
[[ "${TWIRL_FM0_ARTIFACT_GIT_SHA}" == "${EXPECTED_TRAINING_GIT_SHA}" ]] || die "training revision is not the frozen FM0.2.1 revision"
[[ "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" =~ ^[0-9a-f]{64}$ ]] || die "invalid input receipt hash"
[[ "${TWIRL_FM0_STEP2000_POST_VALIDATION_SHA256}" == "${EXPECTED_STEP2000_POST_SHA}" ]] || die "post-validation hash is not the frozen receipt"
[[ "${TWIRL_FM0_STEP2000_REPRESENTATION_RECEIPT_SHA256}" == "${EXPECTED_STEP2000_REPRESENTATION_RECEIPT_SHA}" ]] || die "reference hash is not the frozen receipt"
[[ "${ARTIFACT_ROOT}" == "${PROJECT_ROOT}/reports/stage5_validation/twirl_fm0_2_s56_s64_objective_canary/${EXPECTED_TRAINING_GIT_SHA:0:12}" ]] || die "artifact root/revision mismatch"
[[ "${TRAIN}" == "${ARTIFACT_ROOT}/model_runs/fm0_2_1_canary/seed560067" ]] || die "unexpected training path"
[[ "${CHECKPOINT}" == "${TRAIN}/checkpoint_step_00000000.pt" ]] || die "checkpoint is not the exact initialization"
[[ "${POST}" == "${ARTIFACT_ROOT}/validations/seed560067-real_canary-step00002000/post_validation.receipt.json" ]] || die "unexpected step-2000 post-validation path"
[[ "${REFERENCE_RECEIPT}" == "${ARTIFACT_ROOT}/evaluations/step_00002000/representation_health.receipt.json" ]] || die "unexpected step-2000 representation receipt path"
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" && -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "orchestration repository is not exact and clean"
[[ "$(git -C "${EVALUATOR_REPO}" rev-parse HEAD)" == "${EXPECTED_EVALUATOR_GIT_SHA}" && -z "$(git -C "${EVALUATOR_REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "immutable evaluator repository is absent, dirty, or at the wrong revision"
[[ "$(sha256sum "${REPO}/configs/models/twirl_fm0_2_s56_s64_objective_canary.yaml" | awk '{print $1}')" == "${CONFIG_SHA}" ]] || die "scientific config drift"
[[ "$(sha256sum "${EVALUATOR_CONTRACT}" | awk '{print $1}')" == "${EVALUATOR_CONTRACT_SHA}" ]] || die "evaluator contract drift"
[[ "$(sha256sum "${EVALUATOR_REPO}/scripts/stage5_validation/evaluate_twirl_fm0_representation_health.py" | awk '{print $1}')" == "${EVALUATOR_ENTRY_SHA}" ]] || die "evaluator entrypoint drift"
[[ "$(sha256sum "${EVALUATOR_REPO}/src/twirl/models/fm0/representation_health.py" | awk '{print $1}')" == "${EVALUATOR_MODULE_SHA}" ]] || die "evaluator module drift"
[[ "$(sha256sum "${INPUT_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" ]] || die "input receipt drift"
[[ "$(sha256sum "${POST}" | awk '{print $1}')" == "${TWIRL_FM0_STEP2000_POST_VALIDATION_SHA256}" ]] || die "post-validation receipt drift"
[[ "$(sha256sum "${REFERENCE_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_STEP2000_REPRESENTATION_RECEIPT_SHA256}" ]] || die "step-2000 representation receipt drift"
[[ -s "${REFERENCE_OUTPUT}" && "$(sha256sum "${REFERENCE_OUTPUT}" | awk '{print $1}')" == "${EXPECTED_STEP2000_EVALUATION_SHA}" ]] || die "step-2000 representation output drift"
[[ -s "${CHECKPOINT}" && -s "${CHECKPOINT}.sha256" ]] || die "initialization checkpoint is incomplete"
[[ "$(sha256sum "${CHECKPOINT}" | awk '{print $1}')" == "${EXPECTED_INITIALIZATION_SHA}" ]] || die "initialization checkpoint hash drift"
require_sidecar "${CHECKPOINT}" "${EXPECTED_INITIALIZATION_SHA}"
"${PYTHON}" -c 'import sys; raise SystemExit(0 if sys.flags.optimize == 0 else 2)' || die "optimized Python would disable fail-closed assertions"

"${PYTHON}" - "${POST}" "${REFERENCE_RECEIPT}" "${CHECKPOINT}" "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" <<'PY'
import hashlib, json, sys
from pathlib import Path

post_path = Path(sys.argv[1])
reference_path = Path(sys.argv[2])
checkpoint = Path(sys.argv[3])
post = json.loads(post_path.read_text())
reference = json.loads(reference_path.read_text())
checkpoint_sha = hashlib.sha256(checkpoint.read_bytes()).hexdigest()
assert post.get("schema_version") == "twirl_fm0_2_milestone_post_validation_v1"
assert post.get("passed") is True and post.get("run_kind") == "real_canary"
assert post.get("stop_after_step") == 2000
assert post.get("run_git_sha") == "ddf442aafb8f62966e549e2287abad3474dd556a"
assert post.get("input_receipt_sha256") == sys.argv[4]
assert post.get("seed") == 560067 and post.get("variant") == "TWIRL-FM0.2.1"
assert post.get("sealed_test_access_count") == 0
assert post.get("artifacts_sha256", {}).get(checkpoint.name) == checkpoint_sha
assert reference.get("schema_version") == "twirl_fm0_2_representation_evaluation_receipt_v1"
assert reference.get("passed") is True and reference.get("stop_after_step") == 2000
assert reference.get("artifact_git_sha") == "ddf442aafb8f62966e549e2287abad3474dd556a"
assert reference.get("evaluator_git_sha") == "83816d07975eebe3825d76dfe7096d22b70376f5"
assert reference.get("evaluation_sha256") == "3c802741744999fe553e71bc2dccfbef6c309a29a28c30dfaac693a798bb3718"
assert reference.get("sealed_test_access_count") == 0
assert reference.get("scientific_go_no_go_applied") is False
PY

for path in "${AUTHORITY_CSV}" "${AUTHORITY_CSV}.sha256" "${AUTHORITY_RECEIPT}" "${AUTHORITY_RECEIPT}.sha256"; do
  [[ -s "${path}" ]] || die "the prior evaluator authority is incomplete"
done
readonly AUTHORITY_SHA="$(sha256sum "${AUTHORITY_CSV}" | awk '{print $1}')"
[[ "${AUTHORITY_SHA}" == "${EXPECTED_AUTHORITY_SHA}" ]] || die "authority differs from the step-2000 evaluation"
require_sidecar "${AUTHORITY_CSV}" "${AUTHORITY_SHA}"
require_sidecar "${AUTHORITY_RECEIPT}" "$(sha256sum "${AUTHORITY_RECEIPT}" | awk '{print $1}')"
"${PYTHON}" - "${AUTHORITY_RECEIPT}" "${AUTHORITY_SHA}" "${REFERENCE_RECEIPT}" <<'PY'
import json, sys
from pathlib import Path

authority = json.loads(Path(sys.argv[1]).read_text())
reference = json.loads(Path(sys.argv[3]).read_text())
assert authority.get("schema_version") == "twirl_fm0_2_observation_sector_authority_receipt_v1"
assert authority.get("passed") is True and authority.get("sealed_rows_emitted") == 0
assert authority.get("model_visible") is False
assert authority.get("output_authority", {}).get("sha256") == sys.argv[2]
assert reference.get("observation_sector_authority_sha256") == sys.argv[2]
PY

mkdir -- "${OUTPUT_DIR}" || die "refusing to overwrite, reuse, or race an initialization evaluation directory"
cd "${EVALUATOR_REPO}"
export PYTHONPATH="${EVALUATOR_REPO}/src"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
"${PYTHON}" scripts/stage5_validation/evaluate_twirl_fm0_representation_health.py \
  --run-dir "${TRAIN}" \
  --checkpoint "${CHECKPOINT}" \
  --output "${OUTPUT}" \
  --observation-sector-authority "${AUTHORITY_CSV}" \
  --observation-sector-authority-sha256 "${AUTHORITY_SHA}" \
  --device cpu \
  --max-components 256 \
  --max-pca-train-components 512 \
  --pca-components 32 \
  --batch-size 8 \
  --random-control-seed 0 \
  --expected-git-sha "${EXPECTED_EVALUATOR_GIT_SHA}"

"${PYTHON}" - "${OUTPUT}" "${OUTPUT_RECEIPT}" "${EVALUATOR_CONTRACT}" "${AUTHORITY_SHA}" "${POST}" "${REFERENCE_RECEIPT}" "${CHECKPOINT}" "${TWIRL_EXPECTED_GIT_SHA}" <<'PY'
import hashlib, json, os, sys
from pathlib import Path

output = Path(sys.argv[1])
receipt_path = Path(sys.argv[2])
contract = Path(sys.argv[3])
post = Path(sys.argv[5])
reference = Path(sys.argv[6])
checkpoint = Path(sys.argv[7])
payload = json.loads(output.read_text())
population = payload.get("evaluation_population", {})
assert payload.get("schema_version") == "twirl_fm0_representation_health_v2"
assert payload.get("passed") is True
assert payload.get("evaluator_git_sha") == "83816d07975eebe3825d76dfe7096d22b70376f5"
assert payload.get("run", {}).get("global_step") == 0
assert payload.get("run", {}).get("run_git_sha") == "ddf442aafb8f62966e549e2287abad3474dd556a"
assert payload.get("run", {}).get("variant") == "TWIRL-FM0.2.1"
assert payload.get("run", {}).get("architecture") == "tcn"
assert payload.get("run", {}).get("checkpoint_is_immutable_milestone") is True
assert payload.get("run", {}).get("checkpoint_path") == str(checkpoint)
assert payload.get("run", {}).get("checkpoint_sha256") == "92463070381486f0c6053c190d62da8d0c5c0d31be8072e2bdcd677329ac792c"
assert population.get("selected_leakage_components") == 256
assert population.get("selected_observation_visits") == 383
assert population.get("selected_leakage_components_sha256") == "50e0b59e1bc80ded61e30c1a6c07a56eee6306b03c3bd52e1f30ebdc43f5c72f"
assert population.get("selected_observation_keys_sha256") == "9ffdc49e042cf4e8dccbe53e70fd849723abf96625c36a00adea8f4fde20546a"
assert population.get("observation_sector_authority", {}).get("sha256") == sys.argv[4]
access = payload.get("data_access", {})
assert access.get("sealed_test_access_count") == 0
assert access.get("poc_sealed_test_window_reads") == 0
assert access.get("poc_sealed_test_unique_shards") == 0
assert "trained_encoder" in payload and "random_encoder_control" in payload
assert payload.get("random_encoder_control", {}).get("seed") == 0
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
receipt = {
    "schema_version": "twirl_fm0_2_initialization_representation_evaluation_receipt_v1",
    "passed": True,
    "orchestration_git_sha": sys.argv[8],
    "evaluator_git_sha": "83816d07975eebe3825d76dfe7096d22b70376f5",
    "artifact_git_sha": "ddf442aafb8f62966e549e2287abad3474dd556a",
    "evaluation_sha256": sha(output),
    "evaluator_contract_sha256": sha(contract),
    "observation_sector_authority_sha256": sys.argv[4],
    "step2000_post_validation_sha256": sha(post),
    "step2000_representation_receipt_sha256": sha(reference),
    "checkpoint_sha256": sha(checkpoint),
    "checkpoint_step": 0,
    "checkpoint_seed": 560067,
    "evaluation_role": "exact_same_seed_initialization",
    "output_field_role_mapping": {
        "trained_encoder": "exact_same_seed_initialization",
        "random_encoder_control": "separate_matched_random_control_seed_0",
        "source_paired_trained_minus_random": "exact_initialization_minus_matched_random_control",
    },
    "development_population": {
        "selected_leakage_components": 256,
        "selected_observation_visits": 383,
        "selected_leakage_components_sha256": "50e0b59e1bc80ded61e30c1a6c07a56eee6306b03c3bd52e1f30ebdc43f5c72f",
        "selected_observation_keys_sha256": "9ffdc49e042cf4e8dccbe53e70fd849723abf96625c36a00adea8f4fde20546a",
    },
    "trajectory_evidence": {
        "eligible_evidence_steps_after_receipt_validation": [0, 500, 1000, 2000],
        "performed_by_this_job": False,
    },
    "sealed_test_access_count": 0,
    "development_event_retention_run": False,
    "scientific_go_no_go_applied": False,
    "production_authorized": False,
    "foundation_model_claim_authorized": False,
}
temporary = receipt_path.with_name(receipt_path.name + ".tmp")
temporary.write_text(json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n")
os.chmod(temporary, 0o440)
os.replace(temporary, receipt_path)
digest = sha(receipt_path)
sidecar = receipt_path.with_name(receipt_path.name + ".sha256")
sidecar.write_text(f"{digest}  {receipt_path.name}\n")
os.chmod(sidecar, 0o440)
PY
echo "[fm0.2-initialization] complete=${OUTPUT}; exact initialization is separate from the matched-random control; no sealed access or go/no-go claim"
