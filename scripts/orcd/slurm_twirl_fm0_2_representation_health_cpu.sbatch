#!/bin/bash
# Development-only FM0.2 representation evaluation; never opens the sealed test.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-2-eval
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
readonly CONFIG_SHA="394a0340f9f8518423e39a1c83bc32acde4e2ebff95c4837b0edb5399bd50801"
readonly EVALUATOR_CONTRACT_SHA="ed8004180b28c7e2bd85911d0128f9180eddd3c35bdaae299e2f7edba6e5baac"
readonly MERGE_SUMMARY_SHA="ffb02bc308fd5be6931ce2fed3a54bc99936b85bf93c127aaba539121084bec6"
: "${TWIRL_ORCD_REPO:?missing exact evaluator repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?missing exact evaluator commit}"
: "${TWIRL_FM0_ARTIFACT_GIT_SHA:?missing exact training-artifact commit}"
: "${TWIRL_FM0_ARTIFACT_RUN_ROOT:?missing training-artifact run root}"
: "${TWIRL_FM0_STOP_AFTER_STEP:?missing representation milestone}"
: "${TWIRL_FM0_RUN_OUTPUT:?missing exact run output}"
: "${TWIRL_FM0_RESUME_CHECKPOINT:?missing immutable milestone checkpoint}"
: "${TWIRL_FM0_INPUT_RECEIPT:?missing exact input receipt}"
: "${TWIRL_FM0_INPUT_RECEIPT_SHA256:?missing input receipt hash}"
: "${TWIRL_FM0_PREVIOUS_POST_VALIDATION:?missing milestone post-validation receipt}"
: "${TWIRL_FM0_PREVIOUS_POST_VALIDATION_SHA256:?missing post-validation hash}"
readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly ARTIFACT_ROOT="$(realpath -m -- "${TWIRL_FM0_ARTIFACT_RUN_ROOT}")"
readonly TRAIN="$(realpath -m -- "${TWIRL_FM0_RUN_OUTPUT}")"
readonly CHECKPOINT="$(realpath -m -- "${TWIRL_FM0_RESUME_CHECKPOINT}")"
readonly INPUT_RECEIPT="$(realpath -m -- "${TWIRL_FM0_INPUT_RECEIPT}")"
readonly POST="$(realpath -m -- "${TWIRL_FM0_PREVIOUS_POST_VALIDATION}")"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
readonly RELEASE_ROOT="${PROJECT_ROOT}/exports/fm0/input_releases/twirl_fm0_1_s56_s64_v2"
readonly STAGING_ROOT="${PROJECT_ROOT}/exports/fm0/input_staging/a2v1_s56_s64_fm0_1_v2"
readonly AUTHORITY_DIR="${ARTIFACT_ROOT}/evaluations/observation_sector_authority_v1"
readonly AUTHORITY_CSV="${AUTHORITY_DIR}/observation_sector_authority.csv"
readonly AUTHORITY_RECEIPT="${AUTHORITY_DIR}/observation_sector_authority.receipt.json"
printf -v PADDED_STEP '%08d' "${TWIRL_FM0_STOP_AFTER_STEP}"
readonly OUTPUT_DIR="${ARTIFACT_ROOT}/evaluations/step_${PADDED_STEP}"
readonly OUTPUT="${OUTPUT_DIR}/representation_health.json"
readonly OUTPUT_RECEIPT="${OUTPUT_DIR}/representation_health.receipt.json"
readonly EVALUATOR_CONTRACT="${REPO}/reports/stage5_validation/twirl_fm0_2_design_freeze_v1/development_evaluator_v2.contract.json"
die() { echo "[fm0.2-representation] ERROR: $*" >&2; exit 2; }

[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ && "${TWIRL_FM0_ARTIFACT_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || die "invalid Git revision"
[[ "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" =~ ^[0-9a-f]{64}$ && "${TWIRL_FM0_PREVIOUS_POST_VALIDATION_SHA256}" =~ ^[0-9a-f]{64}$ ]] || die "invalid receipt hash"
case "${TWIRL_FM0_STOP_AFTER_STEP}" in 500|1000|2000) ;; *) die "unauthorized representation milestone" ;; esac
[[ "${ARTIFACT_ROOT}" == "${PROJECT_ROOT}/reports/stage5_validation/twirl_fm0_2_s56_s64_objective_canary/${TWIRL_FM0_ARTIFACT_GIT_SHA:0:12}" ]] || die "artifact root/revision mismatch"
[[ "${TRAIN}" == "${ARTIFACT_ROOT}/model_runs/fm0_2_1_canary/seed560067" ]] || die "unexpected training path"
[[ "${CHECKPOINT}" == "${TRAIN}/checkpoint_step_${PADDED_STEP}.pt" ]] || die "checkpoint/milestone mismatch"
[[ "${POST}" == "${ARTIFACT_ROOT}/validations/seed560067-real_canary-step${PADDED_STEP}/post_validation.receipt.json" ]] || die "post-validation/milestone mismatch"
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" && -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "evaluator repository is not exact and clean"
[[ "$(sha256sum "${REPO}/configs/models/twirl_fm0_2_s56_s64_objective_canary.yaml" | awk '{print $1}')" == "${CONFIG_SHA}" ]] || die "scientific config drift"
[[ "$(sha256sum "${EVALUATOR_CONTRACT}" | awk '{print $1}')" == "${EVALUATOR_CONTRACT_SHA}" ]] || die "evaluator contract drift"
[[ "$(sha256sum "${INPUT_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" ]] || die "input receipt drift"
[[ "$(sha256sum "${POST}" | awk '{print $1}')" == "${TWIRL_FM0_PREVIOUS_POST_VALIDATION_SHA256}" ]] || die "post-validation receipt drift"
[[ -s "${CHECKPOINT}" && -s "${CHECKPOINT}.sha256" ]] || die "milestone checkpoint is incomplete"
[[ "$(sha256sum "${CHECKPOINT}" | awk '{print $1}')" == "$(awk 'NR==1 {print $1}' "${CHECKPOINT}.sha256")" ]] || die "milestone checkpoint hash drift"

"${PYTHON}" - "${POST}" "${TWIRL_FM0_STOP_AFTER_STEP}" "${TWIRL_FM0_ARTIFACT_GIT_SHA}" "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" "${CHECKPOINT}" <<'PY'
import hashlib, json, sys
from pathlib import Path
post = json.loads(Path(sys.argv[1]).read_text())
step = int(sys.argv[2])
checkpoint = Path(sys.argv[5])
observed_checkpoint = hashlib.sha256(checkpoint.read_bytes()).hexdigest()
expected_checkpoint = post.get("artifacts_sha256", {}).get(checkpoint.name)
assert post.get("schema_version") == "twirl_fm0_2_milestone_post_validation_v1"
assert post.get("passed") is True and post.get("run_kind") == "real_canary"
assert post.get("stop_after_step") == step and post.get("run_git_sha") == sys.argv[3]
assert post.get("input_receipt_sha256") == sys.argv[4]
assert post.get("sealed_test_access_count") == 0
assert expected_checkpoint == observed_checkpoint
PY

existing=0
for path in "${AUTHORITY_CSV}" "${AUTHORITY_CSV}.sha256" "${AUTHORITY_RECEIPT}" "${AUTHORITY_RECEIPT}.sha256"; do
  [[ -e "${path}" ]] && existing=$((existing + 1))
done
if [[ "${existing}" == 0 ]]; then
  mkdir -p "${AUTHORITY_DIR}"
  "${PYTHON}" "${REPO}/scripts/stage5_validation/build_twirl_fm0_observation_sector_authority.py" \
    --merged-manifest "${RELEASE_ROOT}/manifest.csv" \
    --merge-summary "${RELEASE_ROOT}/merge_summary.json" \
    --merge-summary-sha256 "${MERGE_SUMMARY_SHA}" \
    --sector-manifest "56=${STAGING_ROOT}/s0056/manifest.csv" \
    --sector-manifest "57=${STAGING_ROOT}/s0057/manifest.csv" \
    --sector-manifest "58=${STAGING_ROOT}/s0058/manifest.csv" \
    --sector-manifest "59=${STAGING_ROOT}/s0059/manifest.csv" \
    --sector-manifest "60=${STAGING_ROOT}/s0060/manifest.csv" \
    --sector-manifest "61=${STAGING_ROOT}/s0061/manifest.csv" \
    --sector-manifest "62=${STAGING_ROOT}/s0062/manifest.csv" \
    --sector-manifest "63=${STAGING_ROOT}/s0063/manifest.csv" \
    --sector-manifest "64=${STAGING_ROOT}/s0064/manifest.csv" \
    --output-csv "${AUTHORITY_CSV}" \
    --output-receipt "${AUTHORITY_RECEIPT}"
  chmod 0440 "${AUTHORITY_CSV}" "${AUTHORITY_CSV}.sha256" "${AUTHORITY_RECEIPT}" "${AUTHORITY_RECEIPT}.sha256"
elif [[ "${existing}" != 4 ]]; then
  die "partial observation-sector authority exists"
fi
readonly AUTHORITY_SHA="$(sha256sum "${AUTHORITY_CSV}" | awk '{print $1}')"
"${PYTHON}" - "${AUTHORITY_RECEIPT}" "${AUTHORITY_SHA}" <<'PY'
import json, sys
from pathlib import Path
receipt = json.loads(Path(sys.argv[1]).read_text())
assert receipt.get("schema_version") == "twirl_fm0_2_observation_sector_authority_receipt_v1"
assert receipt.get("passed") is True and receipt.get("sealed_rows_emitted") == 0
assert receipt.get("model_visible") is False
assert receipt.get("output_authority", {}).get("sha256") == sys.argv[2]
PY

mkdir -p "${OUTPUT_DIR}"
[[ ! -e "${OUTPUT}" && ! -e "${OUTPUT}.sha256" && ! -e "${OUTPUT_RECEIPT}" ]] || die "refusing to overwrite evaluation"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
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
  --expected-git-sha "${TWIRL_EXPECTED_GIT_SHA}"

"${PYTHON}" - "${OUTPUT}" "${OUTPUT_RECEIPT}" "${EVALUATOR_CONTRACT}" "${AUTHORITY_SHA}" "${POST}" "${CHECKPOINT}" "${TWIRL_EXPECTED_GIT_SHA}" "${TWIRL_FM0_ARTIFACT_GIT_SHA}" "${TWIRL_FM0_STOP_AFTER_STEP}" <<'PY'
import hashlib, json, os, sys
from pathlib import Path
output = Path(sys.argv[1])
receipt_path = Path(sys.argv[2])
contract = Path(sys.argv[3])
authority_sha = sys.argv[4]
post = Path(sys.argv[5])
checkpoint = Path(sys.argv[6])
payload = json.loads(output.read_text())
assert payload.get("schema_version") == "twirl_fm0_representation_health_v2"
assert payload.get("passed") is True
assert payload.get("run", {}).get("global_step") == int(sys.argv[9])
assert payload.get("run", {}).get("checkpoint_path") == str(checkpoint)
assert payload.get("data_access", {}).get("sealed_test_access_count") == 0
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
receipt = {
    "schema_version": "twirl_fm0_2_representation_evaluation_receipt_v1",
    "passed": True,
    "evaluator_git_sha": sys.argv[7],
    "artifact_git_sha": sys.argv[8],
    "stop_after_step": int(sys.argv[9]),
    "evaluation_sha256": sha(output),
    "evaluator_contract_sha256": sha(contract),
    "observation_sector_authority_sha256": authority_sha,
    "post_validation_sha256": sha(post),
    "checkpoint_sha256": sha(checkpoint),
    "sealed_test_access_count": 0,
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
echo "[fm0.2-representation] complete=${OUTPUT}; no sealed access or go/no-go claim"
