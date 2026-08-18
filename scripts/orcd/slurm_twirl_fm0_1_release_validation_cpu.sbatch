#!/bin/bash
# Independently scan every immutable shard before FM0.1 loader/H200 admission.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-input-validate
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 47:30:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 027

readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly DESIGN_SHA="94c8672a087884bf8a2c70f5d15315e05de602134af0c6c2073ca1c36232d6f7"
readonly CONFIG_SHA="7de4e8c9e98c0ce27648a21241acc51766e436a537ba932b54f529fbdbf26d8a"
readonly FREEZE_SHA="75092235978c0b582f569be55770ad2d63368e079a6777da0b1c44547f074c25"
: "${TWIRL_ORCD_REPO:?controller must export detached repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?controller must export exact clean SHA}"
: "${TWIRL_FM0_RUN_ROOT:?controller must export run root}"
: "${TWIRL_FM0_MERGED_INPUT_RELEASE:?controller must export merged input release}"

readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly RUN_ROOT="$(realpath -m -- "${TWIRL_FM0_RUN_ROOT}")"
readonly RELEASE="$(realpath -m -- "${TWIRL_FM0_MERGED_INPUT_RELEASE}")"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
readonly OUT_DIR="${RUN_ROOT}/frozen/input_release_validation"
readonly VALIDATION="${OUT_DIR}/validation.json"
readonly RECEIPT="${OUT_DIR}/input_release.receipt.json"

die() { echo "[fm0-input-validate] ERROR: $*" >&2; exit 2; }
sha256() { sha256sum "$1" | awk '{print $1}'; }
for path in "${RUN_ROOT}" "${RELEASE}"; do
  [[ "${path}" == "${PROJECT_ROOT}"/* ]] || die "path escapes project root: ${path}"
done
[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || die "invalid expected Git SHA"
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" ]] || die "wrong Git SHA"
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is not clean"
[[ "$(sha256 "${REPO}/doc/foundation_model_design.md")" == "${DESIGN_SHA}" ]] || die "design drift"
[[ "$(sha256 "${REPO}/configs/models/twirl_fm0_1_s56_s67_poc.yaml")" == "${CONFIG_SHA}" ]] || die "science config drift"
[[ "$(sha256 "${REPO}/reports/stage5_validation/twirl_fm0_1_design_freeze_v1/freeze.json")" == "${FREEZE_SHA}" ]] || die "freeze receipt drift"
[[ -x "${PYTHON}" && -s "${RELEASE}/manifest.csv" ]] || die "runtime or merged release missing"
[[ ! -e "${OUT_DIR}" ]] || die "refusing to overwrite input validation evidence"

mkdir -p "${OUT_DIR}"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
"${PYTHON}" scripts/stage5_validation/validate_twirl_fm0_input_release.py \
  --release-dir "${RELEASE}" --workers 4 \
  --design "${REPO}/doc/foundation_model_design.md" \
  --config "${REPO}/configs/models/twirl_fm0_1_s56_s67_poc.yaml" \
  --freeze-receipt "${REPO}/reports/stage5_validation/twirl_fm0_1_design_freeze_v1/freeze.json" \
  > "${VALIDATION}"

TWIRL_FM0_RELEASE="${RELEASE}" \
TWIRL_FM0_VALIDATION="${VALIDATION}" \
TWIRL_FM0_RECEIPT="${RECEIPT}" \
TWIRL_FM0_EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}" \
TWIRL_FM0_DESIGN_SHA="${DESIGN_SHA}" \
TWIRL_FM0_CONFIG_SHA="${CONFIG_SHA}" \
TWIRL_FM0_FREEZE_SHA="${FREEZE_SHA}" \
"${PYTHON}" - <<'PY'
from __future__ import annotations
from datetime import datetime, timezone
import hashlib, json, os
from pathlib import Path

release = Path(os.environ["TWIRL_FM0_RELEASE"])
validation_path = Path(os.environ["TWIRL_FM0_VALIDATION"])
validation = json.loads(validation_path.read_text())
if validation.get("schema_version") != "twirl_fm0_1_input_release_validation_v1" or validation.get("passed") is not True:
    raise SystemExit("input release validation did not pass")
manifest = release / "manifest.csv"
summary = release / "summary.json"
merge = release / "merge_summary.json"
payload = {
    "schema_version": "twirl_fm0_1_input_release_receipt_v1",
    "created_at_utc": datetime.now(timezone.utc).isoformat(),
    "passed": True,
    "campaign_id": "twirl_fm0_1_s56_s67_poc",
    "scope": "validated_S56_S64_multisector_release_only",
    "expected_git_sha": os.environ["TWIRL_FM0_EXPECTED_SHA"],
    "science_bindings": {
        "design_sha256": os.environ["TWIRL_FM0_DESIGN_SHA"],
        "config_sha256": os.environ["TWIRL_FM0_CONFIG_SHA"],
        "freeze_receipt_sha256": os.environ["TWIRL_FM0_FREEZE_SHA"],
    },
    "input_adapter": "a2v1_hdf5_quality_aware_v1",
    "scientific_training_eligible": True,
    "certifies_full_campaign": False,
    "production_authorized": False,
    "foundation_model_claim_authorized": False,
    "release": {
        "manifest_path": str(manifest),
        "manifest_sha256": hashlib.sha256(manifest.read_bytes()).hexdigest(),
        "summary_path": str(summary),
        "summary_sha256": hashlib.sha256(summary.read_bytes()).hexdigest(),
        "merge_summary_path": str(merge),
        "merge_summary_sha256": hashlib.sha256(merge.read_bytes()).hexdigest(),
        "n_shards": validation["n_observations"],
        "n_cadences": validation["n_cadences"],
    },
    "independent_validation": {
        "path": str(validation_path),
        "sha256": hashlib.sha256(validation_path.read_bytes()).hexdigest(),
    },
    "claim_limit": "input construction and validation only; not a model result",
}
Path(os.environ["TWIRL_FM0_RECEIPT"]).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
for path in "${VALIDATION}" "${RECEIPT}"; do
  sha256sum "${path}" > "${path}.sha256"
done
chmod -R a-w "${OUT_DIR}"
echo "[fm0-input-validate] complete=${RECEIPT}"
