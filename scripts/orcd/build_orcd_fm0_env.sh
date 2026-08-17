#!/usr/bin/env bash
# Build one hash-bound FM0.1 CPU/H200 environment on ORCD project storage.
# Existing prefixes are never updated in place: a completed manifest is
# validated and reused, while an incomplete prefix requires a new version.
set -euo pipefail
umask 027

readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly ENV_PREFIX="${TWIRL_FM0_ENV_PREFIX:-${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1}"
readonly REPO="${TWIRL_ORCD_REPO:-${PROJECT_ROOT}/code/TWIRL}"
readonly ENV_YML="${REPO}/configs/orcd/twirl-fm0-poc-env.yml"
readonly MANIFEST="${ENV_PREFIX}/build-manifest.json"
readonly CONDA_EXPLICIT="${ENV_PREFIX}/conda-explicit.txt"
readonly PIP_FREEZE="${ENV_PREFIX}/pip-freeze.txt"
readonly TORCH_INDEX_URL="https://download.pytorch.org/whl/cu128"
readonly TORCH_VERSION="2.11.0"
readonly EXPECTED_ENV_YML_SHA256="1cfe39eb1946fc4a9bc72ad3130490db9da994cae1e4083a911ed1247a3fa532"
readonly EXPECTED_ORCD_CONFIG_SHA256="733a5c65118d6331979a61e7f63b1dfb9b3e63eff0aa2052c98e6c16290f938f"
readonly EXPECTED_DESIGN_SHA256="94c8672a087884bf8a2c70f5d15315e05de602134af0c6c2073ca1c36232d6f7"
readonly EXPECTED_SCIENCE_CONFIG_SHA256="7de4e8c9e98c0ce27648a21241acc51766e436a537ba932b54f529fbdbf26d8a"
readonly EXPECTED_FREEZE_SHA256="75092235978c0b582f569be55770ad2d63368e079a6777da0b1c44547f074c25"

die() { echo "[fm0-env] ERROR: $*" >&2; exit 2; }
sha256() { sha256sum "$1" | awk '{print $1}'; }

: "${TWIRL_EXPECTED_GIT_SHA:?set the exact reviewed 40-hex commit}"
readonly EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}"
[[ "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]] || die "TWIRL_EXPECTED_GIT_SHA must be full lowercase hex"
[[ "$(realpath -m -- "${ENV_PREFIX}")" == "${PROJECT_ROOT}"/* ]] || die "environment must stay below ${PROJECT_ROOT}"
[[ -d "${REPO}/.git" || -f "${REPO}/.git" ]] || die "repository is not a Git worktree: ${REPO}"
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${EXPECTED_SHA}" ]] || die "repository commit differs from TWIRL_EXPECTED_GIT_SHA"
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is not clean"

declare -A bindings=(
  ["${ENV_YML}"]="${EXPECTED_ENV_YML_SHA256}"
  ["${REPO}/configs/orcd/twirl_fm0_1_s56_s67_poc.json"]="${EXPECTED_ORCD_CONFIG_SHA256}"
  ["${REPO}/doc/foundation_model_design.md"]="${EXPECTED_DESIGN_SHA256}"
  ["${REPO}/configs/models/twirl_fm0_1_s56_s67_poc.yaml"]="${EXPECTED_SCIENCE_CONFIG_SHA256}"
  ["${REPO}/reports/stage5_validation/twirl_fm0_1_design_freeze_v1/freeze.json"]="${EXPECTED_FREEZE_SHA256}"
)
for path in "${!bindings[@]}"; do
  [[ -f "${path}" ]] || die "missing bound file: ${path}"
  [[ "$(sha256 "${path}")" == "${bindings[${path}]}" ]] || die "bound file drift: ${path}"
done

validate_complete() {
  [[ -x "${ENV_PREFIX}/bin/python" ]] || die "completed environment lacks Python"
  [[ -s "${CONDA_EXPLICIT}" && -s "${PIP_FREEZE}" ]] || die "dependency evidence is incomplete"
  "${ENV_PREFIX}/bin/python" - "${MANIFEST}" "${EXPECTED_SHA}" \
    "${EXPECTED_ENV_YML_SHA256}" "${EXPECTED_ORCD_CONFIG_SHA256}" <<'PY'
import hashlib, json, sys
from pathlib import Path
import numpy, h5py, pandas, pyarrow, torch, yaml

manifest = json.loads(Path(sys.argv[1]).read_text())
assert manifest["schema_version"] == "twirl_fm0_1_orcd_environment_v1"
assert manifest["expected_git_sha"] == sys.argv[2]
assert manifest["environment_yml_sha256"] == sys.argv[3]
assert manifest["orcd_config_sha256"] == sys.argv[4]
assert manifest["torch_requested_version"] == "2.11.0"
assert torch.__version__.split("+")[0] == "2.11.0"
assert torch.version.cuda == "12.8"
for key in ("conda_explicit", "pip_freeze"):
    path = Path(manifest[key]["path"])
    assert path.is_file()
    assert hashlib.sha256(path.read_bytes()).hexdigest() == manifest[key]["sha256"]
print("[fm0-env] immutable environment validation passed")
PY
}

if [[ -f "${MANIFEST}" ]]; then
  validate_complete
  exit 0
fi
[[ ! -e "${ENV_PREFIX}" ]] || die "prefix exists without a completed manifest; choose a fresh versioned prefix"

if command -v module >/dev/null 2>&1; then
  module load miniforge/25.11.0-0 || module load miniforge/24.3.0-0 || true
fi
if command -v mamba >/dev/null 2>&1; then
  CONDA_FRONTEND=mamba
elif command -v conda >/dev/null 2>&1; then
  CONDA_FRONTEND=conda
else
  die "mamba/conda is unavailable"
fi

"${CONDA_FRONTEND}" env create --yes --prefix "${ENV_PREFIX}" --file "${ENV_YML}"
"${ENV_PREFIX}/bin/python" -m pip install --no-cache-dir \
  --index-url "${TORCH_INDEX_URL}" "torch==${TORCH_VERSION}"
"${ENV_PREFIX}/bin/python" -m pip install --no-deps --no-build-isolation "${REPO}"
"${CONDA_FRONTEND}" list --explicit --prefix "${ENV_PREFIX}" > "${CONDA_EXPLICIT}"
"${ENV_PREFIX}/bin/python" -m pip freeze > "${PIP_FREEZE}"

TWIRL_FM0_MANIFEST="${MANIFEST}" \
TWIRL_FM0_ENV_PREFIX="${ENV_PREFIX}" \
TWIRL_FM0_EXPECTED_SHA="${EXPECTED_SHA}" \
TWIRL_FM0_ENV_YML_SHA="${EXPECTED_ENV_YML_SHA256}" \
TWIRL_FM0_ORCD_CONFIG_SHA="${EXPECTED_ORCD_CONFIG_SHA256}" \
TWIRL_FM0_CONDA_EXPLICIT="${CONDA_EXPLICIT}" \
TWIRL_FM0_PIP_FREEZE="${PIP_FREEZE}" \
TWIRL_FM0_TORCH_VERSION="${TORCH_VERSION}" \
"${ENV_PREFIX}/bin/python" - <<'PY'
from datetime import datetime, timezone
import hashlib, json, os
from pathlib import Path
import torch

def item(path):
    p = Path(path)
    return {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}

target = Path(os.environ["TWIRL_FM0_MANIFEST"])
payload = {
    "schema_version": "twirl_fm0_1_orcd_environment_v1",
    "created_at_utc": datetime.now(timezone.utc).isoformat(),
    "environment_prefix": os.environ["TWIRL_FM0_ENV_PREFIX"],
    "expected_git_sha": os.environ["TWIRL_FM0_EXPECTED_SHA"],
    "environment_yml_sha256": os.environ["TWIRL_FM0_ENV_YML_SHA"],
    "orcd_config_sha256": os.environ["TWIRL_FM0_ORCD_CONFIG_SHA"],
    "torch_requested_version": os.environ["TWIRL_FM0_TORCH_VERSION"],
    "torch_observed_version": torch.__version__,
    "torch_cuda_build": torch.version.cuda,
    "conda_explicit": item(os.environ["TWIRL_FM0_CONDA_EXPLICIT"]),
    "pip_freeze": item(os.environ["TWIRL_FM0_PIP_FREEZE"]),
}
tmp = target.with_suffix(".json.tmp")
tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
tmp.replace(target)
PY
chmod a-w "${MANIFEST}" "${CONDA_EXPLICIT}" "${PIP_FREEZE}"
validate_complete
echo "[fm0-env] complete: ${ENV_PREFIX}"
