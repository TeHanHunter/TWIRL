#!/bin/bash
# Build one immutable, BLS-free FM0.1 input sub-release on a CPU node.
# The controller serializes sectors deliberately: this is a low-resource,
# restartable data-preparation stage, never an H200 model job.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-input-sector
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
: "${TWIRL_FM0_SECTOR:?controller must export one sector}"
: "${TWIRL_FM0_REGISTRY_DIR:?controller must export the fixed registry}" 
: "${TWIRL_FM0_A2V1_HDF5_MANIFEST:?controller must export bound HDF5 manifest}"
: "${TWIRL_FM0_SECTOR_INPUT_ROOT:?controller must export sector-input root}"

readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly RUN_ROOT="$(realpath -m -- "${TWIRL_FM0_RUN_ROOT}")"
readonly REGISTRY="$(realpath -m -- "${TWIRL_FM0_REGISTRY_DIR}")"
readonly SOURCE_MANIFEST="$(realpath -m -- "${TWIRL_FM0_A2V1_HDF5_MANIFEST}")"
readonly OUT_ROOT="$(realpath -m -- "${TWIRL_FM0_SECTOR_INPUT_ROOT}")"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"

die() { echo "[fm0-sector-input] ERROR: $*" >&2; exit 2; }
sha256() { sha256sum "$1" | awk '{print $1}'; }
[[ "${TWIRL_FM0_SECTOR}" =~ ^(5[6-9]|6[0-4])$ ]] || die "sector must be S56-S64"
for path in "${RUN_ROOT}" "${REGISTRY}" "${SOURCE_MANIFEST}" "${OUT_ROOT}"; do
  [[ "${path}" == "${PROJECT_ROOT}"/* ]] || die "path escapes project root: ${path}"
done
[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || die "invalid expected Git SHA"
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" ]] || die "wrong Git SHA"
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is not clean"
[[ "$(sha256 "${REPO}/doc/foundation_model_design.md")" == "${DESIGN_SHA}" ]] || die "design drift"
[[ "$(sha256 "${REPO}/configs/models/twirl_fm0_1_s56_s67_poc.yaml")" == "${CONFIG_SHA}" ]] || die "science config drift"
[[ "$(sha256 "${REPO}/reports/stage5_validation/twirl_fm0_1_design_freeze_v1/freeze.json")" == "${FREEZE_SHA}" ]] || die "freeze receipt drift"
[[ -x "${PYTHON}" && -s "${SOURCE_MANIFEST}" && -s "${REGISTRY}/summary.json" ]] || die "runtime or bound authority missing"
[[ ! -e "${OUT_ROOT}/s$(printf '%04d' "${TWIRL_FM0_SECTOR}")" ]] || die "refusing to overwrite sector input"

cd "${REPO}"
export PYTHONPATH="${REPO}/src"
"${PYTHON}" scripts/stage5_validation/build_twirl_fm0_sector_input_release.py \
  --sector "${TWIRL_FM0_SECTOR}" \
  --registry-dir "${REGISTRY}" \
  --a2v1-hdf5-manifest "${SOURCE_MANIFEST}" \
  --out-root "${OUT_ROOT}" \
  --producer-git-sha "${TWIRL_EXPECTED_GIT_SHA}" \
  --design "${REPO}/doc/foundation_model_design.md" \
  --config "${REPO}/configs/models/twirl_fm0_1_s56_s67_poc.yaml" \
  --freeze-receipt "${REPO}/reports/stage5_validation/twirl_fm0_1_design_freeze_v1/freeze.json"

readonly OUT="${OUT_ROOT}/s$(printf '%04d' "${TWIRL_FM0_SECTOR}")"
for path in "${OUT}/READY" "${OUT}/sector_summary.json" "${OUT}/summary.json" "${OUT}/manifest.csv" "${OUT}/visit_timing.csv" "${OUT}/a2v1_hdf5_manifest.csv"; do
  [[ -s "${path}" ]] || die "sector input output missing: ${path}"
done
"${PYTHON}" - "${OUT}" "${TWIRL_FM0_SECTOR}" "${TWIRL_EXPECTED_GIT_SHA}" <<'PY'
from __future__ import annotations
import csv, hashlib, json, sys
from pathlib import Path

root, sector, revision = Path(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
summary = json.loads((root / "sector_summary.json").read_text())
manifest = list(csv.DictReader((root / "manifest.csv").open(newline="")))
timing = list(csv.DictReader((root / "visit_timing.csv").open(newline="")))
assert summary["schema_version"] == "twirl_fm0_1_sector_input_release_v1"
assert summary["producer_git_sha"] == revision and int(summary["sector"]) == sector
assert summary["n_observations"] == len(manifest) == len(timing) > 0
assert summary["scientific_training_eligible"] is True
assert summary["input_adapter"] == "a2v1_hdf5_quality_aware_v1"
assert {row["observation_key"] for row in manifest} == {row["observation_key"] for row in timing}
for row in manifest:
    shard = root / row["relative_path"]
    assert shard.is_file()
    assert hashlib.sha256(shard.read_bytes()).hexdigest() == row["sha256"]
PY
chmod -R a-w "${OUT}"
echo "[fm0-sector-input] complete=${OUT}"
