#!/bin/bash
# Build one immutable, provider-neutral S66-S77 FM cadence-quality reference.
# This is a CPU-only authority-normalization job, never an H200 model job.
# S66 resolves to SPOC mission quality; S67+ resolves to TICA in the Python
# contract and is verified again before the output is made read-only.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-quality-ref
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH --signal=B:TERM@600
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 027

readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
: "${TWIRL_ORCD_REPO:?controller must export the detached repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?controller must export the exact clean SHA}"
: "${TWIRL_FM0_SECTOR:?controller must export one sector from 66 through 77}"
: "${TWIRL_FM0_QUALITY_TRANSFER_ROOT:?controller must export the immutable quality-transfer root}"
: "${TWIRL_FM0_QUALITY_TRANSFER_MANIFEST_SHA256:?controller must export the authorized transfer-manifest SHA-256}"
: "${TWIRL_FM0_QUALITY_REFERENCE_ROOT:?controller must export the quality-reference output root}"

readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly TRANSFER="$(realpath -m -- "${TWIRL_FM0_QUALITY_TRANSFER_ROOT}")"
readonly OUT_ROOT="$(realpath -m -- "${TWIRL_FM0_QUALITY_REFERENCE_ROOT}")"
readonly OUT="${OUT_ROOT}/s$(printf '%04d' "${TWIRL_FM0_SECTOR}")"
readonly PARTIAL="${OUT}.partial"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
CHILD_PID=""
OUTPUTS_ABSENT_AT_ENTRY=0

die() { echo "[fm0-quality-reference] ERROR: $*" >&2; exit 2; }
sha256() { sha256sum "$1" | awk '{print $1}'; }
cleanup_owned_tree() {
  local path="$1"
  case "${path}" in
    "${OUT}"|"${PARTIAL}") ;;
    *) echo "[fm0-quality-reference] refusing unsafe cleanup path: ${path}" >&2; return 70 ;;
  esac
  [[ "${path}" == "${PROJECT_ROOT}"/* ]] || return 70
  if [[ -L "${path}" || ( -e "${path}" && ! -d "${path}" ) ]]; then
    echo "[fm0-quality-reference] refusing unsafe cleanup target: ${path}" >&2
    return 70
  fi
  if [[ -d "${path}" ]]; then
    chmod u+rwx -- "${path}"
    rm -rf -- "${path}"
  fi
}
cleanup_owned_outputs() {
  local status=0
  [[ "${OUTPUTS_ABSENT_AT_ENTRY}" == "1" ]] || return 70
  cleanup_owned_tree "${PARTIAL}" || status=$?
  cleanup_owned_tree "${OUT}" || status=$?
  return "${status}"
}
wait_child() {
  local status
  if wait "${CHILD_PID}"; then status=0; else status=$?; fi
  CHILD_PID=""
  return "${status}"
}
stop_child() {
  local signal="$1"
  if [[ -n "${CHILD_PID}" ]]; then
    kill -s "${signal}" "${CHILD_PID}" 2>/dev/null || true
    wait "${CHILD_PID}" 2>/dev/null || true
    CHILD_PID=""
  fi
}
on_error() {
  local status=$?
  trap - ERR TERM INT
  cleanup_owned_outputs || status=70
  exit "${status}"
}
on_term() {
  local signal="$1"
  local status=143
  [[ "${signal}" == "INT" ]] && status=130
  trap - ERR TERM INT
  stop_child "${signal}"
  cleanup_owned_outputs || status=70
  exit "${status}"
}

[[ "${TWIRL_FM0_SECTOR}" =~ ^(6[6-9]|7[0-7])$ ]] || die "sector must be S66-S77"
[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || die "invalid expected Git SHA"
[[ "${TWIRL_FM0_QUALITY_TRANSFER_MANIFEST_SHA256}" =~ ^[0-9a-f]{64}$ ]] || die "invalid transfer-manifest SHA-256"
for path in "${REPO}" "${TRANSFER}" "${OUT_ROOT}"; do
  [[ "${path}" == "${PROJECT_ROOT}"/* ]] || die "path escapes project root: ${path}"
done
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" ]] || die "wrong Git SHA"
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is not clean"
[[ -x "${PYTHON}" ]] || die "FM Python runtime is missing: ${PYTHON}"
[[ -s "${TRANSFER}/manifest.json" && -s "${TRANSFER}/READY" ]] || die "quality-transfer bundle is incomplete"
[[ "$(sha256 "${TRANSFER}/manifest.json")" == "${TWIRL_FM0_QUALITY_TRANSFER_MANIFEST_SHA256}" ]] || die "quality-transfer manifest hash drift"
[[ ! -e "${OUT}" && ! -L "${OUT}" && ! -e "${PARTIAL}" && ! -L "${PARTIAL}" ]] || die "refusing to overwrite sector output: ${OUT}"
OUTPUTS_ABSENT_AT_ENTRY=1

trap on_error ERR
trap 'on_term TERM' TERM
trap 'on_term INT' INT

mkdir -p "${OUT_ROOT}"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[fm0-quality-reference] host=$(hostname) sector=${TWIRL_FM0_SECTOR} start=$(date -Is)"
"${PYTHON}" scripts/stage5_validation/build_twirl_fm0_mission_quality_reference.py \
  --quality-transfer-root "${TRANSFER}" \
  --quality-transfer-manifest-sha256 "${TWIRL_FM0_QUALITY_TRANSFER_MANIFEST_SHA256}" \
  --sector "${TWIRL_FM0_SECTOR}" \
  --output-dir "${OUT}" \
  --producer-git-sha "${TWIRL_EXPECTED_GIT_SHA}" &
CHILD_PID=$!
wait_child

for path in "${OUT}/cadence_quality.csv" "${OUT}/manifest.json" "${OUT}/READY"; do
  [[ -s "${path}" ]] || die "quality-reference builder omitted ${path}"
done

"${PYTHON}" - "${OUT}" "${TWIRL_FM0_SECTOR}" "${TWIRL_EXPECTED_GIT_SHA}" <<'PY' &
from __future__ import annotations

import json
import sys
from pathlib import Path

from twirl.lightcurves.mission_quality_reference import (
    load_mission_quality_reference,
)

root = Path(sys.argv[1])
sector = int(sys.argv[2])
revision = sys.argv[3]
expected_provider = "spoc" if sector < 67 else "tica"
manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
if manifest.get("producer_git_sha") != revision:
    raise SystemExit("quality-reference producer revision drifted")
if manifest.get("mission_quality_provider") != expected_provider:
    raise SystemExit(
        f"provider mismatch: expected {expected_provider}, "
        f"observed {manifest.get('mission_quality_provider')!r}"
    )
reference = load_mission_quality_reference(
    reference_dir=root,
    sector=sector,
    expected_orbits=(2 * sector + 7, 2 * sector + 8),
)
reference.assert_unchanged()
for path in (root, *root.rglob("*")):
    if path.is_symlink() or path.stat().st_mode & 0o222:
        raise SystemExit(f"quality-reference output is not read-only: {path}")
PY
CHILD_PID=$!
wait_child

trap - ERR TERM INT
echo "[fm0-quality-reference] complete=$(date -Is) provider=$(if (( TWIRL_FM0_SECTOR < 67 )); then echo spoc; else echo tica; fi) output=${OUT}"
