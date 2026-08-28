#!/bin/bash
# Build one immutable, full-visit six-view checkpoint for S66-S77.
# This is CPU-only data preparation.  It excludes sealed-test HDF5 content and
# does not authorize a temporal panel or model training.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-six-view
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=32G
#SBATCH -t 2-00:00:00
#SBATCH --signal=B:TERM@600
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 027

readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
: "${TWIRL_ORCD_REPO:?controller must export the detached repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?controller must export the exact clean SHA}"
: "${TWIRL_FM0_SECTOR:?controller must export one sector from 66 through 77}"
: "${TWIRL_FM0_SOURCE_INVENTORY_DIR:?controller must export the sector source-inventory directory}"
: "${TWIRL_FM0_SOURCE_INVENTORY_SUMMARY_SHA256:?controller must export the exact source-inventory summary SHA-256}"
: "${TWIRL_FM0_QUALITY_REFERENCE_DIR:?controller must export the sector mission-quality reference directory}"
: "${TWIRL_FM0_QUALITY_REFERENCE_MANIFEST_SHA256:?controller must export the exact quality-reference manifest SHA-256}"
: "${TWIRL_FM0_HDF5_QUALITY_RECEIPT:?controller must export the sector HDF5-quality receipt}"
: "${TWIRL_FM0_HDF5_QUALITY_RECEIPT_SHA256:?controller must export the exact HDF5-quality receipt SHA-256}"
: "${TWIRL_FM0_SIX_VIEW_RELEASE_DIR:?controller must export the final sector six-view directory}"

readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly SOURCE_INVENTORY="$(realpath -m -- "${TWIRL_FM0_SOURCE_INVENTORY_DIR}")"
readonly QUALITY_REFERENCE="$(realpath -m -- "${TWIRL_FM0_QUALITY_REFERENCE_DIR}")"
readonly HDF5_RECEIPT="$(realpath -m -- "${TWIRL_FM0_HDF5_QUALITY_RECEIPT}")"
readonly OUT_PARENT="$(realpath -m -- "$(dirname -- "${TWIRL_FM0_SIX_VIEW_RELEASE_DIR}")")"
readonly OUT_LEAF="$(basename -- "${TWIRL_FM0_SIX_VIEW_RELEASE_DIR}")"
readonly OUT="${OUT_PARENT}/${OUT_LEAF}"
readonly PARTIAL="${OUT}.partial"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
CHILD_PID=""
OUTPUTS_ABSENT_AT_ENTRY=0

die() { echo "[fm0-later-six-view] ERROR: $*" >&2; exit 2; }
cleanup_owned_tree() {
  local path="$1"
  case "${path}" in
    "${OUT}"|"${PARTIAL}") ;;
    *) echo "[fm0-later-six-view] refusing unsafe cleanup path: ${path}" >&2; return 70 ;;
  esac
  [[ "${path}" == "${PROJECT_ROOT}"/* ]] || return 70
  if [[ -L "${path}" || ( -e "${path}" && ! -d "${path}" ) ]]; then
    echo "[fm0-later-six-view] refusing unsafe cleanup target: ${path}" >&2
    return 70
  fi
  if [[ -d "${path}" ]]; then
    # The builder owns exactly one nested directory, shards/.  Restoring only
    # directory write bits avoids chmod-ing ~200k files before unlinking them.
    chmod u+rwx -- "${path}"
    if [[ -e "${path}/shards" ]]; then
      [[ -d "${path}/shards" && ! -L "${path}/shards" ]] || return 70
      chmod u+rwx -- "${path}/shards"
    fi
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
[[ -n "${OUT_LEAF}" && "${OUT_LEAF}" != "." && "${OUT_LEAF}" != ".." ]] || die "output must name a concrete directory"
for digest in \
  "${TWIRL_FM0_SOURCE_INVENTORY_SUMMARY_SHA256}" \
  "${TWIRL_FM0_QUALITY_REFERENCE_MANIFEST_SHA256}" \
  "${TWIRL_FM0_HDF5_QUALITY_RECEIPT_SHA256}"; do
  [[ "${digest}" =~ ^[0-9a-f]{64}$ ]] || die "invalid caller-authorized SHA-256"
done
for path in "${REPO}" "${SOURCE_INVENTORY}" "${QUALITY_REFERENCE}" "${HDF5_RECEIPT}" "${OUT}"; do
  [[ "${path}" == "${PROJECT_ROOT}"/* ]] || die "path escapes project root: ${path}"
done
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" ]] || die "wrong Git SHA"
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is not clean"
[[ -x "${PYTHON}" ]] || die "FM Python runtime is missing: ${PYTHON}"
for path in \
  "${SOURCE_INVENTORY}/sources.csv" \
  "${SOURCE_INVENTORY}/summary.json" \
  "${SOURCE_INVENTORY}/READY" \
  "${QUALITY_REFERENCE}/cadence_quality.csv" \
  "${QUALITY_REFERENCE}/manifest.json" \
  "${QUALITY_REFERENCE}/READY" \
  "${HDF5_RECEIPT}"; do
  [[ -s "${path}" ]] || die "required authority is missing: ${path}"
done
[[ ! -e "${OUT}" && ! -L "${OUT}" && ! -e "${PARTIAL}" && ! -L "${PARTIAL}" ]] || die "refusing to overwrite sector output: ${OUT}"
OUTPUTS_ABSENT_AT_ENTRY=1
[[ "$(sha256sum "${SOURCE_INVENTORY}/summary.json" | awk '{print $1}')" == "${TWIRL_FM0_SOURCE_INVENTORY_SUMMARY_SHA256}" ]] || die "source-inventory summary SHA-256 drifted"
[[ "$(sha256sum "${QUALITY_REFERENCE}/manifest.json" | awk '{print $1}')" == "${TWIRL_FM0_QUALITY_REFERENCE_MANIFEST_SHA256}" ]] || die "quality-reference manifest SHA-256 drifted"
[[ "$(sha256sum "${HDF5_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_HDF5_QUALITY_RECEIPT_SHA256}" ]] || die "HDF5-quality receipt SHA-256 drifted"

trap on_error ERR
trap 'on_term TERM' TERM
trap 'on_term INT' INT

mkdir -p "$(dirname "${OUT}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[fm0-later-six-view] host=$(hostname) sector=${TWIRL_FM0_SECTOR} start=$(date -Is)"
"${PYTHON}" scripts/stage5_validation/build_twirl_fm0_later_sector_release.py \
  --sector "${TWIRL_FM0_SECTOR}" \
  --source-inventory-dir "${SOURCE_INVENTORY}" \
  --source-inventory-summary-sha256 "${TWIRL_FM0_SOURCE_INVENTORY_SUMMARY_SHA256}" \
  --mission-quality-reference-dir "${QUALITY_REFERENCE}" \
  --mission-quality-reference-manifest-sha256 "${TWIRL_FM0_QUALITY_REFERENCE_MANIFEST_SHA256}" \
  --hdf5-quality-receipt "${HDF5_RECEIPT}" \
  --hdf5-quality-receipt-sha256 "${TWIRL_FM0_HDF5_QUALITY_RECEIPT_SHA256}" \
  --output-dir "${OUT}" \
  --producer-git-sha "${TWIRL_EXPECTED_GIT_SHA}" &
CHILD_PID=$!
wait_child

for path in \
  "${OUT}/manifest.csv" \
  "${OUT}/source_manifest.csv" \
  "${OUT}/visit_timing.csv" \
  "${OUT}/receipt.json" \
  "${OUT}/READY"; do
  [[ -s "${path}" ]] || die "six-view builder omitted ${path}"
done

"${PYTHON}" - "${OUT}" "${TWIRL_FM0_SECTOR}" "${TWIRL_EXPECTED_GIT_SHA}" <<'PY' &
from __future__ import annotations

import sys
from pathlib import Path

from twirl.models.fm0.later_sector_release import validate_later_sector_release

root = Path(sys.argv[1])
sector = int(sys.argv[2])
revision = sys.argv[3]
_receipt_path, receipt, summary = validate_later_sector_release(
    root, require_read_only=True
)
if receipt.get("producer_git_sha") != revision or int(receipt.get("sector", -1)) != sector:
    raise SystemExit("six-view producer/sector binding drifted")
if summary["n_observations"] != receipt["n_observations"]:
    raise SystemExit("six-view validation summary drifted")
PY
CHILD_PID=$!
wait_child

trap - ERR TERM INT
echo "[fm0-later-six-view] complete=$(date -Is) output=${OUT}"
