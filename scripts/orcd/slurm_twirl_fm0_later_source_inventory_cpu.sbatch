#!/bin/bash
# Build one immutable, label-blind retained-source inventory for S66-S77.
# This CPU-only job verifies metadata and hashes; it never opens HDF5 content.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-source-inventory
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 04:00:00
#SBATCH --signal=B:TERM@600
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 027

readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
: "${TWIRL_ORCD_REPO:?controller must export the detached repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?controller must export the exact clean SHA}"
: "${TWIRL_FM0_SECTOR:?controller must export one sector from 66 through 77}"
: "${TWIRL_FM0_SOURCE_RECEIPT:?controller must export the retained-source receipt}"
: "${TWIRL_FM0_SOURCE_RECEIPT_SHA256:?controller must export the retained-source receipt SHA-256}"
: "${TWIRL_FM0_CORPUS_SELECTION:?controller must export the exact corpus-selection CSV}"
: "${TWIRL_FM0_CORPUS_SELECTION_SHA256:?controller must export the corpus-selection SHA-256}"
: "${TWIRL_FM0_CORPUS_SUMMARY:?controller must export the exact corpus-selection summary}"
: "${TWIRL_FM0_CORPUS_SUMMARY_SHA256:?controller must export the corpus-summary SHA-256}"
: "${TWIRL_FM0_GAIA_TIC_ALIAS_AUTHORITY:?controller must export the exact Gaia--TIC alias table}"
: "${TWIRL_FM0_GAIA_TIC_ALIAS_AUTHORITY_SHA256:?controller must export the Gaia--TIC alias-table SHA-256}"
: "${TWIRL_FM0_TARGET_AUTHORITY_SHA256:?controller must export the observation-authority SHA-256}"
: "${TWIRL_FM0_SOURCE_INVENTORY_ROOT:?controller must export the source-inventory output root}"

readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly SOURCE_RECEIPT="$(realpath -m -- "${TWIRL_FM0_SOURCE_RECEIPT}")"
readonly CORPUS_SELECTION="$(realpath -m -- "${TWIRL_FM0_CORPUS_SELECTION}")"
readonly CORPUS_SUMMARY="$(realpath -m -- "${TWIRL_FM0_CORPUS_SUMMARY}")"
readonly GAIA_TIC_ALIAS_AUTHORITY="$(realpath -m -- "${TWIRL_FM0_GAIA_TIC_ALIAS_AUTHORITY}")"
readonly OUT_ROOT="$(realpath -m -- "${TWIRL_FM0_SOURCE_INVENTORY_ROOT}")"
readonly OUT="${OUT_ROOT}/s$(printf '%04d' "${TWIRL_FM0_SECTOR}")"
readonly PARTIAL="${OUT}.partial"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
CHILD_PID=""
OUTPUTS_ABSENT_AT_ENTRY=0

die() { echo "[fm0-source-inventory] ERROR: $*" >&2; exit 2; }
sha256() { sha256sum "$1" | awk '{print $1}'; }
cleanup_owned_tree() {
  local path="$1"
  case "${path}" in
    "${OUT}"|"${PARTIAL}") ;;
    *) echo "[fm0-source-inventory] refusing unsafe cleanup path: ${path}" >&2; return 70 ;;
  esac
  [[ "${path}" == "${PROJECT_ROOT}"/* ]] || return 70
  if [[ -L "${path}" || ( -e "${path}" && ! -d "${path}" ) ]]; then
    echo "[fm0-source-inventory] refusing unsafe cleanup target: ${path}" >&2
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
[[ "${TWIRL_FM0_SOURCE_RECEIPT_SHA256}" =~ ^[0-9a-f]{64}$ ]] || die "invalid source-receipt SHA-256"
[[ "${TWIRL_FM0_CORPUS_SELECTION_SHA256}" =~ ^[0-9a-f]{64}$ ]] || die "invalid corpus-selection SHA-256"
[[ "${TWIRL_FM0_CORPUS_SUMMARY_SHA256}" =~ ^[0-9a-f]{64}$ ]] || die "invalid corpus-summary SHA-256"
[[ "${TWIRL_FM0_GAIA_TIC_ALIAS_AUTHORITY_SHA256}" =~ ^[0-9a-f]{64}$ ]] || die "invalid Gaia--TIC alias-authority SHA-256"
[[ "${TWIRL_FM0_TARGET_AUTHORITY_SHA256}" =~ ^[0-9a-f]{64}$ ]] || die "invalid target-authority SHA-256"
readonly ORBIT_1="$((2 * TWIRL_FM0_SECTOR + 7))"
readonly ORBIT_2="$((2 * TWIRL_FM0_SECTOR + 8))"
for path in "${REPO}" "${SOURCE_RECEIPT}" "${CORPUS_SELECTION}" "${CORPUS_SUMMARY}" "${GAIA_TIC_ALIAS_AUTHORITY}" "${OUT_ROOT}"; do
  [[ "${path}" == "${PROJECT_ROOT}"/* ]] || die "path escapes project root: ${path}"
done
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" ]] || die "wrong Git SHA"
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is not clean"
[[ -x "${PYTHON}" && -s "${SOURCE_RECEIPT}" && -s "${CORPUS_SELECTION}" && -s "${CORPUS_SUMMARY}" && -s "${GAIA_TIC_ALIAS_AUTHORITY}" ]] || die "runtime or authority is missing"
[[ "$(sha256 "${SOURCE_RECEIPT}")" == "${TWIRL_FM0_SOURCE_RECEIPT_SHA256}" ]] || die "source-receipt hash drift"
[[ "$(sha256 "${CORPUS_SELECTION}")" == "${TWIRL_FM0_CORPUS_SELECTION_SHA256}" ]] || die "corpus-selection hash drift"
[[ "$(sha256 "${CORPUS_SUMMARY}")" == "${TWIRL_FM0_CORPUS_SUMMARY_SHA256}" ]] || die "corpus-summary hash drift"
[[ "$(sha256 "${GAIA_TIC_ALIAS_AUTHORITY}")" == "${TWIRL_FM0_GAIA_TIC_ALIAS_AUTHORITY_SHA256}" ]] || die "Gaia--TIC alias-authority hash drift"
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

echo "[fm0-source-inventory] host=$(hostname) sector=${TWIRL_FM0_SECTOR} start=$(date -Is)"
"${PYTHON}" scripts/stage5_validation/build_twirl_fm0_later_source_inventory.py \
  --source-receipt "${SOURCE_RECEIPT}" \
  --source-receipt-sha256 "${TWIRL_FM0_SOURCE_RECEIPT_SHA256}" \
  --sector "${TWIRL_FM0_SECTOR}" \
  --expected-orbit "${ORBIT_1}" \
  --expected-orbit "${ORBIT_2}" \
  --expected-target-authority-sha256 "${TWIRL_FM0_TARGET_AUTHORITY_SHA256}" \
  --target-inventory "${CORPUS_SELECTION}" \
  --target-inventory-sha256 "${TWIRL_FM0_CORPUS_SELECTION_SHA256}" \
  --corpus-summary "${CORPUS_SUMMARY}" \
  --corpus-summary-sha256 "${TWIRL_FM0_CORPUS_SUMMARY_SHA256}" \
  --gaia-tic-alias-authority "${GAIA_TIC_ALIAS_AUTHORITY}" \
  --gaia-tic-alias-authority-sha256 "${TWIRL_FM0_GAIA_TIC_ALIAS_AUTHORITY_SHA256}" \
  --output-dir "${OUT}" \
  --producer-git-sha "${TWIRL_EXPECTED_GIT_SHA}" &
CHILD_PID=$!
wait_child

for path in "${OUT}/aliases.csv" "${OUT}/sources.csv" "${OUT}/summary.json" "${OUT}/READY"; do
  [[ -s "${path}" ]] || die "source-inventory builder omitted ${path}"
done

readonly SUMMARY_SHA256="$(sha256 "${OUT}/summary.json")"
"${PYTHON}" - "${OUT}" "${TWIRL_FM0_SECTOR}" "${TWIRL_EXPECTED_GIT_SHA}" "${SUMMARY_SHA256}" <<'PY' &
from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

from twirl.models.fm0.later_source_inventory import load_later_source_rows

root = Path(sys.argv[1])
sector = int(sys.argv[2])
revision = sys.argv[3]
summary_sha256 = sys.argv[4]
summary_bytes = (root / "summary.json").read_bytes()
summary = json.loads(summary_bytes)
with (root / "sources.csv").open(newline="", encoding="utf-8") as handle:
    sources = list(csv.DictReader(handle))
if summary.get("producer_git_sha") != revision or int(summary.get("sector", -1)) != sector:
    raise SystemExit("source-inventory producer/sector binding drifted")
if summary.get("n_source_rows") != len(sources) or not sources:
    raise SystemExit("source-inventory row count drifted")
if summary.get("hdf5_content_opened") is not False:
    raise SystemExit("inventory unexpectedly claims HDF5 content access")
if summary.get("sealed_hdf5_content_opened") is not False:
    raise SystemExit("inventory unexpectedly claims sealed HDF5 access")
if summary.get("six_view_shards_verified") is not False:
    raise SystemExit("inventory overclaims six-view verification")
if summary.get("panel_admission_authorized") is not False:
    raise SystemExit("inventory overclaims panel admission")
if summary.get("gaia_tic_alias_authority_file_verified") is not True:
    raise SystemExit("Gaia--TIC alias authority was not independently verified")
if (root / "READY").read_text(encoding="utf-8").strip() != hashlib.sha256(summary_bytes).hexdigest():
    raise SystemExit("source-inventory READY marker drifted")
usable = load_later_source_rows(
    root,
    expected_summary_sha256=summary_sha256,
    allowed_source_partitions=("poc_train", "poc_development"),
)
expected_usable = sum(
    int(summary["source_rows_by_partition"][partition])
    for partition in ("poc_train", "poc_development")
)
if len(usable) != expected_usable:
    raise SystemExit("strict nonsealed source inventory count drifted")
for path in (root, *root.rglob("*")):
    if path.is_symlink() or path.stat().st_mode & 0o222:
        raise SystemExit(f"source-inventory output is not read-only: {path}")
PY
CHILD_PID=$!
wait_child

trap - ERR TERM INT
echo "[fm0-source-inventory] complete=$(date -Is) output=${OUT}"
