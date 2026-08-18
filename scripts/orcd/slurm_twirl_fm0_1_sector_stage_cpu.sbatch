#!/bin/bash
# Expand only the selected Gaia-source members from one checksum-bound A2v1
# sector tar. This is a CPU/I/O preparation step, not a model run.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-sector-stage
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 027

readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly DESIGN_SHA="94c8672a087884bf8a2c70f5d15315e05de602134af0c6c2073ca1c36232d6f7"
readonly CONFIG_SHA="7de4e8c9e98c0ce27648a21241acc51766e436a537ba932b54f529fbdbf26d8a"
readonly FREEZE_SHA="75092235978c0b582f569be55770ad2d63368e079a6777da0b1c44547f074c25"
: "${TWIRL_ORCD_REPO:?controller must export the detached repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?controller must export the exact clean SHA}"
: "${TWIRL_FM0_SECTOR:?controller must export one sector}"
: "${TWIRL_FM0_CORPUS_SELECTION:?controller must export the corpus selection}"
: "${TWIRL_FM0_CORPUS_SELECTION_SHA256:?controller must export its SHA-256}"
: "${TWIRL_FM0_SECTOR_ARCHIVE_DIR:?controller must export the verified sector archive directory}"
: "${TWIRL_FM0_QUALITY_TABLE:?controller must export the cadence-quality table}"
: "${TWIRL_FM0_QUALITY_MANIFEST:?controller must export the cadence-quality manifest}"
: "${TWIRL_FM0_SECTOR_STAGE_ROOT:?controller must export the immutable output root}"

readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly SELECTION="$(realpath -m -- "${TWIRL_FM0_CORPUS_SELECTION}")"
readonly ARCHIVE_DIR="$(realpath -m -- "${TWIRL_FM0_SECTOR_ARCHIVE_DIR}")"
readonly QUALITY_TABLE="$(realpath -m -- "${TWIRL_FM0_QUALITY_TABLE}")"
readonly QUALITY_MANIFEST="$(realpath -m -- "${TWIRL_FM0_QUALITY_MANIFEST}")"
readonly STAGE_ROOT="$(realpath -m -- "${TWIRL_FM0_SECTOR_STAGE_ROOT}")"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"

die() { echo "[fm0-sector-stage] ERROR: $*" >&2; exit 2; }
sha256() { sha256sum "$1" | awk '{print $1}'; }
[[ "${TWIRL_FM0_SECTOR}" =~ ^(5[6-9]|6[0-5])$ ]] || die "sector must be S56-S65"
[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || die "invalid expected Git SHA"
[[ "${TWIRL_FM0_CORPUS_SELECTION_SHA256}" =~ ^[0-9a-f]{64}$ ]] || die "invalid selection SHA-256"
for path in "${REPO}" "${SELECTION}" "${ARCHIVE_DIR}" "${QUALITY_TABLE}" "${QUALITY_MANIFEST}" "${STAGE_ROOT}"; do
  [[ "${path}" == "${PROJECT_ROOT}"/* ]] || die "path escapes project root: ${path}"
done
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" ]] || die "wrong Git SHA"
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is not clean"
[[ "$(sha256 "${REPO}/doc/foundation_model_design.md")" == "${DESIGN_SHA}" ]] || die "design drift"
[[ "$(sha256 "${REPO}/configs/models/twirl_fm0_1_s56_s67_poc.yaml")" == "${CONFIG_SHA}" ]] || die "science config drift"
[[ "$(sha256 "${REPO}/reports/stage5_validation/twirl_fm0_1_design_freeze_v1/freeze.json")" == "${FREEZE_SHA}" ]] || die "freeze receipt drift"
[[ -x "${PYTHON}" && -s "${SELECTION}" && -s "${QUALITY_TABLE}" && -s "${QUALITY_MANIFEST}" ]] || die "runtime or authority is missing"
[[ "$(sha256 "${SELECTION}")" == "${TWIRL_FM0_CORPUS_SELECTION_SHA256}" ]] || die "selection hash mismatch"
[[ -s "${ARCHIVE_DIR}/READY_ORCD" ]] || die "sector archive is not ORCD-ready"
[[ ! -e "${STAGE_ROOT}/s00${TWIRL_FM0_SECTOR}" && ! -e "${STAGE_ROOT}/s00${TWIRL_FM0_SECTOR}.partial" ]] || die "refusing to overwrite sector output"

mkdir -p "${STAGE_ROOT}"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
"${PYTHON}" scripts/stage5_validation/prepare_twirl_fm0_sector_archive.py \
  --sector "${TWIRL_FM0_SECTOR}" \
  --selection "${SELECTION}" \
  --selection-sha256 "${TWIRL_FM0_CORPUS_SELECTION_SHA256}" \
  --archive-dir "${ARCHIVE_DIR}" \
  --quality-table "${QUALITY_TABLE}" \
  --quality-manifest "${QUALITY_MANIFEST}" \
  --output-root "${STAGE_ROOT}" \
  --producer-git-sha "${TWIRL_EXPECTED_GIT_SHA}" \
  --workers "${SLURM_CPUS_PER_TASK:-4}"

chmod -R a-w "${STAGE_ROOT}/s00${TWIRL_FM0_SECTOR}"
echo "[fm0-sector-stage] complete=${STAGE_ROOT}/s00${TWIRL_FM0_SECTOR}/READY"
