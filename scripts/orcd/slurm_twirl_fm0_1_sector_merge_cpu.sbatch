#!/bin/bash
# Merge the completed immutable FM0.1 sector source inventories. This is a
# CPU-only provenance step; it does not construct model tensors or train.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-sector-merge
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
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
: "${TWIRL_FM0_SECTOR_STAGE_ROOT:?controller must export the sector-stage root}"
: "${TWIRL_FM0_SECTOR_MERGE_ROOT:?controller must export the immutable merge root}"
: "${TWIRL_FM0_ALLOWED_STAGE_REVISIONS:?controller must export colon-separated reviewed stage revisions}"

readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly STAGE_ROOT="$(realpath -m -- "${TWIRL_FM0_SECTOR_STAGE_ROOT}")"
readonly MERGE_ROOT="$(realpath -m -- "${TWIRL_FM0_SECTOR_MERGE_ROOT}")"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"

die() { echo "[fm0-sector-merge] ERROR: $*" >&2; exit 2; }
sha256() { sha256sum "$1" | awk '{print $1}'; }
[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || die "invalid expected Git SHA"
for path in "${REPO}" "${STAGE_ROOT}" "${MERGE_ROOT}"; do
  [[ "${path}" == "${PROJECT_ROOT}"/* ]] || die "path escapes project root: ${path}"
done
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" ]] || die "wrong Git SHA"
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is not clean"
[[ "$(sha256 "${REPO}/doc/foundation_model_design.md")" == "${DESIGN_SHA}" ]] || die "design drift"
[[ "$(sha256 "${REPO}/configs/models/twirl_fm0_1_s56_s67_poc.yaml")" == "${CONFIG_SHA}" ]] || die "science config drift"
[[ "$(sha256 "${REPO}/reports/stage5_validation/twirl_fm0_1_design_freeze_v1/freeze.json")" == "${FREEZE_SHA}" ]] || die "freeze receipt drift"
[[ -x "${PYTHON}" ]] || die "pinned FM0.1 runtime is missing"
[[ ! -e "${MERGE_ROOT}" ]] || die "refusing to overwrite merge output"

IFS=: read -r -a stage_revisions <<< "${TWIRL_FM0_ALLOWED_STAGE_REVISIONS}"
(( ${#stage_revisions[@]} > 0 )) || die "no allowed sector-stage revision"
merge_args=()
for revision in "${stage_revisions[@]}"; do
  [[ "${revision}" =~ ^[0-9a-f]{40}$ ]] || die "invalid allowed sector-stage revision"
  merge_args+=(--allowed-sector-producer-git-sha "${revision}")
done
for sector in $(seq 56 64); do
  sector_dir="${STAGE_ROOT}/s00${sector}"
  [[ -s "${sector_dir}/READY" && -s "${sector_dir}/summary.json" && -s "${sector_dir}/source_inventory.csv" ]] || die "S${sector} stage is incomplete"
done

cd "${REPO}"
export PYTHONPATH="${REPO}/src"
"${PYTHON}" scripts/stage5_validation/merge_twirl_fm0_sector_inventories.py \
  --sector-root "${STAGE_ROOT}" \
  --from-sector 56 \
  --through-sector 64 \
  --out-dir "${MERGE_ROOT}" \
  --producer-git-sha "${TWIRL_EXPECTED_GIT_SHA}" \
  "${merge_args[@]}"
for path in "${MERGE_ROOT}/source_inventory.csv" "${MERGE_ROOT}/summary.json"; do
  [[ -s "${path}" ]] || die "merge omitted ${path}"
done
printf '%s\n' "${TWIRL_EXPECTED_GIT_SHA}" > "${MERGE_ROOT}/READY"
sha256sum "${MERGE_ROOT}/source_inventory.csv" "${MERGE_ROOT}/summary.json" "${MERGE_ROOT}/READY" > "${MERGE_ROOT}/files.sha256"
chmod -R a-w "${MERGE_ROOT}"
echo "[fm0-sector-merge] complete=${MERGE_ROOT}/READY"
