#!/bin/bash
# Merge ready FM0.1 sector input releases without rereading raw light curves.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-input-merge
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 08:00:00
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
: "${TWIRL_FM0_REGISTRY_DIR:?controller must export fixed registry}"
: "${TWIRL_FM0_SECTOR_INPUT_ROOT:?controller must export sector-input root}"
: "${TWIRL_FM0_SECTOR_INPUT_MERGE_ROOT:?controller must export merged release root}"
: "${TWIRL_FM0_ALLOWED_INPUT_REVISIONS:?controller must export approved revisions}"

readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly REGISTRY="$(realpath -m -- "${TWIRL_FM0_REGISTRY_DIR}")"
readonly INPUT_ROOT="$(realpath -m -- "${TWIRL_FM0_SECTOR_INPUT_ROOT}")"
readonly MERGE_ROOT="$(realpath -m -- "${TWIRL_FM0_SECTOR_INPUT_MERGE_ROOT}")"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"

die() { echo "[fm0-input-merge] ERROR: $*" >&2; exit 2; }
sha256() { sha256sum "$1" | awk '{print $1}'; }
for path in "${REGISTRY}" "${INPUT_ROOT}" "${MERGE_ROOT}"; do
  [[ "${path}" == "${PROJECT_ROOT}"/* ]] || die "path escapes project root: ${path}"
done
[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || die "invalid expected Git SHA"
[[ "${TWIRL_FM0_ALLOWED_INPUT_REVISIONS}" =~ ^[0-9a-f]{40}(:[0-9a-f]{40})*$ ]] || die "invalid approved input revisions"
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" ]] || die "wrong Git SHA"
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is not clean"
[[ "$(sha256 "${REPO}/doc/foundation_model_design.md")" == "${DESIGN_SHA}" ]] || die "design drift"
[[ "$(sha256 "${REPO}/configs/models/twirl_fm0_1_s56_s67_poc.yaml")" == "${CONFIG_SHA}" ]] || die "science config drift"
[[ "$(sha256 "${REPO}/reports/stage5_validation/twirl_fm0_1_design_freeze_v1/freeze.json")" == "${FREEZE_SHA}" ]] || die "freeze receipt drift"
[[ -x "${PYTHON}" && -s "${REGISTRY}/summary.json" ]] || die "runtime or registry missing"
[[ ! -e "${MERGE_ROOT}" ]] || die "refusing to overwrite merged input release"

args=()
IFS=: read -r -a revisions <<< "${TWIRL_FM0_ALLOWED_INPUT_REVISIONS}"
for revision in "${revisions[@]}"; do args+=(--allowed-sector-producer-git-sha "${revision}"); done
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
"${PYTHON}" scripts/stage5_validation/merge_twirl_fm0_sector_input_releases.py \
  --sector-root "${INPUT_ROOT}" \
  --registry-dir "${REGISTRY}" \
  --from-sector 56 --through-sector 64 \
  --out-dir "${MERGE_ROOT}" \
  --producer-git-sha "${TWIRL_EXPECTED_GIT_SHA}" \
  --design "${REPO}/doc/foundation_model_design.md" \
  --config "${REPO}/configs/models/twirl_fm0_1_s56_s67_poc.yaml" \
  --freeze-receipt "${REPO}/reports/stage5_validation/twirl_fm0_1_design_freeze_v1/freeze.json" \
  "${args[@]}"
for path in "${MERGE_ROOT}/READY" "${MERGE_ROOT}/merge_summary.json" "${MERGE_ROOT}/summary.json" "${MERGE_ROOT}/manifest.csv" "${MERGE_ROOT}/visit_timing.csv"; do
  [[ -s "${path}" ]] || die "merged input release missing: ${path}"
done
chmod -R a-w "${MERGE_ROOT}"
echo "[fm0-input-merge] complete=${MERGE_ROOT}"
