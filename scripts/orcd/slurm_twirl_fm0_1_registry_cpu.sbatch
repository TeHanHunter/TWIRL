#!/bin/bash
# Build the Gaia--TIC leakage-component registry on a CPU node.
# This is a pre-freeze build step and does not certify an input release.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-registry
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 06:00:00
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
: "${TWIRL_FM0_RUN_ROOT:?controller must export the fresh run root}"
: "${TWIRL_FM0_ALIASES_TABLE:?stage the alias table below the project root}"
: "${TWIRL_FM0_SOURCE_ADAPTER:?stage the exact FM0.1 source adapter}"

readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly RUN_ROOT="$(realpath -m -- "${TWIRL_FM0_RUN_ROOT}")"
readonly ALIASES="$(realpath -m -- "${TWIRL_FM0_ALIASES_TABLE}")"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
readonly OUT="${RUN_ROOT}/input_build/registry"

die() { echo "[fm0-registry] ERROR: $*" >&2; exit 2; }
sha256() { sha256sum "$1" | awk '{print $1}'; }
for path in "${RUN_ROOT}" "${ALIASES}"; do
  [[ "${path}" == "${PROJECT_ROOT}"/* ]] || die "path escapes project root: ${path}"
done
[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || die "invalid expected Git SHA"
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" ]] || die "wrong Git SHA"
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is not clean"
[[ "$(sha256 "${REPO}/doc/foundation_model_design.md")" == "${DESIGN_SHA}" ]] || die "design drift"
[[ "$(sha256 "${REPO}/configs/models/twirl_fm0_1_s56_s67_poc.yaml")" == "${CONFIG_SHA}" ]] || die "science config drift"
[[ "$(sha256 "${REPO}/reports/stage5_validation/twirl_fm0_1_design_freeze_v1/freeze.json")" == "${FREEZE_SHA}" ]] || die "freeze receipt drift"
case "${TWIRL_FM0_SOURCE_ADAPTER}" in
  strict_npz_fixture_v1)
    : "${TWIRL_FM0_OBSERVATIONS_TABLE:?stage the admitted-observations table below the project root}"
    : "${TWIRL_FM0_SOURCE_MANIFEST:?stage the strict-NPZ source manifest below the project root}"
    readonly OBSERVATIONS="$(realpath -m -- "${TWIRL_FM0_OBSERVATIONS_TABLE}")"
    readonly SOURCE_MANIFEST="$(realpath -m -- "${TWIRL_FM0_SOURCE_MANIFEST}")"
    for path in "${OBSERVATIONS}" "${SOURCE_MANIFEST}"; do
      [[ "${path}" == "${PROJECT_ROOT}"/* ]] || die "path escapes project root: ${path}"
    done
    [[ -s "${OBSERVATIONS}" && -s "${SOURCE_MANIFEST}" ]] || die "fixture authorities are missing"
    args_source=(--observations "${OBSERVATIONS}")
    ;;
  a2v1_hdf5_quality_aware_v1)
    : "${TWIRL_FM0_A2V1_HDF5_SOURCE_INVENTORY:?stage the A2v1 HDF5 source inventory below the project root}"
    readonly SOURCE_INVENTORY="$(realpath -m -- "${TWIRL_FM0_A2V1_HDF5_SOURCE_INVENTORY}")"
    [[ "${SOURCE_INVENTORY}" == "${PROJECT_ROOT}"/* ]] || die "path escapes project root: ${SOURCE_INVENTORY}"
    [[ -s "${SOURCE_INVENTORY}" ]] || die "A2v1 HDF5 source inventory is missing"
    args_source=(--a2v1-hdf5-source-inventory "${SOURCE_INVENTORY}")
    ;;
  *) die "unsupported source adapter: ${TWIRL_FM0_SOURCE_ADAPTER}" ;;
esac
[[ -x "${PYTHON}" && -s "${ALIASES}" ]] || die "runtime or staged authority is missing"
[[ ! -e "${OUT}" ]] || die "refusing to overwrite registry output"

mkdir -p "${RUN_ROOT}/input_build"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
args=(
  --aliases "${ALIASES}"
  "${args_source[@]}"
  --out-dir "${OUT}"
  --design "${REPO}/doc/foundation_model_design.md"
  --config "${REPO}/configs/models/twirl_fm0_1_s56_s67_poc.yaml"
  --freeze-receipt "${REPO}/reports/stage5_validation/twirl_fm0_1_design_freeze_v1/freeze.json"
)
"${PYTHON}" scripts/stage5_validation/build_twirl_fm0_registry.py "${args[@]}"
for path in "${OUT}/components.csv" "${OUT}/aliases.csv" "${OUT}/observations.csv" "${OUT}/summary.json"; do
  [[ -s "${path}" ]] || die "registry builder omitted ${path}"
done
if [[ "${TWIRL_FM0_SOURCE_ADAPTER}" == "a2v1_hdf5_quality_aware_v1" ]]; then
  [[ -s "${OUT}/a2v1_hdf5_manifest.csv" ]] || die "registry builder omitted bound A2v1 HDF5 manifest"
fi
chmod -R a-w "${OUT}"
echo "[fm0-registry] complete=${OUT}"
