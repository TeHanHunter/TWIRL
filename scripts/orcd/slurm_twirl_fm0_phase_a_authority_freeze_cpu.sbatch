#!/bin/bash
# Fully validate and freeze the exact S66-S77 Phase-A/HDF5 authority record.
# This is CPU-only data preparation; it never opens sealed HDF5 content.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-phase-a-freeze
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 08:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 027

readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
: "${TWIRL_ORCD_REPO:?controller must export the detached repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?controller must export the exact clean SHA}"
: "${TWIRL_FM0_PREPARATION_AUTHORITY_MAP_SHA256:?controller must export the map SHA-256}"
: "${TWIRL_FM0_PHASE_A_AUTHORITY_RECORD:?controller must export the output record}"

readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly MAP="${REPO}/configs/orcd/twirl_fm0_s66_s77_full_visit_preparation_v1.json"
readonly OUT="$(realpath -m -- "${TWIRL_FM0_PHASE_A_AUTHORITY_RECORD}")"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"

die() { echo "[fm0-phase-a-freeze] ERROR: $*" >&2; exit 2; }
sha256() { sha256sum "$1" | awk '{print $1}'; }

for path in "${REPO}" "${OUT}"; do
  [[ "${path}" == "${PROJECT_ROOT}"/* ]] || die "path escapes project root: ${path}"
done
[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || die "invalid expected Git SHA"
[[ "${TWIRL_FM0_PREPARATION_AUTHORITY_MAP_SHA256}" =~ ^[0-9a-f]{64}$ ]] || die "invalid authority-map SHA-256"
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" ]] || die "wrong Git SHA"
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is not clean"
[[ -x "${PYTHON}" && -s "${MAP}" ]] || die "runtime or authority map is missing"
[[ "$(sha256 "${MAP}")" == "${TWIRL_FM0_PREPARATION_AUTHORITY_MAP_SHA256}" ]] || die "authority-map hash drift"
[[ ! -e "${OUT}" ]] || die "refusing to overwrite Phase-A authority record"

cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[fm0-phase-a-freeze] host=$(hostname) start=$(date -Is)"
"${PYTHON}" scripts/stage5_validation/manage_twirl_fm0_s66_s77_preparation_authorities.py \
  freeze-phase-a \
  --authority-map "${MAP}" \
  --authority-map-sha256 "${TWIRL_FM0_PREPARATION_AUTHORITY_MAP_SHA256}" \
  --producer-git-sha "${TWIRL_EXPECTED_GIT_SHA}" \
  --output "${OUT}"

[[ -s "${OUT}" ]] || die "Phase-A authority record was not published"
[[ ! -w "${OUT}" ]] || die "Phase-A authority record is writable"
"${PYTHON}" scripts/stage5_validation/manage_twirl_fm0_s66_s77_preparation_authorities.py \
  validate-phase-a \
  --authority-map "${MAP}" \
  --authority-map-sha256 "${TWIRL_FM0_PREPARATION_AUTHORITY_MAP_SHA256}" \
  --producer-git-sha "${TWIRL_EXPECTED_GIT_SHA}"
echo "[fm0-phase-a-freeze] complete=$(date -Is) sha256=$(sha256 "${OUT}")"
