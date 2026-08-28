#!/bin/bash
# Fully revalidate and admit the deferred S66-S77 full-visit preparation pool.
# Passing this CPU gate does not accept A2v1 sectors or authorize training.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-admission-v2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH -t 2-00:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 027

readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
: "${TWIRL_ORCD_REPO:?controller must export the detached repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?controller must export the exact clean SHA}"
: "${TWIRL_FM0_PREPARATION_AUTHORITY_MAP_SHA256:?controller must export the map SHA-256}"
: "${TWIRL_FM0_ADMISSION_V2_RECEIPT:?controller must export the admission receipt}"

readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly MAP="${REPO}/configs/orcd/twirl_fm0_s66_s77_full_visit_preparation_v1.json"
readonly OUT="$(realpath -m -- "${TWIRL_FM0_ADMISSION_V2_RECEIPT}")"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"

die() { echo "[fm0-admission-v2] ERROR: $*" >&2; exit 2; }
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
[[ ! -e "${OUT}" ]] || die "refusing to overwrite admission-v2 receipt"

cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[fm0-admission-v2] host=$(hostname) start=$(date -Is)"
"${PYTHON}" scripts/stage5_validation/manage_twirl_fm0_s66_s77_preparation_authorities.py \
  admit \
  --authority-map "${MAP}" \
  --authority-map-sha256 "${TWIRL_FM0_PREPARATION_AUTHORITY_MAP_SHA256}" \
  --producer-git-sha "${TWIRL_EXPECTED_GIT_SHA}" \
  --output "${OUT}"

[[ -s "${OUT}" ]] || die "admission-v2 receipt was not published"
"${PYTHON}" - "${OUT}" "${TWIRL_EXPECTED_GIT_SHA}" <<'PY'
from __future__ import annotations

import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
revision = sys.argv[2]
receipt = json.loads(path.read_text(encoding="utf-8"))
if receipt.get("schema_version") != "twirl_fm0_later_sector_preparation_pool_receipt_v2":
    raise SystemExit("admission-v2 schema drifted")
if receipt.get("producer_git_sha") != revision:
    raise SystemExit("admission-v2 producer revision drifted")
if receipt.get("preparation_pool_admitted") is not True:
    raise SystemExit("preparation pool did not pass")
for key in (
    "a2v1_accepted",
    "scientific_training_eligible",
    "temporal_panel_frozen",
    "model_training_authorized",
    "sealed_test_access_authorized",
    "sealed_aperture_photometry_opened",
    "sealed_shards_written",
):
    if receipt.get(key) is not False:
        raise SystemExit(f"admission-v2 overclaims {key}")
if receipt.get("preparation_pool_sectors") != list(range(66, 78)):
    raise SystemExit("admission-v2 sector selection drifted")
if [int(row["sector"]) for row in receipt.get("sector_receipts", ())] != list(range(66, 78)):
    raise SystemExit("admission-v2 sector order drifted")
PY
chmod a-w "${OUT}"
echo "[fm0-admission-v2] complete=$(date -Is) sha256=$(sha256 "${OUT}")"
