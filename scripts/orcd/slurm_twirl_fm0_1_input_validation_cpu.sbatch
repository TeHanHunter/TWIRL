#!/bin/bash
# Build and independently validate one immutable FM0.1 input release.  The
# strict-NPZ fixture is mechanics-only; a checksum-bound A2v1 HDF5 manifest
# may produce a scientific-training-eligible S56 canary.  Neither path
# certifies survey completeness, production use, or a foundation-model claim.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-input
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 04:00:00
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
: "${TWIRL_FM0_SOURCE_MANIFEST:?controller must export the staged source manifest}"
: "${TWIRL_FM0_SOURCE_ADAPTER:?controller must export the exact source adapter}"

readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly RUN_ROOT="$(realpath -m -- "${TWIRL_FM0_RUN_ROOT}")"
readonly SOURCE_MANIFEST="$(realpath -m -- "${TWIRL_FM0_SOURCE_MANIFEST}")"
readonly REGISTRY="${RUN_ROOT}/input_build/registry"
readonly RELEASE="${RUN_ROOT}/frozen/input_release"
readonly STAGING="${RELEASE}.staging-${SLURM_JOB_ID:?missing Slurm job ID}"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"

die() { echo "[fm0-input] ERROR: $*" >&2; exit 2; }
sha256() { sha256sum "$1" | awk '{print $1}'; }
for path in "${RUN_ROOT}" "${SOURCE_MANIFEST}"; do
  [[ "${path}" == "${PROJECT_ROOT}"/* ]] || die "path escapes project root: ${path}"
done
[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || die "invalid expected Git SHA"
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" ]] || die "wrong Git SHA"
[[ -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is not clean"
[[ "$(sha256 "${REPO}/doc/foundation_model_design.md")" == "${DESIGN_SHA}" ]] || die "design drift"
[[ "$(sha256 "${REPO}/configs/models/twirl_fm0_1_s56_s67_poc.yaml")" == "${CONFIG_SHA}" ]] || die "science config drift"
[[ "$(sha256 "${REPO}/reports/stage5_validation/twirl_fm0_1_design_freeze_v1/freeze.json")" == "${FREEZE_SHA}" ]] || die "freeze receipt drift"
case "${TWIRL_FM0_SOURCE_ADAPTER}" in
  strict_npz_fixture_v1|a2v1_hdf5_quality_aware_v1) ;;
  *) die "unsupported source adapter: ${TWIRL_FM0_SOURCE_ADAPTER}" ;;
esac
[[ -x "${PYTHON}" && -s "${SOURCE_MANIFEST}" && -s "${REGISTRY}/summary.json" ]] || die "runtime/input/registry is missing"
[[ ! -e "${RELEASE}" && ! -e "${STAGING}" ]] || die "refusing to overwrite input release"

mkdir -p "$(dirname "${RELEASE}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
args=(
  --registry-dir "${REGISTRY}"
  --out-dir "${STAGING}"
  --design "${REPO}/doc/foundation_model_design.md"
  --config "${REPO}/configs/models/twirl_fm0_1_s56_s67_poc.yaml"
  --freeze-receipt "${REPO}/reports/stage5_validation/twirl_fm0_1_design_freeze_v1/freeze.json"
)
if [[ "${TWIRL_FM0_SOURCE_ADAPTER}" == "strict_npz_fixture_v1" ]]; then
  args+=(--raw-manifest "${SOURCE_MANIFEST}")
else
  args+=(--a2v1-hdf5-manifest "${SOURCE_MANIFEST}")
fi
"${PYTHON}" scripts/stage5_validation/build_twirl_fm0_input_release.py "${args[@]}"

TWIRL_FM0_RELEASE="${STAGING}" \
TWIRL_FM0_FINAL_RELEASE="${RELEASE}" \
TWIRL_FM0_EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}" \
TWIRL_FM0_SOURCE_MANIFEST="${SOURCE_MANIFEST}" \
TWIRL_FM0_SOURCE_ADAPTER="${TWIRL_FM0_SOURCE_ADAPTER}" \
TWIRL_FM0_REGISTRY_SUMMARY="${REGISTRY}/summary.json" \
TWIRL_FM0_DESIGN_SHA="${DESIGN_SHA}" \
TWIRL_FM0_CONFIG_SHA="${CONFIG_SHA}" \
TWIRL_FM0_FREEZE_SHA="${FREEZE_SHA}" \
"${PYTHON}" - <<'PY'
from __future__ import annotations
from datetime import datetime, timezone
import csv, hashlib, json, os, stat
from pathlib import Path, PurePosixPath
import numpy as np

release = Path(os.environ["TWIRL_FM0_RELEASE"])
final_release = Path(os.environ["TWIRL_FM0_FINAL_RELEASE"])
manifest_path = release / "manifest.csv"
summary_path = release / "summary.json"
required_columns = [
    "input_release_schema_version", "observation_key", "product_instance_id",
    "source_sha256", "leakage_component_id", "source_partition",
    "product_state", "relative_path", "sha256", "input_source_sha256",
    "n_cadences", "n_segments",
    "view_present_json", "host_visit_offset_cadences",
    "host_visit_gap_cadences", "host_visit_overlaps_previous", "input_adapter",
    "scientific_training_eligible",
]
required_arrays = {
    "flux", "flux_valid", "flux_error", "error_valid",
    "local_time_cadences", "delta_time_cadences", "time_valid",
    "segment_boundary", "segment_id", "view_present",
}
for path in (manifest_path, summary_path):
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"missing input-release artifact: {path}")
summary = json.loads(summary_path.read_text())
if summary.get("certifies_full_campaign") is not False:
    raise SystemExit("builder summary must explicitly deny full-campaign certification")
adapter = os.environ["TWIRL_FM0_SOURCE_ADAPTER"]
if summary.get("input_adapter") != adapter:
    raise SystemExit("input-release adapter does not match the submitted source adapter")
scientific_eligible = adapter == "a2v1_hdf5_quality_aware_v1"
if summary.get("scientific_training_eligible") is not scientific_eligible:
    raise SystemExit("input-release scientific-training eligibility does not match adapter")

rows = list(csv.DictReader(manifest_path.open(newline="")))
if not rows or list(rows[0]) != required_columns:
    raise SystemExit("input manifest is empty or has drifted columns/order")
seen = set()
identity = hashlib.sha256()
for row in rows:
    if row["input_release_schema_version"] != "twirl_fm0_1_input_release_v1":
        raise SystemExit("input-release schema drift")
    if row["input_adapter"] != adapter:
        raise SystemExit("input-release row adapter does not match submitted source adapter")
    if row["scientific_training_eligible"] != str(scientific_eligible):
        raise SystemExit("input-release row scientific-training eligibility mismatch")
    key = row["observation_key"]
    if not key or key in seen:
        raise SystemExit(f"invalid/duplicate observation_key: {key!r}")
    seen.add(key)
    rel = PurePosixPath(row["relative_path"])
    if rel.is_absolute() or ".." in rel.parts:
        raise SystemExit(f"unsafe shard path: {rel}")
    shard = release / rel
    digest = hashlib.sha256(shard.read_bytes()).hexdigest()
    if digest != row["sha256"]:
        raise SystemExit(f"shard hash mismatch: {rel}")
    with np.load(shard, allow_pickle=False) as data:
        if set(data.files) != required_arrays:
            raise SystemExit(f"shard key drift: {rel}: {sorted(data.files)}")
        n = int(row["n_cadences"])
        expected = {
            "flux": (n, 6), "flux_valid": (n, 6),
            "flux_error": (n, 2), "error_valid": (n, 2),
            "local_time_cadences": (n,), "delta_time_cadences": (n,),
            "time_valid": (n,), "segment_boundary": (n,),
            "segment_id": (n,), "view_present": (6,),
        }
        for name, shape in expected.items():
            if data[name].shape != shape:
                raise SystemExit(f"shard shape drift: {rel}:{name}={data[name].shape}")
        for mask in ("flux_valid", "error_valid", "time_valid", "segment_boundary", "view_present"):
            if data[mask].dtype != np.bool_:
                raise SystemExit(f"mask dtype drift: {rel}:{mask}")
        if not np.isfinite(data["flux"][data["flux_valid"]]).all():
            raise SystemExit(f"non-finite valid flux: {rel}")
        if not np.isfinite(data["flux_error"][data["error_valid"]]).all():
            raise SystemExit(f"non-finite valid error: {rel}")
    identity.update(key.encode()); identity.update(b"\0"); identity.update(digest.encode()); identity.update(b"\n")

def hash_file(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
receipt = {
    "schema_version": "twirl_fm0_1_input_release_receipt_v1",
    "created_at_utc": datetime.now(timezone.utc).isoformat(),
    "passed": True,
    "campaign_id": "twirl_fm0_1_s56_s67_poc",
    "scope": "validated_named_release_artifacts_only",
    "input_adapter": adapter,
    "scientific_training_eligible": scientific_eligible,
    "certifies_full_campaign": False,
    "production_authorized": False,
    "foundation_model_claim_authorized": False,
    "expected_git_sha": os.environ["TWIRL_FM0_EXPECTED_SHA"],
    "science_bindings": {
        "design_sha256": os.environ["TWIRL_FM0_DESIGN_SHA"],
        "config_sha256": os.environ["TWIRL_FM0_CONFIG_SHA"],
        "freeze_receipt_sha256": os.environ["TWIRL_FM0_FREEZE_SHA"],
    },
    "source_bindings": {
        "source_manifest": {"path": os.environ["TWIRL_FM0_SOURCE_MANIFEST"], "sha256": hash_file(os.environ["TWIRL_FM0_SOURCE_MANIFEST"])},
        "registry_summary": {"path": os.environ["TWIRL_FM0_REGISTRY_SUMMARY"], "sha256": hash_file(os.environ["TWIRL_FM0_REGISTRY_SUMMARY"])},
    },
    "release": {
        "manifest_path": str(final_release / "manifest.csv"), "manifest_sha256": hash_file(manifest_path),
        "summary_path": str(final_release / "summary.json"), "summary_sha256": hash_file(summary_path),
        "n_shards": len(rows), "ordered_identity_sha256": identity.hexdigest(),
    },
}
target = release / "input_release.receipt.json"
target.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
PY

mv -- "${STAGING}" "${RELEASE}"
sha256sum "${RELEASE}/input_release.receipt.json" > "${RELEASE}/input_release.receipt.json.sha256"
chmod -R a-w "${RELEASE}"
echo "[fm0-input] complete=${RELEASE}/input_release.receipt.json"
