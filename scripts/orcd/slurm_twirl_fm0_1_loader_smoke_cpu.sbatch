#!/bin/bash
# Exercise the official FM0.1 loader/window boundary on two frozen shards.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-loader
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 01:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 027
readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
: "${TWIRL_ORCD_REPO:?missing detached repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?missing exact clean SHA}"
: "${TWIRL_FM0_RUN_ROOT:?missing run root}"
: "${TWIRL_FM0_INPUT_RECEIPT:?missing exact input receipt}"
: "${TWIRL_FM0_INPUT_RECEIPT_SHA256:?missing exact input receipt hash}"
readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly RUN_ROOT="$(realpath -m -- "${TWIRL_FM0_RUN_ROOT}")"
readonly INPUT_RECEIPT="$(realpath -m -- "${TWIRL_FM0_INPUT_RECEIPT}")"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
readonly OUT="${RUN_ROOT}/gates/loader_smoke.receipt.json"
die() { echo "[fm0-loader] ERROR: $*" >&2; exit 2; }
[[ "${RUN_ROOT}" == "${PROJECT_ROOT}"/* && "${INPUT_RECEIPT}" == "${RUN_ROOT}"/* ]] || die "unsafe run/input path"
[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ && "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" =~ ^[0-9a-f]{64}$ ]] || die "invalid hash"
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" && -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is not exact and clean"
[[ "$(sha256sum "${INPUT_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" ]] || die "input receipt hash mismatch"
[[ ! -e "${OUT}" ]] || die "refusing to overwrite loader receipt"
mkdir -p "$(dirname "${OUT}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
# Run the Torch-dependent model, loss, and interrupted-resume contracts in the
# exact pinned environment before publishing a loader receipt that can admit
# an H200 smoke.  These suites are skipped on the local workstation because it
# intentionally has no Torch installation.
"${PYTHON}" -m pytest -q \
  tests/test_twirl_fm0_model.py \
  tests/test_twirl_fm0_training.py
TWIRL_FM0_INPUT_RECEIPT="${INPUT_RECEIPT}" \
TWIRL_FM0_INPUT_RECEIPT_SHA256="${TWIRL_FM0_INPUT_RECEIPT_SHA256}" \
TWIRL_FM0_EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}" \
TWIRL_FM0_LOADER_OUT="${OUT}" \
"${PYTHON}" - <<'PY'
from datetime import datetime, timezone
import csv, hashlib, json, os, stat
from pathlib import Path
import numpy as np
from twirl.models.fm0.input_release import extract_window, load_input_release

receipt_path = Path(os.environ["TWIRL_FM0_INPUT_RECEIPT"])
if receipt_path.stat().st_mode & (stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH):
    raise SystemExit("input receipt is writable")
receipt = json.loads(receipt_path.read_text())
if receipt.get("schema_version") != "twirl_fm0_1_input_release_receipt_v1" or receipt.get("passed") is not True:
    raise SystemExit("input receipt did not pass the exact schema")
if receipt.get("expected_git_sha") != os.environ["TWIRL_FM0_EXPECTED_SHA"]:
    raise SystemExit("input receipt code binding drift")
manifest = Path(receipt["release"]["manifest_path"])
if hashlib.sha256(manifest.read_bytes()).hexdigest() != receipt["release"]["manifest_sha256"]:
    raise SystemExit("manifest drift after input freeze")
rows = sorted(csv.DictReader(manifest.open(newline="")), key=lambda row: row["observation_key"])
if len(rows) < 2:
    raise SystemExit("loader smoke requires at least two observations")
forbidden = {"physical_source_id", "leakage_component_id", "observation_key", "sector", "camera", "ccd", "absolute_bjd", "host_visit_offset_cadences", "host_visit_gap_cadences", "segment_id"}
allowed = {"flux", "flux_valid", "flux_error", "error_valid", "local_time_cadences", "delta_time_cadences", "time_valid", "segment_boundary", "view_present"}
audits = []
for row in rows[:2]:
    shard = manifest.parent / row["relative_path"]
    if hashlib.sha256(shard.read_bytes()).hexdigest() != row["sha256"]:
        raise SystemExit(f"shard drift: {shard}")
    release = load_input_release(shard)
    segment_ids = np.unique(release.segment_id)
    if segment_ids.size == 0:
        raise SystemExit(f"no segment in {shard}")
    window = extract_window(release, segment_id=int(segment_ids[0]), start_offset=0)
    if set(window) != allowed:
        raise SystemExit(f"forward allowlist drift: {sorted(window)}")
    overlap = forbidden.intersection(window)
    if overlap:
        raise SystemExit(f"forbidden forward keys: {sorted(overlap)}")
    audits.append({
        "observation_key": row["observation_key"],
        "shard_sha256": row["sha256"],
        "forward_keys": sorted(window),
        "forward_shapes": {key: list(np.asarray(value).shape) for key, value in window.items()},
    })
out = Path(os.environ["TWIRL_FM0_LOADER_OUT"])
payload = {
    "schema_version": "twirl_fm0_1_loader_smoke_receipt_v1",
    "created_at_utc": datetime.now(timezone.utc).isoformat(),
    "passed": True,
    "expected_git_sha": os.environ["TWIRL_FM0_EXPECTED_SHA"],
    "input_receipt": {"path": str(receipt_path), "sha256": os.environ["TWIRL_FM0_INPUT_RECEIPT_SHA256"]},
    "observations_checked": audits,
    "claim_limit": "CPU loader and model-visible window boundary only",
}
out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
chmod a-w "${OUT}"
sha256sum "${OUT}" > "${OUT}.sha256"
chmod a-w "${OUT}.sha256"
echo "[fm0-loader] complete=${OUT}"
