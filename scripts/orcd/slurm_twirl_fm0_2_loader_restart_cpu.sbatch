#!/bin/bash
# CPU-only exact-code loader/restart gate for the frozen FM0.2.1 canary.
# No GPU is requested and node4900 is explicitly excluded.
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-fm0-2-loader
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
readonly CONFIG_SHA="394a0340f9f8518423e39a1c83bc32acde4e2ebff95c4837b0edb5399bd50801"
readonly FREEZE_SHA="900264edabe13aa77dfffe3c87f5fe0c4967cbcf693b9a68020ab29ff85e691e"
readonly REUSE_SHA="46faed954e6758dfea71f7754670611fab03af63613dad3be433af71dd38e641"
: "${TWIRL_ORCD_REPO:?missing detached repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?missing exact clean SHA}"
: "${TWIRL_FM0_RUN_ROOT:?missing run root}"
: "${TWIRL_FM0_INPUT_RECEIPT:?missing exact upstream input receipt}"
: "${TWIRL_FM0_INPUT_RECEIPT_SHA256:?missing exact upstream input receipt hash}"
: "${TWIRL_FM0_INPUT_REUSE_RECEIPT:?missing frozen FM0.2 input-reuse receipt}"
: "${TWIRL_FM0_INPUT_REUSE_RECEIPT_SHA256:?missing input-reuse receipt hash}"
readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly RUN_ROOT="$(realpath -m -- "${TWIRL_FM0_RUN_ROOT}")"
readonly INPUT_RECEIPT="$(realpath -e -- "${TWIRL_FM0_INPUT_RECEIPT}")"
readonly REUSE_RECEIPT="$(realpath -e -- "${TWIRL_FM0_INPUT_REUSE_RECEIPT}")"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
readonly CONFIG="${REPO}/configs/models/twirl_fm0_2_s56_s64_objective_canary.yaml"
readonly FREEZE="${REPO}/reports/stage5_validation/twirl_fm0_2_design_freeze_v1/freeze.json"
readonly OUT="${RUN_ROOT}/gates/loader_restart.receipt.json"
die() { echo "[fm0.2-loader] ERROR: $*" >&2; exit 2; }

[[ "${RUN_ROOT}" == "${PROJECT_ROOT}"/* ]] || die "unsafe run root"
[[ "${INPUT_RECEIPT}" == "${PROJECT_ROOT}"/* && "${REUSE_RECEIPT}" == "${REPO}"/* ]] || die "unsafe authority path"
[[ "${TWIRL_EXPECTED_GIT_SHA}" =~ ^[0-9a-f]{40}$ ]] || die "invalid Git SHA"
for value in "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" "${TWIRL_FM0_INPUT_REUSE_RECEIPT_SHA256}"; do
  [[ "${value}" =~ ^[0-9a-f]{64}$ ]] || die "invalid receipt hash"
done
[[ "${TWIRL_FM0_INPUT_REUSE_RECEIPT_SHA256}" == "${REUSE_SHA}" ]] || die "unexpected reuse receipt hash"
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" && -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is not exact and clean"
[[ "$(sha256sum "${CONFIG}" | awk '{print $1}')" == "${CONFIG_SHA}" ]] || die "FM0.2 config drift"
[[ "$(sha256sum "${FREEZE}" | awk '{print $1}')" == "${FREEZE_SHA}" ]] || die "FM0.2 freeze drift"
[[ "$(sha256sum "${INPUT_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" ]] || die "upstream input receipt drift"
[[ "$(sha256sum "${REUSE_RECEIPT}" | awk '{print $1}')" == "${REUSE_SHA}" ]] || die "input-reuse receipt drift"
[[ ! -e "${OUT}" ]] || die "refusing to overwrite loader/restart receipt"
mkdir -p "$(dirname "${OUT}")"

cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"${PYTHON}" -m pytest -q \
  tests/test_twirl_fm0_model.py \
  tests/test_twirl_fm0_training.py \
  tests/test_twirl_fm0_real_dataset.py \
  tests/test_twirl_fm0_representation_health.py \
  tests/test_twirl_fm0_representation_health_v2.py \
  tests/test_twirl_fm0_2_contract.py \
  tests/test_twirl_fm0_observation_sector_authority.py

TWIRL_FM0_INPUT_RECEIPT="${INPUT_RECEIPT}" \
TWIRL_FM0_INPUT_RECEIPT_SHA256="${TWIRL_FM0_INPUT_RECEIPT_SHA256}" \
TWIRL_FM0_REUSE_RECEIPT="${REUSE_RECEIPT}" \
TWIRL_FM0_REUSE_RECEIPT_SHA256="${REUSE_SHA}" \
TWIRL_FM0_EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}" \
TWIRL_FM0_LOADER_OUT="${OUT}" \
"${PYTHON}" - <<'PY'
from datetime import datetime, timezone
import csv, hashlib, json, os, stat
from pathlib import Path
import numpy as np
from twirl.models.fm0.input_release import extract_window, load_input_release

input_path = Path(os.environ["TWIRL_FM0_INPUT_RECEIPT"])
reuse_path = Path(os.environ["TWIRL_FM0_REUSE_RECEIPT"])
if input_path.stat().st_mode & (stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH):
    raise SystemExit(f"external upstream authority is writable: {input_path}")
# The reuse receipt is a tracked file inside the already verified clean,
# exact-SHA worktree.  Its frozen SHA-256 is checked by the shell wrapper; its
# ordinary Git mode is therefore not an authority weakness.
upstream = json.loads(input_path.read_text())
reuse = json.loads(reuse_path.read_text())
if upstream.get("schema_version") != "twirl_fm0_1_input_release_receipt_v1" or upstream.get("passed") is not True:
    raise SystemExit("upstream input receipt did not pass")
if reuse.get("schema_version") != "twirl_fm0_2_input_reuse_receipt_v1" or reuse.get("passed") is not True:
    raise SystemExit("FM0.2 input-reuse receipt did not pass")
if reuse.get("scientific_training_eligible") is not True or reuse.get("sealed_test", {}).get("access_count") != 0:
    raise SystemExit("FM0.2 reuse boundary drifted")
if reuse.get("upstream_fm0_1_input_receipt", {}).get("sha256") != os.environ["TWIRL_FM0_INPUT_RECEIPT_SHA256"]:
    raise SystemExit("FM0.2 reuse receipt does not bind the upstream input")
manifest = Path(upstream["release"]["manifest_path"])
manifest_sha = hashlib.sha256(manifest.read_bytes()).hexdigest()
if manifest_sha != upstream["release"]["manifest_sha256"]:
    raise SystemExit("manifest drift after input freeze")
rows = sorted(csv.DictReader(manifest.open(newline="")), key=lambda row: row["observation_key"])
if len(rows) < 2:
    raise SystemExit("loader gate requires at least two observations")
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
    if set(window) != allowed or forbidden.intersection(window):
        raise SystemExit("model-visible loader allowlist drift")
    audits.append({
        "observation_key": row["observation_key"],
        "shard_sha256": row["sha256"],
        "forward_keys": sorted(window),
        "forward_shapes": {key: list(np.asarray(value).shape) for key, value in window.items()},
    })
payload = {
    "schema_version": "twirl_fm0_2_loader_restart_receipt_v1",
    "created_at_utc": datetime.now(timezone.utc).isoformat(),
    "passed": True,
    "expected_git_sha": os.environ["TWIRL_FM0_EXPECTED_SHA"],
    "input_receipt": {"path": str(input_path), "sha256": os.environ["TWIRL_FM0_INPUT_RECEIPT_SHA256"]},
    "input_reuse_receipt": {"path": str(reuse_path), "sha256": os.environ["TWIRL_FM0_REUSE_RECEIPT_SHA256"]},
    "tests": {
        "torch_model_loss_exact_resume_and_milestones": True,
        "real_loader": True,
        "evaluator_v2": True,
    },
    "observations_checked": audits,
    "sealed_test_access_count": 0,
    "gpu_requested": False,
    "resources": {"cpus": 4, "memory_mib": 16384},
    "claim_limit": "CPU loader, exact restart, milestone, and evaluator mechanics only",
}
Path(os.environ["TWIRL_FM0_LOADER_OUT"]).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
sha256sum "${OUT}" > "${OUT}.sha256"
chmod a-w "${OUT}" "${OUT}.sha256"
echo "[fm0.2-loader] complete=${OUT}; no GPU or sealed-test access"
