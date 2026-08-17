#!/bin/bash
# One-H200, FP32, synthetic-only mechanics smoke for TWIRL-FM0.1.1.
# A validated real input receipt is still required as an external gate, but the
# current training CLI does not consume real shards and this is not a real-data
# POC result.
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-fm0-fp32-smoke
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 02:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail
umask 027
readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly DESIGN_SHA="94c8672a087884bf8a2c70f5d15315e05de602134af0c6c2073ca1c36232d6f7"
readonly CONFIG_SHA="7de4e8c9e98c0ce27648a21241acc51766e436a537ba932b54f529fbdbf26d8a"
: "${TWIRL_ORCD_REPO:?missing detached repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?missing exact clean SHA}"
: "${TWIRL_FM0_RUN_ROOT:?missing run root}"
: "${TWIRL_FM0_INPUT_RECEIPT:?missing exact input receipt}"
: "${TWIRL_FM0_INPUT_RECEIPT_SHA256:?missing input receipt hash}"
: "${TWIRL_FM0_LOADER_RECEIPT:?missing loader smoke receipt}"
: "${TWIRL_FM0_LOADER_RECEIPT_SHA256:?missing loader receipt hash}"
: "${TWIRL_FM0_ADMISSION_RECEIPT:?missing fresh Stage-1 admission receipt}"
: "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256:?missing admission receipt hash}"
: "${TWIRL_FM0_SMOKE_OUTPUT:?missing fresh smoke output path}"
readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly RUN_ROOT="$(realpath -m -- "${TWIRL_FM0_RUN_ROOT}")"
readonly INPUT_RECEIPT="$(realpath -m -- "${TWIRL_FM0_INPUT_RECEIPT}")"
readonly LOADER_RECEIPT="$(realpath -m -- "${TWIRL_FM0_LOADER_RECEIPT}")"
readonly ADMISSION_RECEIPT="$(realpath -m -- "${TWIRL_FM0_ADMISSION_RECEIPT}")"
readonly OUT="$(realpath -m -- "${TWIRL_FM0_SMOKE_OUTPUT}")"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
die() { echo "[fm0-fp32-smoke] ERROR: $*" >&2; exit 2; }
for path in "${RUN_ROOT}" "${INPUT_RECEIPT}" "${LOADER_RECEIPT}" "${ADMISSION_RECEIPT}" "${OUT}"; do
  [[ "${path}" == "${RUN_ROOT}" || "${path}" == "${RUN_ROOT}"/* ]] || die "path escapes run root: ${path}"
done
for value in "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" "${TWIRL_FM0_LOADER_RECEIPT_SHA256}" "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}"; do
  [[ "${value}" =~ ^[0-9a-f]{64}$ ]] || die "invalid receipt hash"
done
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" && -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is not exact and clean"
[[ "$(sha256sum "${REPO}/doc/foundation_model_design.md" | awk '{print $1}')" == "${DESIGN_SHA}" ]] || die "design drift"
[[ "$(sha256sum "${REPO}/configs/models/twirl_fm0_1_s56_s67_poc.yaml" | awk '{print $1}')" == "${CONFIG_SHA}" ]] || die "science config drift"
[[ "$(sha256sum "${INPUT_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" ]] || die "input receipt drift"
[[ "$(sha256sum "${LOADER_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_LOADER_RECEIPT_SHA256}" ]] || die "loader receipt drift"
[[ "$(sha256sum "${ADMISSION_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}" ]] || die "admission receipt drift"
[[ ! -e "${OUT}" ]] || die "refusing to overwrite smoke output"

TWIRL_FM0_ADMISSION_RECEIPT="${ADMISSION_RECEIPT}" \
TWIRL_FM0_EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}" \
TWIRL_FM0_INPUT_RECEIPT="${INPUT_RECEIPT}" \
TWIRL_FM0_INPUT_RECEIPT_SHA256="${TWIRL_FM0_INPUT_RECEIPT_SHA256}" \
TWIRL_FM0_LOADER_RECEIPT="${LOADER_RECEIPT}" \
"${PYTHON}" - <<'PY'
import json, os, stat, time
from pathlib import Path
p = Path(os.environ["TWIRL_FM0_ADMISSION_RECEIPT"])
if p.stat().st_mode & (stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH):
    raise SystemExit("admission receipt is writable")
d = json.loads(p.read_text())
input_path = Path(os.environ["TWIRL_FM0_INPUT_RECEIPT"])
loader_path = Path(os.environ["TWIRL_FM0_LOADER_RECEIPT"])
for authority in (input_path, loader_path):
    if authority.stat().st_mode & (stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH):
        raise SystemExit(f"authority is writable: {authority}")
input_receipt = json.loads(input_path.read_text())
loader_receipt = json.loads(loader_path.read_text())
assert input_receipt["schema_version"] == "twirl_fm0_1_input_release_receipt_v1" and input_receipt["passed"] is True
assert loader_receipt["schema_version"] == "twirl_fm0_1_loader_smoke_receipt_v1" and loader_receipt["passed"] is True
assert loader_receipt["input_receipt"]["sha256"] == os.environ["TWIRL_FM0_INPUT_RECEIPT_SHA256"]
assert d["schema_version"] == "twirl_fm0_1_orcd_admission_receipt_v1"
assert d["passed"] is True
assert d["expected_git_sha"] == os.environ["TWIRL_FM0_EXPECTED_SHA"]
assert d["input_receipt"]["path"] == os.environ["TWIRL_FM0_INPUT_RECEIPT"]
assert d["input_receipt"]["sha256"] == os.environ["TWIRL_FM0_INPUT_RECEIPT_SHA256"]
assert d["requested_fm_resources"] == {"h200": 1, "cpus": 4, "memory_mib": 32768}
now = int(time.time())
if now < d["created_at_epoch"] - 30 or now > d["expires_at_epoch"]:
    raise SystemExit("admission receipt is not fresh at job start")
PY

mkdir -p "$(dirname "${OUT}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
module load cuda/13.1.0 || module load cuda/12.9.1 || true
nvidia-smi
"${PYTHON}" -c 'import torch; assert torch.cuda.is_available(); assert torch.cuda.device_count()==1; assert torch.get_default_dtype()==torch.float32; print(torch.__version__, torch.cuda.get_device_name(0))'
"${PYTHON}" scripts/stage5_validation/train_twirl_fm0.py \
  --config "${REPO}/configs/models/twirl_fm0_1_s56_s67_poc.yaml" \
  --freeze-receipt "${REPO}/reports/stage5_validation/twirl_fm0_1_design_freeze_v1/freeze.json" \
  --variant TWIRL-FM0.1.1 \
  --output-dir "${OUT}" \
  --device cuda \
  --precision fp32 \
  --synthetic-smoke \
  --max-steps 8 \
  --seed 560067 \
  --expected-git-sha "${TWIRL_EXPECTED_GIT_SHA}"
for path in run_contract.json checkpoint.pt summary.json; do
  [[ -s "${OUT}/${path}" && -s "${OUT}/${path}.sha256" ]] || die "missing smoke artifact/sidecar: ${path}"
done
TWIRL_FM0_OUT="${OUT}" TWIRL_FM0_ADMISSION_SHA="${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}" \
TWIRL_FM0_INPUT_SHA="${TWIRL_FM0_INPUT_RECEIPT_SHA256}" \
TWIRL_FM0_LOADER_SHA="${TWIRL_FM0_LOADER_RECEIPT_SHA256}" \
"${PYTHON}" - <<'PY'
from datetime import datetime, timezone
import hashlib, json, os
from pathlib import Path
root = Path(os.environ["TWIRL_FM0_OUT"])
def digest(name): return hashlib.sha256((root / name).read_bytes()).hexdigest()
payload = {
  "schema_version": "twirl_fm0_1_orcd_synthetic_fp32_smoke_v1",
  "created_at_utc": datetime.now(timezone.utc).isoformat(), "passed": True,
  "synthetic_only": True, "real_input_consumed": False,
  "claim_limit": "GPU/model/loss/checkpoint mechanics only; not a real-data FM0.1 result",
  "admission_receipt_sha256": os.environ["TWIRL_FM0_ADMISSION_SHA"],
  "input_receipt_sha256": os.environ["TWIRL_FM0_INPUT_SHA"],
  "loader_receipt_sha256": os.environ["TWIRL_FM0_LOADER_SHA"],
  "artifacts": {name: digest(name) for name in ("run_contract.json", "checkpoint.pt", "summary.json")},
}
(root / "orcd_execution_receipt.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
chmod a-w "${OUT}/orcd_execution_receipt.json"
sha256sum "${OUT}/orcd_execution_receipt.json" > "${OUT}/orcd_execution_receipt.json.sha256"
chmod a-w "${OUT}/orcd_execution_receipt.json.sha256"
echo "[fm0-fp32-smoke] complete=${OUT} synthetic_only=true"
