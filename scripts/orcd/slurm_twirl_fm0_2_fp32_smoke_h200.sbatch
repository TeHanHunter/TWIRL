#!/bin/bash
# One-H200 synthetic mechanics smoke for the frozen FM0.2.1 objective path.
# This consumes no real light-curve rows and cannot support a model claim.
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-fm0-2-1-fp32
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
readonly CONFIG_SHA="394a0340f9f8518423e39a1c83bc32acde4e2ebff95c4837b0edb5399bd50801"
readonly FREEZE_SHA="900264edabe13aa77dfffe3c87f5fe0c4967cbcf693b9a68020ab29ff85e691e"
readonly REUSE_SHA="46faed954e6758dfea71f7754670611fab03af63613dad3be433af71dd38e641"
readonly EVENT_SHA="eb84f27866e42bbf942256abe8fa7a532eaa8584e6800ce012aa0552c08a6630"
readonly EVALUATOR_SHA="ed8004180b28c7e2bd85911d0128f9180eddd3c35bdaae299e2f7edba6e5baac"
: "${TWIRL_ORCD_REPO:?missing detached repository}"
: "${TWIRL_EXPECTED_GIT_SHA:?missing exact clean SHA}"
: "${TWIRL_FM0_RUN_ROOT:?missing run root}"
: "${TWIRL_FM0_INPUT_RECEIPT:?missing exact upstream FM0.1 input receipt}"
: "${TWIRL_FM0_INPUT_RECEIPT_SHA256:?missing upstream input receipt hash}"
: "${TWIRL_FM0_INPUT_REUSE_RECEIPT:?missing frozen FM0.2 input-reuse receipt}"
: "${TWIRL_FM0_INPUT_REUSE_RECEIPT_SHA256:?missing input-reuse receipt hash}"
: "${TWIRL_FM0_LOADER_RESTART_RECEIPT:?missing exact-code loader/restart receipt}"
: "${TWIRL_FM0_LOADER_RESTART_RECEIPT_SHA256:?missing loader/restart receipt hash}"
: "${TWIRL_FM0_ADMISSION_RECEIPT:?missing fresh Stage-1-safe admission receipt}"
: "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256:?missing admission receipt hash}"
: "${TWIRL_FM0_SMOKE_OUTPUT:?missing fresh smoke output}"
readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly RUN_ROOT="$(realpath -m -- "${TWIRL_FM0_RUN_ROOT}")"
readonly INPUT_RECEIPT="$(realpath -m -- "${TWIRL_FM0_INPUT_RECEIPT}")"
readonly REUSE_RECEIPT="$(realpath -m -- "${TWIRL_FM0_INPUT_REUSE_RECEIPT}")"
readonly LOADER_RECEIPT="$(realpath -m -- "${TWIRL_FM0_LOADER_RESTART_RECEIPT}")"
readonly ADMISSION_RECEIPT="$(realpath -m -- "${TWIRL_FM0_ADMISSION_RECEIPT}")"
readonly OUT="$(realpath -m -- "${TWIRL_FM0_SMOKE_OUTPUT}")"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
readonly CONFIG="${REPO}/configs/models/twirl_fm0_2_s56_s64_objective_canary.yaml"
readonly FREEZE="${REPO}/reports/stage5_validation/twirl_fm0_2_design_freeze_v1/freeze.json"
readonly EVENT_CONTRACT="${REPO}/reports/stage5_validation/twirl_fm0_2_design_freeze_v1/development_event_retention.contract.json"
readonly EVALUATOR_CONTRACT="${REPO}/reports/stage5_validation/twirl_fm0_2_design_freeze_v1/development_evaluator_v2.contract.json"
die() { echo "[fm0.2-fp32] ERROR: $*" >&2; exit 2; }

[[ "${RUN_ROOT}" == "${PROJECT_ROOT}"/* && "${OUT}" == "${RUN_ROOT}"/* ]] || die "unsafe output path"
[[ "${INPUT_RECEIPT}" == "${PROJECT_ROOT}"/* && "${REUSE_RECEIPT}" == "${REPO}"/* && "${LOADER_RECEIPT}" == "${RUN_ROOT}"/* && "${ADMISSION_RECEIPT}" == "${RUN_ROOT}"/* ]] || die "unsafe authority path"
for value in "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" "${TWIRL_FM0_INPUT_REUSE_RECEIPT_SHA256}" "${TWIRL_FM0_LOADER_RESTART_RECEIPT_SHA256}" "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}"; do
  [[ "${value}" =~ ^[0-9a-f]{64}$ ]] || die "invalid receipt hash"
done
[[ "${TWIRL_FM0_INPUT_REUSE_RECEIPT_SHA256}" == "${REUSE_SHA}" ]] || die "unexpected reuse receipt hash"
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" && -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is not exact and clean"
[[ "$(sha256sum "${CONFIG}" | awk '{print $1}')" == "${CONFIG_SHA}" ]] || die "FM0.2 config drift"
[[ "$(sha256sum "${FREEZE}" | awk '{print $1}')" == "${FREEZE_SHA}" ]] || die "FM0.2 freeze drift"
[[ "$(sha256sum "${REUSE_RECEIPT}" | awk '{print $1}')" == "${REUSE_SHA}" ]] || die "reuse receipt drift"
[[ "$(sha256sum "${EVENT_CONTRACT}" | awk '{print $1}')" == "${EVENT_SHA}" ]] || die "event contract drift"
[[ "$(sha256sum "${EVALUATOR_CONTRACT}" | awk '{print $1}')" == "${EVALUATOR_SHA}" ]] || die "evaluator contract drift"
[[ "$(sha256sum "${INPUT_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" ]] || die "upstream input receipt drift"
[[ "$(sha256sum "${LOADER_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_LOADER_RESTART_RECEIPT_SHA256}" ]] || die "loader/restart receipt drift"
[[ "$(sha256sum "${ADMISSION_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}" ]] || die "admission receipt drift"
[[ ! -e "${OUT}" ]] || die "refusing to overwrite smoke output"

"${PYTHON}" - "${REUSE_RECEIPT}" "${INPUT_RECEIPT}" "${LOADER_RECEIPT}" "${ADMISSION_RECEIPT}" "${TWIRL_EXPECTED_GIT_SHA}" "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" "${TWIRL_FM0_INPUT_REUSE_RECEIPT_SHA256}" "${TWIRL_FM0_LOADER_RESTART_RECEIPT_SHA256}" <<'PY'
import json, stat, sys, time
from pathlib import Path

reuse_path, input_path, loader_path, admission_path = map(Path, sys.argv[1:5])
expected_sha, upstream_sha, reuse_sha, loader_sha = sys.argv[5:9]
for path in (input_path, loader_path, admission_path):
    if path.stat().st_mode & (stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH):
        raise SystemExit(f"authority is writable: {path}")
reuse = json.loads(reuse_path.read_text())
upstream = json.loads(input_path.read_text())
loader = json.loads(loader_path.read_text())
admission = json.loads(admission_path.read_text())
assert reuse["schema_version"] == "twirl_fm0_2_input_reuse_receipt_v1"
assert reuse["passed"] is True and reuse["scientific_training_eligible"] is True
assert reuse["upstream_fm0_1_input_receipt"]["sha256"] == upstream_sha
assert reuse["sealed_test"]["access_count"] == 0
assert upstream["schema_version"] == "twirl_fm0_1_input_release_receipt_v1"
assert upstream["passed"] is True and upstream["scientific_training_eligible"] is True
assert loader["schema_version"] == "twirl_fm0_2_loader_restart_receipt_v1"
assert loader["passed"] is True and loader["expected_git_sha"] == expected_sha
assert loader["input_receipt"]["sha256"] == upstream_sha
assert loader["input_reuse_receipt"]["sha256"] == reuse_sha
assert loader["sealed_test_access_count"] == 0 and loader["gpu_requested"] is False
assert admission["schema_version"] == "twirl_fm0_2_orcd_admission_receipt_v1"
assert admission["passed"] is True and admission["run_kind"] == "synthetic_fp32"
assert admission["expected_git_sha"] == expected_sha
assert admission["upstream_input_receipt"]["sha256"] == upstream_sha
assert admission["input_reuse_receipt"]["sha256"] == reuse_sha
assert admission["loader_restart_receipt"]["sha256"] == loader_sha
assert admission["requested_fm_resources"] == {"h200": 1, "cpus": 4, "memory_mib": 32768}
now = int(time.time())
if now < admission["created_at_epoch"] - 30 or now > admission["expires_at_epoch"]:
    raise SystemExit("admission receipt is not fresh at job start")
PY

mkdir -p "$(dirname "${OUT}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
module load cuda/13.1.0 || module load cuda/12.9.1 || true
nvidia-smi
"${PYTHON}" -c 'import torch; assert torch.cuda.is_available(); assert torch.cuda.device_count() == 1; print(torch.__version__, torch.cuda.get_device_name(0))'
"${PYTHON}" scripts/stage5_validation/train_twirl_fm0.py \
  --config "${CONFIG}" \
  --freeze-receipt "${FREEZE}" \
  --variant TWIRL-FM0.2.1 \
  --synthetic-smoke \
  --input-release-receipt "${INPUT_RECEIPT}" \
  --output-dir "${OUT}" \
  --device cuda \
  --precision fp32 \
  --stop-after-step 8 \
  --micro-batch-windows 8 \
  --seed 560067 \
  --expected-git-sha "${TWIRL_EXPECTED_GIT_SHA}"
for path in run_contract.json checkpoint.pt summary.json; do
  [[ -s "${OUT}/${path}" && -s "${OUT}/${path}.sha256" ]] || die "missing smoke artifact/sidecar: ${path}"
done
chmod -R a-w "${OUT}"
echo "[fm0.2-fp32] complete=${OUT}; synthetic only; no real data or model claim"
