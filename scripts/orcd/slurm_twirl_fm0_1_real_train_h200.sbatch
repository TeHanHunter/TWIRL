#!/bin/bash
# One-H200 real-release training for one controlled TWIRL-FM0.1 candidate.
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-fm0-1-1-real
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 47:30:00
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
: "${TWIRL_FM0_LOADER_RECEIPT:?missing loader receipt}"
: "${TWIRL_FM0_LOADER_RECEIPT_SHA256:?missing loader receipt hash}"
: "${TWIRL_FM0_ADMISSION_RECEIPT:?missing fresh Stage-1 admission receipt}"
: "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256:?missing admission receipt hash}"
: "${TWIRL_FM0_TRAIN_OUTPUT:?missing fresh training output path}"
readonly REPO="$(realpath -m -- "${TWIRL_ORCD_REPO}")"
readonly RUN_ROOT="$(realpath -m -- "${TWIRL_FM0_RUN_ROOT}")"
readonly INPUT_RECEIPT="$(realpath -m -- "${TWIRL_FM0_INPUT_RECEIPT}")"
readonly LOADER_RECEIPT="$(realpath -m -- "${TWIRL_FM0_LOADER_RECEIPT}")"
readonly ADMISSION_RECEIPT="$(realpath -m -- "${TWIRL_FM0_ADMISSION_RECEIPT}")"
readonly OUT="$(realpath -m -- "${TWIRL_FM0_TRAIN_OUTPUT}")"
readonly PYTHON="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1/bin/python"
readonly TARGET_STEP="${TWIRL_FM0_TARGET_STEP:-20000}"
readonly MICRO_BATCH="${TWIRL_FM0_MICRO_BATCH_WINDOWS:-8}"
readonly PRECISION="${TWIRL_FM0_PRECISION:-bf16}"
readonly VARIANT="${TWIRL_FM0_VARIANT:-TWIRL-FM0.1.1}"
readonly SEED="${TWIRL_FM0_SEED:-560067}"
die() { echo "[fm0-real-train] ERROR: $*" >&2; exit 2; }
for path in "${RUN_ROOT}" "${INPUT_RECEIPT}" "${LOADER_RECEIPT}" "${ADMISSION_RECEIPT}" "${OUT}"; do
  [[ "${path}" == "${RUN_ROOT}" || "${path}" == "${RUN_ROOT}"/* ]] || die "path escapes run root: ${path}"
done
[[ "${TARGET_STEP}" =~ ^[0-9]+$ && "${TARGET_STEP}" -ge 1 && "${TARGET_STEP}" -le 20000 ]] || die "invalid target step"
[[ "${MICRO_BATCH}" =~ ^[0-9]+$ ]] || die "micro batch must be numeric"
(( MICRO_BATCH > 0 && 64 % MICRO_BATCH == 0 )) || die "micro batch must divide 64"
[[ "${PRECISION}" == "fp32" || "${PRECISION}" == "bf16" ]] || die "invalid precision"
[[ "${VARIANT}" == "TWIRL-FM0.1.1" || "${VARIANT}" == "TWIRL-FM0.1.2" ]] || die "invalid FM0.1 variant"
[[ "${SEED}" == "560067" || "${SEED}" == "560068" ]] || die "invalid frozen FM0.1 seed"
for value in "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" "${TWIRL_FM0_LOADER_RECEIPT_SHA256}" "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}"; do
  [[ "${value}" =~ ^[0-9a-f]{64}$ ]] || die "invalid receipt hash"
done
[[ "$(git -C "${REPO}" rev-parse HEAD)" == "${TWIRL_EXPECTED_GIT_SHA}" && -z "$(git -C "${REPO}" status --porcelain=v1 --untracked-files=all)" ]] || die "repository is not exact and clean"
[[ "$(sha256sum "${REPO}/doc/foundation_model_design.md" | awk '{print $1}')" == "${DESIGN_SHA}" ]] || die "design drift"
[[ "$(sha256sum "${REPO}/configs/models/twirl_fm0_1_s56_s67_poc.yaml" | awk '{print $1}')" == "${CONFIG_SHA}" ]] || die "science config drift"
[[ "$(sha256sum "${INPUT_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" ]] || die "input receipt drift"
[[ "$(sha256sum "${LOADER_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_LOADER_RECEIPT_SHA256}" ]] || die "loader receipt drift"
[[ "$(sha256sum "${ADMISSION_RECEIPT}" | awk '{print $1}')" == "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}" ]] || die "admission receipt drift"
[[ ! -e "${OUT}" ]] || die "refusing to overwrite training output"

readonly RELEASE_ROOT="$(${PYTHON} - "${INPUT_RECEIPT}" "${LOADER_RECEIPT}" "${ADMISSION_RECEIPT}" "${TWIRL_EXPECTED_GIT_SHA}" "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" <<'PY'
import json, stat, sys, time
from pathlib import Path
input_path, loader_path, admission_path = map(Path, sys.argv[1:4])
expected_sha, input_sha = sys.argv[4:6]
for path in (input_path, loader_path, admission_path):
    if path.stat().st_mode & (stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH):
        raise SystemExit(f"authority is writable: {path}")
input_receipt = json.loads(input_path.read_text())
loader = json.loads(loader_path.read_text())
admission = json.loads(admission_path.read_text())
assert input_receipt["schema_version"] == "twirl_fm0_1_input_release_receipt_v1"
assert input_receipt["passed"] is True and input_receipt["scientific_training_eligible"] is True
assert loader["schema_version"] == "twirl_fm0_1_loader_smoke_receipt_v1" and loader["passed"] is True
assert loader["expected_git_sha"] == expected_sha
assert loader["input_receipt"]["sha256"] == input_sha
assert admission["schema_version"] == "twirl_fm0_1_orcd_admission_receipt_v1" and admission["passed"] is True
assert admission["expected_git_sha"] == expected_sha
assert admission["input_receipt"]["sha256"] == input_sha
assert admission["requested_fm_resources"] == {"h200": 1, "cpus": 4, "memory_mib": 32768}
now = int(time.time())
if now < admission["created_at_epoch"] - 30 or now > admission["expires_at_epoch"]:
    raise SystemExit("admission receipt is not fresh at job start")
manifest = Path(input_receipt["release"]["manifest_path"]).resolve(strict=True)
print(manifest.parent)
PY
)"

mkdir -p "$(dirname "${OUT}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
module load cuda/13.1.0 || module load cuda/12.9.1 || true
nvidia-smi
"${PYTHON}" -c 'import torch; assert torch.cuda.is_available(); assert torch.cuda.device_count()==1; print(torch.__version__, torch.cuda.get_device_name(0))'
"${PYTHON}" scripts/stage5_validation/train_twirl_fm0.py \
  --config "${REPO}/configs/models/twirl_fm0_1_s56_s67_poc.yaml" \
  --freeze-receipt "${REPO}/reports/stage5_validation/twirl_fm0_1_design_freeze_v1/freeze.json" \
  --variant "${VARIANT}" \
  --input-release "${RELEASE_ROOT}" \
  --input-release-receipt "${INPUT_RECEIPT}" \
  --output-dir "${OUT}" \
  --device cuda \
  --precision "${PRECISION}" \
  --max-steps "${TARGET_STEP}" \
  --micro-batch-windows "${MICRO_BATCH}" \
  --seed "${SEED}" \
  --expected-git-sha "${TWIRL_EXPECTED_GIT_SHA}"
for path in run_contract.json checkpoint.pt summary.json; do
  [[ -s "${OUT}/${path}" && -s "${OUT}/${path}.sha256" ]] || die "missing training artifact/sidecar: ${path}"
done
chmod -R a-w "${OUT}"
echo "[fm0-real-train] complete=${OUT} variant=${VARIANT} seed=${SEED} target_step=${TARGET_STEP} precision=${PRECISION}"
