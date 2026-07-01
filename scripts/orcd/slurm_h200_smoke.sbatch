#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-h200-smoke
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --gres=gpu:h200:1
#SBATCH -t 00:10:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs

echo "[h200-smoke] host=$(hostname)"
echo "[h200-smoke] date=$(date -Is)"

module load cuda/13.1.0 || module load cuda/12.9.1 || true

if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi
else
  echo "[h200-smoke] nvidia-smi not found" >&2
  exit 3
fi

PYTHON="${TWIRL_ORCD_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56/bin/python}"
"${PYTHON}" - <<'PY'
import importlib.util
import os
import sys

print("[h200-smoke] CUDA_VISIBLE_DEVICES=", os.environ.get("CUDA_VISIBLE_DEVICES", ""))
print("[h200-smoke] python_executable=", sys.executable)
print("[h200-smoke] python_version=", sys.version.replace("\n", " "))
for name in ("torch", "cupy"):
    spec = importlib.util.find_spec(name)
    print(f"[h200-smoke] {name}_available={spec is not None}")
if importlib.util.find_spec("torch") is not None:
    import torch
    print("[h200-smoke] torch_cuda_available=", torch.cuda.is_available())
    print("[h200-smoke] torch_device_count=", torch.cuda.device_count())
if importlib.util.find_spec("cupy") is not None:
    import cupy
    print("[h200-smoke] cupy_device_count=", cupy.cuda.runtime.getDeviceCount())
PY

echo "[h200-smoke] complete $(date -Is)"
