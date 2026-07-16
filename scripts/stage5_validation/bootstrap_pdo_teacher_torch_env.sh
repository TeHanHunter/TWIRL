#!/usr/bin/env bash
# Build a user-owned PyTorch environment for teacher inference on pdogpu6.
set -euo pipefail

BASE_PYTHON="${TWIRL_PDO_BASE_PYTHON:-/sw/qlp-environment/.venv/bin/python}"
ENV_ROOT="${TWIRL_PDO_TORCH_ENV:-/pdo/users/tehan/envs/twirl-teacher-pdo-v2}"
# PDO is glibc 2.17. PyTorch 2.7+ CUDA wheels require manylinux_2_28,
# while the official 2.6 CUDA 11.8 wheel remains compatible with this host.
TORCH_INDEX="${TWIRL_PDO_TORCH_INDEX_URL:-https://download.pytorch.org/whl/cu118}"
TORCH_SPEC="${TWIRL_PDO_TORCH_SPEC:-torch==2.6.0}"
export LD_LIBRARY_PATH="/sw/python-versions/python-3.11.9/lib:${LD_LIBRARY_PATH:-}"
if [[ -r /etc/pki/tls/certs/ca-bundle.crt ]]; then
  export REQUESTS_CA_BUNDLE=/etc/pki/tls/certs/ca-bundle.crt
  export SSL_CERT_FILE=/etc/pki/tls/certs/ca-bundle.crt
  export PIP_CERT=/etc/pki/tls/certs/ca-bundle.crt
fi
mkdir -p "$(dirname "${ENV_ROOT}")"
exec 9>"/tmp/twirl-teacher-pdo-v2.lock"
if ! flock -n 9; then
  echo "Another PDO teacher environment build is already running" >&2
  exit 9
fi

if [[ ! -x "${ENV_ROOT}/bin/python" ]]; then
  "${BASE_PYTHON}" -m venv --system-site-packages "${ENV_ROOT}"
fi
if ! "${ENV_ROOT}/bin/python" -m pip --version >/dev/null 2>&1; then
  "${ENV_ROOT}/bin/python" -m ensurepip --upgrade
fi
"${ENV_ROOT}/bin/python" -m pip install --upgrade pip
"${ENV_ROOT}/bin/python" -m pip install --index-url "${TORCH_INDEX}" "${TORCH_SPEC}"
"${ENV_ROOT}/bin/python" - <<'PY'
import torch
print("torch", torch.__version__)
print("cuda_available", torch.cuda.is_available())
print("device_count", torch.cuda.device_count())
if not torch.cuda.is_available():
    raise SystemExit("CUDA is not available in the PDO teacher environment")
PY
