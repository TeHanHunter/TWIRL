#!/usr/bin/env bash
# Build the reusable TWIRL S56 Python environment on ORCD project storage.
#
# Run this on ORCD from the staged checkout:
#   cd /orcd/data/mki_aryeh/001/twirl/code/TWIRL
#   bash scripts/orcd/build_orcd_env.sh
set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
ENV_PREFIX="${TWIRL_ORCD_ENV:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56}"
ENV_YML="${TWIRL_ORCD_ENV_YML:-${REPO}/configs/orcd/twirl-s56-env.yml}"
LOG_DIR="${TWIRL_ORCD_LOG_DIR:-/orcd/data/mki_aryeh/001/twirl/logs}"

mkdir -p "$(dirname "${ENV_PREFIX}")" "${LOG_DIR}"

if command -v module >/dev/null 2>&1; then
  module load miniforge/25.11.0-0 || module load miniforge/24.3.0-0 || module load miniforge/23.11.0-0 || true
fi

if command -v mamba >/dev/null 2>&1; then
  CONDA_FRONTEND=mamba
elif command -v conda >/dev/null 2>&1; then
  CONDA_FRONTEND=conda
else
  echo "[orcd-env] neither mamba nor conda is available after loading miniforge" >&2
  exit 3
fi

echo "[orcd-env] host=$(hostname)"
echo "[orcd-env] repo=${REPO}"
echo "[orcd-env] env_prefix=${ENV_PREFIX}"
echo "[orcd-env] env_yml=${ENV_YML}"
echo "[orcd-env] frontend=${CONDA_FRONTEND}"

if [[ -d "${ENV_PREFIX}/conda-meta" ]]; then
  "${CONDA_FRONTEND}" env update --yes --prefix "${ENV_PREFIX}" --file "${ENV_YML}" --prune
else
  "${CONDA_FRONTEND}" env create --yes --prefix "${ENV_PREFIX}" --file "${ENV_YML}"
fi

"${ENV_PREFIX}/bin/python" -m pip install --no-deps -e "${REPO}"

"${ENV_PREFIX}/bin/python" - <<'PY'
import importlib
import sys

required = (
    "numpy",
    "scipy",
    "pandas",
    "astropy",
    "h5py",
    "matplotlib",
    "seaborn",
    "pyarrow",
    "yaml",
    "batman",
    "tess_stars2px",
    "twirl",
)
missing = []
for name in required:
    try:
        importlib.import_module(name)
    except Exception as exc:
        missing.append(f"{name}: {type(exc).__name__}: {exc}")
if missing:
    print("[orcd-env] import smoke failed:", file=sys.stderr)
    for item in missing:
        print(f"  {item}", file=sys.stderr)
    sys.exit(4)
print("[orcd-env] import smoke passed")
PY

{
  echo "created_at=$(date -Is)"
  echo "host=$(hostname)"
  echo "repo=${REPO}"
  echo "env_prefix=${ENV_PREFIX}"
  echo "env_yml=${ENV_YML}"
  echo "python=$("${ENV_PREFIX}/bin/python" --version)"
  echo
  "${ENV_PREFIX}/bin/python" -m pip freeze
} > "${ENV_PREFIX}/twirl-s56-build-info.txt"

echo "[orcd-env] complete: ${ENV_PREFIX}"
echo "[orcd-env] python: ${ENV_PREFIX}/bin/python"
