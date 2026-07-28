#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-ssl-fullpool
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=192G
#SBATCH --gres=gpu:h200:1
#SBATCH --array=0-4%4
#SBATCH -t 2-00:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%A_%a.err

set -euo pipefail

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:-${TWIRL_ROOT}/code/TWIRL}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python}"
REGISTRY="${TWIRL_SSL_FULLPOOL_REGISTRY:?set TWIRL_SSL_FULLPOOL_REGISTRY}"
REGISTRY_SUMMARY="${TWIRL_SSL_FULLPOOL_REGISTRY_SUMMARY:?set TWIRL_SSL_FULLPOOL_REGISTRY_SUMMARY}"
OUT_ROOT="${TWIRL_SSL_FULLPOOL_OUT_ROOT:?set TWIRL_SSL_FULLPOOL_OUT_ROOT}"
EPOCHS="${TWIRL_SSL_FULLPOOL_EPOCHS:-20}"
BATCH_SIZE="${TWIRL_SSL_FULLPOOL_BATCH_SIZE:-64}"
WORKERS="${TWIRL_SSL_FULLPOOL_WORKERS:-8}"
SEED="${TWIRL_SSL_FULLPOOL_SEED:-560064}"
MAX_ROWS="${TWIRL_SSL_FULLPOOL_MAX_ROWS:-}"
SMOKE_OUT_ROOT="${TWIRL_SSL_FULLPOOL_SMOKE_OUT_ROOT:-}"
TRAIN_OUT_ROOT="${OUT_ROOT}"
smoke_args=()

if [[ -n "${MAX_ROWS}" ]]; then
  if [[ ! "${MAX_ROWS}" =~ ^[0-9]+$ ]] || (( 10#${MAX_ROWS} < 2 || 10#${MAX_ROWS} > 10000 )); then
    echo "TWIRL_SSL_FULLPOOL_MAX_ROWS must be an integer in [2,10000]." >&2
    exit 3
  fi
  if [[ "${EPOCHS}" != "1" ]]; then
    echo "A bounded full-pool smoke must set TWIRL_SSL_FULLPOOL_EPOCHS=1." >&2
    exit 3
  fi
  if [[ -z "${SMOKE_OUT_ROOT}" ]]; then
    echo "A bounded smoke requires TWIRL_SSL_FULLPOOL_SMOKE_OUT_ROOT." >&2
    exit 3
  fi
  OUT_ROOT="$(realpath -m -- "${OUT_ROOT}")"
  SMOKE_OUT_ROOT="$(realpath -m -- "${SMOKE_OUT_ROOT}")"
  if [[ "${SMOKE_OUT_ROOT}" == "${OUT_ROOT}" ]]; then
    echo "Smoke and production output roots must be distinct." >&2
    exit 3
  fi
  TRAIN_OUT_ROOT="${SMOKE_OUT_ROOT}"
  smoke_args+=(--max-rows "${MAX_ROWS}")
elif [[ -n "${SMOKE_OUT_ROOT}" ]]; then
  echo "TWIRL_SSL_FULLPOOL_SMOKE_OUT_ROOT requires a bounded MAX_ROWS smoke." >&2
  exit 3
fi

for input in "${REGISTRY}" "${REGISTRY_SUMMARY}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing full-pool SSL input: ${input}" >&2
    exit 3
  fi
done
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "This wrapper must run as the declared fold array." >&2
  exit 4
fi

mkdir -p "${TWIRL_ROOT}/logs" "${TRAIN_OUT_ROOT}"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

module load cuda/13.1.0 || module load cuda/12.9.1 || true
nvidia-smi

echo "[teacher-ssl-fullpool] start=$(date -Is) fold=${SLURM_ARRAY_TASK_ID} epochs=${EPOCHS} max_rows=${MAX_ROWS:-full} out_root=${TRAIN_OUT_ROOT}"
"${PYTHON}" scripts/stage5_validation/train_teacher_ssl_fullpool_fold.py \
  --registry "${REGISTRY}" \
  --registry-summary "${REGISTRY_SUMMARY}" \
  --out-root "${TRAIN_OUT_ROOT}" \
  --fold "${SLURM_ARRAY_TASK_ID}" \
  --epochs "${EPOCHS}" \
  --batch-size "${BATCH_SIZE}" \
  --workers "${WORKERS}" \
  --seed "${SEED}" \
  "${smoke_args[@]}" \
  --resume
echo "[teacher-ssl-fullpool] complete=$(date -Is) fold=${SLURM_ARRAY_TASK_ID}"
