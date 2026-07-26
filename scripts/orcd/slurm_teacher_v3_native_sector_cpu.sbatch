#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-tv3-native
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

SECTOR="${TWIRL_TEACHER_V3_SECTOR:?set TWIRL_TEACHER_V3_SECTOR to one of 56-62}"
if (( SECTOR < 56 || SECTOR > 62 )); then
  echo "Teacher-v3 native preparation is bounded to sectors 56-62." >&2
  exit 2
fi

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:-${TWIRL_ROOT}/code/TWIRL}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python}"
RUN_REL="reports/stage5_validation/teacher_v3_s56_s62_a2v1_current_adp"
TRAINING="${TWIRL_TEACHER_V3_SECTOR_TABLE:-${TWIRL_ROOT}/${RUN_REL}/native_preparation/s${SECTOR}_native_export_table.csv}"
RAW_SOURCE="${TWIRL_TEACHER_V3_RAW_SOURCE_H5:-${TWIRL_ROOT}/exports/teacher_v3/native_sources/s${SECTOR}_teacher_v3_tglc_raw_sources.h5}"
COMPACT_ADP="${TWIRL_TEACHER_V3_COMPACT_ADP_H5:?set the exact S${SECTOR} A2v1 compact ADP HDF5}"
CADENCE_REFERENCE_TABLE="${TWIRL_TEACHER_V3_CADENCE_REFERENCE_TABLE:?set the exact S${SECTOR} cadence-reference table}"
CADENCE_REFERENCE_MANIFEST="${TWIRL_TEACHER_V3_CADENCE_REFERENCE_MANIFEST:?set the exact S${SECTOR} cadence-reference manifest}"
INJECTION_PAIR="${TWIRL_TEACHER_V3_INJECTION_PAIR_H5:-}"
FINAL="${TWIRL_TEACHER_V3_NATIVE_H5:-${TWIRL_ROOT}/${RUN_REL}/native/s${SECTOR}_teacher_v3_native_v2.h5}"
N_SHARDS="${TWIRL_TEACHER_V3_NATIVE_SHARDS:-1}"
RUN_DIR="${TWIRL_ROOT}/runtime/teacher_v3/native_s${SECTOR}_${SLURM_JOB_ID}"

if (( N_SHARDS < 1 || N_SHARDS > SLURM_CPUS_PER_TASK )); then
  echo "TWIRL_TEACHER_V3_NATIVE_SHARDS must be in 1..${SLURM_CPUS_PER_TASK}." >&2
  exit 2
fi
if [[ "${SECTOR}" == "56" && ! -s "${INJECTION_PAIR}" ]]; then
  echo "S56 contains Teacher-v3 injections; set TWIRL_TEACHER_V3_INJECTION_PAIR_H5." >&2
  exit 3
fi
if [[ "${SECTOR}" != "56" && -n "${INJECTION_PAIR}" ]]; then
  echo "S${SECTOR} is real-only; unset TWIRL_TEACHER_V3_INJECTION_PAIR_H5." >&2
  exit 3
fi
for input in \
  "${TRAINING}" \
  "${RAW_SOURCE}" \
  "${COMPACT_ADP}" \
  "${CADENCE_REFERENCE_TABLE}" \
  "${CADENCE_REFERENCE_MANIFEST}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing Teacher-v3 native input: ${input}" >&2
    exit 3
  fi
done
if [[ -e "${FINAL}" || -e "${FINAL%.h5}.summary.json" ]]; then
  echo "Refusing to overwrite an existing Teacher-v3 native artifact." >&2
  echo "output=${FINAL} summary=${FINAL%.h5}.summary.json" >&2
  exit 4
fi

mkdir -p "${TWIRL_ROOT}/logs" "${RUN_DIR}" "$(dirname "${FINAL}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

common_args=(
  --sector "${SECTOR}"
  --training-table "${TRAINING}"
  --raw-source-h5 "${RAW_SOURCE}"
  --compact-adp-h5 "${COMPACT_ADP}"
  --cadence-reference-table "${CADENCE_REFERENCE_TABLE}"
  --cadence-reference-manifest "${CADENCE_REFERENCE_MANIFEST}"
  --repo-root "${REPO}"
  --n-periods 4096
  --n-shards "${N_SHARDS}"
)
if [[ -n "${INJECTION_PAIR}" ]]; then
  common_args+=(--injection-pair-h5 "${INJECTION_PAIR}")
fi

echo "[teacher-v3-native] sector=${SECTOR} shards=${N_SHARDS} start=$(date -Is)"
pids=()
for ((shard=0; shard<N_SHARDS; shard++)); do
  shard_path="${RUN_DIR}/native_${shard}.h5"
  srun --exclusive --exact -N1 -n1 -c1 --mem=12G \
    "${PYTHON}" scripts/stage5_validation/build_teacher_v3_native_inputs.py \
    "${common_args[@]}" \
    --out-h5 "${shard_path}" \
    --shard-index "${shard}" &
  pids+=("$!")
done

shard_status=0
for pid in "${pids[@]}"; do
  if ! wait "${pid}"; then
    shard_status=1
  fi
done
if [[ "${shard_status}" != "0" ]]; then
  echo "At least one S${SECTOR} native shard failed; refusing to merge." >&2
  exit 5
fi

shards=("${RUN_DIR}"/native_*.h5)
"${PYTHON}" scripts/stage5_validation/merge_teacher_v3_native_shards.py \
  --sector "${SECTOR}" \
  --training-table "${TRAINING}" \
  --shards "${shards[@]}" \
  --out-h5 "${FINAL}"

sha256sum "${FINAL}" > "${FINAL}.sha256"
echo "[teacher-v3-native] sector=${SECTOR} complete=$(date -Is) output=${FINAL}"
