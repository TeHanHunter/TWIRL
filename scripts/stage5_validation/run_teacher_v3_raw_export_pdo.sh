#!/usr/bin/env bash
# Export one sector's compact TGLC raw/error sources for Teacher v3.
set -euo pipefail

SECTOR="${1:?usage: $0 <sector> <orbit-one> <orbit-two> <compact-adp-h5> [raw-root] [out-h5]}"
ORBIT_ONE="${2:?usage: $0 <sector> <orbit-one> <orbit-two> <compact-adp-h5> [raw-root] [out-h5]}"
ORBIT_TWO="${3:?usage: $0 <sector> <orbit-one> <orbit-two> <compact-adp-h5> [raw-root] [out-h5]}"
COMPACT_ADP="${4:?usage: $0 <sector> <orbit-one> <orbit-two> <compact-adp-h5> [raw-root] [out-h5]}"
RAW_ROOT="${5:-/pdo/users/tehan/tglc-gpu-production-A2v1}"
OUT_H5="${6:-/pdo/users/tehan/twirl-teacher-v3/native_sources/s${SECTOR}_teacher_v3_tglc_raw_sources.h5}"

if (( SECTOR < 56 || SECTOR > 62 )); then
  echo "Teacher-v3 raw export is bounded to sectors 56-62." >&2
  exit 2
fi
if (( ORBIT_ONE <= 0 || ORBIT_TWO <= 0 || ORBIT_ONE == ORBIT_TWO )); then
  echo "Orbit IDs must be distinct positive integers." >&2
  exit 2
fi
if [[ "${OUT_H5}" != /pdo/users/tehan/* ]]; then
  echo "Refusing a PDO output outside /pdo/users/tehan/: ${OUT_H5}" >&2
  exit 2
fi

REPO="${TWIRL_ROOT:-/pdo/users/tehan/TWIRL}"
PYTHON_BIN="${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}"
RUN_REL="reports/stage5_validation/teacher_v3_s56_s62_a2v1_current_adp"
TRAINING_TABLE="${TWIRL_TEACHER_V3_SECTOR_TABLE:-${REPO}/${RUN_REL}/native_preparation/s${SECTOR}_native_export_table.csv}"
OVERWRITE="${TWIRL_OVERWRITE_TEACHER_V3_RAW_SOURCE:-0}"

for input in "${TRAINING_TABLE}" "${COMPACT_ADP}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing Teacher-v3 raw-export input: ${input}" >&2
    exit 3
  fi
done
for orbit in "${ORBIT_ONE}" "${ORBIT_TWO}"; do
  if [[ ! -d "${RAW_ROOT}/orbit-${orbit}/ffi" ]]; then
    echo "Missing A2v1 raw orbit root: ${RAW_ROOT}/orbit-${orbit}/ffi" >&2
    exit 3
  fi
done
if [[ -e "${OUT_H5}" && "${OVERWRITE}" != "1" ]]; then
  echo "Refusing to overwrite ${OUT_H5}; set TWIRL_OVERWRITE_TEACHER_V3_RAW_SOURCE=1." >&2
  exit 4
fi

mkdir -p "$(dirname "${OUT_H5}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export LD_LIBRARY_PATH="/sw/python-versions/python-3.11.9/lib:/pdo/app/anaconda/anaconda2-4.4.0/lib:${LD_LIBRARY_PATH:-}"
export HDF5_USE_FILE_LOCKING="${HDF5_USE_FILE_LOCKING:-FALSE}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

"${PYTHON_BIN}" scripts/stage5_validation/export_teacher_v3_raw_sources.py \
  --sector "${SECTOR}" \
  --training-table "${TRAINING_TABLE}" \
  --raw-root "${RAW_ROOT}" \
  --compact-adp-h5 "${COMPACT_ADP}" \
  --out-h5 "${OUT_H5}" \
  --orbits "${ORBIT_ONE}" "${ORBIT_TWO}"

sha256sum "${OUT_H5}" > "${OUT_H5}.sha256"
echo "[teacher-v3-raw-export] sector=${SECTOR} complete: ${OUT_H5}"
