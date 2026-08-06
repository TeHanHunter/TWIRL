#!/usr/bin/env bash
# Export one sector's compact TGLC raw/error sources for Teacher v3.
set -euo pipefail

SECTOR="${1:?usage: $0 <sector> <orbit-one> <orbit-two> <compact-adp-h5> [raw-root] [out-h5]}"
ORBIT_ONE="${2:?usage: $0 <sector> <orbit-one> <orbit-two> <compact-adp-h5> [raw-root] [out-h5]}"
ORBIT_TWO="${3:?usage: $0 <sector> <orbit-one> <orbit-two> <compact-adp-h5> [raw-root] [out-h5]}"
COMPACT_ADP="${4:?usage: $0 <sector> <orbit-one> <orbit-two> <compact-adp-h5> [raw-root] [out-h5]}"
RAW_ROOT="${5:-/pdo/users/tehan/tglc-gpu-production-A2v1}"
OUT_H5="${6:-/pdo/users/tehan/twirl-teacher-v3/native_sources/s${SECTOR}_teacher_v3_tglc_raw_sources.h5}"

S63_PROSPECTIVE_CONTRACT="s63_teacher_v3_prospective_v1"
PROSPECTIVE_CONTRACT="${TWIRL_TEACHER_V3_PROSPECTIVE_CONTRACT:-}"
if (( SECTOR >= 56 && SECTOR <= 62 )); then
  if [[ -n "${PROSPECTIVE_CONTRACT}" ]]; then
    echo "Legacy Teacher-v3 raw export must not set a prospective contract." >&2
    exit 2
  fi
elif (( SECTOR == 63 )); then
  if [[ "${PROSPECTIVE_CONTRACT}" != "${S63_PROSPECTIVE_CONTRACT}" ]]; then
    echo "S63 raw export requires TWIRL_TEACHER_V3_PROSPECTIVE_CONTRACT=${S63_PROSPECTIVE_CONTRACT}." >&2
    exit 64
  fi
  if [[ "${ORBIT_ONE},${ORBIT_TWO}" != "133,134" && "${ORBIT_ONE},${ORBIT_TWO}" != "134,133" ]]; then
    echo "The prospective S63 raw export is locked to orbits 133 and 134." >&2
    exit 64
  fi
else
  echo "Teacher-v3 raw export is bounded to legacy S56-S62 plus explicitly authorized prospective S63." >&2
  exit 2
fi
if (( ORBIT_ONE <= 0 || ORBIT_TWO <= 0 || ORBIT_ONE == ORBIT_TWO )); then
  echo "Orbit IDs must be distinct positive integers." >&2
  exit 2
fi
REPO="${TWIRL_ROOT:-/pdo/users/tehan/TWIRL}"
PDO_USER_ROOT="$(realpath -e -- /pdo/users/tehan)"
REPO="$(realpath -e -- "${REPO}")"
OUT_H5="$(realpath -m -- "${OUT_H5}")"
if [[ "${OUT_H5}" != "${PDO_USER_ROOT}"/* \
   || "${OUT_H5}" == "${REPO}" \
   || "${OUT_H5}" == "${REPO}"/* ]]; then
  echo "PDO raw-export output must be canonical, under ${PDO_USER_ROOT}, and outside the checkout: ${OUT_H5}" >&2
  exit 2
fi
PYTHON_BIN="${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}"
RUN_REL="reports/stage5_validation/teacher_v3_s56_s62_a2v1_current_adp"
if [[ "${SECTOR}" == "63" ]]; then
  TRAINING_TABLE="${TWIRL_TEACHER_V3_SECTOR_TABLE:?set the exact prospective S63 candidate table}"
else
  TRAINING_TABLE="${TWIRL_TEACHER_V3_SECTOR_TABLE:-${REPO}/${RUN_REL}/native_preparation/s${SECTOR}_native_export_table.csv}"
fi
OVERWRITE="${TWIRL_OVERWRITE_TEACHER_V3_RAW_SOURCE:-0}"
if [[ "${SECTOR}" == "63" && "${OVERWRITE}" != "0" ]]; then
  echo "Prospective S63 raw export never permits overwrite." >&2
  exit 4
fi

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
if [[ "${SECTOR}" == "63" ]] && {
  [[ -e "${OUT_H5}" ]] ||
  [[ -e "${OUT_H5%.h5}.summary.json" ]] ||
  [[ -e "${OUT_H5}.sha256" ]];
}; then
  echo "Refusing stale/existing prospective S63 raw-export output." >&2
  exit 4
elif [[ -e "${OUT_H5}" && "${OVERWRITE}" != "1" ]]; then
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
if [[ "${SECTOR}" == "63" ]]; then
  if [[ -z "${TWIRL_EXPECTED_GIT_SHA:-}" ]]; then
    echo "Prospective S63 raw export requires TWIRL_EXPECTED_GIT_SHA." >&2
    exit 2
  fi
  "${REPO}/scripts/assert_clean_checkout.sh" "${REPO}"
  export PYTHONNOUSERSITE=1
  "${PYTHON_BIN}" -c 'import astropy, h5py, numpy, pandas; print("Teacher-v3 S63 PDO Python import preflight OK")'
fi

raw_export_args=(
  scripts/stage5_validation/export_teacher_v3_raw_sources.py
  --sector "${SECTOR}"
  --training-table "${TRAINING_TABLE}"
  --raw-root "${RAW_ROOT}"
  --compact-adp-h5 "${COMPACT_ADP}"
  --out-h5 "${OUT_H5}"
  --orbits "${ORBIT_ONE}" "${ORBIT_TWO}"
)
if [[ "${SECTOR}" == "63" ]]; then
  raw_export_args+=(
    --prospective-contract "${PROSPECTIVE_CONTRACT}"
    --producer-git-sha "${TWIRL_EXPECTED_GIT_SHA}"
  )
fi
"${PYTHON_BIN}" "${raw_export_args[@]}"

sha256sum "${OUT_H5}" > "${OUT_H5}.sha256"
echo "[teacher-v3-raw-export] sector=${SECTOR} complete: ${OUT_H5}"
