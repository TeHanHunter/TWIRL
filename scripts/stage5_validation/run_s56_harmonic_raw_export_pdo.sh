#!/usr/bin/env bash
# Export only the native raw/error host light curves needed by the S56 teacher.
set -euo pipefail

REPO="${TWIRL_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
PYTHON_BIN="${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}"
ROOT="${TWIRL_ADJUDICATION_ROOT:-reports/stage5_validation/s56_label_adjudication_real343}"
TRAINING_TABLE="${TWIRL_HARMONIC_TRAINING_TABLE:-${ROOT}/adjudicated_training_table/human_vetting_training_table_adjudicated.csv}"
RAW_ROOT="${TWIRL_HARMONIC_PDO_RAW_ROOT:-/pdo/users/tehan/tglc-gpu-production}"
COMPACT_ADP="${TWIRL_HARMONIC_COMPACT_ADP_H5:-data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp_lc_export_pdo.h5}"
OUT_H5="${TWIRL_HARMONIC_RAW_SOURCE_H5:-${ROOT}/native_inputs/s56_tglc_raw_sources.h5}"
OVERWRITE="${TWIRL_OVERWRITE_HARMONIC_RAW_SOURCE:-0}"

cd "${REPO}"
export PYTHONPATH="${PYTHONPATH:-src}"
export LD_LIBRARY_PATH="/sw/python-versions/python-3.11.9/lib:/pdo/app/anaconda/anaconda2-4.4.0/lib:${LD_LIBRARY_PATH:-}"
export HDF5_USE_FILE_LOCKING="${HDF5_USE_FILE_LOCKING:-FALSE}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

if [[ ! -s "${TRAINING_TABLE}" ]]; then
  printf 'Missing adjudicated training table: %s\n' "${TRAINING_TABLE}" >&2
  exit 2
fi
if [[ ! -d "${RAW_ROOT}/orbit-119/ffi" || ! -d "${RAW_ROOT}/orbit-120/ffi" ]]; then
  printf 'Missing matching S56 raw TGLC roots under: %s\n' "${RAW_ROOT}" >&2
  exit 3
fi
if [[ ! -s "${COMPACT_ADP}" ]]; then
  printf 'Missing authoritative compact ADP product: %s\n' "${COMPACT_ADP}" >&2
  exit 5
fi
if [[ -e "${OUT_H5}" && "${OVERWRITE}" != "1" ]]; then
  printf 'Refusing to overwrite %s; set TWIRL_OVERWRITE_HARMONIC_RAW_SOURCE=1 to rebuild.\n' "${OUT_H5}" >&2
  exit 4
fi

mkdir -p "$(dirname "${OUT_H5}")"
printf '[harmonic-raw-export] training_table=%s\n' "${TRAINING_TABLE}"
printf '[harmonic-raw-export] raw_root=%s\n' "${RAW_ROOT}"
printf '[harmonic-raw-export] compact_adp=%s\n' "${COMPACT_ADP}"
printf '[harmonic-raw-export] out_h5=%s\n' "${OUT_H5}"

"${PYTHON_BIN}" scripts/stage5_validation/export_s56_harmonic_raw_sources.py \
  --training-table "${TRAINING_TABLE}" \
  --raw-root "${RAW_ROOT}" \
  --compact-adp-h5 "${COMPACT_ADP}" \
  --out-h5 "${OUT_H5}" \
  --orbits 119,120

sha256sum "${OUT_H5}" > "${OUT_H5}.sha256"
printf '[harmonic-raw-export] complete: %s\n' "${OUT_H5}"
