#!/usr/bin/env bash
# Build one sector's authoritative A2v1 cadence/quality evidence on PDO.
set -euo pipefail

SECTOR="${1:?usage: $0 <sector> <orbit-one> <orbit-two> [output-dir]}"
ORBIT_ONE="${2:?usage: $0 <sector> <orbit-one> <orbit-two> [output-dir]}"
ORBIT_TWO="${3:?usage: $0 <sector> <orbit-one> <orbit-two> [output-dir]}"
OUTPUT_DIR="${4:-/pdo/users/tehan/twirl-teacher-v3/a2v1_cadence_reference/s${SECTOR}}"

if (( SECTOR < 56 || SECTOR > 62 )); then
  echo "Teacher-v3 cadence preparation is bounded to sectors 56-62." >&2
  exit 2
fi
if (( ORBIT_ONE <= 0 || ORBIT_TWO <= 0 || ORBIT_ONE == ORBIT_TWO )); then
  echo "Orbit IDs must be distinct positive integers." >&2
  exit 2
fi
if [[ "${OUTPUT_DIR}" != /pdo/users/tehan/* ]]; then
  echo "Refusing a PDO output outside /pdo/users/tehan/: ${OUTPUT_DIR}" >&2
  exit 2
fi

REPO="${TWIRL_ROOT:-/pdo/users/tehan/TWIRL}"
PYTHON_BIN="${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}"
QLP_ROOT="${TWIRL_QLP_AUTHORITY_ROOT:-/pdo/qlp-data}"
SPOC_ROOT="${TWIRL_SPOC_FLAG_ROOT:-${QLP_ROOT}/spocflags}"
OVERWRITE="${TWIRL_OVERWRITE_CADENCE_REFERENCE:-0}"

SPOC_TABLE="${OUTPUT_DIR}/spoc_quality.csv"
SPOC_PROVENANCE="${OUTPUT_DIR}/spoc_quality.provenance.json"
REFERENCE_TABLE="${OUTPUT_DIR}/cadence_reference.csv"
REFERENCE_MANIFEST="${OUTPUT_DIR}/cadence_reference.manifest.json"

mkdir -p "${OUTPUT_DIR}"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

quat_args=()
spoc_args=()
qflag_args=()
for orbit in "${ORBIT_ONE}" "${ORBIT_TWO}"; do
  for camera in 1 2 3 4; do
    quat="${QLP_ROOT}/orbit-${orbit}/ffi/run/cam${camera}_quat.txt"
    if [[ ! -s "${quat}" ]]; then
      echo "Missing QLP quaternion authority: ${quat}" >&2
      exit 3
    fi
    quat_args+=(--quat "${orbit},${camera},${quat}")
    for ccd in 1 2 3 4; do
      qflag="${QLP_ROOT}/orbit-${orbit}/ffi/run/cam${camera}ccd${ccd}_qflag.txt"
      if [[ ! -s "${qflag}" ]]; then
        echo "Missing QLP detector qflag authority: ${qflag}" >&2
        exit 3
      fi
      qflag_args+=(--qlp-qflag "${orbit},${camera},${ccd},${qflag}")
    done
  done
done
for camera in 1 2 3 4; do
  for ccd in 1 2 3 4; do
    spoc="${SPOC_ROOT}/spocffiflag_s${SECTOR}_cam${camera}_ccd${ccd}.txt"
    if [[ ! -s "${spoc}" ]]; then
      echo "Missing SPOC detector authority: ${spoc}" >&2
      exit 3
    fi
    spoc_args+=(--spoc-flag "${camera},${ccd},${spoc}")
  done
done

overwrite_args=()
if [[ "${OVERWRITE}" == "1" ]]; then
  overwrite_args+=(--overwrite)
fi

"${PYTHON_BIN}" scripts/stage1_lightcurves/build_a2v1_spoc_quality_table.py \
  --sector "${SECTOR}" \
  --expected-orbit "${ORBIT_ONE}" \
  --expected-orbit "${ORBIT_TWO}" \
  "${quat_args[@]}" \
  "${spoc_args[@]}" \
  --output-table "${SPOC_TABLE}" \
  --output-provenance "${SPOC_PROVENANCE}" \
  "${overwrite_args[@]}"

"${PYTHON_BIN}" scripts/stage1_lightcurves/build_a2v1_cadence_reference.py \
  --sector "${SECTOR}" \
  --expected-orbit "${ORBIT_ONE}" \
  --expected-orbit "${ORBIT_TWO}" \
  "${quat_args[@]}" \
  "${qflag_args[@]}" \
  --spoc-quality-table "${SPOC_TABLE}" \
  --spoc-quality-provenance "${SPOC_PROVENANCE}" \
  --output-table "${REFERENCE_TABLE}" \
  --output-manifest "${REFERENCE_MANIFEST}" \
  "${overwrite_args[@]}"

sha256sum "${REFERENCE_TABLE}" "${REFERENCE_MANIFEST}" > \
  "${OUTPUT_DIR}/cadence_reference.sha256"
echo "[a2v1-cadence-reference] sector=${SECTOR} complete: ${OUTPUT_DIR}"
