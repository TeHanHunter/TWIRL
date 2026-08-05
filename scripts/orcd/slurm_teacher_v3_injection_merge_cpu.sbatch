#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --exclude=node4900
#SBATCH --job-name=twirl-tv3-injmerge
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 06:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:-${TWIRL_ROOT}/code/TWIRL}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python}"
RUN_REL="reports/stage5_validation/teacher_v3_s56_s62_a2v1_current_adp"
SELECTION_TABLE="${TWIRL_TEACHER_V3_S56_TABLE:-${TWIRL_ROOT}/${RUN_REL}/native_preparation/s56_native_export_table.csv}"
SOURCE_ROOT="${TWIRL_TEACHER_V3_INJECTION_SOURCE_ROOT:?set the stable root containing the frozen data_local injection-pair paths}"
PAIR_300="${SOURCE_ROOT}/data_local/stage3_injections/s56_twirlfs_v2_injection_training/recovery50_adp_pair_subset/injected_lightcurves.h5"
PAIR_156="${SOURCE_ROOT}/data_local/stage3_injections/s56_twirlfs_v2_injection_training/recovery50_adp_pair_subset_next4k/injected_lightcurves.h5"
CANONICAL="${TWIRL_TEACHER_V3_CANONICAL_INJECTION_H5:?set the stable canonical S56 raw-injection HDF5}"
OUTPUT="${TWIRL_TEACHER_V3_MERGED_INJECTION_PAIR_H5:-${TWIRL_ROOT}/${RUN_REL}/native_preparation/s56_teacher_v3_injection_pairs_456.h5}"
SUMMARY="${OUTPUT%.h5}.summary.json"

for input in \
  "${SELECTION_TABLE}" \
  "${PAIR_300}" \
  "${PAIR_156}" \
  "${CANONICAL}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing Teacher-v3 S56 injection input: ${input}" >&2
    exit 3
  fi
done
if [[ -e "${OUTPUT}" || -e "${SUMMARY}" ]]; then
  echo "Refusing to overwrite an existing merged injection artifact." >&2
  echo "output=${OUTPUT} summary=${SUMMARY}" >&2
  exit 4
fi

mkdir -p "${TWIRL_ROOT}/logs" "$(dirname "${OUTPUT}")"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

echo "[teacher-v3-injection-merge] start=$(date -Is)"
"${PYTHON}" scripts/stage5_validation/promote_teacher_v3_injection_pairs.py \
  --selection-table "${SELECTION_TABLE}" \
  --source-pair-h5 "${PAIR_300}" "${PAIR_156}" \
  --canonical-injection-h5 "${CANONICAL}" \
  --out-h5 "${OUTPUT}" \
  --summary-json "${SUMMARY}" \
  --repo-root "${SOURCE_ROOT}"
sha256sum "${OUTPUT}" "${SUMMARY}"
echo "[teacher-v3-injection-merge] complete=$(date -Is) output=${OUTPUT}"
