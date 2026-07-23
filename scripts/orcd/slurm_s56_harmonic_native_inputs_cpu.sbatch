#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-harmonic-inputs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH --mem=512G
#SBATCH -t 12:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

REPO="${TWIRL_ORCD_REPO:-/orcd/data/mki_aryeh/001/twirl/code/TWIRL}"
PYTHON="${TWIRL_ORCD_TORCH_PYTHON:-/orcd/data/mki_aryeh/001/twirl/envs/twirl-s56-torch/bin/python}"
ROOT="${TWIRL_ADJUDICATION_ROOT:-reports/stage5_validation/s56_label_adjudication_real343}"
TRAINING="${ROOT}/adjudicated_training_table/human_vetting_training_table_adjudicated.csv"
RAW_SOURCE="${TWIRL_HARMONIC_RAW_SOURCE_H5:-${ROOT}/native_inputs/s56_tglc_raw_sources.h5}"
COMPACT_ADP="${TWIRL_HARMONIC_COMPACT_ADP_H5:-data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp_lc_export_pdo.h5}"
INJECTION_PAIR="${TWIRL_HARMONIC_INJECTION_PAIR_H5:-${ROOT}/native_inputs/s56_adjudicated_injected_adp_pair_v1.h5}"
CADENCE_REFERENCE_TABLE="${TWIRL_A2V1_CADENCE_REFERENCE:?set TWIRL_A2V1_CADENCE_REFERENCE}"
CADENCE_REFERENCE_MANIFEST="${TWIRL_A2V1_CADENCE_REFERENCE_MANIFEST:?set TWIRL_A2V1_CADENCE_REFERENCE_MANIFEST}"
FINAL="${TWIRL_HARMONIC_NATIVE_H5:-${ROOT}/native_inputs/s56_adp_raw_pair_v2.h5}"
N_SHARDS="${TWIRL_HARMONIC_INPUT_SHARDS:-32}"
RUN_DIR="${ROOT}/native_inputs/shards_${SLURM_JOB_ID}"

mkdir -p /orcd/data/mki_aryeh/001/twirl/logs "${RUN_DIR}" "$(dirname "${FINAL}")"
cd "${REPO}"
export PYTHONPATH=src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

if [[ ! -s "${RAW_SOURCE}" ]]; then
  echo "Missing compact raw/error source: ${RAW_SOURCE}" >&2
  echo "Create it from the PDO TGLC files with export_s56_harmonic_raw_sources.py, then stage only that compact HDF5 here." >&2
  exit 3
fi
if [[ ! -s "${INJECTION_PAIR}" ]]; then
  echo "Missing exact adjudicated injection pair product: ${INJECTION_PAIR}" >&2
  exit 6
fi
for evidence in "${CADENCE_REFERENCE_TABLE}" "${CADENCE_REFERENCE_MANIFEST}"; do
  if [[ ! -s "${evidence}" ]]; then
    echo "Missing cadence-reference evidence: ${evidence}" >&2
    exit 7
  fi
done
if [[ -e "${FINAL}" && "${TWIRL_OVERWRITE_HARMONIC_NATIVE:-0}" != "1" ]]; then
  echo "Refusing to overwrite ${FINAL}; set TWIRL_OVERWRITE_HARMONIC_NATIVE=1 for an intentional rebuild." >&2
  exit 4
fi

"${PYTHON}" scripts/stage5_validation/audit_s56_real_label_adjudication.py \
  --root "${ROOT}" \
  --require-complete

pids=()
for ((shard=0; shard<N_SHARDS; shard++)); do
  shard_path="${RUN_DIR}/native_${shard}.h5"
  srun --exclusive --exact -N1 -n1 -c2 --mem=12G "${PYTHON}" \
    scripts/stage5_validation/build_s56_harmonic_native_inputs.py \
    --training-table "${TRAINING}" \
    --raw-source-h5 "${RAW_SOURCE}" \
    --compact-adp-h5 "${COMPACT_ADP}" \
    --injection-pair-h5 "${INJECTION_PAIR}" \
    --cadence-reference-table "${CADENCE_REFERENCE_TABLE}" \
    --cadence-reference-manifest "${CADENCE_REFERENCE_MANIFEST}" \
    --out-h5 "${shard_path}" \
    --repo-root "${REPO}" \
    --n-periods 4096 \
    --shard-index "${shard}" \
    --n-shards "${N_SHARDS}" &
  pids+=("$!")
done
shard_status=0
for pid in "${pids[@]}"; do
  if ! wait "${pid}"; then
    shard_status=1
  fi
done
if [[ "${shard_status}" != "0" ]]; then
  echo "At least one native-input shard failed; refusing to merge partial output." >&2
  exit 5
fi

shards=("${RUN_DIR}"/native_*.h5)
"${PYTHON}" scripts/stage5_validation/merge_s56_harmonic_native_shards.py \
  --training-table "${TRAINING}" \
  --shards "${shards[@]}" \
  --out-h5 "${FINAL}"

echo "[harmonic-inputs] complete $(date -Is): ${FINAL}"
