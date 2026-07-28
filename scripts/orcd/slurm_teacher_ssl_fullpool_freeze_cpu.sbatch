#!/bin/bash
#SBATCH -p pg_mki_aryeh
#SBATCH --job-name=twirl-ssl-pool-freeze
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH -o /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.out
#SBATCH -e /orcd/data/mki_aryeh/001/twirl/logs/%x-%j.err

set -euo pipefail

TWIRL_ROOT="${TWIRL_ORCD_ROOT:-/orcd/data/mki_aryeh/001/twirl}"
REPO="${TWIRL_ORCD_REPO:-${TWIRL_ROOT}/code/TWIRL}"
PYTHON="${TWIRL_ORCD_PYTHON:-${TWIRL_ROOT}/envs/twirl-s56-torch/bin/python}"
RUN_ROOT="${TWIRL_SSL_FULLPOOL_RUN_ROOT:-${TWIRL_ROOT}/reports/stage5_validation/teacher_ssl_fullpool_v1_s56_s62_a2v1_current_adp}"
OUT_DIR="${TWIRL_SSL_FULLPOOL_POOL_DIR:-${RUN_ROOT}/frozen/full_pool}"
SPLIT_REGISTRY="${TWIRL_SSL_FULLPOOL_SPLIT_REGISTRY:-${TWIRL_ROOT}/reports/stage5_validation/teacher_v3_s56_s62_a2v1_current_adp/frozen/tic_split_registry.csv}"
S63_TICS="${TWIRL_SSL_FULLPOOL_S63_TICS:-${RUN_ROOT}/preregistered/s63_reserved_tics.txt}"
S63_SUMMARY="${TWIRL_SSL_FULLPOOL_S63_SUMMARY:-${RUN_ROOT}/preregistered/s63_reserved_tics.summary.json}"

S56_MANIFEST="${TWIRL_SSL_FULLPOOL_S56_MANIFEST:-${TWIRL_ROOT}/data_local/stage1_lightcurves/s56_A2v1/s56_A2v1_adp_pair_rebuilt.manifest.json}"
S57_MANIFEST="${TWIRL_SSL_FULLPOOL_S57_MANIFEST:-${TWIRL_ROOT}/data_local/stage1_lightcurves/s57_A2v1/s57_A2v1_adp_pair.manifest.json}"
S58_MANIFEST="${TWIRL_SSL_FULLPOOL_S58_MANIFEST:-${TWIRL_ROOT}/data_local/stage1_lightcurves/s58_A2v1/s58_A2v1_adp_pair.manifest.json}"
S59_MANIFEST="${TWIRL_SSL_FULLPOOL_S59_MANIFEST:-${TWIRL_ROOT}/data_local/stage1_lightcurves/s59_A2v1/s59_A2v1_adp_pair.manifest.json}"
S60_MANIFEST="${TWIRL_SSL_FULLPOOL_S60_MANIFEST:-${TWIRL_ROOT}/reports/stage5_validation/s60_s62_franklin_enrichment_20260718/inputs/s60_A2v1_adp_pair.manifest.json}"
S61_MANIFEST="${TWIRL_SSL_FULLPOOL_S61_MANIFEST:-${TWIRL_ROOT}/reports/stage5_validation/s60_s62_franklin_enrichment_20260718/inputs/s61_A2v1_adp_pair.manifest.json}"
S62_MANIFEST="${TWIRL_SSL_FULLPOOL_S62_MANIFEST:-${TWIRL_ROOT}/reports/stage5_validation/s60_s62_franklin_enrichment_20260718/inputs/s62_A2v1_adp_pair.manifest.json}"
S56_H5="${TWIRL_SSL_FULLPOOL_S56_H5:-${TWIRL_ROOT}/data_local/stage1_lightcurves/s56_A2v1/s56_A2v1_adp_pair_rebuilt.h5}"

inputs=(
  "${S56_MANIFEST}"
  "${S57_MANIFEST}"
  "${S58_MANIFEST}"
  "${S59_MANIFEST}"
  "${S60_MANIFEST}"
  "${S61_MANIFEST}"
  "${S62_MANIFEST}"
  "${S56_H5}"
  "${SPLIT_REGISTRY}"
  "${S63_TICS}"
  "${S63_SUMMARY}"
)
for input in "${inputs[@]}"; do
  if [[ ! -s "${input}" ]]; then
    echo "Missing full-pool freeze input: ${input}" >&2
    exit 3
  fi
done

mkdir -p "${TWIRL_ROOT}/logs"
cd "${REPO}"
export PYTHONPATH="${REPO}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

echo "[ssl-full-pool-freeze] host=$(hostname) start=$(date -Is)"
"${PYTHON}" scripts/stage5_validation/build_teacher_ssl_full_pool.py \
  --compact-manifest "${S56_MANIFEST}" \
  --compact-manifest "${S57_MANIFEST}" \
  --compact-manifest "${S58_MANIFEST}" \
  --compact-manifest "${S59_MANIFEST}" \
  --compact-manifest "${S60_MANIFEST}" \
  --compact-manifest "${S61_MANIFEST}" \
  --compact-manifest "${S62_MANIFEST}" \
  --compact-h5 "56=${S56_H5}" \
  --split-registry "${SPLIT_REGISTRY}" \
  --s63-reserved-tics "${S63_TICS}" \
  --s63-reserved-summary "${S63_SUMMARY}" \
  --out-dir "${OUT_DIR}"
echo "[ssl-full-pool-freeze] complete=$(date -Is) output=${OUT_DIR}"
