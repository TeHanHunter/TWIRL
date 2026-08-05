#!/bin/bash
# Run a serial, gated A2v1 HDF5-to-FITS sector queue on PDO.
#
# Usage:
#   run_a2v1_sector_queue_pdo.sh <queue-file>
#
# Queue-file format:
#   <sector> <orbit-one>:<tag> <orbit-two>:<tag>
#
# Each sector must pass prepared-input preflight, HDF5 coverage/openability
# validation, FITS production, and full schema validation before the next sector starts. Existing
# HDF5 or FITS products that already pass their gate are retained and skipped.

set -euo pipefail

QUEUE_FILE=${1:?usage: $0 <queue-file>}
REPO=${TWIRL_REPO:-/pdo/users/tehan/TWIRL}
export TWIRL_REPO="$REPO"
SOURCE_ROOT=${TWIRL_SOURCE_ROOT:-/pdo/users/tehan/tglc-gpu-production}
A2V1_ROOT=${TWIRL_A2V1_ROOT:-/pdo/users/tehan/tglc-gpu-production-A2v1}
QLP_PY=${TWIRL_QLP_PY:-/sw/qlp-environment/.venv/bin/python}
SCRIPT_DIR=${TWIRL_SCRIPT_DIR:-$REPO/scripts/stage1_lightcurves}
H5_RUNNER=$SCRIPT_DIR/run_a2v1_reproduction_pdo.sh
FITS_RUNNER=$SCRIPT_DIR/run_a2v1_hlsp_pdo.sh
VALIDATOR=$SCRIPT_DIR/validate_a2v1_product.py
LOG_DIR=${TWIRL_A2V1_LOG_DIR:-$A2V1_ROOT/twirl_logs}
QUEUE_LABEL=${TWIRL_A2V1_QUEUE_LABEL:-a2v1-sector-queue}
QUEUE_LOG=$LOG_DIR/${QUEUE_LABEL}.log
OBSERVATIONS=${TWIRL_OBSERVATIONS:-$REPO/data_local/catalogs/twirl_master_catalog/twirl_wd_tess_observations_v0.fits}

"$REPO/scripts/assert_clean_checkout.sh" "$REPO"

if [ ! -r "$QUEUE_FILE" ]; then
  echo "queue file is not readable: $QUEUE_FILE" >&2
  exit 1
fi
if [ ! -r "$OBSERVATIONS" ]; then
  echo "observation table is not readable: $OBSERVATIONS" >&2
  exit 1
fi

for path in "$H5_RUNNER" "$FITS_RUNNER" "$VALIDATOR"; do
  if [ ! -f "$path" ]; then
    echo "required script is missing: $path" >&2
    exit 1
  fi
done

mkdir -p "$LOG_DIR"

export HDF5_USE_FILE_LOCKING=${HDF5_USE_FILE_LOCKING:-FALSE}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export VECLIB_MAXIMUM_THREADS=${VECLIB_MAXIMUM_THREADS:-1}
export NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS:-1}
export LD_LIBRARY_PATH=/sw/python-versions/python-3.11.9/lib:/pdo/app/anaconda/anaconda2-4.4.0/lib:${LD_LIBRARY_PATH:-}

stamp() { date '+%Y-%m-%d %H:%M:%S %Z'; }
# Detached tmux panes do not consume unbounded stdout. Keep verbose stage output
# in the persistent queue log so a JSON report cannot block the queue on tty I/O.
log() { printf '[%s] %s\n' "$(stamp)" "$*" >> "$QUEUE_LOG"; }
die() { log "ABORT: $*"; exit 1; }

parse_orbit_spec() {
  local spec=$1
  if [[ "$spec" != *:* ]]; then
    die "orbit spec must be '<orbit>:<tag>'; got '$spec'"
  fi
  PARSED_ORBIT=${spec%%:*}
  PARSED_TAG=${spec#*:}
  if ! [[ "$PARSED_ORBIT" =~ ^[0-9]+$ ]] || [ -z "$PARSED_TAG" ]; then
    die "invalid orbit spec '$spec'"
  fi
}

prepared_count() {
  local orbit=$1
  local kind=$2
  local pattern=$3
  local root=$SOURCE_ROOT/orbit-$orbit/ffi
  local count=0
  local directory
  local n
  local directories=()

  shopt -s nullglob
  directories=("$root"/cam*/ccd*/"$kind")
  shopt -u nullglob

  # Source-only prepared sectors legitimately have no legacy ePSF directories.
  # A partial legacy tree remains unsafe because it silently changes reuse mode.
  if [ "$kind" = "epsf" ] && [ "${#directories[@]}" -eq 0 ]; then
    printf '0\n'
    return
  fi
  if [ "${#directories[@]}" -ne 16 ]; then
    die "orbit-$orbit has ${#directories[@]} $kind directories; expected 16"
  fi
  for directory in "${directories[@]}"; do
    [ -d "$directory" ] || die "missing prepared directory: $directory"
    n=$(find "$directory" -maxdepth 1 -type f -name "$pattern" -printf . | wc -c)
    count=$((count + n))
  done
  printf '%s\n' "$count"
}

preflight_orbit() {
  local orbit=$1
  local source_count
  local epsf_count
  local epsf_mode
  source_count=$(prepared_count "$orbit" source 'source_*.pkl')
  epsf_count=$(prepared_count "$orbit" epsf 'epsf_*.npy')
  [ "$source_count" -eq 3136 ] || die "orbit-$orbit has incomplete source preparation"
  if [ "$epsf_count" -eq 3136 ]; then
    epsf_mode=full-reuse-candidate
  elif [ "$epsf_count" -eq 0 ]; then
    epsf_mode=refit-all
  else
    die "orbit-$orbit has partial old-ePSF preparation ($epsf_count/3136)"
  fi
  log "preflight orbit-$orbit source=$source_count epsf=$epsf_count epsf_mode=$epsf_mode"
}

validate_h5() {
  local sector=$1
  local orbit_one=$2
  local orbit_two=$3
  local report=$LOG_DIR/s$(printf '%04d' "$sector")_A2v1_h5_validation.json
  log "validating sector $sector HDF5 checkpoint"
  HDF5_USE_FILE_LOCKING=FALSE "$QLP_PY" "$VALIDATOR" \
    --a2v1-root "$A2V1_ROOT" \
    --observations "$OBSERVATIONS" \
    --sector "$sector" --orbits "$orbit_one" "$orbit_two" \
    --allow-edge-warn-missing --skip-fits --check-h5-open \
    --summary-json "$report" | tee -a "$QUEUE_LOG" >/dev/null
}

validate_product() {
  local sector=$1
  local orbit_one=$2
  local orbit_two=$3
  local report=$LOG_DIR/s${sector}_A2v1_validation_full_schema.json
  log "validating sector $sector HDF5-plus-FITS product"
  HDF5_USE_FILE_LOCKING=FALSE "$QLP_PY" "$VALIDATOR" \
    --a2v1-root "$A2V1_ROOT" \
    --observations "$OBSERVATIONS" \
    --sector "$sector" --orbits "$orbit_one" "$orbit_two" \
    --allow-edge-warn-missing --schema-only --check-h5-open \
    --fits-workers "${TWIRL_A2V1_VALIDATION_WORKERS:-8}" \
    --summary-json "$report" | tee -a "$QUEUE_LOG" >/dev/null
}

run_h5() {
  local sector=$1
  local spec_one=$2
  local spec_two=$3
  log "starting sector $sector A2v1 HDF5 production"
  HDF5_USE_FILE_LOCKING=FALSE \
  TWIRL_A2V1_PREFILL_EMPTY_MASK_EPSFS="${TWIRL_A2V1_PREFILL_EMPTY_MASK_EPSFS:-1}" \
  TWIRL_A2V1_GPU_LIST="${TWIRL_A2V1_GPU_LIST:-0,1,2,3,6,7}" \
  TWIRL_A2V1_GPU_MAX_PARALLEL="${TWIRL_A2V1_GPU_MAX_PARALLEL:-4}" \
  TWIRL_A2V1_EPSFS_NPROCS="${TWIRL_A2V1_EPSFS_NPROCS:-2}" \
  TWIRL_A2V1_RUN_LABEL="${QUEUE_LABEL}-s${sector}" \
  bash "$H5_RUNNER" "$sector" "$spec_one" "$spec_two" | tee -a "$QUEUE_LOG" >/dev/null
}

run_fits() {
  local sector=$1
  local orbit_one=$2
  local orbit_two=$3
  log "starting sector $sector A2v1 FITS production"
  HDF5_USE_FILE_LOCKING=FALSE \
  TWIRL_A2V1_HLSP_WORKERS="${TWIRL_A2V1_HLSP_WORKERS:-8}" \
  bash "$FITS_RUNNER" "$sector" "$orbit_one" "$orbit_two" | tee -a "$QUEUE_LOG" >/dev/null
}

process_sector() {
  local sector=$1
  local spec_one=$2
  local spec_two=$3
  parse_orbit_spec "$spec_one"
  local orbit_one=$PARSED_ORBIT
  parse_orbit_spec "$spec_two"
  local orbit_two=$PARSED_ORBIT

  log "sector $sector queued with $spec_one $spec_two"
  preflight_orbit "$orbit_one"
  preflight_orbit "$orbit_two"

  if validate_h5 "$sector" "$orbit_one" "$orbit_two"; then
    log "sector $sector HDF5 checkpoint already passes"
  else
    log "sector $sector HDF5 checkpoint incomplete; running extraction"
    run_h5 "$sector" "$spec_one" "$spec_two"
    validate_h5 "$sector" "$orbit_one" "$orbit_two" || die "sector $sector HDF5 validation failed"
  fi

  if validate_product "$sector" "$orbit_one" "$orbit_two"; then
    log "sector $sector full A2v1 product already passes"
  else
    log "sector $sector FITS product incomplete; running FITS export"
    run_fits "$sector" "$orbit_one" "$orbit_two"
    validate_product "$sector" "$orbit_one" "$orbit_two" || die "sector $sector full validation failed"
  fi
  log "sector $sector complete"
}

log "A2v1 queue start file=$QUEUE_FILE"
line_number=0
while read -r sector spec_one spec_two extra; do
  line_number=$((line_number + 1))
  if [ -z "${sector:-}" ] || [[ "$sector" == \#* ]]; then
    continue
  fi
  [ -z "${spec_one:-}" ] && die "line $line_number is missing first orbit spec"
  [ -z "${spec_two:-}" ] && die "line $line_number is missing second orbit spec"
  [ -z "${extra:-}" ] || die "line $line_number has unexpected extra fields"
  [[ "$sector" =~ ^[0-9]+$ ]] || die "line $line_number has non-numeric sector '$sector'"
  process_sector "$sector" "$spec_one" "$spec_two"
done < "$QUEUE_FILE"
log "A2v1 queue complete"
