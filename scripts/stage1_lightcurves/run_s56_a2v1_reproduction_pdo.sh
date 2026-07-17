#!/bin/bash
# S56 convenience wrapper for the generic A2v1 PDO reproduction launcher.

set -uo pipefail

REPO=${TWIRL_REPO:-/pdo/users/tehan/TWIRL}
SCRIPT_DIR=${TWIRL_SCRIPT_DIR:-$REPO/scripts/stage1_lightcurves}
SECTOR=${TWIRL_A2V1_SECTOR:-56}

if [ "$#" -gt 0 ]; then
  exec bash "$SCRIPT_DIR/run_a2v1_reproduction_pdo.sh" "$@"
fi

if [ -n "${TWIRL_A2V1_ORBIT_SPECS:-}" ]; then
  read -r -a ORBIT_SPECS <<< "$TWIRL_A2V1_ORBIT_SPECS"
elif [ "$SECTOR" = "56" ]; then
  ORBIT_SPECS=("119:o1" "120:o2")
else
  echo "ABORT: TWIRL_A2V1_ORBIT_SPECS is required when overriding SECTOR=$SECTOR" >&2
  exit 1
fi

exec bash "$SCRIPT_DIR/run_a2v1_reproduction_pdo.sh" "$SECTOR" "${ORBIT_SPECS[@]}"
