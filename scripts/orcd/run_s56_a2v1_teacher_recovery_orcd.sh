#!/usr/bin/env bash
# Backward-compatible Sector 56 entry point for the generic A2v1 workflow.
set -euo pipefail

export TWIRL_RECOVERY_SECTOR="${TWIRL_RECOVERY_SECTOR:-56}"
exec "$(dirname "$0")/run_a2v1_teacher_recovery_orcd.sh" "$@"
