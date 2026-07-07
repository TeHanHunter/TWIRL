#!/usr/bin/env bash
# Serve the S56 recovery-boundary teacher-check queue locally.
set -euo pipefail

cd "${TWIRL_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"

export QUEUE_DIR="${QUEUE_DIR:-reports/stage5_validation/s56_recovery50_teacher_queue}"
export CANDIDATES="${CANDIDATES:-${QUEUE_DIR}/review_queue_1k.csv}"
export LABELS_OUT="${LABELS_OUT:-${QUEUE_DIR}/human_labels_vetted.csv}"
export TWIRL_VET_ROOT="${TWIRL_VET_ROOT:-${QUEUE_DIR}/twirl_vet_sheets}"
export HLSP_ROOT="${HLSP_ROOT:-data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare}"
export APERTURE="${APERTURE:-DET_FLUX_ADP_SML}"

exec scripts/stage5_validation/run_s56_mixed_teacher_vetting_app_local.sh "$@"
