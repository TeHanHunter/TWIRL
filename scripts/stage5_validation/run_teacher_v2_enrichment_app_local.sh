#!/bin/bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
QUEUE_ROOT="${1:-${TWIRL_ENRICH_LOCAL_ROOT:-${REPO_ROOT}/reports/stage5_validation/s56_s64_existing_teacher_enrichment/sector_0056/batch_00}}"
PORT="${2:-${TWIRL_ENRICH_PORT:-5005}}"
PYTHON_BIN="${TWIRL_PYTHON:-python}"

exec "${PYTHON_BIN}" "${REPO_ROOT}/scripts/stage5_validation/run_lightcurve_vetting_app.py" \
  --candidates "${QUEUE_ROOT}/review_queue_1k.csv" \
  --labels-out "${QUEUE_ROOT}/human_labels_vetted.csv" \
  --hlsp-root "${QUEUE_ROOT}/unused_hlsp_fallback" \
  --twirl-vet-root "${QUEUE_ROOT}/twirl_vet_sheets" \
  --host 127.0.0.1 \
  --port "${PORT}" \
  --aperture DET_FLUX_ADP_SML
