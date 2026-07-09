#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
exec python3 franklin_vetting_app.py \
  --queue franklin_review_queue_5k_real.csv \
  --labels-out franklin_labels_vetted.csv \
  --sheet-root vet_sheets \
  --labeler franklin \
  --host 127.0.0.1 \
  --port 5003
