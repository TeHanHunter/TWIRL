# Franklin S57--S59 accepted morphology return

This directory freezes the exact `3,000`-row queue reviewed by Franklin and
the returned labels accepted by TeHan on 2026-07-21.

## Contract

- The training unit is one sector observation; `source_uid` is therefore
  sector/candidate specific.
- `tic` is the grouping key. All observations of one TIC must remain in the
  same train/validation/test split.
- `human_labeler=franklin` records the original decision.
  `morphology_adjudicator=tehan` and `morphology_review_status` record the
  separate batch-level acceptance.
- The rank-one BLS `period_d`, `t0_bjd`, and `duration_min` are unchanged.
- `reported_period_factor` and `reported_period_status` are raw app metadata,
  not verified period truth. `harmonic_supervision_verified=False` for every
  row, so this return contributes zero harmonic targets while retaining all
  seven folds as model inputs.
- The labels are morphology/enrichment decisions, not confirmed planets, EBs,
  or survey candidates.

## Files

- `frozen_review_queue.csv`: exact public handoff queue; SHA-256
  `6431ceefecb565ba8878e0ea2fec0528cff817f1801dbbac4b4fdc78cdeac7b9`.
- `franklin_labels_returned.csv`: exact returned labels; SHA-256
  `2bd4d86870c70091eb7291ced067c63bc908118fd730083bfb2d12d52c5a09bf`.
- `accepted_morphology_labels.csv`: strict exact-key join with acceptance,
  grouping, target-mask, and model-target columns.
- `handoff_summary.json`: original handoff metadata.
- `summary.json`: normalized return counts, hashes, and target policy.

Regenerate the normalized table from the repository root with:

```bash
python scripts/stage5_validation/ingest_franklin_multisector_labels.py \
  --queue reports/stage5_validation/franklin_s57_s59_label_return_20260721/frozen_review_queue.csv \
  --labels reports/stage5_validation/franklin_s57_s59_label_return_20260721/franklin_labels_returned.csv \
  --out-dir reports/stage5_validation/franklin_s57_s59_label_return_20260721 \
  --source-batch-id s57_s59_franklin_hv2_3k_20260717 \
  --morphology-adjudicator tehan \
  --morphology-accepted-utc 2026-07-21T21:10:38.262082+00:00 \
  --expected-sector-count 57=1000 \
  --expected-sector-count 58=1000 \
  --expected-sector-count 59=1000
```

Native HDF5 paths are intentionally blank until the verified per-sector
seven-fold inputs are assembled on ORCD.
