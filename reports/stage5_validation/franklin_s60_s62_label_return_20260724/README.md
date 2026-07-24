# Franklin S60--S62 morphology return

This directory freezes the exact `3,000`-row queue reviewed by Franklin and
the complete label return received on 2026-07-24. The batch-level acceptance
means that the return is structurally valid and may be staged for final
review; it is not row-level scientific adjudication. TeHan's final review of
every planet-like and eclipsing-binary row supersedes these provisional
morphology labels.

## Contract

- The training unit is one sector observation; `source_uid` is therefore
  sector/candidate specific.
- `tic` is the grouping key. All observations of one TIC must remain in the
  same train/validation/test split.
- `human_labeler=franklin` records the original decision.
  `morphology_adjudicator=tehan` and `morphology_review_status` record the
  separate batch-level staging decision.
- The rank-one BLS `period_d`, `t0_bjd`, and `duration_min` are unchanged.
- `reported_period_factor` and `reported_period_status` are raw app metadata,
  not verified period truth. `harmonic_supervision_verified=False` for every
  row, so this return contributes zero harmonic targets while retaining all
  seven folds as model inputs.
- The labels are morphology/enrichment decisions, not confirmed planets, EBs,
  or survey candidates.

## Return summary

- `3,000` complete labels: `1,000` in each of S60, S61, and S62.
- `27` planet-like and `259` eclipsing-binary/PCEB rows advance to TeHan's
  final row-level signal review.
- The remaining labels are `215` wide/broad transit-like, `1,449` uncertain,
  `602` stellar variability, `442` instrumental/systematic, and `6` skipped.

## Files

- `frozen_review_queue.csv`: exact public handoff queue; SHA-256
  `01f1e627b4ba8f38fd30b4c1ef4b0172fa2a8577f975276aed9876a91f0a007f`.
- `franklin_labels_returned.csv`: exact returned labels; SHA-256
  `852a313544b9a30edc8cd131c07a5977f7027f4c026ec92e91362a5df604a785`.
- `accepted_morphology_labels.csv`: strict exact-key join with acceptance,
  grouping, target-mask, and model-target columns; SHA-256
  `3a5ddbb428ca59855bd7de76afd2540e7800b090f80c7d72ec726ea052e88fad`.
- `handoff_summary.json`: original PDO handoff metadata.
- `summary.json`: normalized return counts, hashes, and target policy.

Regenerate the normalized table from the repository root with:

```bash
python scripts/stage5_validation/ingest_franklin_multisector_labels.py \
  --queue reports/stage5_validation/franklin_s60_s62_label_return_20260724/frozen_review_queue.csv \
  --labels reports/stage5_validation/franklin_s60_s62_label_return_20260724/franklin_labels_returned.csv \
  --out-dir reports/stage5_validation/franklin_s60_s62_label_return_20260724 \
  --source-batch-id s60_s62_tehan_hv2_3k_20260721 \
  --morphology-adjudicator tehan \
  --morphology-accepted-utc 2026-07-24T18:15:04.837695+00:00 \
  --expected-sector-count 60=1000 \
  --expected-sector-count 61=1000 \
  --expected-sector-count 62=1000
```

Native HDF5 paths are intentionally blank until the verified per-sector
seven-fold inputs are assembled on ORCD.
