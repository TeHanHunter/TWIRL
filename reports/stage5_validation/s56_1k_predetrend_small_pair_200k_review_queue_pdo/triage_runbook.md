# S56 1k Human Triage Runbook

This is the active quick human-check queue for the S56 raw-flux pre-detrend injection path. It uses the small-aperture pair `DET_FLUX_ADP_SML + DET_FLUX_SML`, the dense `200k` BLS grid, and full WD-tuned LEO-Vetter PDFs only.

The queue is injection-only by design: `1000` BATMAN injections, `1000` LEO reports, `0` LEO metric errors, and `0` LEO plot errors. Use it to answer whether the injected signals look visually recoverable in the LEO sheets before scaling human labels to a mixed real-plus-injected queue.

## Start Locally

```bash
cd /Users/tehan/PycharmProjects/TWIRL
bash scripts/stage5_validation/run_s56_1k_small_pair_vetting_app_local.sh
```

Open:

```text
http://127.0.0.1:5006/
```

If `5006` is already in use, choose another local port:

```bash
PORT=5007 bash scripts/stage5_validation/run_s56_1k_small_pair_vetting_app_local.sh
```

## What To Label

The LEO-Vetter PDF is the primary evidence. The TWIRL light-curve panel is fallback/debug context and should not drive the label when the LEO sheet is present.

Quick labels auto-save and advance:

| Key | Label |
|---|---|
| `1` or `p` | planet_like |
| `2` or `e` | eclipsing_binary_or_pceb |
| `3` or `v` | stellar_variability |
| `4` or `i` | instrumental_or_systematic |
| `5` or `u` | uncertain |
| `0`, `s`, or `x` | skip |

Use the note field only for useful scientific comments, for example "visible only in one aperture" or "LEO fold dominated by long trend".

## Outputs

Human labels are written here:

```text
reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_pdo/human_labels_vetted.csv
```

Do not edit `review_queue.csv`; it contains hidden injection truth and recovery metadata for later scoring. The app hides that provenance during triage.

## After A Labeling Pass

Run the post-labeling audit loop:

```bash
cd /Users/tehan/PycharmProjects/TWIRL
bash scripts/stage5_validation/run_s56_human_label_audit_local.sh
```

This writes the label summary, model-ready training table, next-label priority
list, and training-readiness audit without editing `review_queue.csv` or
`human_labels_vetted.csv`.

For debugging, the individual steps are:

Summarize the labels and join them back to hidden truth metadata:

```bash
PYTHONPATH=src .venv/bin/python scripts/stage5_validation/summarize_human_vetting_labels.py \
  --queue-csv reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_pdo/review_queue.csv \
  --labels-csv reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_pdo/human_labels_vetted.csv \
  --out-dir reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_pdo/human_label_summary
```

The next scientific gate is the label comparison against injected truth, BLS recovery mode, empirical SNR, period, radius, depth, and Tmag bins.

Build the model-ready human-label table without editing the label CSV:

```bash
PYTHONPATH=src .venv/bin/python scripts/stage5_validation/build_human_vetting_training_table.py \
  --queue-csv reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_pdo/review_queue.csv \
  --labels-csv reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_pdo/human_labels_vetted.csv \
  --out-dir reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_pdo/human_training_table
```

This writes a full joined table, a strong-label teacher subset, an audit subset
that keeps `uncertain`, and deterministic train/validation/test split markers.

Select a high-value next batch of unlabeled rows:

```bash
PYTHONPATH=src .venv/bin/python scripts/stage5_validation/select_next_human_labels.py \
  --training-table reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_pdo/human_training_table/human_vetting_training_table.csv \
  --out-dir reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_pdo/human_label_priority_next \
  --n-rows 200 \
  --target-per-cell 25
```

This writes `next_label_priority.csv` and `next_label_row_ids.txt`. The priority
table is a planning artifact; keep using the browser app for the actual labels.

Audit whether the current labels are ready for a model:

```bash
PYTHONPATH=src .venv/bin/python scripts/stage5_validation/audit_human_training_readiness.py \
  --training-table reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_pdo/human_training_table/human_vetting_training_table.csv \
  --priority-table reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_pdo/human_label_priority_next/next_label_priority.csv \
  --out-dir reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_pdo/human_training_readiness
```

The audit intentionally separates injected-signal visibility from object-class
training. The current `101`-label quick pass is enough for
`ready_for_injection_visibility_smoke_only`, but it is not a teacher-model
training set because the queue is injection-only and lacks real false-positive
classes. Use the next priority rows to finish the 1k injected visibility check,
then build a mixed real-plus-injected queue after the peak-ranker selects better
real-candidate ephemerides.
