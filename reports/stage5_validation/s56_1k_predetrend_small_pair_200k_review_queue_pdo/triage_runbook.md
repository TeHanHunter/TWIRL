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

Summarize the labels and join them back to hidden truth metadata:

```bash
PYTHONPATH=src .venv/bin/python scripts/stage5_validation/summarize_human_vetting_labels.py \
  --queue-csv reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_pdo/review_queue.csv \
  --labels-csv reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_pdo/human_labels_vetted.csv \
  --out-dir reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_pdo/human_label_summary
```

The next scientific gate is the label comparison against injected truth, BLS recovery mode, empirical SNR, period, radius, depth, and Tmag bins.
