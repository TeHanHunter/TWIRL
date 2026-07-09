# TWIRL S56 Real-Candidate Vetting Handoff

This package is a real-only S56 review queue for Franklin.

## What Is Included

- `franklin_review_queue_5k_real.csv`: 5,000 real S56 candidate rows.
- `vet_sheets/`: pre-rendered TWIRL two-aperture PNG sheets for the browser app.
- `franklin_vetting_app.py`: standalone local browser app.
- `run_franklin_vetting.sh`: one-command launcher.
- `reference_examples/`: examples from TeHan's real-data labels.
- `reference_examples.csv`: table linking examples to labels and TICs.
- `franklin_labels_vetted.csv`: created by the app as labels are saved.

The main queue intentionally contains no injected rows. Injection truth is not
part of Franklin's label task.

## Start Labeling

From the unzipped package directory:

```bash
python3 franklin_vetting_app.py --check-only
./run_franklin_vetting.sh
```

Then open:

```text
http://127.0.0.1:5003/
```

Labels are saved immediately to `franklin_labels_vetted.csv`.

## Label Meanings

- `planet_like`: preserve a compact transit-like signal.
- `wide_transit_like`: preserve a broad/long-duration transit-like signal, but
  do not merge it blindly with compact planet morphology.
- `eclipsing_binary_or_pceb`: EB/PCEB-like shape, odd/even mismatch, secondary
  eclipse, or stellar-companion morphology.
- `stellar_variability`: astrophysical variability or repeating non-transit
  structure.
- `instrumental_or_systematic`: BLS is catching window-edge leakage, cadence
  artifacts, extraction artifacts, or other non-astrophysical structure.
- `uncertain`: flat/no obvious useful event. TeHan used this heavily for flats.
- `skip`: unusable row or app/sheet problem.

If a row looks period-folded at the wrong harmonic, keep the preserve/reject
label you think is scientifically correct and write a short note such as
`half period`, `refold at P/2`, or `possible harmonic`.

## Keyboard Shortcuts

- `1`: planet_like
- `6`: wide_transit_like
- `2`: eclipsing_binary_or_pceb
- `3`: stellar_variability
- `4`: instrumental_or_systematic
- `5`: uncertain
- `0`: skip
- left/right arrows: previous/next row

## Provenance

Created UTC: `2026-07-08T21:04:44.114310+00:00`

Summary:

```json
{
  "created_utc": "2026-07-08T21:04:44.113381+00:00",
  "example_summary": {
    "example_counts": {
      "eclipsing_binary_or_pceb": 2,
      "instrumental_or_systematic": 2,
      "planet_like": 2,
      "skip": 2,
      "stellar_variability": 2,
      "uncertain": 2,
      "wide_transit_like": 2
    },
    "examples_per_label_requested": 2,
    "missing_example_sheets": []
  },
  "labels_csv": "/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_franklin_real5k_handoff/franklin_labels_vetted.csv",
  "n_available_after_excluding_tehan_labeled_real": 7762,
  "n_excluded_tehan_labeled_real_review_ids": 1543,
  "n_reference_examples": 14,
  "n_review": 5000,
  "n_source_pool": 9000,
  "queue_csv": "/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_franklin_real5k_handoff/franklin_review_queue_5k_real.csv",
  "real_summary": {
    "cadence_alias_selected": 568,
    "duration_q80_min": 30.0,
    "input_rows": 14545,
    "period_cluster_q95": 24.0,
    "sde_quantiles": {
      "q40": 8.87688500234097,
      "q70": 11.199771136272103,
      "q90": 19.272529469949657
    },
    "selected_rows": 9000,
    "selection_bucket_available_counts": {
      "real_aperture_disagreement": 12896,
      "real_cadence_alias_systematic": 1386,
      "real_eb_pceb_like": 10152,
      "real_high_sde_planet_like": 3007,
      "real_low_sde_control": 5535,
      "real_mid_sde_control": 3962,
      "real_variability_broad_duration": 7622,
      "real_wd1856_benchmark": 1
    },
    "selection_bucket_counts": {
      "real_aperture_disagreement": 1000,
      "real_cadence_alias_systematic": 700,
      "real_eb_pceb_like": 2200,
      "real_high_sde_planet_like": 1800,
      "real_low_sde_control": 799,
      "real_mid_sde_control": 1400,
      "real_variability_broad_duration": 1100,
      "real_wd1856_benchmark": 1
    }
  },
  "selection_bucket_counts": {
    "real_aperture_disagreement": 602,
    "real_cadence_alias_systematic": 374,
    "real_eb_pceb_like": 1180,
    "real_high_sde_planet_like": 888,
    "real_low_sde_control": 488,
    "real_mid_sde_control": 819,
    "real_variability_broad_duration": 649
  },
  "source_kind_counts": {
    "real_candidate": 5000
  },
  "verification_failures": [],
  "verification_passed": true
}
```
