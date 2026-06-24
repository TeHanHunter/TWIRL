# LEO WD Threshold Smoke Test

Queue: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_small_pair_200k_review_queue_pdo/review_queue.csv`
Metrics: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_small_pair_200k_review_queue_pdo/leo_metrics.csv`

## Result

| preset                | predicted_pass_n | tp_vs_bls | fp_vs_bls | precision_vs_bls | recall_vs_bls | added_tp_vs_current | added_fp_vs_current |
| --------------------- | ---------------- | --------- | --------- | ---------------- | ------------- | ------------------- | ------------------- |
| current_leo_pc_or_fp  | 77               | 76        | 1         | 98.7%            | 42.9%         | 0                   | 0                   |
| wd_review_high_purity | 97               | 89        | 8         | 91.8%            | 50.3%         | 13                  | 7                   |
| wd_review_balanced    | 114              | 98        | 16        | 86.0%            | 55.4%         | 22                  | 15                  |
| wd_review_aggressive  | 168              | 109       | 59        | 64.9%            | 61.6%         | 33                  | 58                  |

## Interpretation

The current WD-tuned LEO PC/FP label is very high precision but low recall relative to the BLS exact/top-N/harmonic recovery label. A physically motivated WD review pass-through can recover additional BLS-tracked injections, but the purity cost rises quickly if the inherited shape tests are loosened too far.

The most defensible near-term tune is `wd_review_high_purity`: keep `new_N_transit >= 2`, keep the large-object ceiling, require the primary modshift event to dominate secondary/tertiary/positive events, and only relax the under-resolved shape/MES behavior enough to route cases to human review. Do not promote these directly to PC without human labels.

Physical rationale: WD transits last only minutes at 200 s cadence, so trapezoid/planet-shape diagnostics and pruned MES can be unstable when only a few cadences sample each event. Secondary/tertiary modshift peaks, however, remain physically meaningful for rejecting wrong ephemerides, so the tune should use primary dominance rather than dropping uniqueness checks.

## Artifacts

- `scores_csv`: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/leo_wd_tuning_smoke/leo_wd_threshold_smoke_scores.csv`
- `added_review_rows_csv`: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/leo_wd_tuning_smoke/leo_wd_threshold_smoke_added_review_rows.csv`
- `precision_recall_png`: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/leo_wd_tuning_smoke/leo_wd_threshold_smoke_precision_recall.png`
- `precision_recall_pdf`: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/leo_wd_tuning_smoke/leo_wd_threshold_smoke_precision_recall.pdf`
