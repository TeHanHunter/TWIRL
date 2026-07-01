# S56 Peak Handoff Audit

Created UTC: `2026-06-27T00:38:11.579158+00:00`

## Chunked Peak Table

- completed chunks: `0` / `?`
- chunk dirs: `0`
- incomplete preview: `none`

## Gates

| Gate | Ready | Artifact |
|---|---:|---|
| `standard_injected_peak_table` | no | `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training/s56_20k_injection_bls_peaks_chunked_verification.json` |
| `post_bls_peak_gate` | no | `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training_gate_pdo/summary.json` |
| `failure_mode_audit` | no | `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training_gate_pdo/failure_modes/summary.json` |
| `injected_truth_peak_ranker` | no | `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_ranker_pdo/summary.json` |
| `real_bls_peak_table_verified` | no | `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_ranker_selected_real_candidates_pdo/real_peak_table_verification.json` |
| `ranker_selected_real_ephemerides` | no | `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_ranker_selected_real_candidates_pdo/selected_ephemerides.csv` |
| `skip_leo_review_queue_verified` | no | `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_ranker_selected_real_review_queue_pdo/verification.json` |
| `leo_review_queue_verified` | no | `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_ranker_selected_real_leo_queue_pdo/verification.json` |

## Next Action

`standard_injected_peak_table`: Wait for the active chunked build to merge and self-verify.

- ready for branch comparison: `False`
- ready for skip-LEO human queue: `False`
- ready for LEO human queue: `False`
