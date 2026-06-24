# S56 Duration-Aware Recovery Visualizations

BLS input: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/small_pair_200k/injection_bls_recoveries.csv`
LEO comparison input: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_small_pair_200k_review_queue_pdo/review_queue.csv`

## Logistic Boundary Model

- Rows: `10000`
- Recovered exact/top-N/harmonic: `1683`
- AUC: `0.848`
- Log loss: `0.327`
- Converged: `True` in `7` iterations

| feature | raw coefficient |
| --- | ---: |
| `log10_radius2_sqrt_duration_over_period` | `1.181` |
| `tmag` | `-1.427` |

The fitted 50% boundary uses a physically constrained BLS proxy, `R_p^2 * sqrt(duration / period)`, plus `Tmag`. This gives a monotonic radius cutoff in period-duration space and avoids over-interpreting the correlated period-duration-radius sampling as independent physics.

## BLS vs LEO

- Injected rows in LEO queue: `1000`
- BLS exact/top-N/harmonic recovered: `177/1000` = `17.7%`
- LEO PC/FP precision relative to BLS recovery: `98.7%`
- LEO PC/FP recall relative to BLS recovery: `42.9%`

LEO is therefore a high-confidence vetting subset of the BLS recoveries, not yet a near-complete proxy for recovery. BLS recovery is the better statement for whether the injected signal is tracked by the search; LEO recovery is a stricter downstream pass/fail.

The metric diagnostic supports targeted WD retuning, not a blind loosening of every LEO cut. BLS-recovered LEO-FA rows sit between the BLS misses and the LEO PC/FP rows in MES/new-MES, in-transit cadence count, and shape diagnostics.

## Artifacts

- `period_radius_boundary_png`: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/period_radius_50pct_boundary_by_duration.png`
- `period_radius_boundary_pdf`: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/period_radius_50pct_boundary_by_duration.pdf`
- `period_radius_recovery_map_png`: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/period_radius_recovery_fraction_by_duration_tmag.png`
- `period_radius_recovery_map_pdf`: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/period_radius_recovery_fraction_by_duration_tmag.pdf`
- `period_radius_recovery_map_csv`: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/period_radius_recovery_fraction_by_duration_tmag.csv`
- `radius_snr_proxy_png`: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/radius_snr_proxy_recovery_by_tmag.png`
- `radius_snr_proxy_pdf`: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/radius_snr_proxy_recovery_by_tmag.pdf`
- `duration_radius_floor_png`: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/duration_radius_50pct_boundary_floor_maps.png`
- `duration_radius_floor_pdf`: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/duration_radius_50pct_boundary_floor_maps.pdf`
- `bls_vs_leo_png`: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/bls_vs_leo_recovery_contingency.png`
- `bls_vs_leo_pdf`: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/bls_vs_leo_recovery_contingency.pdf`
- `leo_metric_diagnostics_png`: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/leo_bls_metric_diagnostics.png`
- `leo_metric_diagnostics_pdf`: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/leo_bls_metric_diagnostics.pdf`
- `leo_metric_summary_csv`: `/Users/tehan/PycharmProjects/TWIRL/reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/leo_bls_metric_summary.csv`
