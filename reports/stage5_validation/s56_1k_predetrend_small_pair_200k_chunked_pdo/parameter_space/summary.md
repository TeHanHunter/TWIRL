# S56 Injection Recovery Parameter Space

Input recovery table: `reports/stage5_validation/s56_1k_predetrend_small_pair_200k_chunked_pdo/small_pair_200k/injection_bls_recoveries.csv`

- Rows plotted: `1000`
- Strict top-rank BLS recoveries: `137/1000` = `13.7%`
- Exact/top-N/harmonic matches: `177/1000` = `17.7%`

## Tmag Bins

| label        | n   | strict_top1_n | strict_top1_frac | any_exact_or_harmonic_n | any_exact_or_harmonic_frac | median_period_d | median_depth_pct |
| ------------ | --- | ------------- | ---------------- | ----------------------- | -------------------------- | --------------- | ---------------- |
| Tmag [0,16)  | 11  | 10            | 0.909            | 10                      | 0.909                      | 0.493           | 23.6             |
| Tmag [16,17) | 18  | 12            | 0.667            | 12                      | 0.667                      | 1.17            | 29.6             |
| Tmag [17,18) | 68  | 29            | 0.426            | 30                      | 0.441                      | 1.7             | 26.9             |
| Tmag [18,19) | 219 | 54            | 0.247            | 67                      | 0.306                      | 1.32            | 34.4             |
| Tmag [19,20) | 659 | 32            | 0.0486           | 57                      | 0.0865                     | 1.24            | 27.8             |
| Tmag [20,30) | 25  | 0             | 0                | 1                       | 0.04                       | 0.987           | 27.2             |

## Interpretation

The recovery boundary should be read conditionally on Tmag and cadence coverage. A single period-depth panel mixes bright and faint WDs, so the Tmag-sliced recovery-fraction plot is the better diagnostic for the gradual BLS boundary.

## Artifacts

- `period_depth_png`: `reports/stage5_validation/s56_1k_predetrend_small_pair_200k_chunked_pdo/parameter_space/period_depth_recovery_by_tmag.png`
- `period_depth_pdf`: `reports/stage5_validation/s56_1k_predetrend_small_pair_200k_chunked_pdo/parameter_space/period_depth_recovery_by_tmag.pdf`
- `period_depth_fraction_png`: `reports/stage5_validation/s56_1k_predetrend_small_pair_200k_chunked_pdo/parameter_space/period_depth_recovery_fraction_grid.png`
- `period_depth_fraction_pdf`: `reports/stage5_validation/s56_1k_predetrend_small_pair_200k_chunked_pdo/parameter_space/period_depth_recovery_fraction_grid.pdf`
- `period_depth_fraction_tmag_png`: `reports/stage5_validation/s56_1k_predetrend_small_pair_200k_chunked_pdo/parameter_space/period_depth_recovery_fraction_by_tmag.png`
- `period_depth_fraction_tmag_pdf`: `reports/stage5_validation/s56_1k_predetrend_small_pair_200k_chunked_pdo/parameter_space/period_depth_recovery_fraction_by_tmag.pdf`
- `period_depth_50pct_boundary_csv`: `reports/stage5_validation/s56_1k_predetrend_small_pair_200k_chunked_pdo/parameter_space/period_depth_50pct_boundary_by_tmag.csv`
- `period_depth_50pct_boundary_png`: `reports/stage5_validation/s56_1k_predetrend_small_pair_200k_chunked_pdo/parameter_space/period_depth_50pct_boundary_by_tmag.png`
- `period_depth_50pct_boundary_pdf`: `reports/stage5_validation/s56_1k_predetrend_small_pair_200k_chunked_pdo/parameter_space/period_depth_50pct_boundary_by_tmag.pdf`
- `period_radius_png`: `reports/stage5_validation/s56_1k_predetrend_small_pair_200k_chunked_pdo/parameter_space/period_radius_recovery_by_tmag.png`
- `period_radius_pdf`: `reports/stage5_validation/s56_1k_predetrend_small_pair_200k_chunked_pdo/parameter_space/period_radius_recovery_by_tmag.pdf`
