# S56 Injection Recovery Parameter Space

Input recovery table: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/small_pair_200k/injection_bls_recoveries.csv`

- Rows plotted: `10000`
- Strict top-rank BLS recoveries: `1273/10000` = `12.7%`
- Exact/top-N/harmonic matches: `1683/10000` = `16.8%`

## Tmag Bins

| label        | n    | strict_top1_n | strict_top1_frac | any_exact_or_harmonic_n | any_exact_or_harmonic_frac | median_period_d | median_depth_pct |
| ------------ | ---- | ------------- | ---------------- | ----------------------- | -------------------------- | --------------- | ---------------- |
| Tmag [0,16)  | 90   | 70            | 0.778            | 75                      | 0.833                      | 1.38            | 25               |
| Tmag [16,17) | 185  | 130           | 0.703            | 143                     | 0.773                      | 1.28            | 29.6             |
| Tmag [17,18) | 589  | 305           | 0.518            | 337                     | 0.572                      | 1.11            | 28.8             |
| Tmag [18,19) | 1796 | 400           | 0.223            | 485                     | 0.27                       | 1.3             | 28.3             |
| Tmag [19,20) | 6863 | 358           | 0.0522           | 618                     | 0.09                       | 1.25            | 29.9             |
| Tmag [20,30) | 477  | 10            | 0.021            | 25                      | 0.0524                     | 1.12            | 30.7             |

## Interpretation

The recovery boundary should be read conditionally on Tmag and cadence coverage. A single period-depth panel mixes bright and faint WDs, so the Tmag-sliced recovery-fraction plot is the better diagnostic for the gradual BLS boundary.

## Artifacts

- `period_depth_png`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/parameter_space/period_depth_recovery_by_tmag.png`
- `period_depth_pdf`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/parameter_space/period_depth_recovery_by_tmag.pdf`
- `period_depth_fraction_png`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/parameter_space/period_depth_recovery_fraction_grid.png`
- `period_depth_fraction_pdf`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/parameter_space/period_depth_recovery_fraction_grid.pdf`
- `period_depth_fraction_tmag_png`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/parameter_space/period_depth_recovery_fraction_by_tmag.png`
- `period_depth_fraction_tmag_pdf`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/parameter_space/period_depth_recovery_fraction_by_tmag.pdf`
- `period_depth_50pct_boundary_csv`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/parameter_space/period_depth_50pct_boundary_by_tmag.csv`
- `period_depth_50pct_boundary_png`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/parameter_space/period_depth_50pct_boundary_by_tmag.png`
- `period_depth_50pct_boundary_pdf`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/parameter_space/period_depth_50pct_boundary_by_tmag.pdf`
- `period_depth_smoothed_tmag_png`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/parameter_space/period_depth_smoothed_recovery_by_tmag.png`
- `period_depth_smoothed_tmag_pdf`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/parameter_space/period_depth_smoothed_recovery_by_tmag.pdf`
- `period_depth_smoothed_boundary_csv`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/parameter_space/period_depth_smoothed_boundary_by_tmag.csv`
- `period_depth_smoothed_50pct_boundary_png`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/parameter_space/period_depth_smoothed_50pct_boundary_by_tmag.png`
- `period_depth_smoothed_50pct_boundary_pdf`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/parameter_space/period_depth_smoothed_50pct_boundary_by_tmag.pdf`
- `period_radius_png`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/parameter_space/period_radius_recovery_by_tmag.png`
- `period_radius_pdf`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/parameter_space/period_radius_recovery_by_tmag.pdf`
