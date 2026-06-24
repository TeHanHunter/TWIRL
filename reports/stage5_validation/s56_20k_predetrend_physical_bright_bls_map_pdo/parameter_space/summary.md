# S56 Injection Recovery Parameter Space

Input recovery table: `reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/small_pair_200k/injection_bls_recoveries.csv`

- Rows plotted: `20000`
- Strict top-rank BLS recoveries: `7227/20000` = `36.1%`
- Exact/top-N/harmonic matches: `7910/20000` = `39.6%`

## Tmag Bins

| label        | n    | strict_top1_n | strict_top1_frac | any_exact_or_harmonic_n | any_exact_or_harmonic_frac | median_period_d | median_depth_pct |
| ------------ | ---- | ------------- | ---------------- | ----------------------- | -------------------------- | --------------- | ---------------- |
| Tmag [0,16)  | 1690 | 1184          | 0.701            | 1238                    | 0.733                      | 1.25            | 16.7             |
| Tmag [16,17) | 3330 | 1821          | 0.547            | 1916                    | 0.575                      | 1.3             | 20.2             |
| Tmag [17,18) | 4994 | 2118          | 0.424            | 2273                    | 0.455                      | 1.28            | 18.4             |
| Tmag [18,19) | 5026 | 1429          | 0.284            | 1625                    | 0.323                      | 1.2             | 17.6             |
| Tmag [19,20) | 4646 | 647           | 0.139            | 817                     | 0.176                      | 1.25            | 18               |
| Tmag [20,30) | 314  | 28            | 0.0892           | 41                      | 0.131                      | 0.965           | 24.7             |

## Interpretation

The recovery boundary should be read conditionally on Tmag and cadence coverage. A single period-depth panel mixes bright and faint WDs, so the Tmag-sliced recovery-fraction plot is the better diagnostic for the gradual BLS boundary.

## Artifacts

- `period_depth_png`: `reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/parameter_space/period_depth_recovery_by_tmag.png`
- `period_depth_pdf`: `reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/parameter_space/period_depth_recovery_by_tmag.pdf`
- `period_depth_fraction_png`: `reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/parameter_space/period_depth_recovery_fraction_grid.png`
- `period_depth_fraction_pdf`: `reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/parameter_space/period_depth_recovery_fraction_grid.pdf`
- `period_depth_fraction_tmag_png`: `reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/parameter_space/period_depth_recovery_fraction_by_tmag.png`
- `period_depth_fraction_tmag_pdf`: `reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/parameter_space/period_depth_recovery_fraction_by_tmag.pdf`
- `period_depth_50pct_boundary_csv`: `reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/parameter_space/period_depth_50pct_boundary_by_tmag.csv`
- `period_depth_50pct_boundary_png`: `reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/parameter_space/period_depth_50pct_boundary_by_tmag.png`
- `period_depth_50pct_boundary_pdf`: `reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/parameter_space/period_depth_50pct_boundary_by_tmag.pdf`
- `period_depth_smoothed_tmag_png`: `reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/parameter_space/period_depth_smoothed_recovery_by_tmag.png`
- `period_depth_smoothed_tmag_pdf`: `reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/parameter_space/period_depth_smoothed_recovery_by_tmag.pdf`
- `period_depth_smoothed_boundary_csv`: `reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/parameter_space/period_depth_smoothed_boundary_by_tmag.csv`
- `period_depth_smoothed_50pct_boundary_png`: `reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/parameter_space/period_depth_smoothed_50pct_boundary_by_tmag.png`
- `period_depth_smoothed_50pct_boundary_pdf`: `reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/parameter_space/period_depth_smoothed_50pct_boundary_by_tmag.pdf`
- `period_radius_png`: `reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/parameter_space/period_radius_recovery_by_tmag.png`
- `period_radius_pdf`: `reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/parameter_space/period_radius_recovery_by_tmag.pdf`
