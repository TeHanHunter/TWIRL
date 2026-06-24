# S56 Duration-Aware Recovery Visualizations

BLS input: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/small_pair_200k/injection_bls_recoveries.csv`
LEO comparison input: `reports/stage5_validation/s56_10k_predetrend_small_pair_200k_review_queue_pdo/review_queue.csv`

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

The empirical publication map now uses four Tmag panels (`<17`, `17-18`, `18-19`, `>19`) and marginalizes over duration/impact parameter within each slice. Grey cells mean the kernel has too little local injection support; the grey dashed curve is the support boundary; the cyan dotted curve is an illustrative fluid Roche limit; the black 50% contour is the empirical recovery boundary. Thin white contours are the local mean total transit duration in minutes for the accepted BATMAN injections, not an independent analytic edge-on duration curve.

The Roche curve is model-dependent because companion radius does not uniquely define companion density. The plotted curve uses a simple bulk-density envelope: rocky bodies at small radii, lower-density gas giants near Jupiter radius, and a mild density increase into the brown-dwarf-radius regime. It should be read as physical context, not a hard vetting threshold.

The earlier small grey point overlay was only a random subsample of injected examples (`faint = missed`, `brighter = recovered`) drawn on top of the smoothed surface. It is removed from the publication-facing map because the support boundary and grid CSV carry the sample-support information more cleanly.

Counts in the paper-style Tmag panels show whether the map is target-distribution-limited or bright-balanced:

| Tmag bin | injections | BLS recovered | fraction |
| --- | ---: | ---: | ---: |
| Tmag < 17 | 275 | 218 | 79.3% |
| 17 <= Tmag < 18 | 589 | 337 | 57.2% |
| 18 <= Tmag < 19 | 1796 | 485 | 27.0% |
| Tmag > 19 | 7340 | 643 | 8.8% |

For bright-balanced physical runs, the `P, R_p` recovery surface should be interpreted as recovery marginalized over the physically allowed duration/depth scatter at fixed period and radius.

## BLS vs LEO

- Injected rows in LEO queue: `1000`
- BLS exact/top-N/harmonic recovered: `177/1000` = `17.7%`
- LEO PC/FP precision relative to BLS recovery: `98.7%`
- LEO PC/FP recall relative to BLS recovery: `42.9%`

LEO is therefore a high-confidence vetting subset of the BLS recoveries, not yet a near-complete proxy for recovery. BLS recovery is the better statement for whether the injected signal is tracked by the search; LEO recovery is a stricter downstream pass/fail.

The metric diagnostic supports targeted WD retuning, not a blind loosening of every LEO cut. BLS-recovered LEO-FA rows sit between the BLS misses and the LEO PC/FP rows in MES/new-MES, in-transit cadence count, and shape diagnostics.

## Artifacts

- `period_radius_boundary_png`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/period_radius_50pct_boundary_by_duration.png`
- `period_radius_boundary_pdf`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/period_radius_50pct_boundary_by_duration.pdf`
- `period_radius_recovery_map_png`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/period_radius_recovery_fraction_by_duration_tmag.png`
- `period_radius_recovery_map_pdf`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/period_radius_recovery_fraction_by_duration_tmag.pdf`
- `period_radius_recovery_map_csv`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/period_radius_recovery_fraction_by_duration_tmag.csv`
- `period_radius_empirical_publication_png`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/period_radius_empirical_recovery_publication.png`
- `period_radius_empirical_publication_pdf`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/period_radius_empirical_recovery_publication.pdf`
- `period_radius_empirical_publication_csv`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/period_radius_empirical_recovery_publication_grid.csv`
- `radius_snr_proxy_png`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/radius_snr_proxy_recovery_by_tmag.png`
- `radius_snr_proxy_pdf`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/radius_snr_proxy_recovery_by_tmag.pdf`
- `duration_radius_floor_png`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/duration_radius_50pct_boundary_floor_maps.png`
- `duration_radius_floor_pdf`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/duration_radius_50pct_boundary_floor_maps.pdf`
- `bls_vs_leo_png`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/bls_vs_leo_recovery_contingency.png`
- `bls_vs_leo_pdf`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/bls_vs_leo_recovery_contingency.pdf`
- `leo_metric_diagnostics_png`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/leo_bls_metric_diagnostics.png`
- `leo_metric_diagnostics_pdf`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/leo_bls_metric_diagnostics.pdf`
- `leo_metric_summary_csv`: `reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/leo_bls_metric_summary.csv`
