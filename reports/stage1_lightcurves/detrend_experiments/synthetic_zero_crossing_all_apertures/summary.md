# Flux Detrend Model Comparison

## Aggregate Metrics

| variant        | cases | failures | median_mad | median_rms | median_finite_q0 | raw_neg_q0 | retained_neg_q0 | retained_neg_q0_frac |
| -------------- | ----- | -------- | ---------- | ---------- | ---------------- | ---------- | --------------- | -------------------- |
| divisive       | 36    | 15       | 5.2867     | 69.377     | 10842            | 119418     | 119418          | 1                    |
| sub_auto       | 36    | 0        | 0.96011    | 0.96364    | 10842            | 119418     | 119418          | 1                    |
| sub_median     | 36    | 15       | 4.8648     | 4.9142     | 10842            | 119418     | 119418          | 1                    |
| sub_median_abs | 36    | 15       | 4.8648     | 4.9142     | 10842            | 119418     | 119418          | 1                    |

## Injection Preservation

| variant        | injection_cases | median_depth_retention | p16_depth_retention | p84_depth_retention |
| -------------- | --------------- | ---------------------- | ------------------- | ------------------- |
| divisive       | 36              | 1.4648                 | 0.70161             | 5.9135              |
| sub_auto       | 36              | 0.97916                | 0.97028             | 1.0092              |
| sub_median     | 36              | 2.9376                 | 1.1028              | 25.507              |
| sub_median_abs | 36              | 4.887                  | 1.1552              | 47.711              |

## Highest-Risk Case Rows

| case_id              | variant        | tmag   | scale_source | scale     | raw_negative_quality0 | n_det_finite_quality0 | det_mad_quality0 | failed | failure_reason        |
| -------------------- | -------------- | ------ | ------------ | --------- | --------------------- | --------------------- | ---------------- | ------ | --------------------- |
| synthetic001_Primary | sub_median     | 22.553 | median       | 0.0098966 | 5426                  | 10854                 | 5845.1           | True   | unstable_relative_mad |
| synthetic001_Primary | sub_median_abs | 22.553 | median_abs   | 0.0098966 | 5426                  | 10854                 | 5845.1           | True   | unstable_relative_mad |
| synthetic001_Small   | sub_median     | 22.553 | median       | 0.067801  | 5425                  | 10854                 | 1135.5           | True   | unstable_relative_mad |
| synthetic001_Small   | sub_median_abs | 22.553 | median_abs   | 0.067801  | 5425                  | 10854                 | 1135.5           | True   | unstable_relative_mad |
| synthetic001_Large   | sub_median     | 22.553 | median       | -0.17904  | 5431                  | 10854                 | 536.86           | True   | unstable_relative_mad |
| synthetic001_Large   | sub_median_abs | 22.553 | median_abs   | 0.17904   | 5431                  | 10854                 | 536.86           | True   | unstable_relative_mad |
| synthetic000_Primary | sub_median     | 23.253 | median       | 0.38885   | 5399                  | 10853                 | 145.22           | True   | unstable_relative_mad |
| synthetic000_Primary | sub_median_abs | 23.253 | median_abs   | 0.38885   | 5399                  | 10853                 | 145.22           | True   | unstable_relative_mad |
| synthetic000_Large   | sub_median     | 23.253 | median       | 0.84625   | 5388                  | 10853                 | 110.23           | True   | unstable_relative_mad |
| synthetic000_Large   | sub_median_abs | 23.253 | median_abs   | 0.84625   | 5388                  | 10853                 | 110.23           | True   | unstable_relative_mad |
| synthetic000_Small   | sub_median     | 23.253 | median       | 0.98268   | 5361                  | 10853                 | 75.83            | True   | unstable_relative_mad |
| synthetic000_Small   | sub_median_abs | 23.253 | median_abs   | 0.98268   | 5361                  | 10853                 | 75.83            | True   | unstable_relative_mad |
