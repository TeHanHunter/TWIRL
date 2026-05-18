# Flux Detrend Model Comparison

## Aggregate Metrics

| variant        | cases | failures | median_mad | median_rms | median_finite_q0 | raw_neg_q0 | retained_neg_q0 | retained_neg_q0_frac |
| -------------- | ----- | -------- | ---------- | ---------- | ---------------- | ---------- | --------------- | -------------------- |
| divisive       | 8     | 3        | 4.7602     | 247.47     | 10848            | 26312      | 26312           | 1                    |
| sub_auto       | 8     | 0        | 0.95588    | 0.96281    | 10848            | 26312      | 26312           | 1                    |
| sub_median     | 8     | 3        | 4.5737     | 4.6033     | 10848            | 26312      | 26312           | 1                    |
| sub_median_abs | 8     | 3        | 4.5737     | 4.6033     | 10848            | 26312      | 26312           | 1                    |

## Injection Preservation

| variant        | injection_cases | median_depth_retention | p16_depth_retention | p84_depth_retention |
| -------------- | --------------- | ---------------------- | ------------------- | ------------------- |
| divisive       | 8               | 3.1955                 | 1.1891              | 8.1113              |
| sub_auto       | 8               | 0.97895                | 0.97171             | 0.9818              |
| sub_median     | 8               | 4.7154                 | 1.1789              | 43.948              |
| sub_median_abs | 8               | 4.7154                 | 1.1789              | 43.948              |

## Highest-Risk Case Rows

| case_id              | variant        | tmag   | scale_source  | scale   | raw_negative_quality0 | n_det_finite_quality0 | det_mad_quality0 | failed | failure_reason        |
| -------------------- | -------------- | ------ | ------------- | ------- | --------------------- | --------------------- | ---------------- | ------ | --------------------- |
| synthetic000_Primary | sub_median     | 23.253 | median        | 0.38885 | 5399                  | 10853                 | 145.22           | True   | unstable_relative_mad |
| synthetic000_Primary | sub_median_abs | 23.253 | median_abs    | 0.38885 | 5399                  | 10853                 | 145.22           | True   | unstable_relative_mad |
| synthetic001_Primary | sub_median     | 22.153 | median        | 1.2539  | 5352                  | 10859                 | 48.915           | True   | unstable_relative_mad |
| synthetic001_Primary | sub_median_abs | 22.153 | median_abs    | 1.2539  | 5352                  | 10859                 | 48.915           | True   | unstable_relative_mad |
| synthetic000_Primary | divisive       | 23.253 | local_cotrend |         | 5399                  | 10853                 | 18.579           | True   | unstable_relative_mad |
| synthetic002_Primary | sub_median     | 21.053 | median        | 2.9814  | 5180                  | 10858                 | 18.56            | True   | unstable_relative_mad |
| synthetic002_Primary | sub_median_abs | 21.053 | median_abs    | 2.9814  | 5180                  | 10858                 | 18.56            | True   | unstable_relative_mad |
| synthetic001_Primary | divisive       | 22.153 | local_cotrend |         | 5352                  | 10859                 | 15.415           | True   | unstable_relative_mad |
| synthetic002_Primary | divisive       | 21.053 | local_cotrend |         | 5180                  | 10858                 | 12.474           | True   | unstable_relative_mad |
| synthetic003_Primary | divisive       | 19.953 | local_cotrend |         | 4703                  | 10852                 | 6.6434           | False  |                       |
| synthetic003_Primary | sub_median     | 19.953 | median        | 8.7241  | 4703                  | 10852                 | 6.2474           | False  |                       |
| synthetic003_Primary | sub_median_abs | 19.953 | median_abs    | 8.7241  | 4703                  | 10852                 | 6.2474           | False  |                       |
