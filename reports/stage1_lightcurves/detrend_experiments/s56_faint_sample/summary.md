# Flux Detrend Model Comparison

## Aggregate Metrics

| variant        | cases | failures | median_mad | median_rms | median_finite_q0 | raw_neg_q0 | retained_neg_q0 | retained_neg_q0_frac |
| -------------- | ----- | -------- | ---------- | ---------- | ---------------- | ---------- | --------------- | -------------------- |
| divisive       | 11    | 0        | 0.18996    | 9.0698     | 6137             | 41839      | 41839           | 1                    |
| sub_auto       | 11    | 0        | 0.11964    | 0.13489    | 6137             | 41839      | 41839           | 1                    |
| sub_median     | 11    | 0        | 0.16379    | 0.22272    | 6137             | 41839      | 41839           | 1                    |
| sub_median_abs | 11    | 0        | 0.16379    | 0.22272    | 6137             | 41839      | 41839           | 1                    |

## Injection Preservation

| variant        | injection_cases | median_depth_retention | p16_depth_retention | p84_depth_retention |
| -------------- | --------------- | ---------------------- | ------------------- | ------------------- |
| divisive       | 11              | -0.86365               | -0.9203             | -0.78565            |
| sub_auto       | 11              | 0.97957                | 0.97853             | 0.98002             |
| sub_median     | 11              | -1.4133                | -2.4561             | -1.2456             |
| sub_median_abs | 11              | 1.4133                 | 1.2456              | 2.4561              |

## Highest-Risk Case Rows

| case_id                               | variant        | tmag   | scale_source    | scale   | raw_negative_quality0 | n_det_finite_quality0 | det_mad_quality0 | failed | failure_reason |
| ------------------------------------- | -------------- | ------ | --------------- | ------- | --------------------- | --------------------- | ---------------- | ------ | -------------- |
| tic2053333283_s0056_o120_c1-4_Primary | sub_median     | 19.758 | median          | -1774.2 | 4150                  | 6137                  | 0.70777          | False  |                |
| tic2053333283_s0056_o120_c1-4_Primary | sub_median_abs | 19.758 | median_abs      | 1774.2  | 4150                  | 6137                  | 0.70777          | False  |                |
| tic2053292543_s0056_o120_c1-4_Primary | sub_median     | 19.709 | median          | -1329   | 4120                  | 6137                  | 0.40285          | False  |                |
| tic2053292543_s0056_o120_c1-4_Primary | sub_median_abs | 19.709 | median_abs      | 1329    | 4120                  | 6137                  | 0.40285          | False  |                |
| tic2053333283_s0056_o120_c1-4_Primary | divisive       | 19.758 | local_cotrend   |         | 4150                  | 6137                  | 0.36077          | False  |                |
| tic2053292543_s0056_o120_c1-4_Primary | divisive       | 19.709 | local_cotrend   |         | 4120                  | 6137                  | 0.35953          | False  |                |
| tic2053290277_s0056_o120_c1-4_Primary | divisive       | 19.872 | local_cotrend   |         | 4183                  | 6137                  | 0.29713          | False  |                |
| tic2053290277_s0056_o120_c1-4_Primary | sub_median     | 19.872 | median          | -1863.8 | 4183                  | 6137                  | 0.29251          | False  |                |
| tic2053290277_s0056_o120_c1-4_Primary | sub_median_abs | 19.872 | median_abs      | 1863.8  | 4183                  | 6137                  | 0.29251          | False  |                |
| tic2053292543_s0056_o120_c1-4_Primary | sub_auto       | 19.709 | auto_robust_abs | 1919.5  | 4120                  | 6137                  | 0.27892          | False  |                |
| tic2041357363_s0056_o120_c2-4_Primary | sub_median_abs | 19.801 | median_abs      | 6474.7  | 4256                  | 6137                  | 0.25322          | False  |                |
| tic2041357363_s0056_o120_c2-4_Primary | sub_median     | 19.801 | median          | -6474.7 | 4256                  | 6137                  | 0.25322          | False  |                |
