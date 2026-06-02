# Flux Detrend Model Comparison

## Aggregate Metrics

| variant        | cases | failures | median_mad | median_rms | median_finite_q0 | raw_neg_q0 | retained_neg_q0 | retained_neg_q0_frac |
| -------------- | ----- | -------- | ---------- | ---------- | ---------------- | ---------- | --------------- | -------------------- |
| divisive       | 33    | 0        | 0.27078    | 14.62      | 6137             | 117915     | 117915          | 1                    |
| sub_auto       | 33    | 0        | 0.15811    | 0.2447     | 6137             | 117915     | 117915          | 1                    |
| sub_auto_adp03 | 33    | 0        | 0.12911    | 0.15688    | 6137             | 117915     | 117915          | 1                    |
| sub_median     | 33    | 1        | 0.23878    | 0.39186    | 6137             | 117915     | 117915          | 1                    |
| sub_median_abs | 33    | 1        | 0.23878    | 0.39186    | 6137             | 117915     | 117915          | 1                    |

## Injection Preservation

| variant        | injection_cases | median_depth_retention | p16_depth_retention | p84_depth_retention |
| -------------- | --------------- | ---------------------- | ------------------- | ------------------- |
| divisive       | 33              | -0.89458               | -1.1354             | 0.21279             |
| sub_auto       | 33              | 0.99923                | 0.99765             | 0.99962             |
| sub_auto_adp03 | 33              | 0.9795                 | 0.97786             | 0.98002             |
| sub_median     | 33              | -1.8824                | -3.0425             | -1.198              |
| sub_median_abs | 33              | 1.8824                 | 1.2358              | 3.0425              |

## Highest-Risk Case Rows

| case_id                             | variant        | tmag   | scale_source  | scale   | raw_negative_quality0 | n_det_finite_quality0 | det_mad_quality0 | failed | failure_reason        |
| ----------------------------------- | -------------- | ------ | ------------- | ------- | --------------------- | --------------------- | ---------------- | ------ | --------------------- |
| tic2041357363_s0056_o120_c2-4_Large | sub_median     | 19.801 | median        | -1059.9 | 3134                  | 6137                  | 14.635           | True   | unstable_relative_mad |
| tic2041357363_s0056_o120_c2-4_Large | sub_median_abs | 19.801 | median_abs    | 1059.9  | 3134                  | 6137                  | 14.635           | True   | unstable_relative_mad |
| tic2053290277_s0056_o120_c1-4_Small | sub_median     | 19.872 | median        | -28.874 | 3266                  | 6137                  | 6.1463           | False  |                       |
| tic2053290277_s0056_o120_c1-4_Small | sub_median_abs | 19.872 | median_abs    | 28.874  | 3266                  | 6137                  | 6.1463           | False  |                       |
| tic2053333283_s0056_o120_c1-4_Small | sub_median     | 19.758 | median        | -44.642 | 3363                  | 6137                  | 5.875            | False  |                       |
| tic2053333283_s0056_o120_c1-4_Small | sub_median_abs | 19.758 | median_abs    | 44.642  | 3363                  | 6137                  | 5.875            | False  |                       |
| tic1884514391_s0056_o120_c4-2_Small | sub_median_abs | 19.783 | median_abs    | 37.875  | 3291                  | 6137                  | 4.247            | False  |                       |
| tic1884514391_s0056_o120_c4-2_Small | sub_median     | 19.783 | median        | -37.875 | 3291                  | 6137                  | 4.247            | False  |                       |
| tic2014550155_s0056_o120_c3-2_Small | sub_median_abs | 19.302 | median_abs    | 161.68  | 3590                  | 6137                  | 1.553            | False  |                       |
| tic2014550155_s0056_o120_c3-2_Small | sub_median     | 19.302 | median        | -161.68 | 3590                  | 6137                  | 1.553            | False  |                       |
| tic2053333283_s0056_o120_c1-4_Small | divisive       | 19.758 | local_cotrend |         | 3363                  | 6137                  | 1.4156           | False  |                       |
| tic2053290277_s0056_o120_c1-4_Small | divisive       | 19.872 | local_cotrend |         | 3266                  | 6137                  | 1.0044           | False  |                       |
