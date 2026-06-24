# Pre-Detrend Detrending Method Sweep

Empirical SNR threshold: `7`.

## Best Rows

| method                  | raw_aperture | n  | median_depth_retention | p16_depth_retention | p84_depth_retention | median_multi_snr | n_snr_ge_threshold | frac_snr_ge_threshold | median_det_trend_ptp | median_trend_reduction_frac | median_tmag |
| ----------------------- | ------------ | -- | ---------------------- | ------------------- | ------------------- | ---------------- | ------------------ | --------------------- | -------------------- | --------------------------- | ----------- |
| current_adp03q          | Small        | 94 | 0.6485                 | 0.4541              | 0.8139              | 9.194            | 83                 | 0.883                 | 0.03969              | 0.9409                      | 18.54       |
| spline_q05_gap02        | Small        | 94 | 0.6495                 | 0.4541              | 0.8175              | 9.114            | 82                 | 0.8723                | 0.05071              | 0.9221                      | 18.54       |
| spline_q08_gap02        | Small        | 94 | 0.6498                 | 0.4535              | 0.8214              | 8.976            | 80                 | 0.8511                | 0.06954              | 0.898                       | 18.54       |
| current_uniform08_gap05 | Small        | 94 | 0.6529                 | 0.4533              | 0.8339              | 8.946            | 76                 | 0.8085                | 0.08046              | 0.878                       | 18.54       |
| spline_q12_gap05        | Small        | 94 | 0.6535                 | 0.4535              | 0.8342              | 8.841            | 75                 | 0.7979                | 0.1054               | 0.8411                      | 18.54       |
| spline_uniform12_gap05  | Small        | 94 | 0.6539                 | 0.4535              | 0.834               | 8.751            | 75                 | 0.7979                | 0.09889              | 0.8373                      | 18.54       |
| spline_uniform20_gap05  | Small        | 94 | 0.6543                 | 0.4528              | 0.8325              | 8.525            | 72                 | 0.766                 | 0.1467               | 0.7629                      | 18.54       |
| poly3_gap05             | Small        | 94 | 0.6559                 | 0.4512              | 0.8344              | 7.619            | 56                 | 0.5957                | 0.4597               | 0.2084                      | 18.54       |
| poly2_gap05             | Small        | 94 | 0.6559                 | 0.4493              | 0.8345              | 7.569            | 54                 | 0.5745                | 0.5293               | 0.05808                     | 18.54       |
| constant                | Small        | 94 | 0.656                  | 0.447               | 0.8346              | 7.184            | 49                 | 0.5213                | 0.58                 | 0                           | 18.54       |
| current_adp03q          | Primary      | 94 | 0.4059                 | 0.2167              | 0.6823              | 5.2              | 28                 | 0.2979                | 0.04945              | 0.9671                      | 18.54       |
| spline_q05_gap02        | Primary      | 94 | 0.4087                 | 0.2165              | 0.6823              | 5.145            | 26                 | 0.2766                | 0.07748              | 0.9475                      | 18.54       |
| current_uniform08_gap05 | Primary      | 94 | 0.4101                 | 0.2165              | 0.6839              | 4.688            | 26                 | 0.2766                | 0.1482               | 0.91                        | 18.54       |
| spline_q08_gap02        | Primary      | 94 | 0.4088                 | 0.2165              | 0.6825              | 4.851            | 25                 | 0.266                 | 0.1152               | 0.9203                      | 18.54       |
| spline_q12_gap05        | Primary      | 94 | 0.4098                 | 0.2166              | 0.6844              | 4.367            | 25                 | 0.266                 | 0.2106               | 0.8683                      | 18.54       |
| spline_uniform12_gap05  | Primary      | 94 | 0.4108                 | 0.2165              | 0.6838              | 4.471            | 25                 | 0.266                 | 0.2136               | 0.8645                      | 18.54       |
| spline_uniform20_gap05  | Primary      | 94 | 0.4095                 | 0.2165              | 0.6859              | 4.16             | 23                 | 0.2447                | 0.3328               | 0.7952                      | 18.54       |
| current_adp03q          | Large        | 94 | 0.2709                 | 0.09019             | 0.497               | 3.253            | 18                 | 0.1915                | 0.05656              | 0.9743                      | 18.54       |
| spline_q05_gap02        | Large        | 94 | 0.2727                 | 0.09066             | 0.497               | 3.204            | 14                 | 0.1489                | 0.09624              | 0.9523                      | 18.54       |
| spline_q08_gap02        | Large        | 94 | 0.2726                 | 0.09059             | 0.4966              | 3.045            | 12                 | 0.1277                | 0.1501               | 0.9345                      | 18.54       |

## Interpretation Aid

- `n_snr_ge_threshold` is the number of rows whose injected-minus-original truth-window signal reaches the SNR threshold after detrending.
- `median_depth_retention` near `1` means the method preserves the injected BATMAN signal depth.
- `median_det_trend_ptp` is the robust 90-10% range of binned original detrended flux; lower is better.
- `median_trend_reduction_frac` compares that trend range to the raw normalized flux; higher is better.
