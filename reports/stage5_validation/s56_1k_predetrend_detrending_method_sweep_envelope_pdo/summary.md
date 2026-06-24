# Pre-Detrend Detrending Method Sweep

Empirical SNR threshold: `7`.

## Best Rows

| method                  | raw_aperture | n    | median_depth_retention | p16_depth_retention | p84_depth_retention | median_multi_snr | n_snr_ge_threshold | frac_snr_ge_threshold | median_det_trend_ptp | median_trend_reduction_frac | median_tmag |
| ----------------------- | ------------ | ---- | ---------------------- | ------------------- | ------------------- | ---------------- | ------------------ | --------------------- | -------------------- | --------------------------- | ----------- |
| oracle_adp03q           | Small        | 1000 | 0.504                  | 0.3136              | 0.7179              | 1.609            | 126                | 0.126                 | 0.07387              | 0.9331                      | 19.44       |
| median03_gap05          | Small        | 1000 | 0.4974                 | 0.3062              | 0.7122              | 1.591            | 125                | 0.125                 | 0.0486               | 0.9551                      | 19.44       |
| current_adp03q          | Small        | 1000 | 0.496                  | 0.308               | 0.7128              | 1.583            | 124                | 0.124                 | 0.07391              | 0.9327                      | 19.44       |
| savgol03_p2_gap05       | Small        | 1000 | 0.4863                 | 0.2999              | 0.7033              | 1.572            | 123                | 0.123                 | 0.06592              | 0.9387                      | 19.44       |
| spline_q05_gap02        | Small        | 1000 | 0.498                  | 0.3103              | 0.7128              | 1.586            | 123                | 0.123                 | 0.09264              | 0.9134                      | 19.44       |
| pctl60_06_gap05         | Small        | 1000 | 0.5006                 | 0.3114              | 0.7154              | 1.578            | 123                | 0.123                 | 0.1227               | 0.8895                      | 19.44       |
| savgol06_p2_gap05       | Small        | 1000 | 0.4956                 | 0.3065              | 0.7093              | 1.593            | 122                | 0.122                 | 0.07425              | 0.9328                      | 19.44       |
| median06_gap05          | Small        | 1000 | 0.4983                 | 0.3058              | 0.713               | 1.575            | 122                | 0.122                 | 0.1004               | 0.9097                      | 19.44       |
| pctl75_03_gap05         | Small        | 1000 | 0.4977                 | 0.3089              | 0.7126              | 1.58             | 122                | 0.122                 | 0.1249               | 0.8862                      | 19.44       |
| savgol10_p2_gap05       | Small        | 1000 | 0.4982                 | 0.3096              | 0.7128              | 1.582            | 121                | 0.121                 | 0.09195              | 0.9194                      | 19.44       |
| median10_gap05          | Small        | 1000 | 0.5005                 | 0.31                | 0.7129              | 1.549            | 120                | 0.12                  | 0.1573               | 0.8488                      | 19.44       |
| pctl75_06_gap05         | Small        | 1000 | 0.4993                 | 0.3109              | 0.7126              | 1.553            | 119                | 0.119                 | 0.1523               | 0.8616                      | 19.44       |
| current_uniform08_gap05 | Small        | 1000 | 0.5002                 | 0.312               | 0.7146              | 1.563            | 117                | 0.117                 | 0.1339               | 0.8738                      | 19.44       |
| pctl90_06_gap05         | Small        | 1000 | 0.4996                 | 0.3099              | 0.715               | 1.539            | 117                | 0.117                 | 0.2016               | 0.811                       | 19.44       |
| pctl75_10_gap05         | Small        | 1000 | 0.4999                 | 0.311               | 0.7126              | 1.524            | 116                | 0.116                 | 0.1903               | 0.8193                      | 19.44       |
| pctl75_20_gap05         | Small        | 1000 | 0.5015                 | 0.3127              | 0.7131              | 1.474            | 112                | 0.112                 | 0.2663               | 0.7493                      | 19.44       |

## Interpretation Aid

- `n_snr_ge_threshold` is the number of rows whose injected-minus-original truth-window signal reaches the SNR threshold after detrending.
- `median_depth_retention` near `1` means the method preserves the injected BATMAN signal depth.
- `median_det_trend_ptp` is the robust 90-10% range of binned original detrended flux; lower is better.
- `median_trend_reduction_frac` compares that trend range to the raw normalized flux; higher is better.
