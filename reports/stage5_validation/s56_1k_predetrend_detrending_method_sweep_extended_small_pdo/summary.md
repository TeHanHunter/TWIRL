# Pre-Detrend Detrending Method Sweep

Empirical SNR threshold: `7`.

## Best Rows

| method                  | raw_aperture | n    | median_depth_retention | p16_depth_retention | p84_depth_retention | median_multi_snr | n_snr_ge_threshold | frac_snr_ge_threshold | median_det_trend_ptp | median_trend_reduction_frac | median_tmag |
| ----------------------- | ------------ | ---- | ---------------------- | ------------------- | ------------------- | ---------------- | ------------------ | --------------------- | -------------------- | --------------------------- | ----------- |
| oracle_adp03q           | Small        | 1000 | 0.504                  | 0.3136              | 0.7179              | 1.609            | 126                | 0.126                 | 0.07387              | 0.9331                      | 19.44       |
| oracle_q05_gap02        | Small        | 1000 | 0.5043                 | 0.3133              | 0.7162              | 1.601            | 125                | 0.125                 | 0.09286              | 0.9138                      | 19.44       |
| current_adp03q          | Small        | 1000 | 0.496                  | 0.308               | 0.7128              | 1.583            | 124                | 0.124                 | 0.07391              | 0.9327                      | 19.44       |
| spline_q05_gap02        | Small        | 1000 | 0.498                  | 0.3103              | 0.7128              | 1.586            | 123                | 0.123                 | 0.09264              | 0.9134                      | 19.44       |
| savgol06_p2_gap05       | Small        | 1000 | 0.4956                 | 0.3065              | 0.7093              | 1.593            | 122                | 0.122                 | 0.07425              | 0.9328                      | 19.44       |
| median06_gap05          | Small        | 1000 | 0.4983                 | 0.3058              | 0.713               | 1.575            | 122                | 0.122                 | 0.1004               | 0.9097                      | 19.44       |
| savgol10_p2_gap05       | Small        | 1000 | 0.4982                 | 0.3096              | 0.7128              | 1.582            | 121                | 0.121                 | 0.09195              | 0.9194                      | 19.44       |
| median10_gap05          | Small        | 1000 | 0.5005                 | 0.31                | 0.7129              | 1.549            | 120                | 0.12                  | 0.1573               | 0.8488                      | 19.44       |
| current_uniform08_gap05 | Small        | 1000 | 0.5002                 | 0.312               | 0.7146              | 1.563            | 117                | 0.117                 | 0.1339               | 0.8738                      | 19.44       |
| savgol20_p2_gap05       | Small        | 1000 | 0.5005                 | 0.3122              | 0.716               | 1.533            | 117                | 0.117                 | 0.1645               | 0.8482                      | 19.44       |
| median20_gap05          | Small        | 1000 | 0.5018                 | 0.3137              | 0.713               | 1.511            | 116                | 0.116                 | 0.2303               | 0.7853                      | 19.44       |
| poly3_gap05             | Small        | 1000 | 0.5036                 | 0.3144              | 0.7156              | 1.184            | 89                 | 0.089                 | 0.8396               | 0.2185                      | 19.44       |
| poly2_gap05             | Small        | 1000 | 0.5035                 | 0.3146              | 0.7155              | 1.128            | 87                 | 0.087                 | 0.9953               | 0.0653                      | 19.44       |
| constant                | Small        | 1000 | 0.5036                 | 0.3145              | 0.7153              | 1.067            | 81                 | 0.081                 | 1.07                 | 0                           | 19.44       |

## Interpretation Aid

- `n_snr_ge_threshold` is the number of rows whose injected-minus-original truth-window signal reaches the SNR threshold after detrending.
- `median_depth_retention` near `1` means the method preserves the injected BATMAN signal depth.
- `median_det_trend_ptp` is the robust 90-10% range of binned original detrended flux; lower is better.
- `median_trend_reduction_frac` compares that trend range to the raw normalized flux; higher is better.
