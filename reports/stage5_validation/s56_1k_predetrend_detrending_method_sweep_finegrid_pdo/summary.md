# Pre-Detrend Detrending Method Sweep

Empirical SNR threshold: `7`.

## Best Rows

| method                  | raw_aperture | n    | median_depth_retention | p16_depth_retention | p84_depth_retention | median_multi_snr | n_snr_ge_threshold | frac_snr_ge_threshold | median_det_trend_ptp | median_trend_reduction_frac | median_tmag |
| ----------------------- | ------------ | ---- | ---------------------- | ------------------- | ------------------- | ---------------- | ------------------ | --------------------- | -------------------- | --------------------------- | ----------- |
| median015_gap05         | Small        | 1000 | 0.4892                 | 0.3003              | 0.7072              | 1.574            | 128                | 0.128                 | 0.01989              | 0.9815                      | 19.44       |
| median020_gap05         | Small        | 1000 | 0.495                  | 0.3041              | 0.7073              | 1.576            | 127                | 0.127                 | 0.03051              | 0.9721                      | 19.44       |
| median025_gap05         | Small        | 1000 | 0.4955                 | 0.3059              | 0.712               | 1.583            | 127                | 0.127                 | 0.03949              | 0.9638                      | 19.44       |
| median04_gap05          | Small        | 1000 | 0.4975                 | 0.3075              | 0.7126              | 1.583            | 126                | 0.126                 | 0.06641              | 0.9406                      | 19.44       |
| oracle_adp03q           | Small        | 1000 | 0.504                  | 0.3136              | 0.7179              | 1.609            | 126                | 0.126                 | 0.07387              | 0.9331                      | 19.44       |
| median03_gap05          | Small        | 1000 | 0.4974                 | 0.3062              | 0.7122              | 1.591            | 125                | 0.125                 | 0.0486               | 0.9551                      | 19.44       |
| spline_q025_gap02       | Small        | 1000 | 0.4946                 | 0.3076              | 0.7119              | 1.582            | 125                | 0.125                 | 0.07077              | 0.9348                      | 19.44       |
| spline_q02_gap02        | Small        | 1000 | 0.4917                 | 0.306               | 0.7091              | 1.583            | 124                | 0.124                 | 0.06894              | 0.937                       | 19.44       |
| current_adp03q          | Small        | 1000 | 0.496                  | 0.308               | 0.7128              | 1.583            | 124                | 0.124                 | 0.07391              | 0.9327                      | 19.44       |
| spline_q015_gap02       | Small        | 1000 | 0.4894                 | 0.3023              | 0.7079              | 1.583            | 123                | 0.123                 | 0.066                | 0.939                       | 19.44       |
| savgol03_p2_gap05       | Small        | 1000 | 0.4863                 | 0.2999              | 0.7033              | 1.572            | 123                | 0.123                 | 0.06592              | 0.9387                      | 19.44       |
| spline_q05_gap02        | Small        | 1000 | 0.498                  | 0.3103              | 0.7128              | 1.586            | 123                | 0.123                 | 0.09264              | 0.9134                      | 19.44       |
| savgol06_p2_gap05       | Small        | 1000 | 0.4956                 | 0.3065              | 0.7093              | 1.593            | 122                | 0.122                 | 0.07425              | 0.9328                      | 19.44       |
| median05_gap05          | Small        | 1000 | 0.4985                 | 0.3073              | 0.7127              | 1.589            | 122                | 0.122                 | 0.08322              | 0.9238                      | 19.44       |
| median06_gap05          | Small        | 1000 | 0.4983                 | 0.3058              | 0.713               | 1.575            | 122                | 0.122                 | 0.1004               | 0.9097                      | 19.44       |
| current_uniform08_gap05 | Small        | 1000 | 0.5002                 | 0.312               | 0.7146              | 1.563            | 117                | 0.117                 | 0.1339               | 0.8738                      | 19.44       |

## Interpretation Aid

- `n_snr_ge_threshold` is the number of rows whose injected-minus-original truth-window signal reaches the SNR threshold after detrending.
- `median_depth_retention` near `1` means the method preserves the injected BATMAN signal depth.
- `median_det_trend_ptp` is the robust 90-10% range of binned original detrended flux; lower is better.
- `median_trend_reduction_frac` compares that trend range to the raw normalized flux; higher is better.
